#pragma once
// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025
 * Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file TimeSeriesObserver.hpp
 * @brief Lightweight, plotting-oriented time series recorder for Evolution runs.
 *
 * The TimeSeriesObserver is intentionally distinct from DiagnosticsObserver:
 *
 * - DiagnosticsObserver is designed for auditing, unit contracts, regression testing,
 *   and structured per-driver snapshots (JSONL with rich metadata).
 *
 * - TimeSeriesObserver is designed for *plotting and analysis*:
 *   it records a narrow set of scalar quantities (e.g. evolved DOF like temperature
 *   and spin frequency) as a function of time, with minimal overhead and minimal
 *   repetition.
 *
 * Design goals:
 *  - Deterministic column ordering.
 *  - Simple, tooling-friendly output (CSV by default).
 *  - No "contract", warnings, or verbose metadata in each row.
 *  - Observer must never mutate state or context (read-only).
 *
 * Typical output row:
 *   t, sample_index, Tinf_K, Omega_rad_s, dLnTinf_dt_1_s, dOmega_dt_rad_s2, ...
 *
 * Column selection:
 *  - The observer supports an explicit list of columns (keys) to record.
 *  - Columns can be sourced from:
 *      (a) state-only accessors (e.g. ThermalState::Tinf()),
 *      (b) registered driver diagnostics providers (IDriverDiagnostics),
 *      (c) future: derived/computed columns (via callback hooks).
 *
 * NOTE:
 *  - This header defines the interface and configuration only.
 *  - The implementation (cpp) is expected to:
 *      - open/close files in OnStart/OnFinish,
 *      - write a single header line,
 *      - write one row per recorded sample,
 *      - enforce stable ordering and minimal allocations.
 *
 * @ingroup Evolution
 */

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CompactStar/Physics/Evolution/Observers/IObserver.hpp"

// Optional: allow pulling from driver diagnostics providers (same concept as DiagnosticsObserver).
#include "CompactStar/Physics/Driver/Diagnostics/DriverDiagnostics.hpp"
#include "CompactStar/Physics/Evolution/Diagnostics/DiagnosticCatalog.hpp"

#include "Zaki/String/Directory.hpp"

namespace CompactStar::Physics::Evolution
{
class StateVector;
class DriverContext;
} // namespace CompactStar::Physics::Evolution

namespace CompactStar::Physics::Evolution::Observers
{

/**
 * @enum RecordCadence
 * @brief Scheduling policy for time-based recording in TimeSeriesObserver.
 *
 * This enum controls how the observer advances its internal
 * `next_time_trigger_` after a record event is emitted.
 *
 * Important:
 * - RecordCadence affects *only* the observer’s write decision.
 * - It does NOT influence the integrator stepping or sampling cadence.
 * - The observer receives samples whenever EvolutionSystem calls
 *   `NotifySample(t, ...)` and decides independently whether to record.
 */
enum class RecordCadence : std::uint8_t
{
	/**
	 * @brief Linear time cadence.
	 *
	 * A record is triggered whenever:
	 *   \f[
	 *     t \ge t_{\mathrm{next}}
	 *   \f]
	 *
	 * After recording, the next trigger time is advanced as:
	 *   \f[
	 *     t_{\mathrm{next}} \leftarrow t_{\mathrm{next}} + \Delta t
	 *   \f]
	 *
	 * where \f$\Delta t = \texttt{record\_every\_dt}\f$.
	 *
	 * This corresponds to uniformly spaced output in physical time.
	 */
	LinearDt = 0,

	/**
	 * @brief Logarithmic (geometric) time cadence.
	 *
	 * A record is triggered whenever:
	 *   \f[
	 *     t \ge t_{\mathrm{next}}
	 *   \f]
	 *
	 * After recording, the next trigger time is advanced multiplicatively:
	 *   \f[
	 *     t_{\mathrm{next}} \leftarrow t_{\mathrm{next}} \times q
	 *   \f]
	 *
	 * where the geometric ratio \f$q\f$ is determined from:
	 *   - `log_q` if explicitly provided and > 1, otherwise
	 *   - `samples_per_decade`, via
	 *     \f[
	 *       q = 10^{1 / N_{\mathrm{decade}}}
	 *     \f]
	 *
	 * This mode is appropriate for long-term evolutionary problems
	 * (e.g., thermal cooling curves) where dynamics span many decades
	 * in time and logarithmic spacing is desired.
	 */
	LogTime = 1
};

/**
 * @brief Minimal time-series recorder for selected scalar quantities.
 *
 * The TimeSeriesObserver writes a single table-like stream where each row
 * corresponds to a recorded sample time.
 *
 * It is meant to produce a compact dataset suitable for:
 *  - quick plotting (Python, MATLAB, gnuplot),
 *  - quick inspection,
 *  - numerical regression comparisons at the level of *state variables*.
 *
 * It is *not* meant to replace DiagnosticsObserver for deep debugging,
 * driver contracts, or per-driver introspection.
 */
class TimeSeriesObserver final : public IObserver
{
  public:
	/**
	 * @brief Output format for the time series.
	 *
	 * CSV is the default because it is ubiquitous and easy to inspect.
	 * TSV is occasionally useful when units/descriptions include commas.
	 *
	 * JSONL output is intentionally not supported here because JSONL is
	 * already covered by DiagnosticsObserver; the goal here is compact tables.
	 */
	enum class OutputFormat
	{
		CSV, ///< Comma-separated values (default).
		TSV	 ///< Tab-separated values.
	};

	/**
	 * @brief Reference to a scalar in the diagnostics catalog.
	 *
	 * This is the schema-driven version of DriverScalar lookup:
	 * producer = driver DiagnosticsName() (packet producer string)
	 * key      = scalar key within that producer's scalars
	 */
	struct CatalogScalarRef
	{
		std::string producer;
		std::string key;
	};

	/**
	 * @brief Column source kind.
	 *
	 * A "column" is a named scalar recorded in the output table.
	 *
	 * - BuiltinState: value obtained directly from the Evolution::StateVector
	 *   and/or Evolution::DriverContext via known accessors (e.g. Tinf, Omega).
	 *
	 * - DriverScalar: value fetched from a driver diagnostics provider by key.
	 *   This uses IDriverDiagnostics::DiagnoseSnapshot() to populate a packet,
	 *   then extracts the requested scalar key.
	 *
	 * This separation keeps the implementation simple and avoids pushing
	 * state-specific knowledge into drivers.
	 */
	enum class ColumnSource
	{
		BuiltinState,
		DriverScalar
	};

	/**
	 * @brief Column specification for TimeSeriesObserver output.
	 *
	 * A column is defined by:
	 *  - a stable column name (header label),
	 *  - how to obtain the value at each sample,
	 *  - optional provenance metadata (not repeated per row).
	 *
	 * Notes:
	 *  - `key` should be tooling-friendly and stable (snake_case recommended).
	 *  - `unit` and `description` are written once (optional) via a sidecar file
	 *    if enabled; they are not repeated in the table itself.
	 *  - For ColumnSource::DriverScalar, `catalog_ref` must be set.
	 */
	struct Column
	{
		/// Column header label (appears in the output header row).
		std::string key;

		/// Source kind for this column.
		ColumnSource source = ColumnSource::BuiltinState;

		/// Optional unit string for documentation (not repeated per row).
		std::string unit;

		/// Optional description for documentation (not repeated per row).
		std::string description;

		/// ---- DriverScalar fields ----

		// /**
		//  * @brief Driver diagnostics producer name to query (matches DiagnosticsName()).
		//  *
		//  * Only used when source == ColumnSource::DriverScalar.
		//  */
		// std::string driver_name;

		// /**
		//  * @brief Scalar key inside the driver's DiagnosticPacket to extract.
		//  *
		//  * Only used when source == ColumnSource::DriverScalar.
		//  */
		// std::string driver_key;

		/**
		 * @brief Reference to a scalar in the diagnostics catalog.
		 *
		 */
		CatalogScalarRef catalog_ref;

		/// ---- BuiltinState fields ----

		/**
		 * @brief Built-in identifier for common state scalars.
		 *
		 * Only used when source == ColumnSource::BuiltinState.
		 *
		 * This keeps common plots trivial:
		 *   - time, Tinf, Omega, etc.
		 *
		 * For arbitrary expressions, add a callback mechanism in cpp.
		 */
		enum class Builtin
		{
			Time,		 ///< t [s], usually emitted automatically; include only if we want it as a column.
			SampleIndex, ///< sample_index (monotonic output counter).
			StepIndex,	 ///< integrator step_index if provided (0 if unknown).
			Tinf_K,		 ///< redshifted internal temperature [K] from ThermalState.
			Omega_rad_s, ///< spin angular frequency [rad/s] from SpinState (if present).
						 // Future common fields (optional):
						 // dLnTinf_dt_1_s,
						 // dOmega_dt_rad_s2,
		};

		/// Which built-in value to record (only meaningful for BuiltinState columns).
		Builtin builtin = Builtin::Tinf_K;
	};

	// /**
	//  * @brief Configuration for TimeSeriesObserver.
	//  */
	// struct Options
	// {
	// 	/// Output file path for the table (CSV/TSV).
	// 	Zaki::String::Directory output_path = "timeseries.csv";

	// 	/// Output format (CSV by default).
	// 	OutputFormat format = OutputFormat::CSV;

	// 	/// If true, append to existing output file; otherwise truncate.
	// 	bool append = false;

	// 	/// If true, record at OnStart (t=t0) before first OnSample.
	// 	bool record_at_start = true;

	// 	/**
	// 	 * @brief Record every N OnSample calls (step-based scheduling).
	// 	 *
	// 	 * - If 0, step-based triggering is disabled.
	// 	 * - If both step-based and time-based are enabled, either trigger records.
	// 	 */
	// 	std::size_t record_every_n_samples = 1;

	// 	/**
	// 	 * @brief Record every dt in simulation time (time-based scheduling).
	// 	 *
	// 	 * - If <= 0, time-based triggering is disabled.
	// 	 * - If enabled, the observer maintains an internal "next_time_trigger".
	// 	 */
	// 	double record_every_dt = 0.0;

	// 	/**
	// 	 * @brief If true, write a one-line header row containing column keys.
	// 	 *
	// 	 * For CSV/TSV outputs, the header is strongly recommended.
	// 	 */
	// 	bool write_header = true;

	// 	/**
	// 	 * @brief If true, also emit a sidecar metadata file describing columns.
	// 	 *
	// 	 * This keeps the table compact while preserving units and descriptions.
	// 	 * A typical sidecar file could be: "<output_path>.meta.json" or ".meta.txt"
	// 	 * depending on implementation.
	// 	 *
	// 	 * NOTE: This header does not prescribe the exact sidecar format; the .cpp
	// 	 * should implement a deterministic, low-dependency approach.
	// 	 */
	// 	bool write_sidecar_metadata = true;

	// 	/**
	// 	 * @brief Numeric formatting precision for floating point output.
	// 	 *
	// 	 * Use sufficiently high precision for regression comparisons.
	// 	 * Suggested default is 17 (full double precision).
	// 	 */
	// 	int float_precision = 17;

	// 	/**
	// 	 * @brief Columns to record, in the order they will appear in the table.
	// 	 *
	// 	 * The observer will not reorder; this is the authoritative ordering.
	// 	 */
	// 	std::vector<Column> columns;

	// 	/**
	// 	 * @brief If true, build columns automatically from the diagnostics catalog.
	// 	 *
	// 	 * This is the schema-driven path: no hardcoded column list in main.
	// 	 * If enabled, the observer will append (or replace) Options::columns
	// 	 * based on catalog profiles and/or explicit key lists.
	// 	 */
	// 	bool use_catalog = false;

	// 	/**
	// 	 * @brief Optional: path to the catalog JSON (diagnostics.catalog.json).
	// 	 *
	// 	 * Used when use_catalog==true and no in-memory catalog is supplied.
	// 	 */
	// 	Zaki::String::Directory catalog_path = "";

	// 	/**
	// 	 * @brief Optional: which catalog profile(s) to use when auto-building columns.
	// 	 *
	// 	 * Examples:
	// 	 *   "timeseries_default"
	// 	 *   "regression_minimal"
	// 	 *
	// 	 * If empty, you can fall back to "timeseries_default" if present,
	// 	 * otherwise use all scalars or none (cpp policy).
	// 	 */
	// 	std::vector<std::string> catalog_profiles;

	// 	/**
	// 	 * @brief If true, prepend/ensure common builtin columns (time, sample_index).
	// 	 *
	// 	 * Useful when catalog-driven columns only cover driver scalars.
	// 	 */
	// 	bool include_builtin_time = true;
	// 	bool include_builtin_sample_index = true;

	// 	/**
	// 	 * @brief If true, use buffered I/O for output stream.
	// 	 *
	// 	 */
	// 	bool buffered_io = true;

	// 	/**
	// 	 * @brief Buffer flush size in bytes when using buffered I/O.
	// 	 *
	// 	 */
	// 	std::size_t flush_bytes = 1 << 20;
	// };

	/**
	 * @brief Configuration for TimeSeriesObserver.
	 *
	 * This struct configures **what** the observer writes (columns, formatting, output files)
	 * and **when** it writes (sample-based and/or time-based scheduling).
	 *
	 * Key design point:
	 * - The observer receives samples whenever Evolution calls `NotifySample(t, ...)`.
	 * - The observer then applies these Options to decide whether to *record* that sample.
	 *
	 * Scheduling is controlled by two independent triggers:
	 *  1) Sample-count trigger (record_every_n_samples)
	 *  2) Time trigger (record_every_dt + record_cadence + log parameters)
	 *
	 * If both are enabled, **either** trigger may cause a record.
	 *
	 * Notes on cadence:
	 * - `record_every_dt` enables time-based recording. If `record_every_dt <= 0`, time-based
	 *   scheduling is disabled regardless of `record_cadence`.
	 * - `record_cadence` determines how the observer advances `next_time_trigger_` after writing
	 *   a time-triggered record:
	 *     - LinearDt: next_time_trigger_ += record_every_dt
	 *     - LogTime : next_time_trigger_ *= q (geometric progression)
	 *
	 * Important:
	 * - These Options do NOT alter integrator stepping. They only affect what the observer writes.
	 */
	struct Options
	{
		// ---------------------------------------------------------------------
		//  Output destination / format
		// ---------------------------------------------------------------------

		/**
		 * @brief Output file path for the time-series table (CSV/TSV).
		 *
		 * The observer writes a single table-like stream with one row per recorded sample.
		 * Typical usage is a small set of scalars for plotting (e.g., t, Tinf, Lnu, etc.).
		 */
		Zaki::String::Directory output_path = "timeseries.csv";

		/**
		 * @brief Output format for the table.
		 *
		 * CSV is the default because it is ubiquitous. TSV is useful if commas appear in
		 * metadata or if tooling prefers tab-delimited files.
		 */
		OutputFormat format = OutputFormat::CSV;

		/**
		 * @brief If true, append to an existing output file; otherwise truncate.
		 *
		 * - append=false is typical for fresh runs.
		 * - append=true can be useful for continuing a run, but you must ensure
		 *   column definitions remain consistent.
		 *
		 * NOTE:
		 * - If append=true and write_header=true, the implementation typically assumes
		 *   the file already has a header and will not emit it again (cpp policy).
		 */
		bool append = false;

		// ---------------------------------------------------------------------
		//  Scheduling: start behavior
		// ---------------------------------------------------------------------

		/**
		 * @brief Whether to record a sample at the start of the run (t = t0).
		 *
		 * If true:
		 * - OnStart() writes a synthetic row at (t0, sample_index=0).
		 * - OnSample() should then skip the first emitted sample_index==0 to avoid duplication.
		 *
		 * If false:
		 * - The first record will be triggered by sample/time rules.
		 * - For time-based scheduling, the first eligible time is typically:
		 *     t0 + record_every_dt  (LinearDt)
		 *     max(t0, log_t_floor) * q  (LogTime)
		 *   where q is derived from log_q or samples_per_decade.
		 */
		bool record_at_start = true;

		// ---------------------------------------------------------------------
		//  Scheduling: sample-count trigger
		// ---------------------------------------------------------------------

		/**
		 * @brief Record every N OnSample calls (sample-based scheduling).
		 *
		 * If `record_every_n_samples > 0`, the observer records when:
		 *   (sample_index % record_every_n_samples) == 0
		 *
		 * If set to 0, sample-based triggering is disabled.
		 *
		 * This trigger is independent of time-based triggering; if either trigger
		 * fires for a given sample, a row is recorded.
		 */
		std::size_t record_every_n_samples = 1;

		// ---------------------------------------------------------------------
		//  Scheduling: time-based trigger
		// ---------------------------------------------------------------------

		/**
		 * @brief Record based on simulation time (time-based scheduling).
		 *
		 * If `record_every_dt <= 0`, time-based recording is disabled.
		 *
		 * If enabled, the observer maintains an internal `next_time_trigger_` and records when:
		 *   t >= next_time_trigger_
		 *
		 * How `next_time_trigger_` advances after recording is controlled by `record_cadence`.
		 *
		 * Interpretation notes:
		 * - LinearDt:
		 *     Uses record_every_dt as the fixed interval Δt.
		 * - LogTime:
		 *     record_every_dt is used only to initialize the first trigger time if
		 *     record_at_start == false (cpp policy choice; document in implementation).
		 */
		double record_every_dt = 0.0;

		/**
		 * @brief Scheduling policy for time-based recording.
		 *
		 * - LinearDt:  next_time_trigger_ += record_every_dt
		 * - LogTime :  next_time_trigger_ *= q
		 *
		 * This affects only the observer decision-making; the integrator still decides
		 * when it produces samples.
		 */
		RecordCadence record_cadence = RecordCadence::LinearDt;

		// ---------------------------------------------------------------------
		//  LogTime-specific parameters (observer-controlled)
		// ---------------------------------------------------------------------

		/**
		 * @brief Minimum allowed time for logarithmic scheduling (seconds).
		 *
		 * When `record_cadence == RecordCadence::LogTime`, this avoids pathological
		 * behavior near t=0 (e.g., geometric stepping from 0).
		 *
		 * The observer uses an effective time:
		 *   t_eff = max(t_current, log_t_floor)
		 * and advances:
		 *   next_time_trigger_ = max(next_time_trigger_, log_t_floor) * q
		 */
		double log_t_floor = 1.0;

		/**
		 * @brief Desired number of recorded samples per decade in time.
		 *
		 * Used only if `log_q` is not provided (>1).
		 *
		 * q is computed as:
		 *   q = 10^(1 / samples_per_decade)
		 *
		 * Example:
		 * - samples_per_decade = 50 gives ~50 samples between 10^n and 10^(n+1).
		 */
		double samples_per_decade = 50.0;

		/**
		 * @brief Explicit geometric ratio for logarithmic scheduling.
		 *
		 * If `log_q > 1` and finite, it overrides `samples_per_decade` and is used directly.
		 * Otherwise, q is computed from `samples_per_decade`.
		 *
		 * Practical note:
		 * - If q is extremely close to 1, floating-point stagnation can occur
		 *   (tnext == tcurrent). The implementation should guard against this.
		 */
		double log_q = 0.0;

		// ---------------------------------------------------------------------
		//  Output content / documentation
		// ---------------------------------------------------------------------

		/**
		 * @brief If true, write a one-line header row containing column keys.
		 *
		 * Strongly recommended for CSV/TSV outputs.
		 */
		bool write_header = true;

		/**
		 * @brief If true, also emit a sidecar metadata file describing columns.
		 *
		 * This keeps the main table compact while preserving units and descriptions.
		 * The concrete sidecar format is implementation-defined (e.g., JSON, text).
		 */
		bool write_sidecar_metadata = true;

		/**
		 * @brief Numeric formatting precision for floating point output.
		 *
		 * Use sufficiently high precision for regression comparisons.
		 * A good default is 17 (full double precision round-trip in many cases).
		 */
		int float_precision = 17;

		/**
		 * @brief Columns to record, in the order they will appear in the table.
		 *
		 * The observer will not reorder this vector; it defines the authoritative ordering.
		 */
		std::vector<Column> columns;

		/**
		 * @brief If true, build columns automatically from the diagnostics catalog.
		 *
		 * This is the schema-driven path. If enabled, the observer will append or replace
		 * Options::columns based on catalog profiles and/or explicit key lists.
		 */
		bool use_catalog = false;

		/**
		 * @brief Optional path to the diagnostics catalog JSON.
		 *
		 * Used when `use_catalog == true` and no in-memory catalog is provided.
		 */
		Zaki::String::Directory catalog_path = "";

		/**
		 * @brief Optional catalog profile(s) to use when auto-building columns.
		 *
		 * Examples:
		 *   - "timeseries_default"
		 *   - "regression_minimal"
		 *
		 * If empty, the implementation may fall back to a default profile if present,
		 * or apply another deterministic policy.
		 */
		std::vector<std::string> catalog_profiles;

		/**
		 * @brief If true, ensure builtin time column is included (t).
		 *
		 * Useful when catalog-driven columns only cover driver scalars.
		 */
		bool include_builtin_time = true;

		/**
		 * @brief If true, ensure builtin sample_index column is included.
		 */
		bool include_builtin_sample_index = true;

		// ---------------------------------------------------------------------
		//  I/O performance knobs
		// ---------------------------------------------------------------------

		/**
		 * @brief If true, use buffered I/O for the output stream.
		 *
		 * Buffered I/O can substantially reduce overhead for high-frequency sampling.
		 */
		bool buffered_io = true;

		/**
		 * @brief Buffer flush size in bytes when using buffered I/O.
		 *
		 * The implementation should flush when the internal buffer reaches this size,
		 * and always flush on OnFinish().
		 */
		std::size_t flush_bytes = 1 << 20;

		// ---------------------------------------------------------------------
		//  Helpers
		// ---------------------------------------------------------------------

		/**
		 * @brief Compute the geometric ratio q for logarithmic scheduling.
		 *
		 * - If log_q > 1 and finite, returns log_q.
		 * - Else if samples_per_decade > 0 and finite, returns 10^(1/samples_per_decade).
		 * - Else returns a conservative fallback (10.0 => 1 sample/decade).
		 *
		 * This is intended for both:
		 * - scheduling logic (AdvanceTimeTrigger / OnStart initialization),
		 * - configuration logging and debugging.
		 */
		[[nodiscard]] double ComputeLogQ() const noexcept;

		/**
		 * @brief Return a compact, human-readable summary for logging.
		 *
		 * This is intended for "one-line" startup logs.
		 * Keep it stable and tooling-friendly (no localization).
		 */
		[[nodiscard]] std::string ToLogString() const;
	};

	/**
	 * @brief Construct with options only (state-only columns).
	 *
	 * For using DriverScalar columns, we'd use the constructor that accepts drivers.
	 */
	explicit TimeSeriesObserver(const Options &opts);

	/**
	 * @brief Construct with options only (move).
	 */
	explicit TimeSeriesObserver(Options &&opts);

	/**
	 * @brief Construct with options and driver diagnostics providers.
	 *
	 * This enables DriverScalar columns, which are extracted by:
	 *  - asking the driver to populate a DiagnosticPacket,
	 *  - retrieving pkt.GetScalar(driver_key).value.
	 *
	 * Drivers are not owned by this observer; they must outlive the observer.
	 *
	 * @param opts Observer configuration options.
	 * @param drivers Driver diagnostics providers (non-owning).
	 */
	TimeSeriesObserver(Options opts,
					   std::vector<const Driver::Diagnostics::IDriverDiagnostics *> drivers);

	/**
	 * @brief Construct with options, driver diagnostics, and a diagnostics catalog.
	 *
	 * This enables DriverScalar columns, which are extracted by:
	 *  - asking the driver to populate a DiagnosticPacket,
	 *  - retrieving pkt.GetScalar(driver_key).value.
	 * Drivers are not owned by this observer; they must outlive the observer.
	 * @param opts Observer configuration options.
	 * @param drivers Driver diagnostics providers (non-owning).
	 * @param catalog Shared pointer to a diagnostics catalog for schema-driven columns.
	 */
	TimeSeriesObserver(Options opts,
					   std::vector<const Driver::Diagnostics::IDriverDiagnostics *> drivers,
					   std::shared_ptr<const Diagnostics::DiagnosticCatalog> catalog);

	~TimeSeriesObserver() override;

	/**
	 * @brief Called once before integration begins.
	 *
	 * Opens the output stream, writes header/metadata, and resets scheduling state.
	 */
	void OnStart(const RunInfo &run,
				 const Evolution::StateVector &Y0,
				 const Evolution::DriverContext &ctx) override;

	/**
	 * @brief Called whenever the evolution loop emits a sample.
	 *
	 * Writes one row if scheduling conditions are met.
	 */
	void OnSample(const SampleInfo &s,
				  const Evolution::StateVector &Y,
				  const Evolution::DriverContext &ctx) override;

	/**
	 * @brief Called once after integration ends (success or failure).
	 *
	 * Flushes and closes the output stream.
	 */
	void OnFinish(const FinishInfo &fin,
				  const Evolution::StateVector &Yf,
				  const Evolution::DriverContext &ctx) override;

	/**
	 * @brief Observer name for logs.
	 */
	[[nodiscard]] std::string Name() const override { return "TimeSeriesObserver"; }

  private:
	// -----------------------
	//  Scheduling utilities
	// -----------------------

	/**
	 * @brief Return true if a row should be written at time t / sample_index.
	 *
	 * The decision is governed by:
	 *  - record_every_n_samples
	 *  - record_every_dt
	 *
	 * The observer updates internal scheduling state after a successful write.
	 */
	[[nodiscard]] bool ShouldRecord(double t, std::uint64_t sample_index) const;

	/// Advance the next time trigger after recording at time t.
	void AdvanceTimeTrigger(double t);

	// -----------------------
	//  Output utilities
	// -----------------------

	/// Append n bytes from p to the internal output buffer (for buffered I/O).
	inline void AppendBytes_(const char *p, std::size_t n);

	/// Append a single character to the internal output buffer (for buffered I/O).
	inline void AppendChar_(char c);

	/// Append the delimiter string to the internal output buffer (for buffered I/O).
	inline void AppendDelim_();

	/// Append a newline character to the internal output buffer (for buffered I/O).
	inline void AppendNewline_();

	/// Open output stream according to Options.
	void OpenOutput();

	/// Write the table header row (column keys) once.
	void WriteHeader();

	/// Write optional sidecar metadata about columns (units/descriptions).
	void WriteSidecarMetadata(const RunInfo &run) const;

	/// Emit one row (all columns) for the given sample/state/context.
	void WriteRow(const SampleInfo &s,
				  const Evolution::StateVector &Y,
				  const Evolution::DriverContext &ctx);

	/// Delimiter string for current OutputFormat.
	[[nodiscard]] const char *Delim() const;

	// -----------------------
	//  Value extraction
	// -----------------------

	/// Extract a BuiltinState column value.
	double ExtractBuiltin(const Column &col,
						  const SampleInfo &s,
						  const Evolution::StateVector &Y,
						  const Evolution::DriverContext &ctx) const;

	/// Extract a DriverScalar column value (by diagnosing then reading packet scalar).
	double ExtractDriverScalar(const Column &col,
							   const SampleInfo &s,
							   const Evolution::StateVector &Y,
							   const Evolution::DriverContext &ctx) const;

	/// Find a driver diagnostics provider via producer-based lookup (cached lookup recommended).
	const Driver::Diagnostics::IDriverDiagnostics *FindDriverByProducer(const std::string &name) const;

	/// Build/append columns from the catalog according to opts_.catalog_profiles.
	void BuildColumnsFromCatalog();

	/// Load catalog from opts_.catalog_path if catalog_ is not provided.
	void LoadCatalogIfNeeded();

	/// Flush internal buffer to output stream if buffered_io is enabled.
	void FlushBuffer_();

  private:
	// Internal stream buffer for buffered I/O.
	std::vector<char> out_buf_;

	// If true, use internal buffering for output stream.
	// bool buffered_io = true;

	// Buffer size thresholds for flushing.
	// std::size_t flush_bytes = 1 << 20; // 1 MiB flush threshold (reasonable default)

	// Current size of data in out_buf_.
	std::size_t out_buf_bytes = 0;

	// Observer configuration options.
	Options opts_;

	// Non-owning pointers to drivers that support diagnostics snapshots.
	std::vector<const Driver::Diagnostics::IDriverDiagnostics *> drivers_;

	// Output stream for the table.
	std::ofstream out_;

	// Run state
	bool started_ = false;
	bool header_written_ = false;

	// Scheduling state
	double next_time_trigger_ = 0.0;

	// Optional: cache driver-name lookup to avoid linear search each row.
	mutable std::unordered_map<std::string, const Driver::Diagnostics::IDriverDiagnostics *> driver_cache_;

	// Cached schema info for quick validation/lookups:
	mutable std::unordered_map<std::string, const Diagnostics::ProducerCatalog *> producer_catalog_cache_;

	// Diagnostics catalog for schema-driven columns.
	std::shared_ptr<const Diagnostics::DiagnosticCatalog> catalog_;
};

} // namespace CompactStar::Physics::Evolution::Observers