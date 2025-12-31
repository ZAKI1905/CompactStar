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

	/**
	 * @brief Configuration for TimeSeriesObserver.
	 */
	struct Options
	{
		/// Output file path for the table (CSV/TSV).
		Zaki::String::Directory output_path = "timeseries.csv";

		/// Output format (CSV by default).
		OutputFormat format = OutputFormat::CSV;

		/// If true, append to existing output file; otherwise truncate.
		bool append = false;

		/// If true, record at OnStart (t=t0) before first OnSample.
		bool record_at_start = true;

		/**
		 * @brief Record every N OnSample calls (step-based scheduling).
		 *
		 * - If 0, step-based triggering is disabled.
		 * - If both step-based and time-based are enabled, either trigger records.
		 */
		std::size_t record_every_n_samples = 1;

		/**
		 * @brief Record every dt in simulation time (time-based scheduling).
		 *
		 * - If <= 0, time-based triggering is disabled.
		 * - If enabled, the observer maintains an internal "next_time_trigger".
		 */
		double record_every_dt = 0.0;

		/**
		 * @brief If true, write a one-line header row containing column keys.
		 *
		 * For CSV/TSV outputs, the header is strongly recommended.
		 */
		bool write_header = true;

		/**
		 * @brief If true, also emit a sidecar metadata file describing columns.
		 *
		 * This keeps the table compact while preserving units and descriptions.
		 * A typical sidecar file could be: "<output_path>.meta.json" or ".meta.txt"
		 * depending on implementation.
		 *
		 * NOTE: This header does not prescribe the exact sidecar format; the .cpp
		 * should implement a deterministic, low-dependency approach.
		 */
		bool write_sidecar_metadata = true;

		/**
		 * @brief Numeric formatting precision for floating point output.
		 *
		 * Use sufficiently high precision for regression comparisons.
		 * Suggested default is 17 (full double precision).
		 */
		int float_precision = 17;

		/**
		 * @brief Columns to record, in the order they will appear in the table.
		 *
		 * The observer will not reorder; this is the authoritative ordering.
		 */
		std::vector<Column> columns;

		/**
		 * @brief If true, build columns automatically from the diagnostics catalog.
		 *
		 * This is the schema-driven path: no hardcoded column list in main.
		 * If enabled, the observer will append (or replace) Options::columns
		 * based on catalog profiles and/or explicit key lists.
		 */
		bool use_catalog = false;

		/**
		 * @brief Optional: path to the catalog JSON (diagnostics.catalog.json).
		 *
		 * Used when use_catalog==true and no in-memory catalog is supplied.
		 */
		Zaki::String::Directory catalog_path = "";

		/**
		 * @brief Optional: which catalog profile(s) to use when auto-building columns.
		 *
		 * Examples:
		 *   "timeseries_default"
		 *   "regression_minimal"
		 *
		 * If empty, you can fall back to "timeseries_default" if present,
		 * otherwise use all scalars or none (cpp policy).
		 */
		std::vector<std::string> catalog_profiles;

		/**
		 * @brief If true, prepend/ensure common builtin columns (time, sample_index).
		 *
		 * Useful when catalog-driven columns only cover driver scalars.
		 */
		bool include_builtin_time = true;
		bool include_builtin_sample_index = true;

		/**
		 * @brief If true, use buffered I/O for output stream.
		 *
		 */
		bool buffered_io = true;

		/**
		 * @brief Buffer flush size in bytes when using buffered I/O.
		 *
		 */
		std::size_t flush_bytes = 1 << 20;
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