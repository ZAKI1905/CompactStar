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
 * @file RunObservers.cpp
 * @brief Observer factories and option presets for Evolution runs.
 *
 * @ingroup Evolution
 */

#include "CompactStar/Physics/Evolution/Run/RunObservers.hpp"

#include <Zaki/Physics/Constants.hpp>
namespace CompactStar::Physics::Evolution::Run
{
//--------------------------------------------------------------
Observers::DiagnosticsObserver::Options
MakeDefaultDiagnosticsOptions(const RunPaths &p)
{
	Observers::DiagnosticsObserver::Options o;
	o.output_path = p.diagnostics_jsonl;
	o.record_every_n_steps = 0;
	o.record_every_dt = 200.0 * Zaki::Physics::YR_2_SEC; // every 200 years
	o.record_at_start = true;
	o.write_catalog = true;
	o.catalog_output_path = p.diagnostics_catalog_json;
	return o;
}
//--------------------------------------------------------------
Observers::TimeSeriesObserver::Options
MakeDefaultTimeSeriesOptions(const RunPaths &p)
{
	using TS = Observers::TimeSeriesObserver;

	TS::Options o;
	o.output_path = p.timeseries_table;
	o.format = TS::OutputFormat::CSV;
	o.append = false;
	o.record_at_start = true;
	o.record_every_n_samples = 1;
	o.record_every_dt = 0.0;

	o.write_header = true;
	o.write_sidecar_metadata = true;
	o.float_precision = 17;

	// Schema-driven columns
	o.use_catalog = true;
	o.catalog_path = p.diagnostics_catalog_json;
	o.catalog_profiles = {"timeseries_default"};

	// Builtins (ensures t_s and sample_index exist even if catalog profiles are driver-only).
	o.include_builtin_time = true;
	o.include_builtin_sample_index = true;

	return o;
}
//--------------------------------------------------------------
std::shared_ptr<Observers::DiagnosticsObserver>
MakeDiagnosticsObserver(const RunPaths &p,
						const std::vector<const Driver::Diagnostics::IDriverDiagnostics *> &diag_drivers,
						const Observers::DiagnosticsObserver::Options *overrides)
{
	auto o = MakeDefaultDiagnosticsOptions(p);
	if (overrides)
		o = *overrides;

	return std::make_shared<Observers::DiagnosticsObserver>(o, diag_drivers);
}
//--------------------------------------------------------------
std::shared_ptr<Observers::TimeSeriesObserver>
MakeTimeSeriesObserver(const RunPaths &p,
					   const std::vector<const Driver::Diagnostics::IDriverDiagnostics *> &diag_drivers,
					   const Observers::TimeSeriesObserver::Options *overrides)
{
	auto o = MakeDefaultTimeSeriesOptions(p);
	if (overrides)
		o = *overrides;

	return std::make_shared<Observers::TimeSeriesObserver>(o, diag_drivers);
}
//--------------------------------------------------------------
} // namespace CompactStar::Physics::Evolution::Run