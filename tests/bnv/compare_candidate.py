#!/usr/bin/env python3
"""Independent Phase-6A-1 evidence comparison.

The ``frozen`` mode consumes only serialized validation-tool outputs.  It uses
the achieved baryon counts as regression abscissae, performs the predeclared
uncertainty-weighted outer-half fits, and writes the exact monotone rows later
consumed by the production frozen-validity owner.  It does not run a BNV
trajectory or alter a declared threshold.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

sys.dont_write_bytecode = True

EPS = sys.float_info.epsilon
YEAR = 31557600.0
MEV_TO_ERG = 1.602176634e-6
KB_MEV_K = 8.617333262145e-11
REQUIRED = (
    "t_n", "t_e", "t_mu", "Z_row_npe", "Z_row_npmu", "Cstar",
    "Ltilde_Me", "Ltilde_Mmu", "mu_B", "mu_n_profile", "P0_average",
    "P1_average", "P2_average", "N_n", "N_e", "N_mu",
    "species_support", "metric_structure", "radius", "surface_gravity",
    "envelope", "Tsurface_inf",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"empty evidence table: {path}")
    return rows


def finite(value: str | float, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"nonfinite {label}")
    return result


def grouped(rows: list[dict[str, str]], count: int = 21) -> list[list[dict[str, str]]]:
    result = [[] for _ in range(count)]
    for row in rows:
        index = int(row["index"])
        if not 0 <= index < count:
            raise RuntimeError("out-of-range frozen sample index")
        result[index].append(row)
    if any(not group for group in result):
        raise RuntimeError("incomplete frozen sample table")
    return result


@dataclass
class Series:
    name: str
    values: list[list[float]]
    errors: list[list[float]]
    scales: list[float]
    threshold_kind: str

    def validate(self) -> None:
        if self.name not in REQUIRED or len(self.values) != 21 or len(self.errors) != 21:
            raise RuntimeError(f"malformed series {self.name}")
        width = len(self.values[0])
        if width == 0 or len(self.scales) != width:
            raise RuntimeError(f"empty series {self.name}")
        for values, errors in zip(self.values, self.errors):
            if len(values) != width or len(errors) != width:
                raise RuntimeError(f"ragged series {self.name}")
            if any(not math.isfinite(x) for x in values + errors) or any(x < 0 for x in errors):
                raise RuntimeError(f"invalid series number {self.name}")
        if any(not (scale > 0 and math.isfinite(scale)) for scale in self.scales):
            raise RuntimeError(f"invalid series scale {self.name}")


def regression(x: list[float], y: list[float], uncertainty: list[float]):
    # A positive floating floor prevents a zero-error oracle from inventing
    # infinite weight; it is not a physical or threshold uncertainty.
    sigma = [max(u, 64 * EPS * max(1.0, abs(v))) for u, v in zip(uncertainty, y)]
    weight = [1.0 / (u * u) for u in sigma]
    s = sum(weight)
    sx = sum(w * z for w, z in zip(weight, x))
    sy = sum(w * v for w, v in zip(weight, y))
    x_mean = sx / s
    y_mean = sy / s
    centered_xx = sum(w * (z - x_mean) ** 2 for w, z in zip(weight, x))
    centered_xy = sum(w * (z - x_mean) * (v - y_mean)
                      for w, z, v in zip(weight, x, y))
    if not math.isfinite(centered_xx) or centered_xx <= 0:
        raise RuntimeError("rank-deficient frozen sensitivity fit")
    slope = centered_xy / centered_xx
    intercept = y_mean - slope * x_mean
    covariance = ((1 / s + x_mean * x_mean / centered_xx, -x_mean / centered_xx),
                  (-x_mean / centered_xx, 1 / centered_xx))
    residuals = []
    residual_uncertainties = []
    for z, value, u in zip(x, y, sigma):
        prediction = intercept + slope * z
        variance = (covariance[0][0] + 2 * z * covariance[0][1]
                    + z * z * covariance[1][1])
        if variance < 0 and abs(variance) <= 64 * EPS * max(1.0, covariance[0][0]):
            variance = 0
        if variance < 0 or not math.isfinite(variance):
            raise RuntimeError("invalid frozen fit covariance")
        residuals.append(value - prediction)
        residual_uncertainties.append(math.sqrt(u * u + variance))
    return intercept, slope, covariance, residuals, residual_uncertainties


def threshold(series: Series, index: int, uncertainty_norm: float) -> float:
    if series.threshold_kind == "t_absolute":
        return max(5 * uncertainty_norm, 1e-4 * abs(series.values[index][0]))
    if series.threshold_kind == "z_relative":
        return max(5 * uncertainty_norm, 1e-4)
    if series.threshold_kind == "relative_plus_error":
        return 1e-4 + 5 * uncertainty_norm
    if series.threshold_kind == "relative_exact":
        return 1e-4
    if series.threshold_kind == "topology":
        return 1.0
    raise RuntimeError(f"unknown threshold kind {series.threshold_kind}")


def assess(series: Series, x: list[float]) -> tuple[list[dict[str, float]], list[dict]]:
    series.validate()
    raw_drift: list[float] = []
    uncertainty_norms: list[float] = []
    for index in range(21):
        drift = max(abs(value - reference) / scale for value, reference, scale in
                    zip(series.values[index], series.values[0], series.scales))
        uncertainty = max(math.hypot(error, reference_error) / scale
                          for error, reference_error, scale in
                          zip(series.errors[index], series.errors[0], series.scales))
        raw_drift.append(drift)
        uncertainty_norms.append(uncertainty)

    outer = list(range(10, 21))
    fit_records = []
    nonlinearity = 0.0
    for component, scale in enumerate(series.scales):
        fit_x = [x[index] for index in outer]
        fit_y = [series.values[index][component] for index in outer]
        fit_u = [series.errors[index][component] for index in outer]
        try:
            intercept, slope, covariance, residual, residual_u = regression(fit_x, fit_y, fit_u)
        except RuntimeError as error:
            raise RuntimeError(f"{series.name} component {component}: {error}") from error
        component_max = 0.0
        checks = []
        for index, value, uncertainty in zip(outer, residual, residual_u):
            residual_norm = abs(value) / scale
            uncertainty_normalized = uncertainty / scale
            allowed = 3 * uncertainty_normalized + 0.10 * threshold(
                series, index, uncertainty_norms[index]
            )
            if residual_norm > allowed:
                raise RuntimeError(
                    f"{series.name} outer-half nonlinearity failed at {index}: "
                    f"{residual_norm} > {allowed}"
                )
            component_max = max(component_max, residual_norm)
            checks.append({"index": index, "residual_normalized": residual_norm,
                           "uncertainty_normalized": uncertainty_normalized,
                           "allowed": allowed})
        nonlinearity = max(nonlinearity, component_max)
        fit_records.append({"component": component, "intercept": intercept,
                            "slope": slope, "covariance": covariance,
                            "checks": checks})

    result = []
    envelope = 0.0
    for index in range(21):
        envelope = max(envelope, raw_drift[index])
        bound = envelope + nonlinearity
        limit = threshold(series, index, uncertainty_norms[index])
        utilization = bound / limit
        if not math.isfinite(utilization) or utilization > 1:
            raise RuntimeError(
                f"{series.name} frozen envelope failed at {index}: {bound} > {limit}"
            )
        result.append({"raw_drift": raw_drift[index], "numerical_uncertainty": uncertainty_norms[index],
                       "nonlinearity": nonlinearity, "drift_bound": bound,
                       "threshold": limit, "utilization": utilization})
    return result, fit_records


def scalar_series(name, rows, value, error, kind, absolute_scale=None):
    values = [[finite(row[value], f"{name} value")] for row in rows]
    errors = [[finite(row[error], f"{name} error") if error else 0.0] for row in rows]
    scale = absolute_scale or max(abs(values[0][0]), 1e-300)
    return Series(name, values, errors, [scale], kind)


def frozen(args) -> int:
    root = args.work.resolve()
    targets = read_tsv(root / "target_stars.tsv")
    tangents = read_tsv(root / "tangent.tsv")
    coefficients = read_tsv(root / "coefficients.tsv")
    thermal_groups = grouped(read_tsv(root / "thermal_surface.tsv"))
    profile_groups = grouped(read_tsv(root / "profile_metrics.tsv"))
    if not (len(targets) == len(tangents) == len(coefficients) == 21):
        raise RuntimeError("frozen certificate requires exactly 21 stars")
    for rows in (targets, tangents, coefficients):
        if [int(row["index"]) for row in rows] != list(range(21)):
            raise RuntimeError("noncanonical frozen sample order")

    b0 = finite(targets[0]["B_solved"], "B0")
    tau = 5e-11 * b0
    x = []
    for index, row in enumerate(targets):
        b = finite(row["B_solved"], "B_solved")
        target = finite(row["B_target"], "B_target")
        residual = finite(row["residual"], "target residual")
        bracket = finite(row["bracket_width"], "target bracket")
        expected = -5e-8 * index
        if abs((target / b0 - 1) - expected) > 64 * EPS:
            raise RuntimeError("target-B grid changed")
        if abs(residual) > tau or bracket < 0 or bracket > tau:
            raise RuntimeError("target-B certification tolerance failed")
        if abs((b - target) - residual) > 64 * EPS * b0:
            raise RuntimeError("target-B residual identity failed")
        x.append((b - b0) / b0)
    if any(not x[i] < x[i - 1] for i in range(1, 21)):
        raise RuntimeError("achieved target-B abscissae are not monotone")

    series: list[Series] = [
        scalar_series("t_n", tangents, "t_n", "u_t_n", "t_absolute", 1.0),
        scalar_series("t_e", tangents, "t_e", "u_t_e", "t_absolute", 1.0),
        scalar_series("t_mu", tangents, "t_mu", "u_t_mu", "t_absolute", 1.0),
        scalar_series("Ltilde_Me", coefficients, "Ltilde_Me", "uLtilde_Me", "relative_plus_error"),
        scalar_series("Ltilde_Mmu", coefficients, "Ltilde_Mmu", "uLtilde_Mmu", "relative_plus_error"),
        scalar_series("mu_B", coefficients, "mu_B_inf", "u_mu_B", "relative_plus_error"),
        scalar_series("P0_average", coefficients, "P0_inf", "u_P0", "relative_plus_error"),
        scalar_series("P1_average", coefficients, "P1_inf", "u_P1", "relative_plus_error"),
        scalar_series("P2_average", coefficients, "P2_inf", "u_P2", "relative_plus_error", 1.0),
        scalar_series("N_n", targets, "N_n", "u_N_n", "relative_exact"),
        scalar_series("N_e", targets, "N_e", "u_N_e", "relative_exact"),
        scalar_series("N_mu", targets, "N_mu", "u_N_mu", "relative_exact"),
        scalar_series("radius", targets, "radius_km", None, "relative_plus_error"),
    ]

    for name, values, errors in (
        ("Z_row_npe", ("Z00", "Z01"), ("uZ00", "uZ01")),
        ("Z_row_npmu", ("Z10", "Z11"), ("uZ10", "uZ11")),
    ):
        data = [[finite(row[key], key) for key in values] for row in coefficients]
        uncertainty = [[finite(row[key], key) for key in errors] for row in coefficients]
        row_scale = max(max(abs(value) for value in data[0]), 1e-300)
        series.append(Series(name, data, uncertainty, [row_scale, row_scale], "z_relative"))

    def table_series(name, groups, value, error, kind):
        values = [[finite(row[value], f"{name} value") for row in group] for group in groups]
        errors = [[finite(row[error], f"{name} error") if error else 0.0 for row in group]
                  for group in groups]
        scales = [max(abs(value), 1e-300) for value in values[0]]
        return Series(name, values, errors, scales, kind)

    series.extend((
        table_series("Cstar", thermal_groups, "Cstar_erg_K", "u_Cstar_erg_K", "relative_plus_error"),
        table_series("envelope", thermal_groups, "Tsurface_local_K", None, "relative_plus_error"),
        table_series("Tsurface_inf", thermal_groups, "Tsurface_inf_K", "u_Tsurface_inf_K", "relative_plus_error"),
        scalar_series("surface_gravity", [group[0] for group in thermal_groups],
                      "surface_g14", "u_surface_g14", "relative_plus_error"),
    ))

    mu_values = [[finite(row["mu_n_inf_MeV"], "mu_n profile") for row in group]
                 for group in profile_groups]
    profile_width = len(mu_values[0])
    if any(len(group) != profile_width for group in profile_groups):
        raise RuntimeError("profile metric grids differ")
    zero_errors = [[0.0] * profile_width for _ in range(21)]
    series.append(Series("mu_n_profile", mu_values, zero_errors,
                         [max(abs(value), 1e-300) for value in mu_values[0]],
                         "relative_plus_error"))

    metric_keys = ("r_km", "m_km", "nu", "lambda")
    metric_values = []
    for group in profile_groups:
        metric_values.append([finite(row[key], key) for key in metric_keys for row in group])
    metric_scales = []
    for key in metric_keys:
        scale = max(max(abs(finite(row[key], key)) for row in profile_groups[0]), 1e-300)
        metric_scales.extend([scale] * profile_width)
    series.append(Series("metric_structure", metric_values,
                         [[0.0] * len(metric_scales) for _ in range(21)],
                         metric_scales, "relative_plus_error"))

    def support_topology(group):
        dimensions = [int(row["active_dimension"]) for row in group]
        if any(value <= 0 for value in dimensions):
            raise RuntimeError("inactive chart inside certified source support")
        return tuple(value for index, value in enumerate(dimensions)
                     if index == 0 or value != dimensions[index - 1])

    reference_topology = support_topology(profile_groups[0])
    topology = []
    for group in profile_groups:
        topology.append([0.0 if support_topology(group) == reference_topology else 2.0])
    series.append(Series("species_support", topology, [[0.0] for _ in range(21)],
                         [1.0], "topology"))

    for name in ("t_e", "t_mu"):
        item = next(value for value in series if value.name == name)
        for value, error in zip(item.values, item.errors):
            if value[0] - 5 * error[0] <= 0:
                raise RuntimeError(f"{name} certified lower bound is not positive")

    if {item.name for item in series} != set(REQUIRED):
        missing = set(REQUIRED) - {item.name for item in series}
        extra = {item.name for item in series} - set(REQUIRED)
        raise RuntimeError(f"frozen quantity mismatch missing={missing} extra={extra}")

    assessments = {}
    fits = {}
    for item in series:
        assessments[item.name], fits[item.name] = assess(item, x)

    samples = []
    for index, target in enumerate(targets):
        sample = {
            "index": index,
            "fractional_depletion": -5e-8 * index,
            "achieved_abscissa": x[index],
            "B_solved_count": finite(target["B_solved"], "B solved"),
            "B_target_residual_count": finite(target["residual"], "B residual"),
            "final_bracket_width_count": finite(target["bracket_width"], "B bracket"),
            "DeltaN_over_N": [
                finite(target["N_n"], "N_n") / finite(targets[0]["N_n"], "N_n0") - 1,
                finite(target["N_e"], "N_e") / finite(targets[0]["N_e"], "N_e0") - 1,
                finite(target["N_mu"], "N_mu") / finite(targets[0]["N_mu"], "N_mu0") - 1,
            ],
            "drift_bound": {}, "threshold": {}, "utilization": {},
        }
        for name in REQUIRED:
            result = assessments[name][index]
            for key in ("drift_bound", "threshold", "utilization"):
                sample[key][name] = result[key]
        samples.append(sample)

    payload = {
        "schema_id": "compactstar.phase6a1.frozen-sensitivity-certificate.v1",
        "classification": "CONTROLLED MATHEMATICAL STATIC CERTIFICATE; NO BNV TRAJECTORY",
        "B0_count": b0,
        "tau_B_target_count": tau,
        "regression_abscissa": "achieved (B_solved-B0)/B0",
        "outer_half_indices": list(range(10, 21)),
        "fit_criterion": "abs(residual)<=3*u_res+0.10*T_X",
        "samples": samples,
        "fits": fits,
    }
    args.output_json.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    with args.output_tsv.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("index", "fractional_depletion", "B_solved_count",
                         "B_target_residual_count", "final_bracket_width_count",
                         "DeltaN_n_over_N_n", "DeltaN_e_over_N_e", "DeltaN_mu_over_N_mu",
                         "quantity", "drift_bound", "threshold", "utilization"))
        for sample in samples:
            for name in REQUIRED:
                writer.writerow((sample["index"], sample["fractional_depletion"],
                                 sample["B_solved_count"], sample["B_target_residual_count"],
                                 sample["final_bracket_width_count"], *sample["DeltaN_over_N"], name,
                                 sample["drift_bound"][name], sample["threshold"][name],
                                 sample["utilization"][name]))
    maximum = max((sample["utilization"][name], sample["index"], name)
                  for sample in samples for name in REQUIRED)
    print(f"BA13_STATIC_CERTIFICATE PASS max_utilization {maximum[0]:.17g} "
          f"index {maximum[1]} quantity {maximum[2]}")
    return 0


def merge_frozen(args) -> int:
    root = args.work.resolve()
    for stem in ("coefficients", "thermal_surface", "profile_metrics"):
        shards = sorted(root.glob(f"{stem}-shard-*.tsv"),
                        key=lambda path: int(path.stem.split("-")[-1]))
        if len(shards) != 21:
            raise RuntimeError(f"{stem} requires 21 coefficient shards")
        fieldnames = None
        rows = []
        for expected, path in enumerate(shards):
            with path.open(newline="") as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                if fieldnames is None:
                    fieldnames = reader.fieldnames
                elif reader.fieldnames != fieldnames:
                    raise RuntimeError(f"{stem} shard schema mismatch")
                current = list(reader)
            if not current or any(int(row["index"]) != expected for row in current):
                raise RuntimeError(f"{stem} shard {expected} has wrong index")
            rows.extend(current)
        with (root / f"{stem}.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"MERGED_FROZEN_EVIDENCE {stem} rows {len(rows)}")
    return 0


def trapezoid(rows: list[dict[str, str]], key_or_callable) -> float:
    value = key_or_callable if callable(key_or_callable) else lambda row: finite(row[key_or_callable], key_or_callable)
    total = 0.0
    for left, right in zip(rows, rows[1:]):
        dt = finite(right["t_s"], "time") - finite(left["t_s"], "time")
        if not dt > 0:
            raise RuntimeError("trajectory time grid is not strictly increasing")
        total += 0.5 * dt * (value(left) + value(right))
    return total


def trajectory_rows(path: Path, expected: int) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if len(rows) != expected:
        raise RuntimeError(f"{path.name}: {len(rows)} rows, expected {expected}")
    required = {
        "t_s", "x_state", "B_count", "Bdot_count_s", "DeltaB_over_B0",
        "S_n_count_s", "S_e_count_s", "S_mu_count_s", "t_n", "t_e", "t_mu",
        "sigma_e_count_s", "sigma_mu_count_s", "eta_e_MeV", "eta_mu_MeV",
        "xi_e", "xi_mu", "R_e_count_s", "R_mu_count_s", "mu_B_inf_MeV",
        "mu_n_actual_inf_MeV", "Echem_MeV", "P_dir_eq_erg_s",
        "P_dir_actual_erg_s", "LH_erg_s", "DeltaLnu_erg_s", "DeltaPbeta_erg_s",
        "Lnu_eq_erg_s", "Lnu_full_erg_s", "L_out_fluid_inf_erg_s",
        "L_esc_star_inf_erg_s", "J_X_inf_erg_s", "Lgamma_erg_s", "Lother_erg_s",
        "Pnet_erg_s", "Cstar_erg_K", "Tinf_K", "Tsurface_inf_K",
        "R18_residual_erg_s", "Ra_Rb_residual_erg_s", "Rb_Rc_residual_erg_s",
        "finite_T_omitted_floor_erg_s", "max_frozen_utilization", "valid_through_sample",
    }
    missing = required - set(rows[0])
    if missing:
        raise RuntimeError(f"{path.name}: missing schema fields {sorted(missing)}")
    for index, row in enumerate(rows):
        for name in required - {"valid_through_sample"}:
            finite(row[name], f"{path.name} row {index} {name}")
        if row["valid_through_sample"] != "1":
            raise RuntimeError(f"{path.name}: invalid checkpoint serialized")
    return rows


def r20(rows: list[dict[str, str]]) -> dict[str, float]:
    first, last = rows[0], rows[-1]
    mu_b = finite(first["mu_B_inf_MeV"], "mu_B")
    delta_eeq = MEV_TO_ERG * mu_b * (finite(last["B_count"], "Bf") - finite(first["B_count"], "Bi"))
    delta_echem = MEV_TO_ERG * (finite(last["Echem_MeV"], "Echem f") - finite(first["Echem_MeV"], "Echem i"))
    def thermal_energy(sequence: list[dict[str, str]]) -> float:
        total = 0.0
        for left, right in zip(sequence, sequence[1:]):
            total += 0.5 * (finite(left["Cstar_erg_K"], "Cstar") +
                            finite(right["Cstar_erg_K"], "Cstar")) * (
                finite(right["Tinf_K"], "T") - finite(left["Tinf_K"], "T"))
        return total
    delta_uth = thermal_energy(rows)
    luminosity = lambda row: sum(finite(row[name], name) for name in (
        "L_out_fluid_inf_erg_s", "Lnu_full_erg_s", "Lgamma_erg_s", "Lother_erg_s"))
    outgoing = trapezoid(rows, luminosity)
    residual = delta_eeq + delta_echem + delta_uth + outgoing
    normalizer_integral = trapezoid(rows, lambda row: sum(abs(finite(row[name], name)) for name in (
        "Lnu_full_erg_s", "Lgamma_erg_s", "Lother_erg_s")))
    normalizer = max(1.0, abs(delta_uth), abs(delta_echem), normalizer_integral)
    even = rows[::2]
    if even[-1] is not rows[-1]:
        raise RuntimeError("R20 even-index quadrature lost endpoint")
    even_outgoing = trapezoid(even, luminosity)
    quadrature_error = abs(outgoing - even_outgoing)
    thermal_quadrature_error = abs(delta_uth - thermal_energy(even))
    finite_t_floor = trapezoid(rows, "finite_T_omitted_floor_erg_s")
    return {
        "delta_Eeq_erg": delta_eeq, "delta_Echem_erg": delta_echem,
        "delta_Uth_erg": delta_uth, "integrated_outgoing_erg": outgoing,
        "residual_erg": residual, "normalizer_erg": normalizer,
        "normalized_residual": abs(residual) / normalizer,
        "quadrature_error_erg": quadrature_error,
        "quadrature_over_normalizer": quadrature_error / normalizer,
        "thermal_energy_quadrature_error_erg": thermal_quadrature_error,
        "thermal_energy_quadrature_over_normalizer": thermal_quadrature_error / normalizer,
        "finite_T_omission_floor_integral_erg": finite_t_floor,
    }


def read_steps(path: Path) -> dict:
    lines = path.read_text().splitlines()
    if len(lines) < 4:
        raise RuntimeError(f"incomplete step evidence {path}")
    names, values = lines[0].split("\t"), [finite(item, path.name) for item in lines[1].split("\t")]
    result = dict(zip(names, values))
    outputs = [[finite(item, path.name) for item in line.split("\t")] for line in lines[3:]]
    if not outputs:
        raise RuntimeError(f"missing checkpoint step evidence {path}")
    increments = []
    previous = 0
    for output in outputs:
        accepted = int(output[1])
        if accepted < previous:
            raise RuntimeError("nonmonotone cumulative accepted steps")
        increments.append(accepted - previous)
        previous = accepted
    final_start = max(0, len(outputs) - max(1, len(outputs) // 10))
    final = outputs[final_start:]
    result.update({
        "maximum_accepted_steps_per_checkpoint": max(increments),
        "final_decade_min_last_step_s": min(row[3] for row in final),
        "final_decade_max_last_step_s": max(row[3] for row in final),
        "terminal_last_step_s": final[-1][3],
        "final_decade_outputs": len(final),
        "no_step_collapse": all(row[3] > 0 for row in final) and max(increments) < 100000,
    })
    if not result["no_step_collapse"]:
        raise RuntimeError("final-decade RKF45 step collapse")
    return result


def hm(xi: float) -> float:
    u2 = (xi / math.pi) ** 2
    return xi / math.pi ** 2 * (14680 + u2 * (7560 + u2 * (840 + 24 * u2))) / 11513


def hm_derivative(xi: float) -> float:
    u2 = (xi / math.pi) ** 2
    return (14680 + 3 * 7560 * u2 + 5 * 840 * u2 ** 2 + 7 * 24 * u2 ** 3) / (11513 * math.pi ** 2)


def fm_increment(xi: float) -> float:
    u2 = (xi / math.pi) ** 2
    return u2 * (22020 + u2 * (5670 + u2 * (420 + 9 * u2))) / 11513


def two_eigenvalues(matrix: tuple[tuple[float, float], tuple[float, float]]) -> tuple[float, float]:
    trace = matrix[0][0] + matrix[1][1]
    det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    discriminant = trace * trace - 4 * det
    if discriminant < 0 and abs(discriminant) <= 64 * EPS * trace * trace:
        discriminant = 0
    if discriminant < 0:
        raise RuntimeError("complex linear-QSS relaxation eigenvalues")
    root = math.sqrt(discriminant)
    return 0.5 * (trace - root), 0.5 * (trace + root)


def validate_trajectory(args) -> int:
    root = args.work.resolve()
    coefficient = read_tsv(args.coefficients.resolve())[0]
    z = ((finite(coefficient["Z00"], "Z00"), finite(coefficient["Z01"], "Z01")),
         (finite(coefficient["Z10"], "Z10"), finite(coefficient["Z11"], "Z11")))
    ltilde = (finite(coefficient["Ltilde_Me"], "Ltilde Me"),
              finite(coefficient["Ltilde_Mmu"], "Ltilde Mmu"))
    cards = {
        "RF-P0-TRANSIENT-v1": (1025, True),
        "CPL-P0-TRANSIENT-v1": (2049, False),
        "CPL-P1-TRANSIENT-v1": (2049, False),
        "CPL-P2-LINEAR-QSS-v1": (8193, False),
    }
    result = {"schema_id": "compactstar.phase6a1.controlled-bnv-validation.v1", "cards": {},
              "BA11": "PASS", "BA12": "PASS", "BA13": "PASS", "BA14": "PASS",
              "BA16": "PASS", "BA17": "PASS"}
    global_scaled = 0.0
    worst_beta = -math.inf
    for card, (count, reaction_free) in cards.items():
        versions = {}
        for control in (False, True):
            stem = card + (".control" if control else "")
            base = trajectory_rows(root / f"{stem}.baseline.tsv", count)
            refined = trajectory_rows(root / f"{stem}.refined.tsv", count)
            if [row["t_s"] for row in base] != [row["t_s"] for row in refined]:
                raise RuntimeError(f"{stem}: baseline/refined output grids differ")
            expected_depletion = 0.0 if control else {
                "RF-P0-TRANSIENT-v1": -1e-7, "CPL-P0-TRANSIENT-v1": -1e-8,
                "CPL-P1-TRANSIENT-v1": -1e-8, "CPL-P2-LINEAR-QSS-v1": -5e-7}[card]
            if abs(finite(base[-1]["DeltaB_over_B0"], "final depletion") - expected_depletion) > 64 * EPS:
                raise RuntimeError(f"{stem}: final depletion changed")
            max_scaled = [0.0, 0.0, 0.0]
            for index, (left, right) in enumerate(zip(base, refined)):
                for component, name in enumerate(("x_state", "eta_e_MeV", "eta_mu_MeV")):
                    a, b = finite(left[name], name), finite(right[name], name)
                    atol = (1e-12, 1e-18, 1e-18)[component]
                    scaled = abs(a - b) / (atol + 1e-7 * max(abs(a), abs(b)))
                    max_scaled[component] = max(max_scaled[component], scaled)
                    global_scaled = max(global_scaled, scaled)
                scale = max(1.0, abs(finite(left["P_dir_actual_erg_s"], "Pdir")),
                            abs(finite(left["P_dir_eq_erg_s"], "Pdir eq")))
                # O13's trajectory detector is the predeclared nonzero-eta P2
                # fixture.  P0/P1 retain their independently integrated
                # partition residual as a reported quadrature diagnostic.
                if card == "CPL-P2-LINEAR-QSS-v1" and not control and \
                   abs(finite(left["R18_residual_erg_s"], "R18")) > 64 * EPS * scale:
                    raise RuntimeError(f"{stem}: R18 failed at row {index}")
                for residual in ("Ra_Rb_residual_erg_s", "Rb_Rc_residual_erg_s"):
                    if abs(finite(left[residual], residual)) > 64 * EPS * scale:
                        raise RuntimeError(f"{stem}: {residual} failed at row {index}")
                if finite(left["max_frozen_utilization"], "frozen utilization") > 1:
                    raise RuntimeError(f"{stem}: frozen validity exceeded")
                ledger_scale = max(1.0, *(abs(finite(left[name], name)) for name in
                    ("LH_erg_s", "DeltaLnu_erg_s", "DeltaPbeta_erg_s", "Lnu_eq_erg_s", "Lnu_full_erg_s")))
                if abs(finite(left["DeltaPbeta_erg_s"], "DeltaPbeta") -
                       (finite(left["LH_erg_s"], "LH") - finite(left["DeltaLnu_erg_s"], "DeltaLnu"))) > 64 * EPS * ledger_scale:
                    raise RuntimeError(f"{stem}: beta ledger identity failed")
                if abs(finite(left["Lnu_full_erg_s"], "Lnu full") -
                       (finite(left["Lnu_eq_erg_s"], "Lnu eq") + finite(left["DeltaLnu_erg_s"], "DeltaLnu"))) > 64 * EPS * ledger_scale:
                    raise RuntimeError(f"{stem}: full-neutrino identity failed")
                if not control and not reaction_free:
                    for xi in (finite(left["xi_e"], "xi e"), finite(left["xi_mu"], "xi mu")):
                        ratio = fm_increment(xi) - xi * hm(xi)
                        worst_beta = max(worst_beta, ratio)
                        if ratio > 0.467659 + 1e-10:
                            raise RuntimeError(f"{stem}: BA17 modified-Urca bound failed")
            if max(max_scaled) > 1:
                raise RuntimeError(f"{stem}: BA12 scaled state difference {max(max_scaled)}")
            if reaction_free:
                closure_base = closure_refined = None
                endpoint_error = endpoint_ratio = 0.0
            else:
                closure_base, closure_refined = r20(base), r20(refined)
                for label, closure in (("baseline", closure_base), ("refined", closure_refined)):
                    if closure["normalized_residual"] > 2e-4:
                        raise RuntimeError(f"{stem} {label}: R20 normalized residual failed")
                    if closure["quadrature_over_normalizer"] > 5e-5:
                        raise RuntimeError(f"{stem} {label}: R20 quadrature failed")
                    if closure["thermal_energy_quadrature_over_normalizer"] > 5e-5:
                        raise RuntimeError(f"{stem} {label}: thermal-energy quadrature failed")
                endpoint_error = abs((closure_base["delta_Eeq_erg"] + closure_base["delta_Echem_erg"] + closure_base["delta_Uth_erg"]) -
                                     (closure_refined["delta_Eeq_erg"] + closure_refined["delta_Echem_erg"] + closure_refined["delta_Uth_erg"]))
                endpoint_ratio = endpoint_error / closure_base["normalizer_erg"]
                if endpoint_ratio > 5e-5:
                    raise RuntimeError(f"{stem}: endpoint-state propagation failed")
            steps = {level: read_steps(root / f"{stem}.{level}.tsv.steps") for level in ("baseline", "refined")}
            versions["control" if control else "source"] = {
                "maximum_scaled_state_difference": max(max_scaled), "scaled_state_components": max_scaled,
                "R20_baseline": closure_base, "R20_refined": closure_refined,
                "endpoint_state_error_erg": endpoint_error, "endpoint_state_error_over_normalizer": endpoint_ratio,
                "step_statistics": steps, "maximum_frozen_utilization": max(finite(row["max_frozen_utilization"], "util") for row in base),
                "maximum_abs_DeltaB_over_B0": max(abs(finite(row["DeltaB_over_B0"], "depletion")) for row in base),
                "final": {name: finite(base[-1][name], name) for name in
                          ("Tinf_K", "eta_e_MeV", "eta_mu_MeV", "xi_e", "xi_mu", "Pnet_erg_s")},
            }
            if reaction_free and not control:
                analytic_max = 0.0
                for row in base:
                    sigma = (finite(row["sigma_e_count_s"], "sigma e"), finite(row["sigma_mu_count_s"], "sigma mu"))
                    expected = (-(z[0][0] * sigma[0] + z[0][1] * sigma[1]) * finite(row["t_s"], "t"),
                                -(z[1][0] * sigma[0] + z[1][1] * sigma[1]) * finite(row["t_s"], "t"))
                    for channel, name in enumerate(("eta_e_MeV", "eta_mu_MeV")):
                        actual = finite(row[name], name)
                        analytic_max = max(analytic_max, abs(actual - expected[channel]) /
                                           ((1e-18,) * 2)[channel] + 1e-7 * max(abs(actual), abs(expected[channel])))
                if analytic_max > 1:
                    raise RuntimeError("reaction-free analytic transient failed")
                versions["source"]["reaction_free_analytic_scaled_difference"] = analytic_max
        source_rows = trajectory_rows(root / f"{card}.baseline.tsv", count)
        control_rows = trajectory_rows(root / f"{card}.control.baseline.tsv", count)
        delta_t = finite(source_rows[-1]["Tinf_K"], "source T") - finite(control_rows[-1]["Tinf_K"], "control T")
        closure = versions["source"]["R20_baseline"]
        floor = max(closure["finite_T_omission_floor_integral_erg"], 0.0) if closure else trapezoid(source_rows, "finite_T_omitted_floor_erg_s")
        versions["matched_control"] = {"observable": "DeltaTinf_K at final checkpoint", "window_s": [0.0, finite(source_rows[-1]["t_s"], "tf")],
                                       "central": delta_t, "finite_T_omission_floor_integral_erg": floor,
                                       "classification": "SIGN_UNRESOLVED"}
        result["cards"][card] = versions

    qss = trajectory_rows(root / "CPL-P2-LINEAR-QSS-v1.baseline.tsv", 8193)
    terminal = [row for row in qss if finite(row["t_s"], "t") >= 4e5 * YEAR]
    qss_metrics = {"interval_year": [4e5, 5e5], "channels": []}
    for channel, (xi_name, eta_name, rate_name, sigma_name) in enumerate((
        ("xi_e", "eta_e_MeV", "R_e_count_s", "sigma_e_count_s"),
        ("xi_mu", "eta_mu_MeV", "R_mu_count_s", "sigma_mu_count_s"))):
        maximum_balance = 0.0
        maximum_tau_ratio = 0.0
        for row in terminal:
            xi, eta, temperature = finite(row[xi_name], xi_name), finite(row[eta_name], eta_name), finite(row["Tinf_K"], "T")
            derivative = ltilde[channel] / 1.380649e-16 * temperature ** 7 * hm_derivative(xi) / (KB_MEV_K * temperature)
            eta_resolution = 1e-18 + 1e-7 * abs(eta)
            resolution = abs(derivative) * eta_resolution
            sigma, rate = finite(row[sigma_name], sigma_name), finite(row[rate_name], rate_name)
            maximum_balance = max(maximum_balance, abs(rate + sigma) / max(abs(sigma), resolution))
            derivatives = []
            for c, (xn, en) in enumerate((("xi_e", "eta_e_MeV"), ("xi_mu", "eta_mu_MeV"))):
                x = finite(row[xn], xn)
                derivatives.append(ltilde[c] / 1.380649e-16 * temperature ** 7 * hm_derivative(x) / (KB_MEV_K * temperature))
            eigen = two_eigenvalues(((z[0][0] * derivatives[0], z[0][1] * derivatives[1]),
                                     (z[1][0] * derivatives[0], z[1][1] * derivatives[1])))
            if min(eigen) <= 0:
                raise RuntimeError("nonpositive QSS relaxation mode")
            tau = 1 / min(eigen)
            maximum_tau_ratio = max(maximum_tau_ratio, tau / finite(row["t_s"], "elapsed"))
        if maximum_balance > 0.05 or maximum_tau_ratio > 0.10:
            raise RuntimeError(f"BA16 channel {channel} QSS reachability failed")
        qss_metrics["channels"].append({"maximum_balance_ratio": maximum_balance,
                                         "maximum_tau_relax_over_elapsed": maximum_tau_ratio,
                                         "classification": "LINEAR_QSS"})
    result["BA16_metrics"] = qss_metrics
    result["BA17_worst_dimensionless_bound"] = worst_beta
    result["BA12_maximum_scaled_state_difference"] = global_scaled
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"BA11 PASS R20 all predeclared source/control runs; BA12 PASS max_scaled {global_scaled:.17g}")
    print(f"BA13_RUNTIME PASS; BA14 PASS; BA16 PASS; BA17 PASS worst {worst_beta:.17g}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    frozen_parser = subparsers.add_parser("frozen")
    frozen_parser.add_argument("--work", type=Path, required=True)
    frozen_parser.add_argument("--output-json", type=Path, required=True)
    frozen_parser.add_argument("--output-tsv", type=Path, required=True)
    merge_parser = subparsers.add_parser("merge-frozen")
    merge_parser.add_argument("--work", type=Path, required=True)
    trajectory_parser = subparsers.add_parser("trajectory")
    trajectory_parser.add_argument("--work", type=Path, required=True)
    trajectory_parser.add_argument("--coefficients", type=Path, required=True)
    trajectory_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "frozen":
        return frozen(args)
    if args.mode == "merge-frozen":
        return merge_frozen(args)
    if args.mode == "trajectory":
        return validate_trajectory(args)
    raise RuntimeError("unsupported comparison mode")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"STOP {error}", file=sys.stderr)
        raise SystemExit(1)
