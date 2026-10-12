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


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    frozen_parser = subparsers.add_parser("frozen")
    frozen_parser.add_argument("--work", type=Path, required=True)
    frozen_parser.add_argument("--output-json", type=Path, required=True)
    frozen_parser.add_argument("--output-tsv", type=Path, required=True)
    merge_parser = subparsers.add_parser("merge-frozen")
    merge_parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "frozen":
        return frozen(args)
    if args.mode == "merge-frozen":
        return merge_frozen(args)
    raise RuntimeError("unsupported comparison mode")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"STOP {error}", file=sys.stderr)
        raise SystemExit(1)
