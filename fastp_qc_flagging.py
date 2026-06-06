# -*- coding: utf-8 -*-
"""
fastp_qc_flagging.py
--------------------
Generates sample-specific QC interpretation JSON files from:
  - fastp_thresholds.csv     : per-variable higher/lower limits and which levels to check
  - fastp_distributions.csv  : cohort-level distribution statistics (count, mean, std, min, 25%, 50%, 75%, max)
  - fastp_qc_interpretations.json : biological/NGS interpretations keyed by variable and flag level
  - model_samples.csv        : one column per sample, one row per variable (rows = variables, cols = samples)

Output:
  - One JSON file per sample: <sample_name>_qc_report.json
    containing per-variable flag status, observed value, thresholds,
    cohort context, and full biological/NGS interpretation.

Usage:
    python fastp_qc_flagging.py \
        --thresholds  fastp_thresholds.csv \
        --distributions fastp_distributions.csv \
        --interpretations fastp_qc_interpretations.json \
        --samples model_samples.csv \
        --outdir ./qc_reports

Author: generated 2026-06-05
"""

import json
import argparse
import sys
from pathlib import Path
from datetime import date

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_thresholds(path: Path) -> dict:
    """
    Returns dict keyed by variable name.
    Each value: { higher_limit, lower_limit, check_levels: list[str] }
    """
    df = pd.read_csv(path)
    required = {"variables", "higher_limit", "lower_limit", "levels_to_check"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"[ERROR] Threshold file missing columns: {missing}")

    thresholds = {}
    for _, row in df.iterrows():
        var = str(row["variables"]).strip()
        levels = [lv.strip() for lv in str(row["levels_to_check"]).split(",")]
        thresholds[var] = {
            "higher_limit": float(row["higher_limit"]),
            "lower_limit": float(row["lower_limit"]),
            "check_levels": levels,
        }
    return thresholds


def load_distributions(path: Path) -> dict:
    """
    Returns dict keyed by variable name.
    Each value: { count, mean, std, min, p25, p50, p75, max }
    The distributions CSV has stat-name in the first column (Unnamed: 0)
    and one column per variable.
    """
    df = pd.read_csv(path)
    stat_col = df.columns[0]          # 'Unnamed: 0' holding count/mean/std...
    df = df.rename(columns={stat_col: "stat"})
    df["stat"] = df["stat"].str.strip()

    stat_map = {
        "count": "count",
        "mean":  "mean",
        "std":   "std",
        "min":   "min",
        "25%":   "p25",
        "50%":   "p50",
        "75%":   "p75",
        "max":   "max",
    }

    distributions = {}
    variable_cols = [c for c in df.columns if c != "stat"]

    for var in variable_cols:
        dist = {}
        for _, row in df.iterrows():
            key = stat_map.get(str(row["stat"]).strip())
            if key:
                try:
                    dist[key] = float(row[var])
                except (ValueError, TypeError):
                    dist[key] = None
        distributions[var] = dist

    return distributions


def load_interpretations(path: Path) -> dict:
    """Loads the QC interpretations JSON."""
    with open(path) as fh:
        data = json.load(fh)
    return data.get("variables", {})


def load_samples(path: Path) -> dict:
    """
    model_samples.csv: first column = variable names, subsequent columns = samples.
    Returns dict: { sample_name: { variable: value } }
    """
    df = pd.read_csv(path)
    var_col = df.columns[0]
    df = df.rename(columns={var_col: "variable"})
    df["variable"] = df["variable"].str.strip()

    sample_cols = [c for c in df.columns if c != "variable"]
    samples = {}
    for col in sample_cols:
        sample_name = str(col).strip()
        samples[sample_name] = {}
        for _, row in df.iterrows():
            try:
                samples[sample_name][row["variable"]] = float(row[col])
            except (ValueError, TypeError):
                samples[sample_name][row["variable"]] = None
    return samples


# ---------------------------------------------------------------------------
# Flag logic
# ---------------------------------------------------------------------------

def determine_flag(value: float, higher_limit: float, lower_limit: float,
                   check_levels: list) -> str:
    """
    Returns 'high', 'low', or 'normal'.
    Only evaluates directions listed in check_levels.
    """
    if value is None:
        return "missing"

    if "high" in check_levels and value > higher_limit:
        return "high"
    if "low" in check_levels and value < lower_limit:
        return "low"
    return "normal"


def cohort_context(value: float, dist: dict) -> dict:
    """
    Builds a small cohort context block: z-score, percentile band, and
    whether the value is above/below the median.
    """
    ctx = {}
    if not dist or value is None:
        return ctx

    ctx["cohort_mean"] = round(dist.get("mean", None), 6) if dist.get("mean") is not None else None
    ctx["cohort_std"]  = round(dist.get("std",  None), 6) if dist.get("std")  is not None else None
    ctx["cohort_median"] = round(dist.get("p50", None), 6) if dist.get("p50") is not None else None
    ctx["cohort_p25"]  = round(dist.get("p25",  None), 6) if dist.get("p25")  is not None else None
    ctx["cohort_p75"]  = round(dist.get("p75",  None), 6) if dist.get("p75")  is not None else None
    ctx["cohort_min"]  = round(dist.get("min",  None), 6) if dist.get("min")  is not None else None
    ctx["cohort_max"]  = round(dist.get("max",  None), 6) if dist.get("max")  is not None else None

    mean = dist.get("mean")
    std  = dist.get("std")
    if mean is not None and std is not None and std > 0:
        ctx["z_score"] = round((value - mean) / std, 3)
    else:
        ctx["z_score"] = None

    # Rough percentile band from the quartile boundaries
    p25 = dist.get("p25")
    p50 = dist.get("p50")
    p75 = dist.get("p75")
    mn  = dist.get("min")
    mx  = dist.get("max")

    if all(x is not None for x in [mn, p25, p50, p75, mx]):
        if value < mn:
            band = "below_min"
        elif value < p25:
            band = "Q1 (0-25%)"
        elif value < p50:
            band = "Q2 (25-50%)"
        elif value < p75:
            band = "Q3 (50-75%)"
        elif value <= mx:
            band = "Q4 (75-100%)"
        else:
            band = "above_max"
        ctx["cohort_quartile_band"] = band
    else:
        ctx["cohort_quartile_band"] = None

    return ctx


# ---------------------------------------------------------------------------
# Per-variable report builder
# ---------------------------------------------------------------------------

def build_variable_report(var: str, value, threshold: dict,
                           dist: dict, interpretation: dict) -> dict:
    """
    Constructs a complete per-variable report dict.
    """
    higher_limit  = threshold["higher_limit"]
    lower_limit   = threshold["lower_limit"]
    check_levels  = threshold["check_levels"]

    flag_direction = determine_flag(value, higher_limit, lower_limit, check_levels)

    report = {
        "variable": var,
        "display_name": interpretation.get("display_name", var),
        "description": interpretation.get("description", ""),
        "unit": interpretation.get("unit", ""),
        "observed_value": round(value, 6) if value is not None else None,
        "thresholds": {
            "higher_limit": higher_limit,
            "lower_limit": lower_limit,
            "levels_checked": check_levels,
        },
        "flag": flag_direction.upper(),
        "cohort_context": cohort_context(value, dist),
    }

    # Attach interpretation details for the flagged level
    interp_levels = interpretation.get("interpretations", {})
    if flag_direction in interp_levels:
        interp = interp_levels[flag_direction]
        report["flag_code"]      = interp.get("flag", "UNKNOWN")
        report["severity"]       = interp.get("severity", "UNKNOWN")
        report["short_summary"]  = interp.get("short_summary", "")
        report["biological_causes"]    = interp.get("biological_causes", [])
        report["library_prep_causes"]  = interp.get("library_prep_causes", [])
        report["sequencing_causes"]    = interp.get("sequencing_causes", [])
        report["recommended_action"]   = interp.get("recommended_action", "")
    elif "normal" in interp_levels:
        interp = interp_levels["normal"]
        report["flag_code"]     = interp.get("flag", "PASS")
        report["severity"]      = interp.get("severity", "OK")
        report["short_summary"] = interp.get("short_summary", "")
        report["biological_causes"]   = []
        report["library_prep_causes"] = []
        report["sequencing_causes"]   = []
        report["recommended_action"]  = ""
    elif flag_direction == "missing":
        report["flag_code"]     = "MISSING"
        report["severity"]      = "UNKNOWN"
        report["short_summary"] = "Value not available for this sample."
        report["biological_causes"]   = []
        report["library_prep_causes"] = []
        report["sequencing_causes"]   = []
        report["recommended_action"]  = "Check source fastp JSON for this variable."
    else:
        report["flag_code"]     = "PASS"
        report["severity"]      = "OK"
        report["short_summary"] = "Within normal range."
        report["biological_causes"]   = []
        report["library_prep_causes"] = []
        report["sequencing_causes"]   = []
        report["recommended_action"]  = ""

    return report


# ---------------------------------------------------------------------------
# Sample-level summary
# ---------------------------------------------------------------------------

SEVERITY_ORDER = {"CRITICAL": 0, "WARNING": 1, "INFO": 2, "OK": 3, "UNKNOWN": 4}


def build_sample_summary(variable_reports: list) -> dict:
    """
    Aggregates all variable reports into a top-level sample summary.
    """
    flagged = [r for r in variable_reports if r["flag"] not in ("NORMAL", "MISSING")]
    critical = [r for r in flagged if r.get("severity") == "CRITICAL"]
    warnings = [r for r in flagged if r.get("severity") == "WARNING"]
    info     = [r for r in flagged if r.get("severity") == "INFO"]

    if critical:
        overall_status = "CRITICAL"
    elif warnings:
        overall_status = "WARNING"
    elif info:
        overall_status = "INFO"
    else:
        overall_status = "PASS"

    summary = {
        "overall_status": overall_status,
        "total_variables_checked": len(variable_reports),
        "total_flagged": len(flagged),
        "flags_by_severity": {
            "CRITICAL": len(critical),
            "WARNING":  len(warnings),
            "INFO":     len(info),
            "PASS":     len(variable_reports) - len(flagged),
        },
        "flagged_variables": [
            {
                "variable":      r["variable"],
                "flag":          r["flag"],
                "flag_code":     r.get("flag_code", ""),
                "severity":      r.get("severity", ""),
                "observed_value": r["observed_value"],
                "short_summary": r.get("short_summary", ""),
            }
            for r in sorted(flagged, key=lambda x: SEVERITY_ORDER.get(x.get("severity", "UNKNOWN"), 99))
        ],
    }
    return summary


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def generate_sample_report(sample_name: str, sample_values: dict,
                            thresholds: dict, distributions: dict,
                            interpretations: dict) -> dict:
    """
    Builds the full report dict for one sample.
    """
    variable_reports = []

    for var, threshold in thresholds.items():
        value = sample_values.get(var)
        dist  = distributions.get(var, {})
        interp = interpretations.get(var, {})

        report = build_variable_report(
            var=var,
            value=value,
            threshold=threshold,
            dist=dist,
            interpretation=interp,
        )
        variable_reports.append(report)

    sample_summary = build_sample_summary(variable_reports)

    full_report = {
        "report_metadata": {
            "sample_name":  sample_name,
            "generated_on": str(date.today()),
            "tool":         "fastp_qc_flagging.py",
            "description":  "Sample-specific fastp QC interpretation report with biological and NGS context.",
        },
        "sample_summary": sample_summary,
        "variable_reports": variable_reports,
    }

    return full_report


def run(args):
    print("[INFO] Loading input files...")
    thresholds     = load_thresholds(Path(args.thresholds))
    distributions  = load_distributions(Path(args.distributions))
    interpretations = load_interpretations(Path(args.interpretations))
    samples        = load_samples(Path(args.samples))

    print(f"[INFO] Threshold variables  : {len(thresholds)}")
    print(f"[INFO] Distribution variables: {len(distributions)}")
    print(f"[INFO] Interpretation entries: {len(interpretations)}")
    print(f"[INFO] Samples to process   : {list(samples.keys())}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for sample_name, sample_values in samples.items():
        print(f"\n[INFO] Processing sample: {sample_name}")

        # Warn about any threshold variables missing from the sample data
        missing_vars = [v for v in thresholds if v not in sample_values or sample_values[v] is None]
        if missing_vars:
            print(f"  [WARN] Variables missing in sample data: {missing_vars}")

        report = generate_sample_report(
            sample_name=sample_name,
            sample_values=sample_values,
            thresholds=thresholds,
            distributions=distributions,
            interpretations=interpretations,
        )

        out_path = outdir / f"{sample_name}_qc_report.json"
        with open(out_path, "w") as fh:
            json.dump(report, fh, indent=2)

        # Console summary
        s = report["sample_summary"]
        print(f"  Overall status : {s['overall_status']}")
        print(f"  Flagged vars   : {s['total_flagged']} / {s['total_variables_checked']}")
        for sv in s["flagged_variables"]:
            print(f"    [{sv['severity']:<8}] {sv['variable']:<40} = {sv['observed_value']}  → {sv['flag_code']}")
        print(f"  Output written : {out_path}")

    print(f"\n[DONE] All reports saved to: {outdir.resolve()}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate sample-specific fastp QC interpretation JSON reports."
    )
    parser.add_argument("--thresholds",      required=True,
                        help="Path to fastp_thresholds.csv")
    parser.add_argument("--distributions",   required=True,
                        help="Path to fastp_distributions.csv")
    parser.add_argument("--interpretations", required=True,
                        help="Path to fastp_qc_interpretations.json")
    parser.add_argument("--samples",         required=True,
                        help="Path to model_samples.csv")
    parser.add_argument("--outdir",          default="./qc_reports",
                        help="Output directory for per-sample JSON reports (default: ./qc_reports)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)
