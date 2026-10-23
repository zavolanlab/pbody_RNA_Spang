"""
Shared helpers of the Sanity_utils-based fracSanity pipeline (Snakefile-fracSanity) and of the
analysis notebook: sample metadata, comparison design, atomic pickling.
"""

import os
import pickle

import pandas as pd

# controls (Mookherjee, Rommel et al. bioRxiv 2026): untreated for -D, H2O2, DTT; DMSO (vehicle) for
# CCCP and tunicamycin; DMSO is itself compared with UT
CONTROL_OF = {"minusD": "UT", "DTT": "UT", "H2O2": "UT", "DMSO": "UT", "CCCP": "DMSO", "Tun": "DMSO"}
STRESSES = ["minusD", "DTT", "H2O2", "CCCP", "Tun"]
CONDITION_ORDER = ["UT", "DMSO", "minusD", "DTT", "H2O2", "CCCP", "Tun"]
CONDITION_COLORS = dict(zip(CONDITION_ORDER, ["grey", "skyblue", "orangered", "lime", "royalblue",
                                              "yellow", "fuchsia"]))
TIME_POINTS = ["10min", "60min"]
TIME_MARKERS = {"10min": "o", "60min": "s"}

_PRESENT_CONDITIONS = {
    "D": "minusD", "min": "minusD", "DTT": "DTT", "U": "UT", "Untreat": "UT", "CCCP": "CCCP",
    "H2O2": "H2O2", "DMSO": "DMSO", "Tun": "Tun", "Tunica": "Tun",
}


def _condition(raw):
    if raw in _PRESENT_CONDITIONS:
        return _PRESENT_CONDITIONS[raw]
    best = ""
    for key in _PRESENT_CONDITIONS:
        if raw.startswith(key) and len(key) > len(best):
            best = key
    return _PRESENT_CONDITIONS[best]


def sample_metadata(start_samples: pd.DataFrame, exclude_cccp_total_60_r123: bool = True) -> pd.DataFrame:
    """
    One row per merged sample: 'sample' (= merged_sample), 'time_point', 'fraction', 'condition',
    'replicate', 'cohort' ('<time>;<condition>'), 'group' ('<fraction>;<time>;<condition>'),
    'sample_label_within_experiment'. Same parsing as get_metadata_with_labels_and_colors in
    pbody_RNA_Spang.current.ipynb, including the exclusion of replicates 1-3 of CCCP Total 60 min.
    """
    df = start_samples[["merged_sample"]].drop_duplicates().reset_index(drop=True)
    tp = df["merged_sample"].str.split("_").str[1].str.startswith("26025").map({True: "60min", False: "10min"})
    df["time_point"] = tp
    df["fraction"] = df["merged_sample"].str.contains("-PD-").map({True: "PD", False: "Total"})
    cond = []
    for s, t in zip(df["merged_sample"], df["time_point"]):
        raw = s.split("-")[5].split("_")[0] if t == "10min" else s.split("-")[2].split("_")[0]
        cond.append(_condition(raw))
    df["condition"] = cond
    reps = []
    for s, t, c, f in zip(df["merged_sample"], df["time_point"], df["condition"], df["fraction"]):
        unadj = s.split("_S")[0]
        if t == "10min":
            reps.append("R" + (unadj.split("-")[-2][-1] if (c in ("DMSO", "UT") and f == "PD") else unadj[-1]))
        else:
            reps.append("R" + unadj.split("-")[-1])
    df["replicate"] = reps
    if exclude_cccp_total_60_r123:
        bad = (df["time_point"] == "60min") & (df["fraction"] == "Total") & (df["condition"] == "CCCP") & \
              df["replicate"].isin(["R1", "R2", "R3"])
        df = df.loc[~bad]
    df["sample"] = df["merged_sample"]
    df["cohort"] = df["time_point"] + ";" + df["condition"]
    df["group"] = df["fraction"] + ";" + df["cohort"]
    df["sample_label_within_experiment"] = df["time_point"] + ";" + df["fraction"] + ";" + \
        df["condition"] + ";" + df["replicate"]
    return df.sort_values(["fraction", "time_point", "condition", "replicate"]).reset_index(drop=True)


def comparisons() -> pd.DataFrame:
    """
    Pairwise fracSanity comparisons: every stress (and DMSO) against its control within each time
    point, and 60 vs 10 min within each condition. Columns 'label', 'test', 'reference', 'kind'.
    """
    rows = []
    for t in TIME_POINTS:
        for cond, ctrl in CONTROL_OF.items():
            rows.append((f"{t};{cond}_vs_{ctrl}", f"{t};{cond}", f"{t};{ctrl}", "stress_vs_control"))
    for cond in CONDITION_ORDER:
        rows.append((f"{cond};60min_vs_10min", f"60min;{cond}", f"10min;{cond}", "time"))
    return pd.DataFrame(rows, columns=["label", "test", "reference", "kind"])


def safe_name(label: str) -> str:
    return label.replace(";", "__")


def dump_pickle(obj, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)


def load_pickle(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)
