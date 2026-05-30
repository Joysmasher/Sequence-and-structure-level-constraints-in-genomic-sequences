import os
import re
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import ttest_ind, t

# ============================================================
# CURRENT WORKING DIRECTORY
# ============================================================

BASE_DIR = Path.cwd()

print("\nWorking directory:")
print(BASE_DIR)

# ============================================================
# FILE PATHS
# ============================================================

input_file = BASE_DIR / "CQWW_frames.csv"

unique_cds_file = BASE_DIR / "unique_cds_CQWW.csv"

cache_file = BASE_DIR / "centroid_cache_CQWW.csv"

output_file = BASE_DIR / "centroid_output_CQWW.csv"

summary_file = BASE_DIR / "centroid_summary_CQWW.csv"

ttest_file = BASE_DIR / "centroid_ttests_95CI_CQWW.csv"

# ============================================================
# BASIC FUNCTIONS
# ============================================================

def normalize_shift(x):

    x = str(x).strip()

    if x in ["0", "0.0"]:
        return "0"

    if x in ["1", "1.0"]:
        return "1"

    if x in ["2", "2.0"]:
        return "2"

    return x


def clean_sequence(seq):

    seq = str(seq).upper()

    seq = seq.replace("T", "U")

    seq = re.sub(r"[^AUGC]", "", seq)

    return seq


def make_cds_id(seq):

    return hashlib.md5(seq.encode()).hexdigest()


def clean_position_plus_one(pos):

    try:

        pos = str(pos).strip()

        if "|" in pos:
            pos = pos.split("|")[0]

        return int(float(pos)) + 1

    except Exception:
        return None


def calculate_structure_metrics(struct):

    struct = str(struct).strip()

    length = len(struct)

    if length == 0:

        return np.nan, np.nan

    loop_count = struct.count(".")

    basepair_count = (
        struct.count("(") +
        struct.count(")")
    )

    loop_percentage = (
        loop_count / length
    ) * 100

    basepair_percentage = (
        basepair_count / length
    ) * 100

    return (
        loop_percentage,
        basepair_percentage
    )

# ============================================================
# STEP 3: MAP STRUCTURES
# ============================================================

def map_structures_to_15mers():

    print(
        "\nSTEP 3: Mapping structures..."
    )

    df = pd.read_csv(
        input_file,
        dtype=str,
        low_memory=False
    )

    # normalize shifts
    df["shift"] = (
        df["shift"]
        .apply(normalize_shift)
    )

    # keep only 0,1,2
    df = df[
        df["shift"].isin(
            ["0", "1", "2"]
        )
    ].copy()

    # clean sequence
    df["cds_clean"] = (
        df["cds_sequence"]
        .apply(clean_sequence)
    )

    # generate cds id
    df["cds_id"] = (
        df["cds_clean"]
        .apply(make_cds_id)
    )

    # load cache
    cache_df = pd.read_csv(
        cache_file,
        dtype=str,
        low_memory=False
    )

    cache_dict = {

        row["cds_id"]: (

            str(
                row["centroid_structure"]
            ),

            row["centroid_energy"]

        )

        for _, row in cache_df.iterrows()
    }

    output_rows = []

    for _, row in tqdm(

        df.iterrows(),

        total=len(df),

        desc="Extracting structures",

        unit="rows"

    ):

        cds_id = row["cds_id"]

        position = clean_position_plus_one(
            row["position"]
        )

        centroid_structure, centroid_energy = (

            cache_dict.get(
                cds_id,
                ("", "")
            )

        )

        local_15mer_structure = ""

        if position is not None:

            start = position - 1

            end = start + 15

            if (

                centroid_structure != ""

                and start >= 0

                and end <= len(
                    centroid_structure
                )

            ):

                local_15mer_structure = (

                    centroid_structure[
                        start:end
                    ]

                )

        loop_pct, basepair_pct = (
            calculate_structure_metrics(
                local_15mer_structure
            )
        )

        out = row.to_dict()

        out[
            "centroid_15mer_structure_from_full_CDS"
        ] = local_15mer_structure

        out[
            "centroid_15mer_loop_percentage"
        ] = loop_pct

        out[
            "centroid_15mer_basepair_percentage"
        ] = basepair_pct

        out[
            "cds_centroid_energy"
        ] = centroid_energy

        output_rows.append(out)

    out_df = pd.DataFrame(output_rows)

    out_df.to_csv(
        output_file,
        index=False
    )

    print(
        "\nMapped output saved:"
    )

    print(output_file)

# ============================================================
# WELCH T-TEST
# ============================================================

def welch_ttest_with_ci(
    x,
    y,
    confidence=0.95
):

    x = (
        pd.Series(x)
        .dropna()
        .astype(float)
    )

    y = (
        pd.Series(y)
        .dropna()
        .astype(float)
    )

    n1 = len(x)

    n2 = len(y)

    if n1 < 2 or n2 < 2:

        return (
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            n1,
            n2
        )

    mean1 = x.mean()

    mean2 = y.mean()

    var1 = x.var(ddof=1)

    var2 = y.var(ddof=1)

    diff = mean1 - mean2

    se = np.sqrt(
        var1 / n1 +
        var2 / n2
    )

    df_num = (
        var1 / n1 +
        var2 / n2
    ) ** 2

    df_den = (

        ((var1 / n1) ** 2 / (n1 - 1))

        +

        ((var2 / n2) ** 2 / (n2 - 1))

    )

    df_welch = df_num / df_den

    t_stat, p_value = ttest_ind(
        x,
        y,
        equal_var=False
    )

    alpha = 1 - confidence

    t_critical = t.ppf(
        1 - alpha / 2,
        df_welch
    )

    ci_low = diff - t_critical * se

    ci_high = diff + t_critical * se

    return (

        t_stat,
        p_value,
        diff,
        ci_low,
        ci_high,
        n1,
        n2
    )

# ============================================================
# STEP 4: SUMMARY + TTESTS
# ============================================================

def run_summary_and_ttests():

    print(
        "\nSTEP 4: Running statistics..."
    )

    df = pd.read_csv(
        output_file,
        dtype=str,
        low_memory=False
    )

    df["shift"] = (
        df["shift"]
        .apply(normalize_shift)
    )

    numeric_cols = [

        "centroid_15mer_loop_percentage",

        "centroid_15mer_basepair_percentage"

    ]

    for col in numeric_cols:

        df[col] = pd.to_numeric(
            df[col],
            errors="coerce"
        )

    # ========================================================
    # SUMMARY
    # ========================================================

    summary = df.groupby("shift").agg(

        centroid_loop_mean=(

            "centroid_15mer_loop_percentage",
            "mean"

        ),

        centroid_loop_median=(

            "centroid_15mer_loop_percentage",
            "median"

        ),

        centroid_loop_std=(

            "centroid_15mer_loop_percentage",
            "std"

        ),

        centroid_basepair_mean=(

            "centroid_15mer_basepair_percentage",
            "mean"

        ),

        centroid_basepair_median=(

            "centroid_15mer_basepair_percentage",
            "median"

        ),

        centroid_basepair_std=(

            "centroid_15mer_basepair_percentage",
            "std"

        )

    ).reset_index()

    summary = summary.round(4)

    summary.to_csv(
        summary_file,
        index=False
    )

    print(
        "\nSummary saved:"
    )

    print(summary_file)

    # ========================================================
    # T-TESTS
    # ========================================================

    groups = sorted(
        df["shift"]
        .dropna()
        .unique()
    )

    metrics = [

        (
            "loop_percentage",
            "centroid_15mer_loop_percentage"
        ),

        (
            "basepair_percentage",
            "centroid_15mer_basepair_percentage"
        )

    ]

    ttest_results = []

    for i in range(len(groups)):

        for j in range(i + 1, len(groups)):

            g1 = groups[i]

            g2 = groups[j]

            for metric_name, col in metrics:

                x = df[
                    df["shift"] == g1
                ][col]

                y = df[
                    df["shift"] == g2
                ][col]

                (
                    t_stat,
                    p_value,
                    mean_diff,
                    ci_low,
                    ci_high,
                    n1,
                    n2

                ) = welch_ttest_with_ci(
                    x,
                    y
                )

                ttest_results.append({

                    "metric":
                        metric_name,

                    "group_1":
                        g1,

                    "group_2":
                        g2,

                    "n_group_1":
                        n1,

                    "n_group_2":
                        n2,

                    "mean_difference":
                        mean_diff,

                    "ci_95_low":
                        ci_low,

                    "ci_95_high":
                        ci_high,

                    "t_statistic":
                        t_stat,

                    "p_value":
                        p_value
                })

    ttest_df = pd.DataFrame(
        ttest_results
    )

    ttest_df = ttest_df.round(6)

    ttest_df.to_csv(
        ttest_file,
        index=False
    )

    print(
        "\nT-tests saved:"
    )

    print(ttest_file)

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    print(
        "\nUsing existing cache file..."
    )

    map_structures_to_15mers()

    run_summary_and_ttests()

    print("\nDONE!")