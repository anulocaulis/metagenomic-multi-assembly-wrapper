#!/usr/bin/env python3
import argparse
import os
import re
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_size_to_mb(text: str) -> Optional[float]:
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(Mbp|Mb|Kbp|Kb|bp)", text)
    if not match:
        return None
    value = float(match.group(1))
    unit = match.group(2)
    if unit in {"Mbp", "Mb"}:
        return value
    if unit in {"Kbp", "Kb"}:
        return value / 1000.0
    if unit == "bp":
        return value / 1_000_000.0
    return None


def parse_int(text: str) -> Optional[int]:
    match = re.search(r"([0-9][0-9,]*)", text)
    if not match:
        return None
    return int(match.group(1).replace(",", ""))


def parse_assembly_meta(path: str) -> Dict[str, Optional[str]]:
    m = re.search(r"/assemblies/(S\d+)/assembly\.([^/]+)/", path)
    if not m:
        return {"sample": None, "assembler": None}
    return {"sample": m.group(1), "assembler": m.group(2)}


def extract_gc(lines: List[str]) -> Optional[float]:
    for idx, line in enumerate(lines):
        if line.strip().startswith("A\tC\tG\tT\tN"):
            if idx + 1 < len(lines):
                parts = lines[idx + 1].strip().split("\t")
                if len(parts) >= 8:
                    try:
                        return float(parts[7])
                    except ValueError:
                        return None
    return None


def parse_block(block: str) -> Optional[Dict[str, object]]:
    lines = [line.rstrip("\n") for line in block.splitlines() if line.strip()]
    if not lines:
        return None

    assembly_line = next((line for line in lines if line.startswith("Assembly:")), None)
    if not assembly_line:
        return None

    assembly_path = assembly_line.split("Assembly:", 1)[1].strip()
    meta = parse_assembly_meta(assembly_path)

    values: Dict[str, object] = {
        "assembly_path": assembly_path,
        "sample": meta["sample"],
        "assembler": meta["assembler"],
        "gc_fraction": extract_gc(lines),
        "scaffold_total": None,
        "sequence_total_mb": None,
        "n50_kbp": None,
        "l50": None,
        "max_scaffold_kbp": None,
        "scaffolds_gt_50kb": None,
    }

    for line in lines:
        if line.startswith("Main genome scaffold total:"):
            values["scaffold_total"] = parse_int(line)
        elif line.startswith("Main genome scaffold sequence total:"):
            values["sequence_total_mb"] = parse_size_to_mb(line)
        elif line.startswith("Main genome scaffold N/L50:"):
            m = re.search(r"([0-9][0-9,]*)\s*/\s*([0-9]+(?:\.[0-9]+)?)\s*(Kbp|Kb|Mbp|Mb|bp)", line)
            if m:
                values["l50"] = int(m.group(1).replace(",", ""))
                n50_val = float(m.group(2))
                n50_unit = m.group(3)
                if n50_unit in {"Kbp", "Kb"}:
                    values["n50_kbp"] = n50_val
                elif n50_unit in {"Mbp", "Mb"}:
                    values["n50_kbp"] = n50_val * 1000.0
                elif n50_unit == "bp":
                    values["n50_kbp"] = n50_val / 1000.0
        elif line.startswith("Max scaffold length:"):
            max_mb = parse_size_to_mb(line)
            if max_mb is not None:
                values["max_scaffold_kbp"] = max_mb * 1000.0
        elif line.startswith("Number of scaffolds > 50 KB:"):
            values["scaffolds_gt_50kb"] = parse_int(line)

    return values


def parse_summary_stats_log(log_path: str) -> pd.DataFrame:
    with open(log_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    blocks: List[str] = []
    current_block: List[str] = []

    for line in lines:
        if line.startswith("Assembly:"):
            if current_block:
                blocks.append("".join(current_block))
                current_block = []
        if current_block or line.startswith("Assembly:"):
            current_block.append(line)

    if current_block:
        blocks.append("".join(current_block))

    parsed = [parse_block(block) for block in blocks]
    parsed = [row for row in parsed if row is not None]
    return pd.DataFrame(parsed)


def make_plots(df: pd.DataFrame, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")

    plot_df = df.dropna(subset=["sample", "assembler"]).copy()
    sample_order = sorted(plot_df["sample"].unique(), key=lambda x: int(x[1:]) if x and x.startswith("S") else x)
    assembler_order = sorted(plot_df["assembler"].unique())
    plot_df["sample"] = pd.Categorical(plot_df["sample"], categories=sample_order, ordered=True)
    plot_df["assembler"] = pd.Categorical(plot_df["assembler"], categories=assembler_order, ordered=True)

    # 1) Total assembly size
    size_df = plot_df.dropna(subset=["sequence_total_mb"])
    if not size_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=size_df, x="sample", y="sequence_total_mb", hue="assembler", errorbar=None)
        plt.title("Assembly Size by Sample and Assembler")
        plt.ylabel("Total Scaffold Sequence (Mb)")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "assembly_size_mb.png"), dpi=180)
        plt.close()

    # 2) N50
    n50_df = plot_df.dropna(subset=["n50_kbp"])
    if not n50_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=n50_df, x="sample", y="n50_kbp", hue="assembler", errorbar=None)
        plt.title("Scaffold N50 by Sample and Assembler")
        plt.ylabel("N50 (Kbp)")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "assembly_n50_kbp.png"), dpi=180)
        plt.close()

    # 3) Max scaffold length
    max_df = plot_df.dropna(subset=["max_scaffold_kbp"])
    if not max_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=max_df, x="sample", y="max_scaffold_kbp", hue="assembler", errorbar=None)
        plt.title("Max Scaffold Length by Sample and Assembler")
        plt.ylabel("Max Scaffold Length (Kbp)")
        plt.xlabel("Sample")
        plt.legend(title="Assembler", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "max_scaffold_kbp.png"), dpi=180)
        plt.close()

    # 4) GC vs size scatter
    plt.figure(figsize=(10, 8))
    gc_df = plot_df.dropna(subset=["gc_fraction", "sequence_total_mb"])
    if not gc_df.empty:
        plt.figure(figsize=(10, 8))
        sns.scatterplot(
            data=gc_df,
            x="gc_fraction",
            y="sequence_total_mb",
            hue="assembler",
            style="sample",
            s=140,
        )
        plt.title("GC Fraction vs Assembly Size")
        plt.xlabel("GC Fraction")
        plt.ylabel("Total Scaffold Sequence (Mb)")
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "gc_vs_size.png"), dpi=180)
        plt.close()

    # 5) Within-sample comparison (size)
    if not size_df.empty:
        g = sns.catplot(
            data=size_df,
            x="assembler",
            y="sequence_total_mb",
            col="sample",
            col_wrap=4,
            kind="bar",
            order=assembler_order,
            sharey=False,
            height=3.4,
            aspect=1.1,
            errorbar=None,
        )
        g.set_axis_labels("Assembler", "Total Scaffold Sequence (Mb)")
        g.set_titles("{col_name}")
        for ax in g.axes.flatten():
            ax.tick_params(axis="x", rotation=35)
        g.fig.suptitle("Within-Sample Assembler Comparison: Assembly Size", y=1.03)
        g.fig.tight_layout()
        g.fig.savefig(os.path.join(outdir, "within_sample_size_facets.png"), dpi=180)
        plt.close(g.fig)

    # 6) Within-sample comparison (N50)
    if not n50_df.empty:
        g = sns.catplot(
            data=n50_df,
            x="assembler",
            y="n50_kbp",
            col="sample",
            col_wrap=4,
            kind="bar",
            order=assembler_order,
            sharey=False,
            height=3.4,
            aspect=1.1,
            errorbar=None,
        )
        g.set_axis_labels("Assembler", "N50 (Kbp)")
        g.set_titles("{col_name}")
        for ax in g.axes.flatten():
            ax.tick_params(axis="x", rotation=35)
        g.fig.suptitle("Within-Sample Assembler Comparison: N50", y=1.03)
        g.fig.tight_layout()
        g.fig.savefig(os.path.join(outdir, "within_sample_n50_facets.png"), dpi=180)
        plt.close(g.fig)

    # 7) Heatmap view for quick within-sample comparisons
    size_matrix = size_df.pivot_table(index="assembler", columns="sample", values="sequence_total_mb", aggfunc="first")
    if not size_matrix.empty and size_matrix.notna().any().any():
        plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
        sns.heatmap(size_matrix, annot=True, fmt=".1f", cmap="viridis")
        plt.title("Assembly Size (Mb): Assembler vs Sample")
        plt.xlabel("Sample")
        plt.ylabel("Assembler")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "within_sample_size_heatmap.png"), dpi=180)
        plt.close()

    n50_matrix = n50_df.pivot_table(index="assembler", columns="sample", values="n50_kbp", aggfunc="first")
    if not n50_matrix.empty and n50_matrix.notna().any().any():
        plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
        sns.heatmap(n50_matrix, annot=True, fmt=".1f", cmap="magma")
        plt.title("N50 (Kbp): Assembler vs Sample")
        plt.xlabel("Sample")
        plt.ylabel("Assembler")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "within_sample_n50_heatmap.png"), dpi=180)
        plt.close()
    plt.close()

    n50_matrix = plot_df.pivot_table(index="assembler", columns="sample", values="n50_kbp", aggfunc="first")
    if not n50_matrix.empty:
        plt.figure(figsize=(max(8, len(sample_order) * 1.2), max(4, len(assembler_order) * 0.7 + 2)))
        sns.heatmap(n50_matrix, annot=True, fmt=".1f", cmap="magma")
        plt.title("N50 (Kbp): Assembler vs Sample")
        plt.xlabel("Sample")
        plt.ylabel("Assembler")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "within_sample_n50_heatmap.png"), dpi=180)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot assembly stats from summary_stats_log.txt")
    parser.add_argument("--log", default="summary_stats_log.txt", help="Path to summary stats log")
    parser.add_argument("--outdir", default="plots/summary_stats", help="Directory to write plots and parsed CSV")
    args = parser.parse_args()

    df = parse_summary_stats_log(args.log)
    if df.empty:
        raise SystemExit("No assembly blocks parsed from log. Check --log path and format.")

    os.makedirs(args.outdir, exist_ok=True)
    df.to_csv(os.path.join(args.outdir, "parsed_summary_stats.csv"), index=False)
    make_plots(df, args.outdir)

    print(f"Parsed {len(df)} assembly entries")
    print(f"Wrote CSV: {os.path.join(args.outdir, 'parsed_summary_stats.csv')}")
    print(f"Wrote plots to: {args.outdir}")


if __name__ == "__main__":
    main()
