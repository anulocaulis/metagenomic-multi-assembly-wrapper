#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


LOG_PATTERN = re.compile(r"^(?P<assembly>.+)_(?P<jobid>\d+)\.(?P<ext>out|err)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export per-assembly SLURM walltime from log filenames and sacct."
    )
    parser.add_argument(
        "--logs-dir",
        type=Path,
        default=Path("logs"),
        help="Directory containing SLURM .out/.err logs (default: logs)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Output CSV path (default: <logs-dir>/slurm_walltime_by_assembly.csv)",
    )
    parser.add_argument(
        "--txt",
        type=Path,
        default=None,
        help="Output TXT path (default: <logs-dir>/slurm_walltime_by_assembly.txt)",
    )
    return parser.parse_args()


def collect_log_pairs(logs_dir: Path) -> dict[tuple[str, str], dict[str, str]]:
    records: dict[tuple[str, str], dict[str, str]] = {}
    for path in logs_dir.iterdir():
        if not path.is_file():
            continue
        match = LOG_PATTERN.match(path.name)
        if not match:
            continue

        assembly = match.group("assembly")
        job_id = match.group("jobid")
        key = (assembly, job_id)

        record = records.setdefault(
            key,
            {
                "assembly": assembly,
                "job_id": job_id,
                "out_file": "",
                "err_file": "",
            },
        )

        if match.group("ext") == "out":
            record["out_file"] = path.name
        else:
            record["err_file"] = path.name

    return records


def query_sacct(job_ids: list[str]) -> dict[str, dict[str, str]]:
    if not job_ids:
        return {}

    command = [
        "sacct",
        "-j",
        ",".join(job_ids),
        "--format=JobIDRaw,State,Elapsed",
        "-P",
        "-n",
    ]

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        return {}

    sacct_by_job_id: dict[str, dict[str, str]] = {}
    for line in result.stdout.splitlines():
        parts = line.strip().split("|")
        if len(parts) < 3:
            continue

        job_id_raw, state, elapsed = parts[0].strip(), parts[1].strip(), parts[2].strip()

        if not job_id_raw.isdigit():
            continue

        if job_id_raw not in sacct_by_job_id:
            sacct_by_job_id[job_id_raw] = {"state": state, "elapsed": elapsed}

    return sacct_by_job_id


def build_rows(records: dict[tuple[str, str], dict[str, str]], sacct_by_job_id: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for assembly, job_id in sorted(records.keys(), key=lambda item: (item[0], int(item[1]))):
        record = records[(assembly, job_id)]
        sacct_record = sacct_by_job_id.get(job_id, {})

        rows.append(
            {
                "assembly": assembly,
                "job_id": job_id,
                "state": sacct_record.get("state", "NA"),
                "walltime_elapsed": sacct_record.get("elapsed", "NA"),
                "out_file": record["out_file"],
                "err_file": record["err_file"],
            }
        )

    return rows


def write_csv(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["assembly", "job_id", "state", "walltime_elapsed", "out_file", "err_file"],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_txt(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        handle.write("assembly\tjob_id\tstate\twalltime_elapsed\n")
        for row in rows:
            handle.write(
                f"{row['assembly']}\t{row['job_id']}\t{row['state']}\t{row['walltime_elapsed']}\n"
            )


def main() -> int:
    args = parse_args()

    logs_dir = args.logs_dir
    if not logs_dir.exists() or not logs_dir.is_dir():
        raise SystemExit(f"logs directory not found: {logs_dir}")

    output_csv = args.csv or (logs_dir / "slurm_walltime_by_assembly.csv")
    output_txt = args.txt or (logs_dir / "slurm_walltime_by_assembly.txt")

    records = collect_log_pairs(logs_dir)
    job_ids = sorted({record["job_id"] for record in records.values()}, key=int)
    sacct_by_job_id = query_sacct(job_ids)
    rows = build_rows(records, sacct_by_job_id)

    write_csv(rows, output_csv)
    write_txt(rows, output_txt)

    print(f"Wrote {len(rows)} rows to {output_csv}")
    print(f"Wrote {len(rows)} rows to {output_txt}")

    if not sacct_by_job_id:
        print("Warning: sacct data unavailable; state/walltime fields set to NA.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
