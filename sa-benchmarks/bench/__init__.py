"""Driver for the suffix-array benchmark suites.

The Rust harness (`sa-benchmarks`) measures exactly one configuration per invocation and knows
nothing about sweeps. Everything that makes a sweep a *measurement* rather than a pile of numbers
lives here: preflight, building one binary per arm, ordering the cells so machine drift cannot land
on one arm, imposing memory ceilings, and refusing to call a delta an effect when it is smaller than
the spread the rig itself produces.

Modules, in dependency order:

* `profile`  — where this machine keeps its index, peptide files, k-mer tables and scratch space.
* `config`   — a suite's TOML: arms, axes, and how they expand into cells.
* `rig`      — the machine: root, cgroup v2, systemd-run, drop_caches, load average, git state.
* `build`    — one binary per arm, built up front and proven distinct.
* `runner`   — running cells: ordering, ceilings, thread counts, resumability.
* `records`  — reading the JSONL back and computing the statistics every suite needs.
* `report`   — rendering tables that never show a delta without its noise floor.
* `suites/`  — per-suite cell definitions and analysis.
* `fullreport` — `run.sh all`: every suite in one session, one report.

Stdlib only (`tomllib` included), except `plot.py`, which needs matplotlib and is optional.
"""

__all__ = ["build", "config", "profile", "records", "report", "rig", "runner"]
