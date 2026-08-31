"""Where this machine keeps its index, its query files and its scratch space.

Every driver script this package replaces hardcoded its paths, and they had already drifted apart:
some pointed at `/mnt/data/sa-improvements/...`, some at `$REPO/uniprot-2025-04`, while the checkout
they ran from held `uniprot-2026-01`. Several carried a comment warning that the dated directory
name differs per machine and to check before running — a warning that only helps someone who reads
it. A profile makes the machine's layout one file that is validated once, up front, instead of a
constant in ten scripts.

A profile is TOML:

    index_dir = "/mnt/data/.../uniprot-2025-04/suffix-array"
    scratch   = "/mnt/data/tmp/bench"

    [peptides]
    small = "../peptides/small.txt"      # relative paths resolve against index_dir
    mixed = "../peptides/peptides_5_50.txt"

    [kmer_tables]
    k5 = "kmer-tables/5mer_table.bin"
    k6 = "kmer-tables/6mer_table.bin"

`profiles/local.toml` is the default and is gitignored; `profiles/example.toml` is the committed
template. Suites refer to peptide files and k-mer tables by *name* (`mixed`, `k6`), never by path,
so the same suite definition runs on any machine that has a profile.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

from .rig import count_lines

#: Files a suffix-array index directory must contain. `warmup.txt` is included because every suite
#: warms up before timing, and a missing warmup file fails deep inside a run rather than here.
INDEX_FILES = ("sa.bin", "proteins.bin", "mapping.bin", "warmup.txt")


class ProfileError(Exception):
    """A profile is missing, malformed, or points at something that is not there."""


@dataclass
class Profile:
    """One machine's layout. Paths are absolute and have been checked to exist."""

    name: str
    repo: Path
    index_dir: Path
    scratch: Path
    peptides: dict[str, Path] = field(default_factory=dict)
    kmer_tables: dict[str, Path] = field(default_factory=dict)

    def peptide_file(self, name: str) -> Path:
        """The peptide file a suite knows as `name`, or a listing of what this profile has."""
        try:
            return self.peptides[name]
        except KeyError:
            have = ", ".join(sorted(self.peptides)) or "(none)"
            raise ProfileError(
                f"profile '{self.name}' has no peptide file named '{name}' (has: {have})"
            ) from None

    def kmer_table(self, name: str | None) -> Path | None:
        """The pre-built table a suite knows as `name`, or None if this profile has no such file.

        None does NOT mean "run without a table" — the caller falls back to building the same table
        in-process (see `kmer_k`). Silently dropping a table a suite asked for would change what is
        measured; paying to build it only costs time.
        """
        return self.kmer_tables.get(name) if name else None

    def describe(self) -> list[tuple[str, str]]:
        """Every input this run read, as (label, value) pairs for the provenance table.

        The whole dataset, not just the index: which peptide files were queried and how long they
        are (a file shorter than `runs * amount` silently shortens the later reps), and which k-mer
        tables were attached and how big they are (the 6-mer is 3 GB resident against the 5-mer's
        127 MB, which is the entire trade under a memory ceiling). A report that names the index but
        not the queries cannot be reproduced from.
        """
        rows = [("profile", self.name), ("index dir", str(self.index_dir))]
        total = 0
        for filename in INDEX_FILES:
            size = (self.index_dir / filename).stat().st_size
            total += size
            rows.append((f"  {filename}", _size(size)))
        rows.append(("  index total", _size(total)))

        for name, path in sorted(self.peptides.items()):
            rows.append((f"  peptides:{name}", f"{_size(path.stat().st_size)}  ({count_lines(path):,} lines)  {path}"))
        if self.kmer_tables:
            for name, path in sorted(self.kmer_tables.items()):
                rows.append((f"  kmer:{name}", f"{_size(path.stat().st_size)}  {path}"))
        else:
            rows.append(("  kmer tables", "none in this profile — suites that ask for one build it in-process"))
        rows.append(("scratch", str(self.scratch)))
        return rows


def _size(byte_count: int) -> str:
    """Bytes at whichever unit reads without a string of zeros or a leading 0.00."""
    for unit, scale in (("GB", 2**30), ("MB", 2**20), ("KB", 2**10)):
        if byte_count >= scale:
            return f"{byte_count / scale:,.2f} {unit}"
    return f"{byte_count} B"


def profiles_dir(repo: Path) -> Path:
    return repo / "sa-benchmarks" / "profiles"


def load(name: str, repo: Path) -> Profile:
    """Loads and fully validates `profiles/<name>.toml`.

    Validation is deliberately eager — every path is resolved and stat'ed here rather than when a
    cell first needs it. A sweep that discovers a missing peptide bucket four hours in has wasted
    four hours, and that is the failure mode these scripts actually had.
    """
    path = profiles_dir(repo) / f"{name}.toml"
    if not path.exists():
        available = sorted(p.stem for p in profiles_dir(repo).glob("*.toml"))
        raise ProfileError(
            f"no profile '{name}' at {path}\n"
            f"  available: {', '.join(available) or '(none)'}\n"
            f"  copy profiles/example.toml to profiles/local.toml and edit it for this machine"
        )

    with path.open("rb") as handle:
        raw = tomllib.load(handle)

    index_dir = _require_dir(raw, "index_dir", path, base=repo)
    for filename in INDEX_FILES:
        target = index_dir / filename
        if not target.is_file() or target.stat().st_size == 0:
            raise ProfileError(f"{path}: index_dir is missing (or has an empty) {filename}: {target}")

    scratch = _abs(raw.get("scratch") or (repo / "target" / "bench-results"), base=repo)
    scratch.mkdir(parents=True, exist_ok=True)

    peptides = {
        key: _require_file(value, f"peptides.{key}", path, base=index_dir)
        for key, value in raw.get("peptides", {}).items()
    }
    if not peptides:
        raise ProfileError(f"{path}: no [peptides] entries — every suite needs at least one query file")

    kmer_tables = {
        key: _require_file(value, f"kmer_tables.{key}", path, base=index_dir)
        for key, value in raw.get("kmer_tables", {}).items()
    }

    return Profile(
        name=name,
        repo=repo,
        index_dir=index_dir,
        scratch=scratch,
        peptides=peptides,
        kmer_tables=kmer_tables,
    )


def _abs(value: str | Path, base: Path) -> Path:
    """Resolves `value` against `base` when relative, so a profile can be written either way."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def _require_dir(raw: dict, key: str, source: Path, base: Path) -> Path:
    if key not in raw:
        raise ProfileError(f"{source}: missing required key '{key}'")
    path = _abs(raw[key], base)
    if not path.is_dir():
        raise ProfileError(f"{source}: {key} is not a directory: {path}")
    return path


def _require_file(value: str, key: str, source: Path, base: Path) -> Path:
    path = _abs(value, base)
    if not path.is_file():
        raise ProfileError(f"{source}: {key} does not exist: {path}")
    if path.stat().st_size == 0:
        raise ProfileError(f"{source}: {key} is empty: {path}")
    return path
