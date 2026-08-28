# Files that belong in the BRANCH checkout, not this one

The comparison is driven from the `feature/preloaded-sa-improvements` tree — that is where the
`preloaded` and `mmap` arms are built and where `run.sh` is invoked. These four files are the suite
definitions it needs, and they live here so that one `git fetch` of this branch carries the whole
comparison to a new machine instead of half of it.

They are carried here rather than committed to the feature branch on purpose: that branch is
finished, and a one-off comparison should not land in it. Copy them in, run, delete.

    bash sa-benchmarks/branch-side/install.sh /path/to/unipept-index

`install.sh` refuses to overwrite a file that already exists and differs, so re-running it after
editing a suite on the branch side will not silently discard the edit.
