# memory-hints

![Test](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/test.yml?logo=github&label=test)
![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&flag=memory-hints&logo=codecov)

Best-effort hints to the memory subsystem. The index is far larger than any cache and larger than
the TLB can map at 4 KB granularity — the suffix array alone is 160 GB of a 223 GB full-UniProt
index — so both storage backends spend most of their time waiting on memory. Nothing here affects
correctness.

## Three hints, three disciplines

| module | what it asks for | issued |
|---|---|---|
| `prefetch` | start a load early, so its ~80-100 ns DRAM latency overlaps useful work | per access, on the hot path |
| `hugepages` | 2 MB pages, so walking a multi-gigabyte buffer costs far fewer TLB misses | once per allocation, at load time |
| `warmup` | fault a mapped section in up front, so the reads that follow find it resident | once per mapped section, before serving or timing |

They are grouped here because they share a rationale, not a usage pattern, and getting any of them
wrong is silent:

* `prefetch_read` is `#[inline(always)]` on purpose — a `call`/`ret` pair around one hint
  instruction costs more than the hint saves.
* `advise_capacity` must be issued **between reserving an allocation and first writing to it**.
  One line too late it is not merely weaker, it is worthless, and it still looks like it is working.
  That ordering bug once left 160 GB of suffix array and 43 GB of text on 4 KB pages.
* `touch_all_pages` launders its read through `black_box`; without that the optimizer deletes the
  loop and the warmup silently does nothing.

The first two are advisory — the CPU and the kernel may ignore them — and are silent no-ops off
their target platform, so callers need no `cfg` guards. `warmup` is the exception: it does the
faulting itself and costs real time, which is why it runs at startup rather than per query.

## Where it sits

Depends on `libc` and `memmap2`. Used by `bitarray`, `protein-text`, `protein-metadata` and
`sa-index`, none of which depend on one another — which is why this is a crate rather than a module.
`prefetch` and `hugepages` were separate crates before, the latter inside `bitarray`; they merged
once it was clear four of `hugepages`' six call sites advise plain `Vec`s in `sa-index` with no bit
array in sight.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p memory-hints --open`
