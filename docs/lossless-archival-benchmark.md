# Lossless archival benchmark for this dataset

Measured 2026-09-20/21 by [Slid Phi Labs](https://slidphilabs.com) on the
76-file corpus behind this repository: 44 `.mat` files, 22 `.csv`, 5 `.txt`,
2 `.ipynb`, 2 `.py`, 1 `.mplstyle` — **273,190,871 bytes** total (the Zenodo
raw data plus the processing code in this repo).

## Method

Every candidate below was run for real on every file. The per-file winner is
the smallest **actual** compressed bytes (exact selection — no estimates, no
skipped candidates). Every compressed output was decoded and SHA-256-verified
against the original: **zero decode failures** across 380 verified
encode/decode roundtrips.

## Results (totals over all 76 files)

| Method | Total bytes | Ratio | Encode wall time |
|---|---|---:|---:|
| Original | 273,190,871 | 1.0000 | — |
| gzip -9 | 230,352,166 | 0.8432 | 20 s |
| zstd -3 | 203,914,264 | 0.7464 | 6 s |
| zstd -19 | 146,059,751 | 0.5346 | 174 s |
| xz -6 | 141,264,984 | 0.5171 | 236 s |
| PCCX 0.3.0¹ | 165,693,847 | 0.6065 | 670 s |
| PCC (classic) | 136,096,642 | 0.4982 | 5,055 s |
| TNSSRC dev | 136,096,642 | 0.4982 | 10,849 s |
| Coaster smallest | 133,397,937 | 0.4883 | 8,040 s |
| Coaster fastpass | 133,397,937 | 0.4883 | 4,563 s |
| **Per-file best (exact selection)** | **133,395,218** | **0.4883** | — |

¹ PCCX 0.3.0 refuses incompressible inputs by design (no raw passthrough), so
its total covers 65 of 76 files; the other 11 (small text/config files) are
better served by xz/zstd.

Per-file winners: Coaster smallest/fastpass won 69 files, PCCX 0.3.0 won 6
(text/CSV), PCC classic won 1. xz -6 still won 18 files outright (mostly small
CSVs and notebooks, by small margins) — the honest summary is that the lab
codecs win the corpus by **7,869,766 bytes (5.57%) over xz -6**, not every
individual file.

## Caveats

- These codecs trade encode time for bytes: Coaster took ~2.2 h on this
  corpus vs 4 min for xz -6. Decode is fast (all methods < 70 s total).
- An archive is **not** a loadable format: decode back to the original files
  before use.
- Codecs used are research builds from Slid Phi Labs (PCC is AGPL-3.0);
  nothing here changes this repository's licensing — this is a docs-only
  contribution with measured numbers.

## Reproduce

Full per-file tables, drivers, and logs:
<https://github.com/ceedot-rock/penn-nv-lossless>

Corpus sources: data via <https://doi.org/10.5281/zenodo.20969058>,
processing code via <https://doi.org/10.5281/zenodo.20970922>.
