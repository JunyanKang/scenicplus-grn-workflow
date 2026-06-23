# SCENIC+ GRN Version Lock

This installer pins the GRN analysis environment to the successful local dry-run combination validated on 2026-06-09. The macOS arm64 recipe keeps the local lock. The Linux recipe is conda-first and uses a CentOS7/glibc 2.17-compatible lock for packages whose newest Linux builds no longer solve cleanly on older servers.

For Linux, `environment-linux-64.yml` pins top-level package versions and `locks/environment-linux-64.solved-lock.yml` records the full dry-run solution with transitive dependency build strings.

GitHub source commits are also bundled as exact `vendor/github/*.tar.gz`
archives. During installation the script tries GitHub first, then falls back to
these bundled archives after the configured retry count. The bundled archives
are listed in `vendor/github/SHA256SUMS` and checksum-validated before use. For
packages that derive their version from git metadata, the installer passes the
pinned versions to setuptools-scm during local archive installation, so the offline
tarballs do not need a `.git` directory.

## Core Python / CLI

```text
Python                       3.10.20 macOS; 3.10.19 Linux glibc2.17 recipe
setuptools                   80.10.2
wheel                        0.47.0
MACS2                        2.2.9.1
Snakemake                    7.32.4
MALLET                       2.0.8, bundled archive installed under $CONDA_PREFIX/opt/mallet-2.0.8 with $CONDA_PREFIX/bin/mallet Python wrapper
scanpy                       1.10.4
anndata                      0.10.8
mudata                       0.2.3
pandas                       1.5.3
numpy                        1.23.5
scipy                        1.14.1
scikit-learn                 1.7.2
pyscenic                     0.12.1
ctxcore                      0.2.0
arboreto                     0.1.6
dask/distributed             2023.4.1
ray-default                  2.53.0 macOS; 2.9.3 Linux glibc2.17 recipe
pyarrow                      23.0.0 macOS; 14.0.2 Linux glibc2.17 recipe
polars                       1.41.2 macOS; 1.35.2 Linux glibc2.17 recipe
flatbuffers                  25.12.19 macOS; 25.9.23 Linux glibc2.17 recipe
python-flatbuffers           25.9.23
leidenalg                    0.12.0 macOS; 0.10.2 Linux glibc2.17 recipe
python-igraph                1.0.0 macOS; 0.11.8 Linux glibc2.17 recipe
lxml                         6.1.1 macOS; 5.3.1 Linux glibc2.17 recipe
bokeh                        3.9.1
tornado                      6.5.7 macOS; 6.5.2 Linux glibc2.17 recipe
```

The Linux `leidenalg`/`python-igraph` versions are intentionally held at the
0.10/0.11 line. Newer `leidenalg 0.12` and `python-igraph 1.0` pulled a newer
`libxml2/icu` stack that conflicted with the Linux R 4.4 layer during solving.

## Aerts / SCENIC+ Source Commits

```text
pycisTopic                   219225df56b32738d82cd14532b187a1483de04f
pycistarget                  5aa517604e4842539a7531c16905825dc7cb80fb
scenicplus                   e82b82f14b76618b850dfe442efc2421bb34f3b4
create_cisTarget_databases   304d5dc1b15e5c923908a50a1ec291c3faaccf9c
Cluster-Buster / cbust       5911cd6201b767a43316ce613afc6c9255dc3511
LoomXpy                      61995ff10940968eac2cee8fe48300ab477a15d0
```

Bundled GitHub archives:

```text
pycisTopic-219225df56b32738d82cd14532b187a1483de04f.tar.gz
pycistarget-5aa517604e4842539a7531c16905825dc7cb80fb.tar.gz
scenicplus-e82b82f14b76618b850dfe442efc2421bb34f3b4.tar.gz
create_cisTarget_databases-304d5dc1b15e5c923908a50a1ec291c3faaccf9c.tar.gz
cluster-buster-5911cd6201b767a43316ce613afc6c9255dc3511.tar.gz
LoomXpy-61995ff10940968eac2cee8fe48300ab477a15d0.tar.gz
hdWGCNA-afa09abb890f5be087b63e510a7346e8e1952ecc.tar.gz
autozyme-35f91f2229eb44d82710470803865d3c15102716.tar.gz
SHA256SUMS
```

Bundled non-GitHub archive:

```text
vendor/mallet/mallet-2.0.8.zip
```

## AutoZyme Acceleration Layer

```text
AutoZyme Python/R            0.3.1, commit 35f91f2229eb44d82710470803865d3c15102716
r-rcppdist                   0.1.1, conda-provided for R AutoZyme LinkingTo
RcppParallel                 5.1.11-1, installed from CRAN source with dependencies=FALSE if absent
```

AutoZyme is installed as a no-dependency overlay. Python installation uses the
bundled `autozyme_py` source with `pip --no-deps --no-build-isolation`. R
installation uses the bundled `autozyme_r` source with `dependencies=FALSE` and
will not upgrade or downgrade existing R packages. macOS dry-runs showed that
conda `r-rcppparallel` would downgrade `tbb/libhwloc`, so RcppParallel is not
added as a conda dependency.

## R / Bioconductor / hdWGCNA

```text
R                            4.5.3 macOS; 4.4.3 Linux glibc2.17 recipe
cmake                        4.3.3 Linux glibc2.17 recipe
pkg-config                   0.29.2 Linux glibc2.17 recipe
libuv                        1.52.1 Linux glibc2.17 recipe
xz / liblzma / liblzma-devel 5.8.3 Linux glibc2.17 recipe
Seurat                       5.5.0
Signac                       1.17.1 macOS/R 4.5; 1.16.0 Linux/R 4.4
WGCNA                        1.74
R igraph                     2.1.4 Linux/R 4.4
hdWGCNA                      0.4.11
hdWGCNA GitHub commit        afa09abb890f5be087b63e510a7346e8e1952ecc
GenomicRanges                1.62.1 macOS/R 4.5; 1.58.0 Linux/R 4.4
GeneOverlap                  1.46.0 macOS/R 4.5; 1.42.0 Linux/R 4.4
UCell                        2.14.0 macOS/R 4.5; 2.10.1 Linux/R 4.4
impute                       1.84.0 macOS/R 4.5; 1.80.0 Linux/R 4.4
preprocessCore               1.72.0 macOS/R 4.5; 1.68.0 Linux/R 4.4
fs                           2.1.0 Linux/R 4.4
Hmisc                        5.2-5 Linux/R 4.4
htmlTable                    2.5.0 Linux/R 4.4
htmlwidgets                  1.6.4 Linux/R 4.4
rmarkdown                    2.31 Linux/R 4.4
```

On Linux, the CMake/libuv stack, xz/liblzma development files, Seurat, Signac,
WGCNA, R igraph, and the heavy Seurat dependency stack are installed through
conda before `install_r.R` runs. This keeps the R layer from source-compiling
R igraph/Seurat and avoids the `ld: cannot find -llzma` failure observed on
the server. R igraph is held at `2.1.4` for compatibility with the pinned
Python `lxml 5.3.1` / `libxml2 2.13` stack.

## Pinned R Source Dependencies

```text
systemfonts                  1.3.2
tweenr                       2.0.3
WriteXLS                     6.8.0
rjson                        0.2.23
RhpcBLASctl                  0.23-42
graphlayouts                 1.2.3
tidygraph                    1.3.1 macOS/R 4.5; 1.3.0 Linux/R 4.4
ggforce                      0.5.0
proxy                        0.4-29
tester                       0.3.0
enrichR                      3.4
harmony                      2.0.4 macOS/R 4.5; 2.0.2 Linux/R 4.4
ggraph                       2.2.2
```

## Pinned Pip Dependencies

The full pip pin set is in `pip-constraints.txt`. The installer uses this constraints file and force-reinstalls the pinned pip layer.
