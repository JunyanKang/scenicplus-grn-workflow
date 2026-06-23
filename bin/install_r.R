#!/usr/bin/env Rscript

cran_candidates <- unique(c(
  Sys.getenv("R_CRAN_REPOS", unset = NA_character_),
  "https://cloud.r-project.org",
  "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
  "https://mirrors.ustc.edu.cn/CRAN",
  "https://mirrors.aliyun.com/CRAN"
))
cran_candidates <- cran_candidates[!is.na(cran_candidates) & nzchar(cran_candidates)]

choose_cran <- function() {
  for (repo in cran_candidates) {
    message("Checking CRAN mirror: ", repo)
    ok <- tryCatch({
      ap <- available.packages(repos = repo)
      nrow(ap) > 1000
    }, warning = function(w) FALSE, error = function(e) FALSE)
    if (isTRUE(ok)) return(repo)
  }
  stop("No reachable CRAN mirror. Set R_CRAN_REPOS=https://your.cran.mirror and rerun.")
}

cran_repo <- choose_cran()
options(repos = c(CRAN = cran_repo), timeout = max(600, getOption("timeout")))
options(enrichR.live = FALSE)
Sys.setenv(
  OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", unset = "1"),
  OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", unset = "1"),
  MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", unset = "1"),
  VECLIB_MAXIMUM_THREADS = Sys.getenv("VECLIB_MAXIMUM_THREADS", unset = "1"),
  NUMEXPR_NUM_THREADS = Sys.getenv("NUMEXPR_NUM_THREADS", unset = "1"),
  MAKEFLAGS = Sys.getenv("MAKEFLAGS", unset = "-j1")
)
message("Using CRAN mirror: ", cran_repo)

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- "--file="
script_path <- sub(script_file_arg, "", script_args[grepl(script_file_arg, script_args)][1])
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
vendor_github_dir <- Sys.getenv("SCENICPLUS_VENDOR_GITHUB_DIR", unset = file.path(script_dir, "vendor", "github"))
github_tries <- suppressWarnings(as.integer(Sys.getenv("GITHUB_TRIES", unset = "3")))
if (is.na(github_tries) || github_tries < 0L) github_tries <- 3L
r_install_ncpus <- suppressWarnings(as.integer(Sys.getenv("R_INSTALL_NCPUS", unset = "2")))
if (is.na(r_install_ncpus) || r_install_ncpus < 1L) r_install_ncpus <- 2L
message("R package install parallelism: ", r_install_ncpus)
options(Ncpus = r_install_ncpus)

install_deps <- c("Depends", "Imports", "LinkingTo")
hdwgcna_commit <- "afa09abb890f5be087b63e510a7346e8e1952ecc"
hdwgcna_archive <- file.path(vendor_github_dir, paste0("hdWGCNA-", hdwgcna_commit, ".tar.gz"))
autozyme_commit <- "35f91f2229eb44d82710470803865d3c15102716"
autozyme_archive <- Sys.getenv("SCENICPLUS_AUTOZYME_ARCHIVE", unset = file.path(vendor_github_dir, paste0("autozyme-", autozyme_commit, ".tar.gz")))
install_autozyme_r <- identical(Sys.getenv("INSTALL_AUTOZYME_R", unset = "1"), "1")
r_version <- getRversion()
is_linux <- identical(Sys.info()[["sysname"]], "Linux")
cran_versions <- c(
  BiocManager = "1.30.27",
  remotes = "2.5.0",
  WGCNA = "1.74",
  Seurat = "5.5.0",
  Signac = if (r_version >= "4.5.0") "1.17.1" else "1.16.0",
  systemfonts = "1.3.2",
  tweenr = "2.0.3",
  WriteXLS = "6.8.0",
  rjson = "0.2.23",
  RhpcBLASctl = "0.23-42",
  graphlayouts = "1.2.3",
  tidygraph = if (is_linux && r_version < "4.5.0") "1.3.0" else "1.3.1",
  ggforce = "0.5.0",
  proxy = "0.4-29",
  tester = "0.3.0",
  enrichR = "3.4",
  harmony = if (is_linux && r_version < "4.5.0") "2.0.2" else "2.0.4",
  ggraph = "2.2.2",
  RcppParallel = "5.1.11-1"
)

cleanup_locks <- function() {
  for (lib in .libPaths()) {
    locks <- Sys.glob(file.path(lib, "00LOCK*"))
    if (length(locks) > 0) {
      message("Removing stale R package locks under: ", lib)
      unlink(locks, recursive = TRUE, force = TRUE)
    }
  }
}

cleanup_locks()

installed_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

version_matches <- function(pkg, version) {
  installed <- installed_version(pkg)
  if (is.na(installed)) return(FALSE)
  normalize <- function(x) gsub("-", ".", x, fixed = TRUE)
  normalize(installed) == normalize(version)
}

install_cran_exact <- function(pkg) {
  version <- unname(cran_versions[[pkg]])
  if (is.na(version)) stop("No pinned CRAN version for: ", pkg)
  if (!version_matches(pkg, version)) {
    message("Installing CRAN package: ", pkg, " ", version)
    remotes::install_version(
      pkg,
      version = version,
      repos = getOption("repos"),
      dependencies = install_deps,
      upgrade = "never",
      Ncpus = r_install_ncpus
    )
  }
  if (!version_matches(pkg, version)) {
    stop("Failed to install/load pinned CRAN package: ", pkg, " ", version)
  }
}

install_cran_exact_no_deps <- function(pkg) {
  version <- unname(cran_versions[[pkg]])
  if (is.na(version)) stop("No pinned CRAN version for: ", pkg)
  if (!version_matches(pkg, version)) {
    message("Installing CRAN package without dependency changes: ", pkg, " ", version)
    remotes::install_version(
      pkg,
      version = version,
      repos = getOption("repos"),
      dependencies = FALSE,
      upgrade = "never",
      Ncpus = r_install_ncpus
    )
  }
  if (!version_matches(pkg, version)) {
    stop("Failed to install/load pinned CRAN package: ", pkg, " ", version)
  }
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  message("Installing CRAN package: remotes")
  install.packages("remotes", repos = getOption("repos"), dependencies = install_deps, Ncpus = r_install_ncpus)
}
install_cran_exact("remotes")
install_cran_exact("BiocManager")

install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing Bioconductor package: ", pkg)
    old_repos <- getOption("repos")
    options(repos = BiocManager::repositories())
    on.exit(options(repos = old_repos), add = TRUE)
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Failed to install/load Bioconductor package: ", pkg)
  }
}

for (pkg in c("GenomicRanges", "GeneOverlap", "UCell", "impute", "preprocessCore")) {
  install_bioc_if_missing(pkg)
}

for (pkg in c(
  "WGCNA", "Seurat", "Signac", "systemfonts", "tweenr", "WriteXLS",
  "rjson", "RhpcBLASctl", "graphlayouts", "tidygraph", "ggforce",
  "proxy", "tester", "enrichR", "harmony", "ggraph"
)) {
  install_cran_exact(pkg)
}

hdwgcna_sha <- function() {
  desc <- system.file("DESCRIPTION", package = "hdWGCNA")
  if (!nzchar(desc) || !file.exists(desc)) return(NA_character_)
  fields <- read.dcf(desc)
  sha_fields <- intersect(c("RemoteSha", "GithubSHA1"), colnames(fields))
  if (length(sha_fields) == 0) return(NA_character_)
  value <- fields[1, sha_fields[1]]
  if (length(value) == 0 || is.na(value) || !nzchar(value)) return(NA_character_)
  unname(value)
}

hdwgcna_vendor_commit <- function() {
  marker <- system.file("SCENICPLUS_VENDOR_COMMIT", package = "hdWGCNA")
  if (!nzchar(marker) || !file.exists(marker)) return(NA_character_)
  value <- readLines(marker, warn = FALSE, n = 1)
  if (length(value) == 0 || !nzchar(value[1])) return(NA_character_)
  value[1]
}

hdwgcna_commit_matches <- function() {
  identical(hdwgcna_sha(), hdwgcna_commit) || identical(hdwgcna_vendor_commit(), hdwgcna_commit)
}

mark_hdwgcna_vendor_commit <- function() {
  desc <- system.file("DESCRIPTION", package = "hdWGCNA")
  if (!nzchar(desc) || !file.exists(desc)) return(invisible(FALSE))
  writeLines(hdwgcna_commit, file.path(dirname(desc), "SCENICPLUS_VENDOR_COMMIT"))
  invisible(TRUE)
}

sha256_file <- function(file) {
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    out <- system2(sha256sum, file, stdout = TRUE, stderr = TRUE)
    status <- attr(out, "status")
    if (is.null(status) || identical(status, 0L)) {
      return(strsplit(out[[1]], "[[:space:]]+")[[1]][[1]])
    }
  }

  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    out <- system2(shasum, c("-a", "256", file), stdout = TRUE, stderr = TRUE)
    status <- attr(out, "status")
    if (is.null(status) || identical(status, 0L)) {
      return(strsplit(out[[1]], "[[:space:]]+")[[1]][[1]])
    }
  }

  if (requireNamespace("openssl", quietly = TRUE)) {
    return(as.character(openssl::sha256(file(file, "rb"))))
  }

  stop("Cannot calculate SHA256: install sha256sum/shasum or the R openssl package.")
}

verify_vendor_archive <- function(archive) {
  if (!file.exists(archive) || file.info(archive)$size <= 0) return(FALSE)
  manifest <- file.path(dirname(archive), "SHA256SUMS")
  if (!file.exists(manifest) || file.info(manifest)$size <= 0) return(TRUE)
  sums <- read.table(manifest, stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(sums) < 2) stop("Malformed vendor checksum manifest: ", manifest)
  hit <- sums[sums[[2]] == basename(archive), , drop = FALSE]
  if (nrow(hit) == 0) stop("No checksum entry for bundled source archive: ", basename(archive))
  observed <- unname(sha256_file(archive))
  if (!identical(observed, hit[[1]][1])) {
    stop("Checksum mismatch for bundled source archive: ", basename(archive))
  }
  TRUE
}

install_hdwgcna <- function() {
  message("Installing hdWGCNA from GitHub commit: ", hdwgcna_commit)
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
  for (attempt in seq_len(github_tries)) {
    message("Installing hdWGCNA from GitHub (attempt ", attempt, "/", github_tries, ")")
    ok <- tryCatch({
      remotes::install_github("smorabit/hdWGCNA", ref = hdwgcna_commit, upgrade = "never", dependencies = install_deps)
      TRUE
    }, error = function(e) {
      message("GitHub install failed for hdWGCNA: ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) return(invisible(TRUE))
    Sys.sleep(attempt)
  }

  if (verify_vendor_archive(hdwgcna_archive)) {
    message("GitHub install failed after ", github_tries, " attempts; installing hdWGCNA from bundled source archive: ", hdwgcna_archive)
    remotes::install_local(hdwgcna_archive, upgrade = "never", dependencies = install_deps)
    mark_hdwgcna_vendor_commit()
    return(invisible(TRUE))
  }

  stop("Failed to install hdWGCNA from GitHub and bundled source archive was not available: ", hdwgcna_archive)
}

if (!requireNamespace("hdWGCNA", quietly = TRUE) || !hdwgcna_commit_matches()) {
  install_hdwgcna()
}

if (!requireNamespace("hdWGCNA", quietly = TRUE) || !hdwgcna_commit_matches()) {
  stop("Failed to install/load pinned hdWGCNA commit: ", hdwgcna_commit)
}

suppressPackageStartupMessages(library(hdWGCNA))
cat("hdWGCNA", as.character(packageVersion("hdWGCNA")), "\n")
cat("MetacellsByGroups", exists("MetacellsByGroups", where = asNamespace("hdWGCNA"), inherits = FALSE), "\n")

install_autozyme_r_layer <- function() {
  if (!isTRUE(install_autozyme_r)) {
    message("INSTALL_AUTOZYME_R=0: skipping R AutoZyme.")
    return(invisible(FALSE))
  }
  message("Installing R AutoZyme from pinned bundled archive with dependency changes disabled.")
  if (!requireNamespace("RcppParallel", quietly = TRUE) || !version_matches("RcppParallel", cran_versions[["RcppParallel"]])) {
    install_cran_exact_no_deps("RcppParallel")
  }
  required <- c("Rcpp", "RcppParallel", "data.table", "RcppArmadillo", "BH", "RcppEigen", "RcppAnnoy", "RcppDist")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Cannot install R AutoZyme without changing existing package versions. Missing required packages: ",
      paste(missing, collapse = ", "),
      ". The environment recipe should provide these; rerun install.sh after updating the conda environment."
    )
  }
  verify_vendor_archive(autozyme_archive)
  tmpdir <- tempfile("autozyme_r_src_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(autozyme_archive, exdir = tmpdir)
  top <- list.dirs(tmpdir, full.names = TRUE, recursive = FALSE)
  if (length(top) != 1L) stop("AutoZyme archive did not contain exactly one top-level directory: ", autozyme_archive)
  autozyme_r_dir <- file.path(top[[1]], "autozyme_r")
  if (!dir.exists(autozyme_r_dir)) stop("AutoZyme archive lacks autozyme_r directory: ", autozyme_archive)
  remotes::install_local(autozyme_r_dir, upgrade = "never", dependencies = FALSE)
  if (!requireNamespace("autozyme", quietly = TRUE)) stop("Failed to install/load R AutoZyme")
  message("R AutoZyme ", as.character(utils::packageVersion("autozyme")), " installed")
  suppressPackageStartupMessages(library(autozyme))
  ok <- tryCatch({
    autozyme::activate("seurat")
    TRUE
  }, error = function(e) {
    message("R AutoZyme installed, but seurat activation check failed: ", conditionMessage(e))
    FALSE
  })
  message("R AutoZyme seurat activation check: ", ok)
  invisible(ok)
}

install_autozyme_r_layer()
