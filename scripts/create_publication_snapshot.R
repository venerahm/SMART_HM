#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/create_publication_snapshot.R <output_dir>",
      "[--include-optional] [--overwrite] [--init-git]",
      "[--branch=<name>] [--remote-url=<url>]",
      "[--commit-message=<message>] [--allow-inside-source]"
    ),
    call. = FALSE
  )
}

output_dir <- args[[1]]
flags <- args[-1]
include_optional <- "--include-optional" %in% flags
overwrite <- "--overwrite" %in% flags
init_git <- "--init-git" %in% flags
allow_inside_source <- "--allow-inside-source" %in% flags

branch_flag <- grep("^--branch=", flags, value = TRUE)
remote_flag <- grep("^--remote-url=", flags, value = TRUE)
commit_flag <- grep("^--commit-message=", flags, value = TRUE)
if (length(branch_flag) > 1) {
  stop("Use at most one --branch=<name> flag.", call. = FALSE)
}
if (length(remote_flag) > 1) {
  stop("Use at most one --remote-url=<url> flag.", call. = FALSE)
}
if (length(commit_flag) > 1) {
  stop("Use at most one --commit-message=<message> flag.", call. = FALSE)
}

branch <- if (length(branch_flag) == 1) {
  sub("^--branch=", "", branch_flag)
} else {
  "main"
}
remote_url <- if (length(remote_flag) == 1) {
  sub("^--remote-url=", "", remote_flag)
} else {
  ""
}
commit_message <- if (length(commit_flag) == 1) {
  sub("^--commit-message=", "", commit_flag)
} else {
  ""
}

allowed_flags <- c(
  "--include-optional",
  "--overwrite",
  "--init-git",
  "--allow-inside-source"
)
flag_like <- flags[
  !grepl("^--branch=", flags) &
    !grepl("^--remote-url=", flags) &
    !grepl("^--commit-message=", flags)
]
unknown_flags <- setdiff(flag_like, allowed_flags)
if (length(unknown_flags) > 0) {
  stop("Unknown flag(s): ", paste(unknown_flags, collapse = ", "), call. = FALSE)
}
if (!nzchar(branch)) {
  stop("--branch must not be empty.", call. = FALSE)
}
if (nzchar(remote_url) && !init_git) {
  stop("--remote-url requires --init-git.", call. = FALSE)
}
if (nzchar(commit_message) && !init_git) {
  stop("--commit-message requires --init-git.", call. = FALSE)
}

normalize_for_compare <- function(path) {
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path <- file.path(getwd(), path)
  }
  if (file.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  parent <- dirname(path)
  parent_normalized <- if (file.exists(parent)) {
    normalizePath(parent, winslash = "/", mustWork = TRUE)
  } else {
    normalizePath(parent, winslash = "/", mustWork = FALSE)
  }
  file.path(parent_normalized, basename(path))
}

path_is_inside <- function(child, parent) {
  child <- normalize_for_compare(child)
  parent <- normalize_for_compare(parent)
  child == parent || startsWith(child, paste0(parent, "/"))
}

git <- Sys.which("git")
if (!allow_inside_source && nzchar(git)) {
  repo_root_output <- system2(
    git,
    c("rev-parse", "--show-toplevel"),
    stdout = TRUE,
    stderr = FALSE
  )
  repo_status <- attr(repo_root_output, "status")
  if (is.null(repo_status) && length(repo_root_output) > 0) {
    source_root <- repo_root_output[[1]]
    if (path_is_inside(output_dir, source_root)) {
      stop(
        "Refusing to create a publication snapshot inside the source repository: ",
        output_dir,
        "\nChoose a sibling/outside directory, or pass --allow-inside-source intentionally.",
        call. = FALSE
      )
    }
  }
}

manifest_path <- "PUBLICATION_MANIFEST.tsv"
if (!file.exists(manifest_path)) {
  stop("Missing manifest: ", manifest_path, call. = FALSE)
}

manifest <- utils::read.delim(
  manifest_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

statuses <- "keep_required"
if (include_optional) {
  statuses <- c(statuses, "keep_optional")
}

to_copy <- manifest[
  manifest$status %in% statuses & manifest$type == "file",
  ,
  drop = FALSE
]

missing <- to_copy[!file.exists(to_copy$path), , drop = FALSE]
if (nrow(missing) > 0) {
  message("Missing manifest files:")
  for (path in missing$path) {
    message("  - ", path)
  }
  stop("Cannot create publication snapshot with missing files.", call. = FALSE)
}

if (dir.exists(output_dir)) {
  existing <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
  if (length(existing) > 0 && !overwrite) {
    stop(
      "Output directory exists and is not empty. Re-run with --overwrite: ",
      output_dir,
      call. = FALSE
    )
  }
} else {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

for (path in to_copy$path) {
  destination <- file.path(output_dir, path)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(path, destination, overwrite = overwrite, copy.date = TRUE)
  if (!ok) {
    stop("Failed to copy: ", path, call. = FALSE)
  }
}

message("Publication snapshot created: ", normalizePath(output_dir, mustWork = FALSE))
message("Files copied: ", nrow(to_copy))
if (include_optional) {
  message("Included optional development-reference files.")
}

if (init_git) {
  if (!nzchar(git)) {
    stop("Cannot initialize git repository because git is not on PATH.", call. = FALSE)
  }
  old_wd <- setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)

  if (!dir.exists(".git")) {
    result <- system2(git, c("init", "-b", shQuote(branch)), stdout = TRUE, stderr = TRUE)
    message(paste(result, collapse = "\n"))
  } else {
    current_branch <- system2(git, c("branch", "--show-current"), stdout = TRUE)
    if (length(current_branch) == 0 || !nzchar(current_branch[[1]])) {
      system2(git, c("checkout", "-B", shQuote(branch)), stdout = TRUE, stderr = TRUE)
    }
  }

  add_output <- system2(git, c("add", "."), stdout = TRUE, stderr = TRUE)
  add_status <- attr(add_output, "status")
  if (!is.null(add_status) && add_status != 0) {
    message(paste(add_output, collapse = "\n"))
    stop("git add failed in publication snapshot.", call. = FALSE)
  }

  if (nzchar(remote_url)) {
    existing_remote <- system2(git, c("remote"), stdout = TRUE, stderr = TRUE)
    if ("origin" %in% existing_remote) {
      remote_output <- system2(
        git,
        c("remote", "set-url", "origin", shQuote(remote_url)),
        stdout = TRUE,
        stderr = TRUE
      )
    } else {
      remote_output <- system2(
        git,
        c("remote", "add", "origin", shQuote(remote_url)),
        stdout = TRUE,
        stderr = TRUE
      )
    }
    remote_status <- attr(remote_output, "status")
    if (!is.null(remote_status) && remote_status != 0) {
      message(paste(remote_output, collapse = "\n"))
      stop("Failed to configure origin remote.", call. = FALSE)
    }
    message("Configured origin remote: ", remote_url)
  }

  committed <- FALSE
  if (nzchar(commit_message)) {
    staged_files <- system2(
      git,
      c("diff", "--cached", "--name-only"),
      stdout = TRUE,
      stderr = TRUE
    )
    if (length(staged_files) == 0) {
      stop("No staged files available for initial commit.", call. = FALSE)
    }

    commit_output <- system2(
      git,
      c("commit", "-m", shQuote(commit_message)),
      stdout = TRUE,
      stderr = TRUE
    )
    commit_status <- attr(commit_output, "status")
    if (!is.null(commit_status) && commit_status != 0) {
      message(paste(commit_output, collapse = "\n"))
      stop("Initial publication commit failed.", call. = FALSE)
    }
    message(paste(commit_output, collapse = "\n"))
    committed <- TRUE
  }

  if (committed) {
    message("Initialized git repository on branch ", branch, " and created initial commit.")
  } else {
    message("Initialized git repository on branch ", branch, " and staged publication files.")
  }
}
