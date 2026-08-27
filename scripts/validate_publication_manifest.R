#!/usr/bin/env Rscript

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

required_columns <- c("status", "type", "path", "notes")
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns) > 0) {
  stop(
    "Manifest is missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

required <- manifest[manifest$status == "keep_required", , drop = FALSE]
missing_required <- required[!file.exists(required$path), , drop = FALSE]

if (nrow(missing_required) > 0) {
  message("Missing required publication files:")
  for (path in missing_required$path) {
    message("  - ", path)
  }
  stop("Publication manifest has missing required files.", call. = FALSE)
}

active_r_files <- c(
  "GEE_Code/HM/Y0_Baseline/HM_Data_Application_Conditional.R",
  "GEE_Code/HM/Y0_Baseline/Derive_Data_Like_Simulation_Values.R",
  "GEE_Code/HM/Y0_Baseline/Manuscript_Simulation.R",
  "GEE_Code/HM/Y0_Baseline/Simulation_Stress_Tests.R",
  "GEE_Code/HM/Y0_Baseline/functions_HM_conditional.R",
  "GEE_Code/HM/Y0_Baseline/generateSMART_GEE_conditional.R",
  "GEE_Code/HM/Y0_Baseline/fit_gcomp_conditional.R",
  "GEE_Code/HM/Y0_Baseline/plot_conditional_results.R",
  "GEE_Code/HM/Y0_Baseline/simulation_wrapper_utils.R"
)

missing_active_r <- active_r_files[!file.exists(active_r_files)]
if (length(missing_active_r) > 0) {
  stop(
    "Missing active workflow R files: ",
    paste(missing_active_r, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(active_r_files, parse))

git_available <- nzchar(Sys.which("git"))
if (git_available) {
  ignored_required <- character()
  for (path in required$path) {
    result <- system2(
      "git",
      c("check-ignore", "--quiet", "--", path),
      stdout = FALSE,
      stderr = FALSE
    )
    if (identical(result, 0L)) {
      ignored_required <- c(ignored_required, path)
    }
  }

  if (length(ignored_required) > 0) {
    message("Required publication files are ignored by .gitignore:")
    for (path in ignored_required) {
      message("  - ", path)
    }
    stop("Required publication files must be git-addable.", call. = FALSE)
  }

  generated_patterns <- c(
    "GEE_Code/HM/Y0_Baseline/New_Surrogate_Results/Sim_Results_Cond/Manuscript_Primary/niter_5000/iteration_estimates_manuscript_primary.csv",
    "GEE_Code/HM/Y0_Baseline/New_Surrogate_Results/Sim_Results_Cond/Manuscript_Primary/niter_5000/truth_manuscript_primary.csv",
    "GEE_Code/HM/Y0_Baseline/New_Surrogate_Results/Sim_Results_Cond/Manuscript_Primary/niter_5000/positive_cell_diagnostics_manuscript_primary.csv"
  )
  generated_existing <- generated_patterns[file.exists(generated_patterns)]
  unignored_generated <- character()
  for (path in generated_existing) {
    result <- system2(
      "git",
      c("check-ignore", "--quiet", "--", path),
      stdout = FALSE,
      stderr = FALSE
    )
    if (!identical(result, 0L)) {
      unignored_generated <- c(unignored_generated, path)
    }
  }

  if (length(unignored_generated) > 0) {
    message("Large generated simulation outputs are not ignored:")
    for (path in unignored_generated) {
      message("  - ", path)
    }
    stop("Large generated simulation outputs should stay outside git.", call. = FALSE)
  }
}

active_text <- unlist(lapply(active_r_files, readLines, warn = FALSE), use.names = FALSE)
bad_output_path <- grepl(
  "New_Surrogate/Sim_Results_Cond",
  active_text,
  fixed = TRUE
)
if (any(bad_output_path)) {
  stop(
    "Found obsolete New_Surrogate/Sim_Results_Cond path in active workflow.",
    call. = FALSE
  )
}

data_app_text <- readLines(
  "GEE_Code/HM/Y0_Baseline/HM_Data_Application_Conditional.R",
  warn = FALSE
)
if (!any(grepl("MCOACH_WEEKLY_PATH", data_app_text, fixed = TRUE)) ||
    !any(grepl("MCOACH_OUTCOME_PATH", data_app_text, fixed = TRUE))) {
  stop(
    "Data application must support MCOACH_WEEKLY_PATH and MCOACH_OUTCOME_PATH.",
    call. = FALSE
  )
}

text_file_pattern <- "(\\.(R|r|md|tsv|csv|txt)$|^\\.gitignore$)"
required_text <- required$path[grepl(text_file_pattern, required$path)]
private_path_patterns <- c(
  paste0("/", "Users", "/"),
  paste("University of Michigan", "Dropbox")
)
private_path_hits <- character()
for (path in required_text) {
  lines <- readLines(path, warn = FALSE)
  hit <- Reduce(`|`, lapply(private_path_patterns, grepl, x = lines, fixed = TRUE))
  if (any(hit)) {
    private_path_hits <- c(
      private_path_hits,
      paste0(path, ":", which(hit), ": ", lines[hit])
    )
  }
}
if (length(private_path_hits) > 0) {
  message("Private/local path references found in required publication files:")
  for (hit in private_path_hits) {
    message("  - ", hit)
  }
  stop("Remove private/local paths from required publication files.", call. = FALSE)
}

public_docs <- c(
  "README.md",
  "REQUIREMENTS.md",
  "GEE_Code/HM/README.md",
  "GEE_Code/HM/Y0_Baseline/README.md",
  "GEE_Code/HM/Y0_Baseline/Data_Results_Cond_New_Surrogate/README.md",
  "GEE_Code/HM/Y0_Baseline/New_Surrogate_Results/README.md"
)
stale_doc_patterns <- c(
  "Y0_Outcome",
  "GLMM_Code",
  "Code_Checks",
  "GEE_Code/HM/Archive",
  "GEE_Code/HM/NB",
  "Batch/",
  "logs/"
)
stale_doc_hits <- character()
for (path in public_docs[file.exists(public_docs)]) {
  lines <- readLines(path, warn = FALSE)
  hit <- Reduce(`|`, lapply(stale_doc_patterns, grepl, x = lines, fixed = TRUE))
  if (any(hit)) {
    stale_doc_hits <- c(
      stale_doc_hits,
      paste0(path, ":", which(hit), ": ", lines[hit])
    )
  }
}
if (length(stale_doc_hits) > 0) {
  message("Public-facing docs mention folders omitted from the publication snapshot:")
  for (hit in stale_doc_hits) {
    message("  - ", hit)
  }
  stop("Remove stale omitted-folder references from public-facing docs.", call. = FALSE)
}

message("Publication manifest validation passed.")
message("Required files present: ", nrow(required))
message("Active workflow R files parse: ", length(active_r_files))
