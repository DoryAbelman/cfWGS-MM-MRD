#!/usr/bin/env Rscript

# =============================================================================
# run_manuscript_workflow.R
#
# Purpose:
#   One command-line entry point for the manuscript workflow, while preserving
#   the original numbered-script organization.
#
# What this script does:
#   1. Dry-runs or executes the numbered source-pipeline scripts.
#   2. Refreshes the stage-ordered script-to-artifact map.
#   3. Validates the direct manuscript-output tree in final_manuscript_objects/.
#   4. Optionally runs the separate reproducible_workflow generation/validation
#      harness when --run-reference-workflow is supplied.
#
# Guardrails:
#   - Dry-run is the default.
#   - The cache-sensitive nested-CV/model-training script is skipped unless
#     --include-cache-sensitive is explicitly supplied.
#   - Scientific logic remains in the numbered analysis scripts. This file only
#     orchestrates command-line execution and validation.
#
# Typical use from the project root:
#   Rscript Scripts_2025/Final_Scripts/run_manuscript_workflow.R
#   Rscript Scripts_2025/Final_Scripts/run_manuscript_workflow.R --execute
#   Rscript Scripts_2025/Final_Scripts/run_manuscript_workflow.R --execute --skip-source
#   Rscript Scripts_2025/Final_Scripts/run_manuscript_workflow.R --execute --run-reference-workflow
#
# Manuscript outputs created/updated:
#   - None directly. This top-level orchestration script coordinates numbered
#     source-stage execution, artifact-map refresh, validation, and manuscript
#     export staging while preserving the numbered scripts as the scientific
#     source of truth.
# =============================================================================

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

message_log <- function(..., log_file = NULL) {
  msg <- paste0(timestamp(), " | ", paste0(..., collapse = ""))
  message(msg)
  if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
}

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  normalizePath(getwd(), mustWork = TRUE)
}

parse_flag_value <- function(args, flag, default = NULL) {
  inline <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(inline)) return(sub(paste0("^", flag, "="), "", inline[1]))
  idx <- which(args == flag)
  if (!length(idx)) return(default)
  if (idx[1] == length(args)) stop(flag, " requires a value.", call. = FALSE)
  args[idx[1] + 1]
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if ("--help" %in% args || "-h" %in% args) {
    cat(
      "Usage:\n",
      "  Rscript Scripts_2025/Final_Scripts/run_manuscript_workflow.R [options]\n\n",
      "Purpose:\n",
      "  One command-line entry point for the cleaned Final_Scripts manuscript pipeline.\n",
      "  By default this is a dry run. Add --execute to run the numbered scripts.\n\n",
      "Common options:\n",
      "  --execute                  Run numbered scripts instead of dry-run.\n",
      "  --skip-source              Skip numbered-script execution; validate existing outputs only.\n",
      "  --run-reference-workflow   Also run run_analysis.R --mode generate and --mode validate.\n",
      "  --include-cache-sensitive  Include cache-sensitive model/nested-CV training stages.\n",
      "  --from ID --to ID          Run a contiguous numbered-script range, e.g. --from 2_0 --to 2_4.\n",
      "  --only ID                  Run one numbered script ID, e.g. --only 4_1.\n",
      "  --check-packages           Check package availability before stage execution.\n",
      "  --keep-going               Continue after a numbered-script failure where possible.\n",
      "  --help, -h                 Show this help text.\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }
  list(
    execute = "--execute" %in% args,
    skip_source = "--skip-source" %in% args,
    refresh_frozen = "--refresh-frozen" %in% args,
    run_reference_workflow = "--run-reference-workflow" %in% args,
    include_cache_sensitive = "--include-cache-sensitive" %in% args,
    keep_going = "--keep-going" %in% args,
    check_packages = "--check-packages" %in% args,
    from = parse_flag_value(args, "--from"),
    to = parse_flag_value(args, "--to"),
    only = parse_flag_value(args, "--only")
  )
}

add_optional_value_arg <- function(args, flag, value) {
  if (!is.null(value) && nzchar(value)) c(args, flag, value) else args
}

run_rscript <- function(project_root, script_path, args = character(), label, log_file) {
  cmd <- file.path(R.home("bin"), "Rscript")
  full_args <- c(script_path, args)
  message_log("Starting ", label, ": ", paste(c("Rscript", full_args), collapse = " "), log_file = log_file)
  status <- system2(cmd, args = full_args)
  message_log("Finished ", label, " with status ", status, log_file = log_file)
  if (!identical(status, 0L)) {
    stop(label, " failed with status ", status, call. = FALSE)
  }
  invisible(status)
}

source_pipeline_args <- function(args) {
  out <- character()
  if (isTRUE(args$execute)) out <- c(out, "--execute")
  if (isTRUE(args$include_cache_sensitive)) out <- c(out, "--include-cache-sensitive")
  if (isTRUE(args$keep_going)) out <- c(out, "--keep-going")
  if (isTRUE(args$check_packages)) out <- c(out, "--check-packages")
  out <- add_optional_value_arg(out, "--from", args$from)
  out <- add_optional_value_arg(out, "--to", args$to)
  out <- add_optional_value_arg(out, "--only", args$only)
  out
}

main <- function() {
  args <- parse_args()
  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  source(file.path(script_dir, "pipeline_metadata.R"))

  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(script_dir, "pipeline_logs", paste0("manuscript_workflow_", run_id))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(run_dir, "run_manuscript_workflow.log")

  message_log("Project root: ", project_root, log_file = log_file)
  message_log("Run directory: ", run_dir, log_file = log_file)
  if (!args$execute) {
    message_log("Dry-run mode. Add --execute to run numbered source scripts.", log_file = log_file)
  }
  if (!args$include_cache_sensitive) {
    message_log(
      "Cache-sensitive nested-CV/model training remains skipped. Add --include-cache-sensitive only for deliberate recomputation; regenerated artifacts may differ unless the exact inputs, seeds, package versions, and execution environment are preserved.",
      log_file = log_file
    )
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  if (!args$skip_source) {
    run_rscript(
      project_root = project_root,
      script_path = file.path("Scripts_2025", "Final_Scripts", "run_pipeline.R"),
      args = source_pipeline_args(args),
      label = "numbered source pipeline",
      log_file = log_file
    )
  } else {
    message_log("Skipping numbered source pipeline because --skip-source was supplied.", log_file = log_file)
  }

  run_rscript(
    project_root = project_root,
    script_path = file.path("Scripts_2025", "Final_Scripts", "build_stage_artifact_map.R"),
    label = "stage-ordered script-to-artifact map",
    log_file = log_file
  )

  if (isTRUE(args$refresh_frozen)) {
    run_rscript(
      project_root = project_root,
      script_path = "run_analysis.R",
      args = c("--mode", "frozen"),
      label = "frozen reference refresh",
      log_file = log_file
    )
  }

  if (isTRUE(args$execute) && isTRUE(args$run_reference_workflow)) {
    run_rscript(
      project_root = project_root,
      script_path = "run_analysis.R",
      args = c("--mode", "generate"),
      label = "reference manuscript generation",
      log_file = log_file
    )
    run_rscript(
      project_root = project_root,
      script_path = "run_analysis.R",
      args = c("--mode", "validate"),
      label = "reference generated artifact validation",
      log_file = log_file
    )
  } else if (isTRUE(args$execute)) {
    message_log(
      "Skipping separate reproducible_workflow generation/validation. ",
      "Add --run-reference-workflow to run it as an additional reference check.",
      log_file = log_file
    )
  } else {
    message_log(
      "Dry-run: reference generation/validation not executed. Add --execute --run-reference-workflow to run that optional check.",
      log_file = log_file
    )
  }

  run_rscript(
    project_root = project_root,
    script_path = file.path("Scripts_2025", "Final_Scripts", "validate_manuscript_outputs.R"),
    label = "direct manuscript-output validation",
    log_file = log_file
  )

  message_log(
    "Direct manuscript outputs are in Scripts_2025/Final_Scripts/final_manuscript_objects/.",
    log_file = log_file
  )

  stage_paths <- fs_write_stage_artifact_map(project_root)
  message_log("Stage-ordered map: ", stage_paths$md, log_file = log_file)
  message_log("Workflow command complete", log_file = log_file)
}

main()
