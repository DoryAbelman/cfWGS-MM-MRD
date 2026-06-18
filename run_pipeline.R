# =============================================================================
# run_pipeline.R
#
# Command-line runner for the original numbered cfWGS-MM-MRD scripts.
#
# This file intentionally lives beside the original scripts and runs those
# scripts directly. It does not replace the original analysis code. The separate
# reproducible_workflow/ directory is used as a provenance/validation harness.
#
# Usage from the project root:
#   Rscript Scripts_2025/Final_Scripts/run_pipeline.R
#   Rscript Scripts_2025/Final_Scripts/run_pipeline.R --execute
#
# Usage from this directory:
#   Rscript run_pipeline.R
#   Rscript run_pipeline.R --execute --from 2_0 --to 2_4
#
# Purpose:
#   Dry-run or execute the original numbered scripts in source-pipeline order,
#   with guardrails for cache-sensitive/model-training stages.
#
# Manuscript outputs created/updated:
#   - None directly. This support runner documents and executes the original
#     numbered script order but does not itself contain scientific analysis.
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
  list(
    execute = "--execute" %in% args,
    include_cache_sensitive = "--include-cache-sensitive" %in% args,
    keep_going = "--keep-going" %in% args,
    check_packages = "--check-packages" %in% args,
    from = parse_flag_value(args, "--from"),
    to = parse_flag_value(args, "--to"),
    only = parse_flag_value(args, "--only")
  )
}

select_plan_rows <- function(plan, args) {
  selected <- plan
  if (!is.null(args$only)) {
    wanted <- trimws(strsplit(args$only, ",", fixed = TRUE)[[1]])
    selected <- selected[selected$script_id %in% wanted | selected$script %in% wanted, , drop = FALSE]
  }
  if (!is.null(args$from)) {
    if (!args$from %in% plan$script_id) stop("--from does not match a script_id: ", args$from, call. = FALSE)
    selected <- selected[selected$order >= plan$order[match(args$from, plan$script_id)], , drop = FALSE]
  }
  if (!is.null(args$to)) {
    if (!args$to %in% plan$script_id) stop("--to does not match a script_id: ", args$to, call. = FALSE)
    selected <- selected[selected$order <= plan$order[match(args$to, plan$script_id)], , drop = FALSE]
  }
  if (!args$include_cache_sensitive) {
    selected <- selected[selected$run_policy != "cache_sensitive", , drop = FALSE]
  }
  selected[order(selected$order), , drop = FALSE]
}

check_packages <- function(script_dir) {
  config_env <- new.env(parent = baseenv())
  sys.source(file.path(script_dir, "config.R"), envir = config_env)
  missing <- config_env$packages[!vapply(config_env$packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R packages:\n",
      paste(missing, collapse = "\n"),
      "\nInstall missing packages in the analysis environment before running --execute.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

run_one_script <- function(project_root, script_dir, script_name, log_path) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)
  script_path <- file.path(script_dir, script_name)
  if (!file.exists(script_path)) stop("Missing script: ", script_path, call. = FALSE)
  system2(file.path(R.home("bin"), "Rscript"), args = shQuote(script_path), stdout = log_path, stderr = log_path)
}

main <- function() {
  args <- parse_args()
  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
  source(file.path(script_dir, "pipeline_metadata.R"))

  plan <- fs_read_source_pipeline(project_root)
  source_map <- fs_read_artifact_source_map(project_root, required = TRUE)
  plan <- fs_annotate_plan_with_outputs(plan, source_map)
  script_index_path <- file.path(script_dir, "script_index.tsv")
  utils::write.table(plan, script_index_path, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  stage_map_paths <- fs_write_stage_artifact_map(project_root)

  selected <- select_plan_rows(plan, args)
  if (!nrow(selected)) stop("No scripts selected.", call. = FALSE)

  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(script_dir, "pipeline_logs", run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  main_log <- file.path(run_dir, "run_pipeline.log")

  message_log("Project root: ", project_root, log_file = main_log)
  message_log("Script directory: ", script_dir, log_file = main_log)
  message_log("Script index: ", script_index_path, log_file = main_log)
  message_log("Stage-ordered artifact map: ", stage_map_paths$tsv, log_file = main_log)
  message_log("Selected scripts: ", nrow(selected), log_file = main_log)
  if (!args$execute) message_log("Dry run only. Add --execute to run scripts.", log_file = main_log)
  if (!args$include_cache_sensitive) {
    message_log(
      "Cache-sensitive nested-CV/model training is skipped by default. Add --include-cache-sensitive only for deliberate recomputation; regenerated artifacts may differ unless the exact inputs, seeds, package versions, and execution environment are preserved.",
      log_file = main_log
    )
  }

  if (args$check_packages) {
    message_log("Checking required R packages from config.R", log_file = main_log)
    check_packages(script_dir)
    message_log("Package check passed", log_file = main_log)
  }

  manifest <- selected
  manifest$status <- if (args$execute) "pending" else "dry_run_not_executed"
  manifest$status_code <- NA_integer_
  manifest$started_at <- ""
  manifest$finished_at <- ""
  manifest$elapsed_seconds <- NA_real_
  manifest$log_path <- ""

  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    log_path <- file.path(run_dir, paste0(sprintf("%02d", row$order), "_", row$script_id, ".log"))
    manifest$log_path[i] <- log_path
    outputs <- row$manuscript_outputs
    if (is.na(outputs) || !nzchar(outputs)) {
      outputs <- "no direct mapped manuscript figure/table; upstream or intermediate dependency"
    }

    if (!args$execute) {
      message_log(
        "DRY RUN: ", row$script_id, " | ", row$script, " | ", row$stage,
        " | manuscript outputs: ", outputs,
        log_file = main_log
      )
      next
    }

    message_log("Running ", row$script_id, ": ", row$script, log_file = main_log)
    message_log("Manuscript outputs: ", outputs, log_file = main_log)
    start <- Sys.time()
    manifest$started_at[i] <- format(start, "%Y-%m-%d %H:%M:%S")
    status <- run_one_script(project_root, script_dir, row$script, log_path)
    finish <- Sys.time()
    manifest$finished_at[i] <- format(finish, "%Y-%m-%d %H:%M:%S")
    manifest$elapsed_seconds[i] <- as.numeric(difftime(finish, start, units = "secs"))
    manifest$status_code[i] <- status
    manifest$status[i] <- if (identical(status, 0L)) "success" else "failed"
    message_log("Finished ", row$script_id, " with status ", status, log_file = main_log)

    if (!identical(status, 0L) && !args$keep_going) {
      utils::write.table(manifest, file.path(run_dir, "run_manifest.tsv"), sep = "\t", row.names = FALSE, quote = TRUE, na = "")
      stop("Pipeline stopped after failed script ", row$script_id, ". See log: ", log_path, call. = FALSE)
    }
  }

  manifest_path <- file.path(run_dir, "run_manifest.tsv")
  utils::write.table(manifest, manifest_path, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  message_log("Run manifest: ", manifest_path, log_file = main_log)
  message_log("Done", log_file = main_log)
}

main()
