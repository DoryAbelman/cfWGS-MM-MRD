# =============================================================================
# run_pipeline.R
#
# Command-line runner for the numbered cfWGS-MM-MRD analysis scripts.
#
# This file intentionally lives beside the numbered analysis scripts and runs
# them directly. Scientific logic remains in those stage scripts; this runner
# only selects the script range, applies cache-sensitive guardrails, records
# logs/manifests, and reports which manuscript outputs each stage affects.
#
# Usage from the project root:
#   Rscript Scripts_2025/Final_Scripts/run_pipeline.R
#   Rscript Scripts_2025/Final_Scripts/run_pipeline.R --execute
#
# Usage from this directory:
#   Rscript run_pipeline.R
#   Rscript run_pipeline.R --execute --from 2_0 --to 2_4
#   Rscript run_pipeline.R --execute --only 5_3 --include-post-cv
#
# Purpose:
#   Dry-run or execute the retained numbered scripts in source-pipeline order.
#
# Workflow scope:
#   The executed scripts come only from config/source_pipeline.tsv. The default
#   plan covers the retained preprocessing, preserved-model scoring, and
#   manuscript figure/table generation scripts. The 50-repeat
#   patient-grouped nested-CV workflow remains separate because it requires
#   three parameterized long-running fits with unique run IDs. After those fits
#   and their plots are complete, add --include-post-cv to rebuild and validate
#   the figure source-data workbooks through the 5_3 runner.
#
# Inputs:
#   * config/source_pipeline.tsv supplies the ordered script list and run policy.
#   * config.R is read only when --check-packages is requested.
#   * Each selected numbered script has its own data inputs and assumptions.
#
# Outputs:
#   In --execute mode, each numbered script writes its normal outputs. This
#   runner additionally writes script_index.tsv and a timestamped
#   pipeline_logs/<run-id>/ directory containing individual logs and
#   run_manifest.tsv. Dry-run mode writes nothing.
#
# Manuscript outputs created/updated:
#   - None directly. This support runner documents and executes the numbered
#     script order but does not itself contain scientific analysis.
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
    include_post_cv = "--include-post-cv" %in% args,
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
  if (!args$include_post_cv) {
    selected <- selected[selected$run_policy != "post_cv", , drop = FALSE]
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
  source(file.path(script_dir, "pipeline_metadata.R"))
  project_root <- fs_project_root_from_script_dir(script_dir)

  plan <- fs_read_source_pipeline(project_root)

  selected <- select_plan_rows(plan, args)
  if (!nrow(selected)) stop("No scripts selected.", call. = FALSE)

  allowed_policies <- c("run", "post_cv")
  unknown_policies <- setdiff(unique(plan$run_policy), allowed_policies)
  if (length(unknown_policies)) {
    stop("Unknown run_policy value(s): ", paste(unknown_policies, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(plan$order)) stop("config/source_pipeline.tsv contains duplicate order values.", call. = FALSE)
  if (anyDuplicated(plan$script_id)) stop("config/source_pipeline.tsv contains duplicate script_id values.", call. = FALSE)
  missing_scripts <- selected$script[!file.exists(file.path(script_dir, selected$script))]
  if (length(missing_scripts)) {
    stop("Selected pipeline script(s) are missing: ", paste(missing_scripts, collapse = ", "), call. = FALSE)
  }

  if (args$check_packages) {
    message_log("Checking required R packages from config.R")
    check_packages(script_dir)
    message_log("Package check passed")
  }

  if (!args$execute) {
    message_log("Project root: ", project_root)
    message_log("Read-only dry run. Add --execute to run scripts.")
    if (!args$include_post_cv) {
      message_log("Post-CV source-workbook rebuilding is excluded. Add --include-post-cv after grouped-CV outputs are complete.")
    }
    for (i in seq_len(nrow(selected))) {
      message_log(
        "DRY RUN: ", selected$script_id[i], " | ", selected$script[i],
        " | ", selected$stage[i], " | purpose: ", selected$notes[i]
      )
    }
    message_log("No files were written.")
    return(invisible(selected))
  }

  script_index_path <- file.path(script_dir, "script_index.tsv")
  utils::write.table(plan, script_index_path, sep = "\t", row.names = FALSE, quote = TRUE, na = "")

  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(script_dir, "pipeline_logs", run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  main_log <- file.path(run_dir, "run_pipeline.log")

  message_log("Project root: ", project_root, log_file = main_log)
  message_log("Script directory: ", script_dir, log_file = main_log)
  message_log("Script index: ", script_index_path, log_file = main_log)
  message_log("Selected scripts: ", nrow(selected), log_file = main_log)
  if (!args$include_post_cv) {
    message_log(
      "Post-CV source-workbook rebuilding is skipped. Add --include-post-cv after grouped-CV outputs are complete.",
      log_file = main_log
    )
  }

  manifest <- selected
  manifest$status <- "pending"
  manifest$status_code <- NA_integer_
  manifest$started_at <- ""
  manifest$finished_at <- ""
  manifest$elapsed_seconds <- NA_real_
  manifest$log_path <- ""

  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    log_path <- file.path(run_dir, paste0(sprintf("%02d", row$order), "_", row$script_id, ".log"))
    manifest$log_path[i] <- log_path

    message_log("Running ", row$script_id, ": ", row$script, log_file = main_log)
    message_log("Purpose: ", row$notes, log_file = main_log)
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
