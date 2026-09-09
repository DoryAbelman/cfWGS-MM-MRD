# Contributing

Thank you for helping improve this research software. Changes should be small,
reviewable, and tied to a scientific or reproducibility need.

## Before opening a change

1. Do not add protected clinical records, raw sequencing data, names, exact
   collection dates, reidentification keys, tokens, workstation-specific
   paths, or identifiers other than the manuscript display IDs approved for
   this repository.
2. State which stage an included derived table represents. Panel-level source
   data do not reproduce the upstream FASTQ/BAM processing steps.
3. Preserve patient-grouped data splitting. Any model-evaluation change must
   demonstrate that all samples from a patient remain in one fold.
4. Set and document random seeds, software versions, input checksums, analysis
   denominators, and the unit of analysis.
5. Never overwrite locked manuscript figures with a new run. Write new outputs
   separately and compare them explicitly.
6. Update `docs/SCRIPT_GUIDE.md` and `docs/FIGURE_TABLE_MAP.md` when a script's
   role or manuscript output changes.

## Checks

Before proposing a documentation-only change, inspect the diff and parse every
R script that changed:

```bash
git diff --check
Rscript -e 'parse(file = "path/to/changed_script.R")'
```

Parse and test every R script that was changed. For analysis changes, rerun the
affected workflow against a fixed input snapshot and compare the new tables,
figures, and row-level outputs with the previous version. Changes to the model
validation require a reduced smoke test followed by the complete 32-model,
50-repeat run before release; label reduced-run results as tests rather than
manuscript estimates.

## Reporting scientific changes

Describe the input version, cohort/denominator, grouping unit, random seed,
changed files, generated artifacts, QC results, and whether manuscript values
or only presentation changed. State known limitations directly.
