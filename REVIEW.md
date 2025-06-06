# Code Review Notes

This repository contains a large collection of R scripts used for processing clinical and genomic data. While the scripts are functional, the following observations may help improve maintainability and reduce redundancy.

## General suggestions

1. **Remove generated files** – `.DS_Store` is included in the repository and should be deleted and added to `.gitignore`.
2. **Consistent file naming** – several scripts include spaces or misspellings (e.g. `...inetegrate...`). Normalising names will make sourcing and execution easier.
3. **Centralise configuration** – many scripts hard-code input paths. Consider creating a single configuration script or using the `here` package so paths only need to be set once.
4. **Library management** – each script loads packages independently, sometimes with duplicate `library()` calls. A shared setup script would reduce repetition.
5. **Factor out common functions** – operations such as cleaning sample IDs or extracting timepoints appear multiple times. Packaging these as helper functions would avoid code duplication.
6. **Error handling** – some scripts assume input columns always exist. Adding checks with informative error messages would make failures easier to diagnose.
7. **Documentation** – although many scripts have comment blocks at the top, consider adding more detailed explanations (inputs/outputs) and example commands for reproducibility.

## Potential errors

- `1_7B_Process_fragment_score_and_inetegrate_fragmentomics_data.R` contains a typo in the word *integrate* within the filename.
- `2_3_Feature_Concordance_And_Mutation_Counts.R` loads `Hmisc` twice. Removing the duplicate call is recommended.
- Several file names contain spaces (e.g. `2_1_Clinical Demographics Table.R`) which can cause issues when running scripts non-interactively.

These points are not exhaustive but should help guide future refactoring.
