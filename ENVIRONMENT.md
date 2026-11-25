# Reproduction environment

The analysis used 64-bit Windows, Python 3.12.10, R 4.6.0 (ucrt), and Arial for the finalized figures. Python modules are installed in the project's `.venv`. The R packages in the table below are in `.venv/R/library`; recommended R packages ship with the R runtime.

The confirmed inputs are MC-MED v1.0.1, MIMIC-IV-ED v2.2 and MIMIC-IV v3.1. The access routes, archive layout and required filenames are listed in [README.md](README.md).

## Python

Create the virtual environment with Python 3.12.10. Install the recorded direct dependencies in that environment.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install numpy==2.5.5 pandas==2.3.3 scipy==1.18.1 duckdb==2.0.1 pyarrow==26.0.1 statsmodels==0.15.0 patsy==1.0.3 lifelines==0.30.3 matplotlib==3.11.2 pillow==13.0.0
```

Use the exact listed versions for comparisons to the study artifacts. Pinning compiled-library versions and preserving platform and language settings helps reproduce numerical optimizer outputs. The P2/P6 R launcher runs with one BLAS/OpenMP thread and a fixed random seed.

## R

Install R 4.6.0 and add `Rscript` to `PATH`. From the project root, create `.venv/R/library`, then install the matching Windows binaries into that library. The package versions below reflect the tested environment. Check the CRAN archives or an authorized binary mirror if a listed binary version is no longer current. The R runtime supplied `Matrix` 1.7-5 and the recommended `MASS` and `nlme` (3.1-169) packages.

| Package | Version | Package | Version |
| --- | --- | --- | --- |
| lcmm | 2.2.2 | lmerTest | 3.2-1 |
| lme4 | 2.1-0 | ggplot2 | 4.0.4 |
| cli | 3.6.6 | cpp11 | 0.5.5 |
| doParallel | 1.0.17 | farver | 2.1.2 |
| foreach | 1.5.2 | glue | 1.8.1 |
| gtable | 0.3.6 | isoband | 0.3.0 |
| iterators | 1.0.14 | labeling | 0.4.3 |
| lifecycle | 1.0.5 | marqLevAlg | 2.0.8 |
| minqa | 1.2.8 | mvtnorm | 1.4-2 |
| nloptr | 2.2.1 | numDeriv | 2016.8-1.1 |
| R6 | 2.6.1 | rbibutils | 2.4.1 |
| RColorBrewer | 1.1-3 | Rcpp | 1.1.2 |
| RcppEigen | 0.3.4.0.2 | Rdpack | 2.6.6 |
| reformulas | 0.4.4 | rlang | 1.3.0 |
| S7 | 1.0.0 | scales | 1.4.0 |
| spacefillr | 0.4.0 | vctrs | 0.7.3 |
| viridisLite | 0.4.3 | withr | 3.0.3 |

R launchers take the local library path from the project layout. Figure generation checks for local Arial and exports PDF with embedded fonts and matched 600 dpi PNG/TIFF images. See [README.md](README.md) for source files and the command order.

## Submission figures and supplementary analyses

The submission step uses the same Python and R environment. Pillow exports RGB TIFFs at 600 dpi with LZW compression, and Matplotlib writes PDFs with embedded Arial fonts. The post hoc hemoglobin extraction uses the existing DuckDB installation.
