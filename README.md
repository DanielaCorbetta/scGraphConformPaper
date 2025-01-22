# Code for conformal prediction for cell type annotation with graph-structured constrains

Analysis for the paper on [conformal prediction for cell type annotation with graph-structured constrains](https://arxiv.org/abs/2410.23786)



## Files
- `1_ConformalvsGraph.Rmd`: code to reproduce results of Section 4
- `2_resampling.Rmd`: code to reproduce results of Section 4.1
- `3_covid.Rmd`: code to reproduce results of Section 5
  - data: `sce_Covid_Bcells.rds` are openly available in `Covid case/objects/sce Covid Bcells.rds` at 
<https://doi.org/10.5281/zenodo.10391097>.
  - code for model: `params_EM_81929.rda`, `cell_type_identification4_clean.R` downloaded from <https://github.com/igrabski/scRNAseq-cell-type>
