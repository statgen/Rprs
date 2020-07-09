# Rprs
R package for computing Polygenic Risk Scores (PRS) from dosages (DS).
Supports VCF and BCF input formats.

## Install
Clone repository using `git clone` and run the following command:
```
R CMD INSTALL Rprs
```

## Run

### P-value Thresholding (PT) method.

- Uses an LD clumped weight file as input
- Weights file must have the following columns (in the same order): CHROM, POS, OA, EA, WEIGHT, PVALUE.
- If no samples file is provided (i.e. NULL), then all samples are used.

```
library(Rprs)
p <- prs_pt("weights.txt", "dosages.vcf.gz", "samples.txt", pvalues = c(1.0, 5e-8))
write.table(p, col.names = TRUE, row.names = FALSE, sep = "\t", quote = F)
```

### Genome-wide PRS method.

- Weights file must have the following columns (in arbitrary order): CHROM, POS, OA, EA, and multiple columns for weights.
- All weight columns must have the same prefix (e.g. "LDpred_", "LDpred", "WEIGHT_", "WEIGHT", or similar). You must specify this prefix using the `weight_col` function argument. Default is `WEIGHT_`.
- If no samples file is provided (i.e. NULL), then all samples are used.

```
library(Rprs)
p <- prs_gw("weights.txt", "dosages.vcf.gz", "samples.txt", weight_col = "LDpred_")
write.table(p, col.names = TRUE, row.names = FALSE, sep = "\t", quote = F)
```
