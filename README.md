# Rprs
R package for computing Polygenic Risk Scores (PRS) from dosages (DS).
Supports VCF, BCF and savvy input formats.

## Install
Clone repository using `git clone` and run the following command:
```
R CMD INSTALL Rprs
```

## Run

### Pruning and Thresholding (PT) method.

- Weights file must have the following columns (in the same order): CHROM, POS, OA, EA, WEIGHT, PVALUE.
- If no samples file is provided (i.e. NULL), then all samples are used.

```
library(Rprs)
p <- prs_pt("weights.txt", "dosages.vcf.gz", "samples.txt", pvalues = c(1.0, 5e-8))
write.table(p, col.names = TRUE, row.names = FALSE, sep = "\t", quote = F)
```

### Genome-wide PRS method.

- Weights file must have the following columns (in arbitrary order): CHROM, POS, OA, EA, PVALUE, and multiple columns for weights.
- The weight columns must have the sampe prefix (e.g. "LDpred_", "WEIGHT_", or similar). You can specify this prefix using the `weight_col` function argument.
- If no samples file is provided (i.e. NULL), then all samples are used.

```
library(Rprs)
p <- prs_gw("weights.txt", "dosages.vcf.gz", "samples.txt", pvalue = 1.0, weight_col = "LDpred_")
write.table(p, col.names = TRUE, row.names = FALSE, sep = "\t", quote = F)
```
