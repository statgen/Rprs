# Rprs
R package for computing Polygenic Risk Scores (PRS).
Supports VCF, BCF and savvy input formats.

## Install
Clone repository using `git clone` and run the following command:
```
R CMD INSTALL Rprs
```

## Run
```
library(Rprs)
p <- prs("weights.txt", "dosages.vcf.gz", "samples.txt")
write.table(p, col.names = TRUE, row.names = FALSE, sep = "\t", quote = F)
```
