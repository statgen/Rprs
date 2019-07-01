prs_gw <- function(weights_file, genotypes_file, samples_file, pvalue = 1.0, weight_col = 'WEIGHT_') {
   if (is.null(weights_file)) {
      stop("Provide weights file.")
   }
   if (is.null(genotypes_file)) {
      stop("Provide VCF/BCV/savvy file with genotypes and dosages.")
   }

   lookahead <- read.table(weights_file, header = TRUE, nrows = 1, stringsAsFactors = FALSE) # read first line to determine columns and setup their types (needed for faster read.table)
   col_classes <- sapply(lookahead, class)
   col_classes[c("CHROM", "OA", "EA")] <- "character"
   col_classes["POS"] <- "integer"
   col_classes["PVALUE"] <- "numeric"
   col_classes[which(startsWith(colnames(lookahead), weight_col))] <- "numeric"
   col_classes[-c(which(colnames(lookahead) %in% c("CHROM", "POS", "OA", "EA", "PVALUE")), which(startsWith(colnames(lookahead), weight_col)))] <- list(NULL)

   weights <- read.table(weights_file, header = TRUE, stringsAsFactors = FALSE, colClasses = col_classes)
   weights <- weights[with(weights, order(CHROM, POS, decreasing = FALSE)),] # ensure that positions are sorted
   weights <- weights[weights$PVALUE < pvalue, ] # subset to weights which satisfy p-value threshold

   # set samples to use
   if (!is.null(samples_file)) {
      samples <- readLines(samples_file)
   } else {
      samples <- as.character(c()) # use all samples
   }

   genotypes <- new(Genotypes, genotypes_file) # create VCF/BCF/savvy reader
   genotypes$open(samples) # open file for reading and set samples
   geno_samples <- genotypes$get_selected_samples() # get samples in the same order as in the genotype file

   weight_cols <- colnames(weights)[which(startsWith(colnames(weights), weight_col))]
   individual_prs <- matrix(nrow = length(geno_samples), ncol = length(weight_cols), data = 0, dimnames = list(samples = geno_samples, scores = paste("PRS_", weight_cols, sep="")))
   
   for (row_idx in 1:nrow(weights)) {
      ea <- weights[row_idx, "EA"]
      oa <- weights[row_idx, "OA"] 
      dosage <- genotypes$read_variant(weights[row_idx, "CHROM"], weights[row_idx, "POS"], ea, oa)
      if (nrow(dosage) == 0) {
         next
      }
      for (i in seq_along(weight_cols)) {
         weight <- weights[row_idx, weight_cols[i]]
         if (weight >= 0) {
            individual_prs[, i] <- individual_prs[, i] + weight * dosage
         } else {
            individual_prs[, i] <- individual_prs[, i] + weight * (dosage - 2.0)
         }
      }
      if (row_idx %% 10000 == 0) {
         message(sprintf("Processed %d/%d variants", row_idx, nrow(weights)))
      }
   }
   message("Done")
   d <- cbind(IID = dimnames(individual_prs)$samples, data.frame(individual_prs, row.names = NULL))
   if (length(samples) > 0) {
      d <- d[samples, ] # order samples as in input samples file
   }
   return(d)  
}

prs_pt <- function(weights_file, genotypes_file, samples_file, pvalues = c(1.0)) {
   if (is.null(weights_file)) {
      stop("Provide weights file.")
   }
   if (is.null(genotypes_file)) {
      stop("Provide VCF/BCV/savvy file with genotypes and dosages.")
   }

   weights <- read.table(weights_file, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "integer", "character", "character", "numeric", "numeric"))
   weights <- weights[with(weights, order(CHROM, POS, decreasing = FALSE)),] # ensure that positions are sorted 
   
   
   weights <- weights[weights$PVALUE <= max(pvalues), ] # subset to weights which satisfy all p-value thresholds
   pvalues <- sort(pvalues, decreasing = TRUE)

   # set samples to use
   if (!is.null(samples_file)) {
      samples <- readLines(samples_file)
   } else {
      samples <- as.character(c()) # use all samples
   }

   genotypes <- new(Genotypes, genotypes_file) # create VCF/BCF/savvy reader
   genotypes$open(samples) # open file for reading and set samples
   geno_samples <- genotypes$get_selected_samples() # get samples in the same order as in the genotype file

   individual_prs <- matrix(nrow = length(geno_samples), ncol = length(pvalues), data = 0, dimnames = list(samples = geno_samples, scores = paste("PRS_", pvalues, sep="")))
   for (row_idx in 1:nrow(weights)) {
      weight <- weights[row_idx, "WEIGHT"]
      if (weight >= 0) {
         ra <- weights[row_idx, "EA"]
         pa <- weights[row_idx, "OA"]
         risk <- weight
      } else {
         ra <- weights[row_idx, "OA"]
         pa <- weights[row_idx, "EA"]
         risk <- -1 * weight
      }
      pvalue <- weights[row_idx, 'PVALUE']
      if (row_idx %% 10000 == 0) {
         message(sprintf("Processed %d/%d variants", row_idx, nrow(weights)))
      }
      dosage <- genotypes$read_variant(weights[row_idx, "CHROM"], weights[row_idx, "POS"], ra, pa)
      if (nrow(dosage) == 0) {
         next
      }
      for (i in seq_along(pvalues)) {
         if (pvalue < pvalues[i]) {
            individual_prs[, i] <- individual_prs[, i] + risk * dosage
         } else {
            break
         }
      }
   }
   message("Done")
   if (length(samples) > 0) {
      individual_prs <- individual_prs[samples, ] # order samples as in input samples file
   }
   
   d <- cbind(IID = dimnames(individual_prs)$samples, data.frame(individual_prs))
   return(d)
}
