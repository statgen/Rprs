prs <- function(weights_file, genotypes_file, samples_file, pvalue = 1) {
   if (is.null(weights_file)) {
      stop("Provide weights file.")
   }
   if (is.null(genotypes_file)) {
      stop("Provide VCF/BCV/savvy file with genotypes and dosages.")
   }

   weights <- read.table(weights_file, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "integer", "character", "character", "numeric", "numeric"))
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

   individual_prs <- matrix(nrow = length(samples), 1, data = 0, dimnames = list(samples = samples, scores = c("PRS")))
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
      if (row_idx %% 10000 == 0) {
         message(sprintf("Processed %d/%d variants", row_idx, nrow(weights)))
      }
      dosage <- genotypes$read_variant(weights[row_idx, "CHROM"], weights[row_idx, "POS"], ra, pa)
      if (nrow(dosage) == 0) {
         next
      }
      individual_prs <- individual_prs + risk * dosage
   }
   message("Done")
   d <- cbind(IID = dimnames(individual_prs)$samples, data.frame(individual_prs, row.names = NULL))
   return(d)
}
