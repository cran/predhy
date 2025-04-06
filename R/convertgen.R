#' @title Convert Genotype
#' @description Convert genotypes in HapMap format or in numeric format for predhy package.
#' @param input_geno  genotype in HapMap format or in numeric format. The names of individuals should be provided. Missing (NA) values are allowed.
#' @param type the type of genotype. There are three options: "hmp1" for genotypes in HapMap format with single bit, "hmp2" for genotypes in HapMap format with double bit, and "num" for genotypes in numeric format.
#' @param missingrate max missing percentage for each SNP, default is 0.2.
#' @param maf minor allele frequency for each SNP, default is 0.05.
#' @param impute logical. If TRUE, imputation. Default is TRUE.
#'
#' @return A matrix of genotypes in numeric format, coded as 1, 0, -1 for AA, Aa, aa. Each row represents an individual and each column represents a marker. The rownames of the matrix are the names of individuals.
#' @export
#'
#' @examples
#' ## load genotype in HapMap format with double bit
#' data(input_geno)
#'
#' ## convert genotype for hypred package
#' inbred_gen <- convertgen(input_geno, type = "hmp2")
#'
#'
#' ## load genotype in numeric format
#' data(input_geno1)
#' head(input_geno1)
#'
#' ## convert genotype for hypred package
#' inbred_gen1 <- convertgen(input_geno1, type = "num")
#'
#'
convertgen<-function(input_geno, type = c("hmp1", "hmp2", "num"), missingrate = 0.2, maf = 0.05, impute = TRUE) {
  {
    if (type == "num") {
      genotype <- as.matrix(input_geno)
	  rownames(genotype)<-input_geno[,1]
      gen_missingrate <- apply(genotype, 1, function(x) {
        rate <- length(which(is.na(x)))/ncol(genotype)
        return(rate)
      })
      num_filter <- which(gen_missingrate <= missingrate)
      gen_filter <- genotype[num_filter, ]
	   row.names(gen_filter)<-genotype[,1]
      gen <- apply(gen_filter, 2, function(x) {
        x1 <- as.numeric(x)
        return(x1)
      })
      if (impute) {
        gene1 <- apply(gen, 1, function(x) {
          x[which(is.na(x))] <- mean(x, na.rm = TRUE)
          return(x)
        })
      } else {
        gene1 <- t(gen)
      }
    }
    if (type != "num") {
      input_geno[input_geno == "N"] = NA
      input_geno[input_geno == "NN"] = NA
      genotype <- input_geno[, -c(1:11)]
	  rownames(genotype)<-input_geno[,1]
      gen_missingrate <- apply(genotype, 1, function(x) {
        rate <- length(which(is.na(x)))/ncol(genotype)
        return(rate)
      })
      num_filter <- which(gen_missingrate <= missingrate)
      gen_filter <- genotype[num_filter, ]
      map <- input_geno[num_filter, c(1:4)]
      x1 <- substr(map[, 2], 1, 1)
      x2 <- substr(map[, 2], 3, 3)
      aa <- paste(x1, x1, sep = "")
      bb <- paste(x2, x2, sep = "")
      cc <- paste(x1, x2, sep = "")
      if (type == "hmp1") {
        gen <- apply(gen_filter, 2, function(x) {
          x[x == x1] <- 1
          x[x == x2] <- -1
          x[x %in% c("R", "Y", "S", "W", "K", "M")] <- 0
          x1 <- as.numeric(x)
          return(x1)
        })
      }
      if (type == "hmp2") {
        gen <- apply(gen_filter, 2, function(x) {
          x[x == aa] <- 1
          x[x == bb] <- -1
          x[x == cc] <- 0
          x1 <- as.numeric(x)
          return(x1)
        })
      }
       rownames(gen)<-rownames(gen_filter )
      if (impute) {
        gene1 <- apply(gen, 1, function(x) {
          x[which(is.na(x))] <- mean(x, na.rm = TRUE)
          return(x)
        })
      } else {
        gene1 <- t(gen)
      }
    }
  }
  gen_maf <- apply(gene1, 2, function(x) {
    #maf_rate <- (length(which(x == -1)) + (1/2 * (length(which(x == 0)))))/nrow(gene1)
    maf_rate <- (length(which(x == -1)) + (1/2 * (length(which(x == 0)))))/ (length(which(x == -1))+length(which(x == 0))+length(which(x == 1)))
    return(maf_rate)
  })
  num_gene1 <- which(gen_maf >= maf)
  gene2 <- gene1[, num_gene1]
 
 
  return(gene2)
}
