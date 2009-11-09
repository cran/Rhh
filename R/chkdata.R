`chkdata` <-
function(genotypes) {

   genotypes <- as.matrix(genotypes)

   if (ncol(genotypes) %% 2 == 1) {
      stop("Odd number of columns in the input.")
   }

   individuals <- nrow(genotypes)
   loci <- ncol(genotypes) / 2
   n_alleles <- array(loci)

   for (l in 1:loci) {
      g <- 2 * l - 1
      h <- 2 * l
      frequencies <- list(table(genotypes[, g:h]))
      n_alleles[l] <- length(frequencies[[1]])
   }

   if (sum(n_alleles < 2) > 0) {
      loci_list <- ""
      for (i in 1:loci) {
         if (n_alleles[i] < 2) {
            if (loci_list == "") {
               loci_list <- i
            }
            else {
               loci_list <- paste(loci_list, i, sep = ", ")
            }
         }
      }
      stop(paste("Only one allele in loci:", loci_list))
   }

   write(n_alleles, file="number_of_alleles.txt", ncolumns=1)

   k <- 0
   ind_loci_list <- "\n"

   for (i in 1:individuals) {

      for (l in 1:loci) {

         g <- 2 * l - 1
         h <- 2 * l

         if (xor(!is.na(genotypes[i, g]), !is.na(genotypes[i, h]))) {
            ind_loci_list <- paste(ind_loci_list, i, l, "\n")
            k <- k + 1
         }
      }
   }

   if (k > 0) {
      warning(paste("One or more individuals are missing one allele in one or more loci:\nIndividual  Locus", ind_loci_list))
   }

}

