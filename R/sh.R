`sh` <-
function(genotypes) {

   genotypes <- as.matrix(genotypes)

   individuals <- nrow(genotypes)
   loci <- ncol(genotypes) / 2
   sh <- array(NA, dim=c(individuals, 1))
   heterozygosity <- array()

   for (l in 1:loci) {
      g <- 2 * l - 1
      h <- 2 * l
      freq <- table((genotypes[, g] == genotypes[, h]))
      if (nrow(freq) == 2) {
         heterozygosity[l] <- sweep(freq, 1, sum(freq), "/")[["FALSE"]]
      }
      else {
         if (labels(freq)[[1]][1] == TRUE) {
            heterozygosity[l] <- 1 - sweep(freq, 1, sum(freq), "/")[["TRUE"]]
         }
         else {
            heterozygosity[l] <- sweep(freq, 1, sum(freq), "/")[["FALSE"]]
         }
      }

   }

   for (i in 1:individuals) {

      H <- 0
      N <- 0
      mh <- 0

      for (l in 1:loci) {

         g <- 2 * l - 1
         h <- 2 * l

         if ((!is.na(genotypes[i, g])) && (!is.na(genotypes[i, h]))) {
            mh <- mh + heterozygosity[l] 
            N <- N + 1
            if (genotypes[i, g] != genotypes[i, h]) {
               H <- H + 1
            }
         }
      }

      sh[i] <- (H / N) / (mh / N)

   }

   sh

}

