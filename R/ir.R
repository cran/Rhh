`ir` <-
function(genotypes) {

   genotypes <- as.matrix(genotypes)

   individuals <- nrow(genotypes)
   loci <- ncol(genotypes) / 2
   ir <- array(NA, dim=c(individuals, 1))
   frequencies <- array()

   for (l in 1:loci) {
      g <- 2 * l - 1
      h <- 2 * l
      frequencies[l] <- list(table(genotypes[, g:h]))
   }

   for (i in 1:individuals) {

      H <- 0
      N <- 0
      f <- 0

      for (l in 1:loci) {

         g <- 2 * l - 1
         h <- 2 * l

         if ((!is.na(genotypes[i, g])) && (!is.na(genotypes[i, h]))) {
            N <- N + 1
            if (genotypes[i, g] == genotypes[i, h]) {
               H <- H + 1
               c <- as.character(genotypes[i, g])
               f <- f + (2 * frequencies[[l]][[c]] - 2) / (sum(frequencies[[l]]) - 2)
            }
            else {
               c <- as.character(genotypes[i, g])
               f <- f + (frequencies[[l]][[c]] - 1) / (sum(frequencies[[l]]) - 2)
               c <- as.character(genotypes[i, h])
               f <- f + (frequencies[[l]][[c]] - 1) / (sum(frequencies[[l]]) - 2)
            }
         }
      }

      ir[i] <- (2 * H - f) / (2 * N - f)

   }

   ir

}

