`hl` <-
function(genotypes) {

   genotypes <- as.matrix(genotypes)

   individuals <- nrow(genotypes)
   loci <- ncol(genotypes) / 2
   hl <- array(NA, dim=c(individuals, 1))
   E <- array(loci)
   frequencies <- array(loci)

   for (l in 1:loci) {
      E[l] <- 1
      g <- 2 * l - 1
      h <- 2 * l
      frequencies[l] <- list(table(genotypes[, g:h]))
      E[l] <- 1 - sum((frequencies[[l]] / sum(frequencies[[l]]))^2)
   }

   for (i in 1:individuals) {

      sum.Eh <- 0
      sum.Ej <- 0

      for (l in 1:loci) {

         g <- 2 * l - 1
         h <- 2 * l

         if ((!is.na(genotypes[i, g])) && (!is.na(genotypes[i, h]))) {
            if (genotypes[i, g] == genotypes[i, h]) {
               sum.Eh <- sum.Eh + E[l]
            }
            else {
               sum.Ej <- sum.Ej + E[l]
            }
         }
      }

      hl[i] <- sum.Eh / (sum.Eh + sum.Ej)

   }

   hl

}

