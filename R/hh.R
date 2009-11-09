`hh` <-
function(genotypes, n, method) {

   if (method == "sh" || method == "ir" || method == "hl") {

      individuals <- nrow(genotypes)
      loci <- ncol(genotypes) / 2
      h1 <- array(0, dim=c(individuals, n))
      h2 <- array(0, dim=c(individuals, n))
      corr <- 1:n
      
      cat("Starting...", "\n")
      flush.console()

      for (i in 1:n) {

         if (i %% 10 == 1) {
            cat("Iteration:", i, "\n")
            flush.console()
         }

         samp <- sample(loci)

         for (l in 1:(floor(loci / 2))) {
            g <- 2 * samp[l] - 1
            h <- 2 * samp[l]
            if (l == 1) {
               gen1 <- genotypes[, g:h]
            }
            else {
               gen1 <- cbind(gen1, genotypes[, g:h])
            }
         }

         for (l in (floor(loci / 2) + 1):loci) {
            g <- 2 * samp[l] - 1
            h <- 2 * samp[l]
            if (l == floor(loci / 2) + 1) {
               gen2 <- genotypes[, g:h]
            }
            else {
               gen2 <- cbind(gen2, genotypes[, g:h])
            }
         }

         if (method == "sh") {
            h1[, i] <- sh(gen1)
            h2[, i] <- sh(gen2)
         }
         if (method == "ir") {
            h1[, i] <- ir(gen1)
            h2[, i] <- ir(gen2)
         }
         if (method == "hl") {
            h1[, i] <- hl(gen1)
            h2[, i] <- hl(gen2)
         }

         corr[i] <- cor(h1[, i], h2[, i], use="complete.obs")
      }

      cat("Finished.", "\n")
      flush.console()
      
   }

   else {
      stop("Invalid method.")
   }

   corr

}

