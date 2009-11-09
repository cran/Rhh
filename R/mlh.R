`mlh` <-
function(input, output, na.string, n.digits) {

   # Read data:
   g <- read.table(input, na.strings = na.string)

   # Check data:
   chkdata(g[, 2:ncol(g)])

   # Reserve array for results:
   h <- as.data.frame(array(dim = c(nrow(g), 4)), stringsAsFactors = FALSE)

   # Calculate heterozygosity estimates:
   h[, 1] <- g[, 1]
   h[, 2] <- ir(g[, 2:ncol(g)])
   h[, 3] <- sh(g[, 2:ncol(g)])
   h[, 4] <- hl(g[, 2:ncol(g)])
   colnames(h) <- c("ID", "IR", "SH", "HL")

   # Write results to the output file:
   write.table(format(h, digits=n.digits), file = output, col.names = c("ID", "IR", "SH", "HL"), row.names = FALSE, quote = FALSE, sep="\t")

   h

}

