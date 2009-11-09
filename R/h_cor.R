`h_cor` <-
function(input, na.string, n, method) {

   # Read data:
   g <- read.table(input, na.strings = na.string)

   # Check data:
   chkdata(g[, 2:ncol(g)])

   # Correlation sampling:
   r <- hh(g[, 2:ncol(g)], n, method)
   
   r

}

