#causalTree.anova <- function(y, offset, parms, wt)
causalTree.anova <- function(y, offset, parms, wt)
{
    if (!is.null(offset)) y <- y - offset
    #list(y = y, parms = NULL, numresp = 1L, numy = 1L,
      list(y = y, parms = NULL, numresp = 1L, numy = 1L,
	 summary = function(yval, dev, wt, ylevel, digits) {
	   # paste0("  mean=", formatg(yval, digits),
	   #                  ", MSE=" , formatg(dev, digits))
	   paste0("  causal effect=", formatg(yval, digits),
	          ", error=" , formatg(dev, digits))
     
	     #paste0("  mean=", formatg(yval, digits),
        #            ", MSE=" , formatg(dev/wt, digits))
         },
	 text = function(yval, dev, wt, ylevel, digits, n, use.n ) {
	     if (use.n) paste0(formatg(yval, digits), "\nn=", n) else
             formatg(yval, digits)
         })
}
