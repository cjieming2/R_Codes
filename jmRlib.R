###############################################################################
## one can get this library via:
## source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")
## source("/Users/jiemingchen/R_codes/jmRlib.R")

#######################################################################
## 

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## a function to round numbers in the conventional way
## round2(2.5)=3 not 2
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n  ## makes it as close to an integer to the number of places you are round to
  z = z + 0.5      ## add half of that "integer"
  z = trunc(z)     ## then you truncate
  z = z/10^n       ## convert back your "number"
  z*posneg
}

#######################################################################
## A function to add transparent color, where alpha=1 is opaque
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#######################################################################
## A function to find the second highest value
max2 <- function(x) 
{
  max( x[x!=max(x)] )
}

#######################################################################
## A function to use identify to select points, and overplot the
## points with another symbol as they are selected
identifyPch <- function(x, y=NULL, n=length(x), pch=19, tag=data$sample.id, sizeofpoint = 3,...)
{
  xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, length(x)); res <- integer(0)
  while(sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel],labels = as.character(tag), col="red", ...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch, cex = sizeofpoint)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  res
}

#######################################################################
## process a nxn symmetric matrix of sequence identities
## and returns the lower triangular as a vector n*(n+1) / 2 rows * 1
processMatrix <- function(data.mod)
{
  ## remove diagonals by first setting them as NA then remove NA
  ## also transforms them into a vector
  diag(data.mod) = NA
  data.ltri=vech(data.mod)
  d <- data.ltri[!is.na(data.ltri)]
  
  return(d)
}

###############################################################################
## this function calculates the 
## min, max, quantiles, variance, sd, and 95% CI
###############################################################################
summary.stats = function(data,filename) {
  cat("<SUMMARY STATISTICS> for", filename, "\n")
  cat("Minimum value: ", min(data), "\n")
  cat("Maximum value: ", max(data), "\n")
  cat("\n---QUANTILES---\n")
  cat("10%\t25%\t50%\t75%\n")
  cat(quantile(data, probs=c(10,25,50,75)/100, type=8, names = TRUE), "\n")
  cat("Median: ", median(data), "\n")
  cat("\n---MEAN & 95% CONFIDENCE INTERVAL OF MEAN---\n")
  cat("Mean: ", mean(data), "\n")
  cat("Variance: ", var(data), "\n")
  cat("Standard deviation: ", sd(data), "\n")
  
  ## this calc the normal CI
  error = qnorm(0.975)*sd(data)/sqrt(length(data))
  cat("95% Confidence Interval: left:", mean(data)-error, "; right:", mean(data)-error, "\n")
}

###############################################################################
### this is a modified function of vioplot from the package vioplot
### the only difference is las=2, so the x labels are vertical
### updated: cexMed to control size of median point
###############################################################################
vioplot.las2 <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                          horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                          lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, cexMed = 1,
                          at, add = FALSE, wex = 1, drawRect = TRUE) 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label, las=2)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed, cex = cexMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed, cex = cexMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}