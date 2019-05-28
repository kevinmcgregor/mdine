## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE---------------------------------------------------------
#  if (!require(devtools)) {
#    install.packages("devtools")
#    library(devtools)
#  }
#  
#  install_github("kevinmcgregor/mdine", dependencies=TRUE)

## ------------------------------------------------------------------------
library(mdine)
data(crohns)

# Covariate data
head(crohns$covars)
# OTU table
head(crohns$otu.counts)

## ------------------------------------------------------------------------
X <- model.matrix(~disease, data=crohns$covars)
head(X)

## ------------------------------------------------------------------------
# Estimated precision matrix for control samples (Z=0):
md.fit$post_mean$invsigma0

# Estimated precision matrix for Crohn's samples (Z=1):
md.fit$post_mean$invsigma1

## ------------------------------------------------------------------------
# Weighted adjacency matrices based on each precision matrix
adj <- ci2adj(md.fit, weighted = TRUE)
adj

## ---- fig.height=4, fig.width=6, fig.align='center'----------------------
# Plotting the two networks
plot_networks(md.fit)

## ---- fig.height=4, fig.width=4, fig.align='center'----------------------
# Weighted adjacency matrices based on each precision matrix
ig0 <- adj2ig(adj$adj0)
igraph::plot.igraph(ig0)

