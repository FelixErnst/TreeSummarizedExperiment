% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/changeTree.R
\name{changeTree}
\alias{changeTree}
\title{Change the row or column tree}
\usage{
changeTree(
  x,
  rowTree = NULL,
  rowNodeLab = NULL,
  colTree = NULL,
  colNodeLab = NULL,
  whichRowTree = 1,
  whichColTree = 1
)
}
\arguments{
\item{x}{A TreeSummarizedExperiment object}

\item{rowTree}{A phylo object. A new row tree.}

\item{rowNodeLab}{A character string. It provides the labels of nodes that
the rows of assays tables corresponding to. If NULL (default), the row
names of the assays tables are used.}

\item{colTree}{A phylo object. A new column tree.}

\item{colNodeLab}{A character string. It provides the labels of nodes that
the columns of assays tables corresponding to. If NULL (default), the
column names of the assays tables are used.}

\item{whichRowTree}{Which row tree to be replaced? Default is 1 (the first
tree in the \code{rowTree} slot).}

\item{whichColTree}{Which column tree to be replaced? Default is 1 (the first
tree in the \code{colTree} slot).}
}
\value{
A TreeSummarizedExperiment object
}
\description{
\code{changeTree} changes a row or column tree in a
\code{TreeSummarizedExperiment} object.
}
\examples{

library(ape)
set.seed(1)
treeR <- ape::rtree(10)

# the count table
count <- matrix(rpois(160, 50), nrow = 20)
rownames(count) <- paste0("entity", 1:20)
colnames(count) <- paste("sample", 1:8, sep = "_")
# The sample information
sampC <- data.frame(condition = rep(c("control", "trt"), 
                                    each = 4),
                    gender = sample(x = 1:2, size = 8, 
                                    replace = TRUE))
rownames(sampC) <- colnames(count)
# build a TreeSummarizedExperiment object
tse <- TreeSummarizedExperiment(assays = list(count),
                                colData = sampC,
                                rowTree = treeR,
                                rowNodeLab = rep(treeR$tip.label, each =2))

treeR2 <- drop.tip(phy = treeR, tip = c("t10", "t9", "t8"))

# if rownames are not used in node labels of the tree, provide rowNodeLab
use <- changeTree(x = tse, rowTree = treeR2, 
                  rowNodeLab = rep(treeR$tip.label, each =2))
use

# if rownames are used in node labels of tree, rowNodeLab is not required.

rownames(tse) <- rep(treeR$tip.label, each =2)
cse <- changeTree(x = tse, rowTree = treeR2)
cse
}
\author{
Ruizhu Huang
}
