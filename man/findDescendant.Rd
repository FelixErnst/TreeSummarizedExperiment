% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecate_Fun.R, R/tree_findDescendant.R
\name{findOS}
\alias{findOS}
\alias{findDescendant}
\title{Find descendants (or offsprings)}
\usage{
findOS(tree, node, only.leaf = TRUE, self.include = FALSE, use.alias = FALSE)

findDescendant(
  tree,
  node,
  only.leaf = TRUE,
  self.include = FALSE,
  use.alias = FALSE
)
}
\arguments{
\item{tree}{A phylo object.}

\item{node}{An internal node. It could be the node number or the node
label.}

\item{only.leaf}{A logical value, TRUE or FALSE. The default is TRUE. If
default, only the leaf nodes in the descendant nodes would be returned.}

\item{self.include}{A logical value, TRUE or FALSE. The default is FALSE. If
TRUE, the node specified in \strong{node} is included and the leaf node
itself is returned as its descendant.}

\item{use.alias}{A logical value, TRUE or FALSE. The default is FALSE, and
the node label would be used to name the output; otherwise, the alias of
node label would be used to name the output. The alias of node label is
created by adding a prefix \code{"alias_"} to the node number.}
}
\value{
A vector of nodes. The numeric value is the node number, and the
  vector name is the corresponding node label. If a node has no label, it
  would have NA as name when \code{use.alias = FALSE}, and have the alias of
  node label as name when \code{use.alias = TRUE}.
}
\description{
\code{findDescendant} finds descendants of a node.
}
\examples{
data(tinyTree)

library(ggtree)
ggtree(tinyTree) +
geom_text2(aes(label = node), color = "darkblue",
               hjust = -0.5, vjust = 0.7) +
geom_hilight(node = 17, fill = 'steelblue', alpha = 0.5) +
geom_text2(aes(label = label), color = "darkorange",
           hjust = -0.1, vjust = -0.7)

(tips <- findDescendant(tree = tinyTree, node = c(17), only.leaf = TRUE))
}
\author{
Ruizhu Huang
}
