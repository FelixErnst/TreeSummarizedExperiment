
# ============================================================================
### The documentation of accessor functions
### The corresponding codes are in the file allGenerics.R
# ============================================================================

# -----------------------------------------------------------------------------
### TreeSummarizedExperiment
# -----------------------------------------------------------------------------
#' TreeSummarizedExperiment-accessors
#'
#' All accessor functions that work on
#' \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}} should work on
#' \strong{TreeSummarizedExperiment}. Additionally, new accessors
#' \code{rowLinks} \code{colLinks}, \code{rowTree} and \code{colTree} accessor
#' function are available for \strong{TreeSummarizedExperiment}.
#'
#' @param x A TreeSummarizedExperiment object
#' @param i,j The row, column index to subset \code{x}. The arguments of the
#'   subset function \code{[]}
#' @param drop A logical value, TRUE or FALSE. The argument from the subset
#'   function \code{[]}
#' @param value the new rownames or colnames as a \code{character} value.
#'   See \code{\link[BiocGenerics:row_colnames]{BiocGenerics}}.
#' @param ... The argument from the subset function \code{[]}
#' @param rowNode A vector of nodes that are used to subset rows. One could use
#'   the node number, the node label or the node alias to specify nodes, but not
#'   a mixture of them.
#' @param colNode A vector of nodes that are used to subset columns. One could
#'   use the node number, the node label or the node alias to specify nodes, but
#'   not a mixture of them.
#' @name TreeSummarizedExperiment-accessor
#' @return Elements from \code{TreeSummarizedExperiment}.
#' @seealso \code{\link{TreeSummarizedExperiment}}
#'   \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}}
#'
#' @author Ruizhu HUANG
#' @examples
#'
#' # the assay table
#' set.seed(1)
#' y <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(y) <- paste(rep(LETTERS[1:3], each = 10), rep(1:10,3), sep = "_")
#' rownames(y) <- tinyTree$tip.label
#'
#' # the row data
#' rowInf <- DataFrame(var1 = sample(letters[1:3], 10, replace = TRUE),
#'                     var2 = sample(c(TRUE, FALSE), 10, replace = TRUE))
#' # the column data
#' colInf <- DataFrame(gg = factor(sample(1:3, 30, replace = TRUE)),
#'                     group = rep(LETTERS[1:3], each = 10))
#'
#' # the tree structure on the rows of assay tables
#' data("tinyTree")
#'
#' # the tree structure on the columns of assay tables
#' sampTree <- ape::rtree(30)
#' sampTree$tip.label <- colnames(y)
#'
#' # create the TreeSummarizedExperiment object
#' toy_tse <- TreeSummarizedExperiment(assays = list(y),
#'                                     rowData = rowInf,
#'                                     colData = colInf,
#'                                     rowTree = tinyTree,
#'                                     colTree = sampTree)
#'
#' ## extract the rowData
#' (rowD <- rowData(x = toy_tse))
#'
#' ## extract the colData
#' (colD <- colData(x = toy_tse))
#'
#' ## extract the linkData
#' # on rows
#' (rowL <- rowLinks(x = toy_tse))
#' # on columns
#' (colL <- colLinks(x = toy_tse))
#'
#'  ## extract the treeData
#' # on rows
#' (rowT <- rowTree(x = toy_tse))
#' # on columns
#' (colT <- colTree(x = toy_tse))
#'
NULL

# -----------------------------------------------------------------------------
### Accessors for TreeSummarizedExperiment
# -----------------------------------------------------------------------------

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("rowLinks", function(x)
    standardGeneric("rowLinks")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowLinks", signature("TreeSummarizedExperiment"),
          function(x) {
              x@rowLinks
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colLinks", function(x)
    standardGeneric("colLinks")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colLinks", signature("TreeSummarizedExperiment"),
          function(x) {
              x@colLinks
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("rowTree", function(x)
    standardGeneric("rowTree")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowTree", signature("TreeSummarizedExperiment"),
          function(x) {
              x@rowTree$phylo
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colTree", function(x)
    standardGeneric("colTree")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colTree", signature("TreeSummarizedExperiment"),
          function(x) {
              x@colTree$phylo
          })


#' @importFrom methods callNextMethod
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @rdname TreeSummarizedExperiment-accessor
#' @export
#'
setMethod("[", signature(x = "TreeSummarizedExperiment"),
          function(x, i, j, ..., drop = TRUE){

              # Subset the rowLinks
              lr <- rowLinks(x)
              rt <- rowTree(x)
              if (!missing(i) & !is.null(rt)) {
                  # match with rownames
                  # multiple rows in assays might have the same name
                  if (is.character(i)) {
                      isRn <- all(i %in% rownames(x))
                      if (isRn) {
                          i <- which(rownames(x) %in% i)
                      } else {
                          stop(i, " can't be found in rownames")
                      }
                  }

                  nlr <- lr[i, , drop = FALSE]
              } else {
                  nlr <- lr
              }

              # Subset the colLinks
              lc <- colLinks(x)
              ct <- colTree(x)
              if (!missing(j) & !is.null(ct)) {
                  # match with colnames
                  # multiple columns in assays might have the same name
                  if (is.character(j)) {
                      isCn <- all(j %in% colnames(x))
                      if (isCn) {
                          j <- which(colnames(x) %in% j)
                      } else {
                          stop(j, " can't be found in colnames")
                      }
                  }
                  nlc <- lc[j, , drop = FALSE]
              } else {
                  nlc <- lc
              }


              # Subset the traditional slots from SummarizedExperiment
              nx <- callNextMethod()

              # update slots
              final <- BiocGenerics:::replaceSlots(nx,
                                                   rowLinks = nlr,
                                                   colLinks = nlc)

              return(final)
          })

#' @importFrom methods callNextMethod
#' @rdname TreeSummarizedExperiment-accessor
#' @export
setReplaceMethod("rownames", signature(x = "TreeSummarizedExperiment"),
                 function(x, value){
                     x <- callNextMethod()
                     if(!is.null(x@rowLinks)){
                         rownames(x@rowLinks) <- value
                     }
                     x
                 }
)

#' @importFrom methods callNextMethod
#' @rdname TreeSummarizedExperiment-accessor
#' @export
setReplaceMethod("colnames", signature(x = "TreeSummarizedExperiment"),
                 function(x, value){
                     x <- callNextMethod()
                     if(!is.null(x@colLinks)){
                         rownames(x@colLinks) <- value
                     }
                     x
                 }
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("subsetByNode", function(x, rowNode, colNode)
    standardGeneric("subsetByNode")
)

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("subsetByNode", signature(x = "TreeSummarizedExperiment"),
          function(x, rowNode, colNode){
              # row link
              rl <- rowLinks(x)
              if (!missing(rowNode)) {
                  if (!is.numeric(rowNode)) {
                      rowNode <- convertNode(tree = rowTree(x), node = rowNode)
                  }
                  x <- x[which(rl$nodeNum %in% rowNode),]
              }

              # column link
              cl <- colLinks(x)
              if (!missing(colNode)) {
                  if (!is.numeric(colNode)) {
                      colNode <- convertNode(tree = colTree(x), node = colNode)
                  }
                  x <- x[, which(cl$nodeNum %in% colNode)]
              }
              return(x)
          }
)

#' @keywords internal
#' @importFrom methods callNextMethod
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "TreeSummarizedExperiment", function(object) {
    callNextMethod()

    rt <- rowTree(object)
    ct <- colTree(object)



    rlk <- rowLinks(object)
    clk <- colLinks(object)

    # on row
    if (is.null(rt)) {
        msg1a <- "rowLinks: NULL"
        msg1b <- "rowTree: NULL"
    } else {
        msg1a <- sprintf("rowLinks: a LinkDataFrame (%d %s)",
                         nrow(rlk), "rows")


        # the number of leaf nodes & internal nodes
        nlr <- countLeaf(rt)
        nnr <- countNode(rt) - countLeaf(rt)
        msg1b <- sprintf("rowTree: a %s (%d leaves)", class(rt), nlr)
    }

    # on column
    if (is.null(ct)) {
        msg2a <- "colLinks: NULL"
        msg2b <- "colTree: NULL"
    } else {
        msg2a <- sprintf("colLinks: a LinkDataFrame (%d %s)", nrow(clk), "rows")

        # the number of leaf nodes & internal nodes
        nlc <- countLeaf(ct)
        nnc <- countNode(ct) - countLeaf(ct)
        msg2b <- sprintf("colTree: a %s (%d leaves)", class(ct), nlc)
    }

    cat(msg1a, "\n", msg1b, "\n",
        msg2a, "\n", msg2b, "\n",
        sep = "")
})








