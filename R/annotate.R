#' Annotate a gene x sample data.frame with biomaRt gene names
#'
#' Looks up each row name in Ensembl via biomaRt and appends
#' \code{external_gene_name}, \code{description}, and \code{gene_biotype}.
#'
#' @param data Data.frame whose row names are Ensembl gene IDs.
#' @param dataset biomaRt dataset name; default
#'   \code{"hsapiens_gene_ensembl"}.
#' @return The input data.frame with three new columns; rows that don't
#'   resolve are dropped.
#' @export
add_geneName <- function(data, dataset = "hsapiens_gene_ensembl") {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required. Install with: ",
         "BiocManager::install('biomaRt').", call. = FALSE)
  }
  ensembl  <- biomaRt::useMart("ensembl")
  ensemblh <- biomaRt::useDataset(dataset, mart = ensembl)
  ensinfos <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name",
                   "description", "gene_biotype"),
    filters = "ensembl_gene_id",
    values  = rownames(data),
    mart    = ensemblh
  )
  rownames(ensinfos) <- ensinfos$ensembl_gene_id
  ensinfos <- ensinfos[, -1, drop = FALSE]
  out <- merge(data, ensinfos, by = "row.names")
  rownames(out) <- out[, 1]
  out[, -1, drop = FALSE]
}
