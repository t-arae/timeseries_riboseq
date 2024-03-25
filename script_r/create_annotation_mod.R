# Original code are from riboWaltz@v1.1.0 and slightly modified.
# Instead of specifying the package name in the `txdb` argument,
# we changed it so that we can pass our own txdb object directly to the `txdbanno` argument.
create_annotation_mod <- function(txdbanno = NULL) {
  exon <-
    suppressWarnings(GenomicFeatures::exonsBy(txdbanno, by = "tx", use.names = T))
  utr5 <-
    suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdbanno, use.names = T))
  cds <-
    suppressWarnings(GenomicFeatures::cdsBy(txdbanno, by = "tx", use.names = T))
  utr3 <-
    suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdbanno, use.names = T))
  exon <- data.table::as.data.table(exon[unique(names(exon))])
  utr5 <- data.table::as.data.table(utr5[unique(names(utr5))])
  cds <- data.table::as.data.table(cds[unique(names(cds))])
  utr3 <- data.table::as.data.table(utr3[unique(names(utr3))])
  anno_df <- exon[, list(l_tr = sum(width)), by = list(transcript = group_name)]
  l_utr5 <- utr5[, list(l_utr5 = sum(width)), by = list(transcript = group_name)]
  l_cds <- cds[, list(l_cds = sum(width)), by = list(transcript = group_name)]
  l_utr3 <- utr3[, list(l_utr3 = sum(width)), by = list(transcript = group_name)]
  merge_allx <- function(x, y) merge(x, y, all.x = TRUE)
  anno_df <- BiocGenerics::Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
  anno_df[is.na(anno_df)] <- 0
  return(anno_df)
}
