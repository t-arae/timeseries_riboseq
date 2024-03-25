
library(magrittr)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)

#' Extract GRanges object of gene(s) of interest
#' @param txdb TxDb-object
#' @param goi character vector of gene(s) of interest or NULL
#' @returns A GRanges
#' @export
#' @examples 
#' AGI_GOI <- "AT1G01010"
#' txdb_arabi <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
#' extract_gr_gene(txdb_arabi)
#' extract_gr_gene(txdb_arabi, AGI_GOI)
#' 
extract_gr_gene <- function(txdb, goi = NULL) {
  if(is.null(goi)) {
    FILTER <- NULL
  } else {
    FILTER <- list("gene_id" = goi)
  }
  GenomicFeatures::genes(txdb, filter = FILTER)
}

#' Extract GRanges object of transcript(s) of gene(s) interest
#' @param txdb TxDb-object
#' @param goi character vector of gene(s) of interest or NULL
#' @returns A GRanges
#' @export
#' @examples 
#' AGI_GOI <- "AT1G01010"
#' txdb_arabi <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
#' extract_gr_transcript(txdb_arabi)
#' extract_gr_transcript(txdb_arabi, AGI_GOI)
#' 
extract_gr_transcript <- function(txdb, goi = NULL) {
  if(is.null(goi)) {
    FILTER <- NULL
  } else {
    FILTER <- list("gene_id" = goi)
  }
  GenomicFeatures::transcripts(
    txdb,
    columns = c("gene_id", "tx_id", "tx_name"),
    filter = FILTER
  )
}

#' Extract GRanges object of exon(s) of gene(s) interest
#' @param txdb TxDb-object
#' @param goi character vector of gene(s) of interest or NULL
#' @returns A GRanges
#' @export
#' @examples 
#' AGI_GOI <- "AT1G01010"
#' txdb_arabi <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
#' extract_gr_exon(txdb_arabi)
#' extract_gr_exon(txdb_arabi, AGI_GOI)
#' 
extract_gr_exon <- function(txdb, goi = NULL) {
  if(is.null(goi)) {
    FILTER <- NULL
  } else {
    FILTER <- list("gene_id" = goi)
  }
  GenomicFeatures::exons(
    txdb,
    columns = c("gene_id", "tx_id", "tx_name", "exon_id", "exon_name", "exon_rank"),
    filter = FILTER
  )
}

#' Extract GRanges object of CDS(s) of gene(s) interest
#' @param txdb TxDb-object
#' @param goi character vector of gene(s) of interest or NULL
#' @returns A GRanges
#' @export
#' @examples 
#' AGI_GOI <- "AT1G01010"
#' txdb_arabi <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
#' extract_gr_cds(txdb_arabi)
#' extract_gr_cds(txdb_arabi, AGI_GOI)
#' 
extract_gr_cds <- function(txdb, goi = NULL) {
  if(is.null(goi)) {
    FILTER <- NULL
  } else {
    FILTER <- list("gene_id" = goi)
  }
  GenomicFeatures::cds(
    txdb,
    columns = c("gene_id", "tx_id", "tx_name", "exon_id", "exon_name", "exon_rank",
                "cds_id", "cds_name", "cds_phase"),
    filter = FILTER
  )
}

#' Extract GRanges object of 5'UTR(s) of gene(s) interest
#' @param txdb TxDb-object
#' @param goi character vector of gene(s) of interest or NULL
#' @returns A GRanges
#' @export
#' @examples 
#' AGI_GOI <- "AT1G01010"
#' txdb_arabi <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
#' extract_gr_futr(txdb_arabi)
#' extract_gr_futr(txdb_arabi, AGI_GOI)
#' 
extract_gr_futr <- function(txdb, goi = NULL) {
  gr_futr <-
    GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE) %>%
    unlist() %>%
    {plyranges::mutate(., tx_name = names(.))}
  if(is.null(goi)) {
    return(gr_futr)
  } else {
    plyranges::filter(gr_futr, grepl(paste(goi, collapse = "|"), tx_name))
  }
}

#' Extract GRanges object of 3'UTR(s) of gene(s) interest
#' @param txdb TxDb-object
#' @param goi character vector of gene(s) of interest or NULL
#' @returns A GRanges
#' @export
#' @examples 
#' AGI_GOI <- "AT1G01010"
#' txdb_arabi <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
#' extract_gr_tutr(txdb_arabi)
#' extract_gr_tutr(txdb_arabi, AGI_GOI)
#' 
extract_gr_tutr <- function(txdb, goi = NULL) {
  gr_tutr <-
    GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE) %>%
    unlist() %>%
    {plyranges::mutate(., tx_name = names(.))}
  if(is.null(goi)) {
    return(gr_tutr)
  } else {
    plyranges::filter(gr_tutr, grepl(paste(goi, collapse = "|"), tx_name))
  }
}

#' Return an extended GRanges of the specified gene region
#' @param txdb TxDb-object
#' @param AGI character vector of gene of interest
#' @returns A GRanges
#' @export
#' @examples 
#' AGI_GOI <- "AT1G01010"
#' txdb_arabi <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
#' make_gr_goi_all(txdb_arabi, AGI_GOI)
#' 
make_gr_goi_all <- function(txdb, AGI) {
  if(length(AGI) != 1L) stop("the length of AGI must be 1.")
  
  gr_gene_goi <- extract_gr_gene(txdb, AGI)
  gr_goi_all <- gr_gene_goi %>% plyranges::select(!gene_id)
  strand(gr_goi_all) <- "*"
  gr_goi_all %>% {plyranges::stretch(., extend = width(.) * 0.02)}
}

#' Extract GRanges objects overlapped specified region and return as a list
#' @param txdb TxDb-object
#' @param goi character vector of gene(s) of interest or NULL
#' @returns a list of GRanges
#' @export
#' @examples 
#' txdb_arabi <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
#' gr_all <- plyranges::as_granges(data.frame(seqnames = "1", start = 1, end = 10000))
#' extract_li_gr(txdb_arabi, gr_all)
#' 
extract_li_gr <- function(txdb, gr) {
  purrr::map(
    list(
      "transcript" = extract_gr_transcript,
      "exon" = extract_gr_exon,
      "cds" = extract_gr_cds,
      "futr" = extract_gr_futr,
      "tutr" = extract_gr_tutr
    ),
    ~ .x(txdb)
  ) %>%
    purrr::map(plyranges::filter_by_overlaps, y = gr)
}

#' Transform a list of GRanges to a list of tibble
#' @param li_gr an output from `etract_li_gr()`
#' @export
#' 
li_gr2li_tbl <- function(li_gr) {
  tbl_plot_start <-
    tibble::as_tibble(li_gr$transcript) %>%
    dplyr::mutate(tss = ifelse(strand == "+", start, end))
  tbl_plot_tail <-
    tibble::as_tibble(li_gr$transcript) %>%
    dplyr::mutate(tts = ifelse(strand == "+", end, start),
                  start = ifelse(strand == "+", tts-1, tts+1))
  li <- purrr::map(li_gr, tibble::as_tibble)
  li[["start"]] <- tbl_plot_start
  li[["tail"]] <- tbl_plot_tail
  purrr::map(li, tidyr::unnest, cols = tidyselect::where(is.list))
}



#' Extract a transcript data and add transcript coordinates
#' @param li_tbl_plot an output from `li_gr2li_tbl()`
#' @param tx_name_goi a character vector to specify the transcript
#' @export
#' 
make_li_tbl_plot_tx <- function(li_tbl_plot, tx_name_goi) {
  li_tbl_plot_tx <-
    li_tbl_plot %>%
    purrr::map(dplyr::filter, tx_name == tx_name_goi) %>%
    purrr::modify_if(
      .p = ~ any(names(.x) == "exon_rank"),
      .f = ~ dplyr::arrange(.x, exon_rank)
    ) %>%
    purrr::map(dplyr::mutate, tx_end = cumsum(width)) %>%
    purrr::map(dplyr::mutate, tx_start = tx_end - width)
  
  len_futr <- sum(li_tbl_plot_tx$futr$width)
  li_tbl_plot_tx$cds$tx_start <- li_tbl_plot_tx$cds$tx_start + len_futr
  li_tbl_plot_tx$cds$tx_end <- li_tbl_plot_tx$cds$tx_end + len_futr
  len_futr_cds <- max(li_tbl_plot_tx$cds$tx_end)
  li_tbl_plot_tx$tutr$tx_start <- li_tbl_plot_tx$tutr$tx_start + len_futr_cds
  li_tbl_plot_tx$tutr$tx_end <- li_tbl_plot_tx$tutr$tx_end + len_futr_cds
  
  li_tbl_plot_tx %>%
    purrr::map(dplyr::select, tx_start, tx_end, dplyr::everything())
}

#' Show the transcript feature length information
#' @param li_tbl_plot an output from `li_gr2li_tbl()`
#' @export
#' 
show_transcript_info <- function(li_tbl_plot) {
  tx_names <- li_tbl_plot$transcript$tx_name
  f <- function(tbl, name) sum(dplyr::filter(tbl, tx_name %in% name)$width)
  len_exons <- purrr::map_int(tx_names, ~ f(li_tbl_plot$exon, .x))
  len_futrs <- purrr::map_int(tx_names, ~ f(li_tbl_plot$futr, .x))
  len_cdss <- purrr::map_int(tx_names, ~ f(li_tbl_plot$cds, .x))
  len_tutrs <- purrr::map_int(tx_names, ~ f(li_tbl_plot$tutr, .x))
  
  purrr::pmap_dfr(
    .l = list(tx_names, len_exons, len_futrs, len_cdss, len_tutrs),
    .f = ~ {
      tibble::tribble(
        ~ tx_name, ~ total, ~ "5'UTR", ~ CDS, ~ "3'UTR",
        ..1, ..2, ..3, ..4, ..5)
    }
  ) %>%
    dplyr::arrange(tx_name) %>%
    knitr::kable()
}

#' Plot gene structure with sequence
#' @param li_tbl_plot_tx an output from `make_li_tbl_plot_tx()`
#' @param zoom_x_min a numeric, default: NA
#' @param zoom_x_max a numeric, default: NA
#' @export
#' 
plot_structure <- function(li_tbl_plot_tx, zoom_x_min = NA, zoom_x_max = NA) {
  TX_X_MIN <- li_tbl_plot_tx$transcript$start
  TX_X_MAX <- li_tbl_plot_tx$transcript$end
  shimashima <-
    annotate("rect",
             xmin = seq(TX_X_MIN, TX_X_MAX, by = 2) -.5,
             xmax = seq(TX_X_MIN, TX_X_MAX, by = 2) +.5,
             ymin = .5, ymax = 2.5,
             fill = "grey", alpha = .5)
  gp_check_transcript <-
    ggplot(mapping = aes(xmin = start, xmax = end, y = tx_name)) +
    annotate("rect", xmin = zoom_x_min, xmax = zoom_x_max,
             ymin = -Inf, ymax = Inf, fill = "grey") +
    geom_linerange(data = li_tbl_plot_tx$transcript, aes(group = tx_id)) +
    geom_linerange(data = li_tbl_plot_tx$cds, aes(group = exon_id), linewidth = 2) +
    geom_linerange(data = li_tbl_plot_tx$futr, aes(group = exon_id), linewidth = 1) +
    geom_linerange(data = li_tbl_plot_tx$tutr, aes(group = exon_id), linewidth = 1) +
    geom_point(data = li_tbl_plot_tx$start, aes(x = tss, group = tx_id)) +
    geom_segment(data = li_tbl_plot_tx$tail,
                 aes(x = start, xend = tts, yend = tx_name, group = tx_id),
                 arrow = arrow(type = "closed", length = unit(5, "points")))
    gp_check_transcript_mini <-
    ggplot(mapping = aes(xmin = start, xmax = end, y = 1)) +
    shimashima +
    geom_linerange(data = li_tbl_plot_tx$transcript, aes(group = tx_id)) +
    geom_linerange(data = li_tbl_plot_tx$cds, aes(group = exon_id), linewidth = 2, lineend = "butt") +
    geom_linerange(data = li_tbl_plot_tx$futr, aes(group = exon_id), linewidth = 1, lineend = "butt") +
    geom_linerange(data = li_tbl_plot_tx$tutr, aes(group = exon_id), linewidth = 1, lineend = "butt")
  gp_check_seq <-
    bsg_tair[[as.character(li_tbl_plot_tx$transcript$seqnames[1])]] %>%
    Biostrings::subseq(start = TX_X_MIN, end = TX_X_MAX) %>%
    as.character() %>%
    stringr::str_extract_all(".") %>%
    {tibble::tibble(coord = TX_X_MIN:TX_X_MAX, nucleotide = .[[1]])} %>%
    ggplot(aes(coord, 1, label = nucleotide, color = after_stat(label))) +
    shimashima +
    geom_text() +
    geom_text(aes(y = 2, label = bstringr::dstr_complement(nucleotide)))
  
  patchwork::wrap_plots(
    gp_check_transcript,
    gp_check_transcript_mini +
      coord_cartesian(xlim = c(zoom_x_min, zoom_x_max)),
    gp_check_seq +
      coord_cartesian(xlim = c(zoom_x_min, zoom_x_max), ylim = c(.5, 2.5)),
    ncol = 1, heights = c(3, 3, 3)) &
    theme(panel.background = element_blank())
}

plot_struct <- function(tbl_plot, linewidth = 1, color = "#191970", fill = "#191970") {
  if(nrow(tbl_plot) == 0) return(geom_linerange())
  geom_rect(
    data = tbl_plot,
    aes(
      xmin = tx_start - 0.5, xmax = tx_end + 0.5,
      ymin = ymin, ymax = ymax
    ),
    linewidth = linewidth, color = color, fill = fill)
}

plot_uorf <- function(tbl_plot, linewidth = 1, color = "#191970", fill = "#191970") {
  if(nrow(tbl_plot) == 0) return(geom_linerange())
  geom_rect(
    data = tbl_plot,
    aes(
      xmin = tx_start - 0.5, xmax = tx_end + 0.5,
      ymin = ymin, ymax = ymax,
      alpha = is_active
    ),
    linewidth = linewidth, color = color, fill = fill)
}

plot_uorf2 <- function(tbl_plot, linewidth = 1, color = "#191970", fill = "#191970") {
  if(nrow(tbl_plot) == 0) return(geom_linerange())
  geom_rect(
    data = tbl_plot,
    aes(
      xmin = tx_start - 0.5, xmax = tx_end + 0.5,
      ymin = ymin, ymax = ymax,
      alpha = is_active,
      linewidth = is_overlap
    ),
    # linewidth = linewidth, color = color, fill = fill)
    color = color, fill = fill, linetype = "dashed")
}

#' Plot gene structure with sequence
#' @param li_tbl_plot_tx an output from `make_li_tbl_plot_tx()`
#' @param tbl_uorf an output from `make_tbl_uorf()`
#' @param zoom_x_min a numeric, default: NA
#' @param zoom_x_max a numeric, default: NA
#' @export
#' 
plot_structure_w_uorf <- function(li_tbl_plot_tx, tbl_uorf, zoom_x_min = NA, zoom_x_max = NA) {
  TX_X_MIN <- min(li_tbl_plot_tx$exon$tx_start) + 1
  TX_X_MAX <- max(li_tbl_plot_tx$exon$tx_end)
  zoom_x_min <- ifelse(is.na(zoom_x_min), TX_X_MIN, zoom_x_min)
  zoom_x_max <- ifelse(is.na(zoom_x_max), TX_X_MAX, zoom_x_max)
  shimashima <-
    annotate("rect",
             xmin = seq(TX_X_MIN, TX_X_MAX, by = 2) -.5,
             xmax = seq(TX_X_MIN, TX_X_MAX, by = 2) +.5,
             ymin = -Inf, ymax = Inf,
             fill = "grey", alpha = .5)
    
  tbl_uorf_plot <-
    tbl_uorf %>%
    dplyr::with_groups(c(ID, is_active), dplyr::summarise,
                       tx_start = min(tx_start)+1, tx_end = max(tx_end)) %>%
    dplyr::mutate(frame = tx_start %% 3)
  
  summarise <- function(tbl) dplyr::summarise(tbl, tx_start = min(tx_start) + 1, tx_end = max(tx_end), .by = tx_name)
  mutate_ymin_ymax <- function(tbl, y, offset) dplyr::mutate(summarise(tbl), ymin = y - offset, ymax = y + offset)
  mutate_y_uorf <- function(tbl, offset) dplyr::mutate(tbl, ymin = frame - offset + 3, ymax = frame + offset + 3)
  
  gp_check_transcript <-
    ggplot() +
    annotate("rect", xmin = zoom_x_min, xmax = zoom_x_max,
             ymin = -Inf, ymax = Inf, fill = "grey") +
    plot_struct(mutate_ymin_ymax(li_tbl_plot_tx$cds, 1, 1), fill = NA) +
    plot_struct(mutate_ymin_ymax(li_tbl_plot_tx$futr, 1, .3), color = NA) +
    plot_struct(mutate_ymin_ymax(li_tbl_plot_tx$tutr, 1, .3), color = NA) +
    plot_uorf(mutate_y_uorf(tbl_uorf_plot, .45), color = NA) +
    scale_alpha_manual(values = c("TRUE" = .8, "FALSE" = .4))
  
  gp_check_transcript_mini <-
    ggplot() +
    shimashima +
    plot_struct(mutate_ymin_ymax(li_tbl_plot_tx$cds, 1, 1), fill = NA) +
    plot_struct(mutate_ymin_ymax(li_tbl_plot_tx$futr, 1, .3), color = NA) +
    plot_struct(mutate_ymin_ymax(li_tbl_plot_tx$tutr, 1, .3), color = NA) +
    plot_uorf(mutate_y_uorf(tbl_uorf_plot, .45), color = NA) +
    scale_alpha_manual(values = c("TRUE" = .8, "FALSE" = .4))
  
  seqname <- as.character(li_tbl_plot_tx$transcript$seqnames[1])
  exon_seq <-
    bsg_tair[[seqname]] %>%
    as.character() %>%
    bstringr::as_dstr(n = seqname) %>%
    bstringr::bstr_sub_all(from = li_tbl_plot_tx$exon$start,
                           to = li_tbl_plot_tx$exon$end) %>%
    .[[1]]
  
  if(li_tbl_plot_tx$transcript$strand == "-")
    exon_seq <- bstringr::dstr_rev_comp(exon_seq)
  
  gp_check_seq <-
    exon_seq %>% 
    paste(collapse = "") %>%
    stringr::str_extract_all(".") %>%
    {tibble::tibble(coord = TX_X_MIN:TX_X_MAX, nucleotide = .[[1]])} %>%
    ggplot(aes(coord, 1, label = nucleotide, color = after_stat(label))) +
    shimashima +
    geom_text() +
    theme(legend.position = "none")
  
  patchwork::wrap_plots(
    gp_check_transcript,
    gp_check_transcript_mini +
      coord_cartesian(xlim = c(zoom_x_min, zoom_x_max)),
    gp_check_seq +
      coord_cartesian(xlim = c(zoom_x_min, zoom_x_max), ylim = c(.5, 1.5)),
    ncol = 1, heights = c(3, 3, 3)) &
    theme(panel.background = element_blank())
}


theme_gene_structure <- function(gr_all) {
  list(
    labs(x = as.character(seqnames(gr_all)), y = "transcripts"),
    theme_minimal(),
    theme(
      line = element_blank(),
      axis.ticks.x = element_line(),
      plot.margin = unit(c(0,0,0,0), "pt")
    ),
    scale_x_continuous(
      limits = c(start(gr_all), end(gr_all)),
      expand = expansion(mult = c(0, 0)),
      breaks = scales::breaks_pretty(n = 6),
      labels = scales::label_number(
        big.mark = ",",
        suffix = "bp",
        scale_cut = c(0, " k" = 1000)
      )
    )
  )
}


#' @param bigwig_files path to the bigwig files
#' @param gr_region
read_bigwig_region_directed <- function(bigwig_files, gr_region, gr_feature) {
  strand <- ifelse(as.character(gr_feature@strand) == "+", "plus", "minus")
  bw_selection <-
    gr_goi_all %>%
    rtracklayer::BigWigSelection()
  
  path_bigwig_files %>%
    {.[stringr::str_detect(., strand)]} %>%
    purrr::map(rtracklayer::import, selection = bw_selection)
}

#' Read coverage on the specified region from bigwig files
#' @param bigwig_files path to the bigwig files
#' @param gr_region
read_bigwig_region <- function(bigwig_files, gr_region) {
  bw_selection <-
    gr_region %>%
    rtracklayer::BigWigSelection()
  
  bigwig_files %>%
    purrr::map(rtracklayer::import, selection = bw_selection)
}

#' Convert GenomicRanges to tibble
#' @param gr GRanges object
gr2tbl <- function(gr) {
  if(!inherits(gr, "GRanges")) stop("gr is not a GRanges object.")
  if(length(gr) == 0) return(NULL)
  tibble::as_tibble(gr) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coord = list(start:end)) %>%
    tidyr::unnest(cols = c(coord)) %>%
    dplyr::select(seqnames, coord, strand, score)
}

# "https://raw.githubusercontent.com/adibender/pammtools/master/R/ggplot-extensions.R" %>%
#   readr::read_lines() %>%
#   readr::write_lines(fs::path(wd, "script_r", "geom_stepribbon.R"))
source(fs::path(wd, "script_r", "geom_stepribbon.R"))

mutate_rel_to_maximum_roi <- function(tbl, start, end) {
  maximum <-
    tbl %>%
    dplyr::filter(start <= coord, coord <= end) %>%
    dplyr::pull(score) %>%
    max()
  tbl %>%
    dplyr::mutate(score2 = score / maximum)
}

theme_coverage <- function() {
  list(
    theme_minimal(),
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      line = element_blank(),
      legend.position = "",
      panel.background = element_blank(),
      plot.margin = unit(c(0, 0, -5, 0), "pt")
    ),
    scale_x_continuous(
      limits = c(start(gr_goi_all), end(gr_goi_all)),
      expand = expansion(mult = c(0, 0))
    ),
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(.1, .1))),
    scale_fill_manual(values = LABEL_PALETTE)
  )
}

