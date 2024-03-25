library(magrittr)

wb_create <- function(sheet_name = "Sheet 1", ...) {
  wb <- openxlsx::createWorkbook(...)
  if(!is.null(sheet_name)) {
    openxlsx::addWorksheet(wb = wb, sheetName = sheet_name)
  }
  wb
}

wb_add_sheet <- function(wb, sheet_name = "Sheet 1", ...) {
  wb2 <- wb$copy(shallow = FALSE)
  openxlsx::addWorksheet(wb = wb2, sheetName = sheet_name, ...)
  wb2
}

sheet2idx <- function(wb, sheet) {
  if(is.character(sheet)) {
    idx <- which(names(wb) == sheet)
  } else {
    idx <- sheet
  }
  idx
}

wb_get_row_max <- function(wb, sheet) {
  idx <- sheet2idx(wb, sheet)
  maximum <- suppressWarnings(suppressMessages(
    max(wb$worksheets[[idx]]$sheet_data$rows)
    ))
  ifelse(is.infinite(maximum), NA, maximum)
}

wb_get_col_max <- function(wb, sheet) {
  idx <- sheet2idx(wb, sheet)
  maximum <- suppressWarnings(suppressMessages(
    max(wb$worksheets[[idx]]$sheet_data$cols)
    ))
  ifelse(is.infinite(maximum), NA, maximum)
}


wb_write <- function(wb, sheet = 1, data, start_row = 2, style_header = NULL, ...) {
  if(!inherits(data, "data.frame")) stop("data must be a data.frame object.")
  
  col_max <- wb_get_col_max(wb, sheet)
  wb2 <- wb$copy(shallow = FALSE)
  openxlsx::writeData(
    wb = wb2, sheet = sheet, x = data,
    startCol = ifelse(is.na(col_max), 0, col_max) + 1,
    startRow = start_row, headerStyle = style_header,
    ...
  )
  wb2
}

wb_add_spacer_col <- function(wb, sheet = 1, width = 3) {
  row_max <- wb_get_row_max(wb, sheet)
  col_max <- wb_get_col_max(wb, sheet)
  
  wb2 <- wb$copy(shallow = FALSE)
  openxlsx::writeData(
    wb2, sheet,
    x = rep("", ifelse(is.na(row_max), 1, row_max)),
    startCol = ifelse(is.na(col_max), 0, col_max) + 1,
    startRow = 1
  )
  openxlsx::setColWidths(
    wb2, sheet,
    cols = ifelse(is.na(col_max), 0, col_max) + 1,
    widths = width
  )
  wb2
}

wb_merge_write <- function(wb, sheet = 1, rows, cols, x, style = NULL) {
  force(wb)
  row_max <- wb_get_row_max(wb, sheet)
  col_max <- wb_get_col_max(wb, sheet)
  
  if(all(rows < 0)) rows <- seq_len(row_max)[row_max + (rows + 1)]
  if(all(cols < 0)) cols <- seq_len(col_max)[col_max + (cols + 1)]
  mes_1 <- paste(range(rows), collapse = '-')
  mes_2 <- paste(range(cols), collapse = '-')
  message(stringr::str_glue("merge {mes_1} x {mes_2}"))
    
  row <- min(rows)
  col <- min(cols)
  wb2 <- wb$copy(shallow = FALSE)
  openxlsx::mergeCells(wb2, sheet, rows = rows, cols = cols)
  openxlsx::writeData(wb2, sheet, x, startRow = row, startCol = col)
  if(!is.null(style))
    openxlsx::addStyle(wb2, sheet, style = style, rows = row, cols = col, stack = TRUE)
  wb2
}

wb_apply_style <- function(wb, sheet = 1, rows, cols, style) {
  force(wb)
  row_max <- wb_get_row_max(wb, sheet)
  col_max <- wb_get_col_max(wb, sheet)
  
  if(all(rows < 0)) rows <- seq_len(row_max)[row_max + (rows + 1)]
  if(all(cols < 0)) cols <- seq_len(col_max)[col_max + (cols + 1)]
  mes_1 <- paste(range(rows), collapse = '-')
  mes_2 <- paste(range(cols), collapse = '-')
  message(stringr::str_glue("apply style to {mes_1} x {mes_2}"))
  
  wb2 <- wb$copy(shallow = FALSE)
  openxlsx::addStyle(wb2, sheet, style = style, rows = rows, cols = cols,
                     gridExpand = TRUE, stack = TRUE)
  wb2
}

wb_set_colwidth <- function(wb, sheet = 1, cols, width) {
  force(wb)
  col_max <- wb_get_col_max(wb, sheet)
  
  if(all(cols < 0)) cols <- seq_len(col_max)[col_max + (cols + 1)]
  mes_1 <- paste(range(cols), collapse = '-')
  message(stringr::str_glue("set column {mes_1} to width {width}"))
  
  wb2 <- wb$copy(shallow = FALSE)
  openxlsx::setColWidths(wb2, sheet, cols = cols, widths = width)
  wb2
}

wb_save <- function(wb, fpath, overwrite = TRUE, returnValue = FALSE) {
  openxlsx::saveWorkbook(wb = wb, file = fpath,
                         overwrite = overwrite, returnValue = returnValue)
}

### セルスタイル設定
# header
hs <- openxlsx::createStyle(
  halign = "center", valign = "center", textDecoration = "Bold")
# centerize text
center <- openxlsx::createStyle(halign = "center", valign = "center")
# cell border line
b <- openxlsx::createStyle(border = "LeftRight")
lb <- openxlsx::createStyle(border = "left", borderStyle = "hair")
rb <- openxlsx::createStyle(border = "right", borderStyle = "hair")

# 条件付き書式設定
up <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
dn <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")

# wb <-
#   wb_create(NULL) %>%
#   wb_add_sheet() %>%
#   wb_add_spacer_col() %>%
#   wb_write(data = iris) %>%
#   wb_add_spacer_col() %>%
#   wb_write(data = iris[1:50,])
# wb_save(wb, "~/Desktop/temp.xlsx")
