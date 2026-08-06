# Shared .xlsx writer/reader for the Extended Data tables.
#
# The six ExtendedDataTable*.R scripts all end by calling write_table_xlsx(),
# and the matching analysis/ExtendedDataTable*.Rmd pages all read the workbook
# back with read_table_xlsx(), so the published files share one format: a bold,
# frozen header row, blank (not "NA") cells for missing values, and column
# widths capped so a 200-character gene list doesn't push the numeric columns
# off screen.
#
# .xlsx rather than .csv because the Comments and lineage names carry non-ASCII
# characters (the Greek gamma/delta of "gd T", "TCRVg3"). A .csv has no place to
# record its encoding, so Excel opens a BOM-less UTF-8 file as legacy 8-bit
# text and those characters come out garbled; a workbook stores the encoding
# itself and round-trips them intact.

write_table_xlsx <- function(df, file, sheet = "Table", max_width = 60) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Please install it with install.packages('openxlsx').")
  }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df, colNames = TRUE, na.string = "")

  openxlsx::addStyle(wb, sheet, openxlsx::createStyle(textDecoration = "bold"),
                     rows = 1, cols = seq_along(df), gridExpand = TRUE)
  openxlsx::freezePane(wb, sheet, firstActiveRow = 2)

  # Width from the widest cell in the column, header included, capped: the
  # signature-gene columns run to several hundred characters.
  widths <- vapply(seq_along(df), function(j) {
    cells <- as.character(df[[j]])
    longest <- suppressWarnings(max(nchar(c(names(df)[j], cells)), na.rm = TRUE))
    min(max(longest, 8) + 2, max_width)
  }, numeric(1))
  openxlsx::setColWidths(wb, sheet, cols = seq_along(df), widths = widths)

  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  invisible(file)
}

# Read a table written by write_table_xlsx() back as a plain data.frame, with
# the column names left exactly as stored ("Top Genes +", "Prop. active cells")
# and the blank cells back as "" rather than the NA readxl would give -- an
# empty Comments cell means "no comment", and should preview as empty.
read_table_xlsx <- function(file, sheet = 1) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required. Please install it with install.packages('readxl').")
  }
  tab <- as.data.frame(readxl::read_excel(file, sheet = sheet, .name_repair = "minimal"),
                       check.names = FALSE, stringsAsFactors = FALSE)
  chr <- vapply(tab, is.character, logical(1))
  tab[chr] <- lapply(tab[chr], function(x) ifelse(is.na(x), "", x))
  tab
}
