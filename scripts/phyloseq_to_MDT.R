# Save a phyloseq object as upload files for the GBIF Metabarcoding Data Toolkit (MDT).
# See the MDT user guide, appendix "Tools for power users".
#
# phyloseq_to_MDT(physeq)                      -> one Excel workbook (needs openxlsx)
# phyloseq_to_MDT(physeq, format = "tsv")      -> OTU_table.tsv, Taxonomy.tsv, Samples.tsv
#
# The OTU table is written with OTUs as rows and samples as columns.
# Column names are not changed: rename them to Darwin Core terms before upload
# where needed, and add missing fields.

library(phyloseq)

phyloseq_to_MDT <- function(physeq, format = c("xlsx", "tsv"),
                            filename = "MDT_upload.xlsx", dir = ".") {
  format <- match.arg(format)

  otus <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) otus <- t(otus)   # OTUs as rows
  otus     <- as.data.frame(otus)
  taxonomy <- as.data.frame(as(tax_table(physeq), "matrix"))
  samples  <- as(sample_data(physeq), "data.frame")
  seqs <- refseq(physeq, errorIfNULL = FALSE)
  if (!is.null(seqs)) taxonomy$DNA_sequence <- as.character(seqs[rownames(taxonomy)])

  tables <- list(
    OTU_table = cbind(id = rownames(otus), otus),
    Taxonomy  = cbind(id = rownames(taxonomy), taxonomy),
    Samples   = cbind(id = rownames(samples), samples)
  )

  if (format == "xlsx") {
    wb <- openxlsx::createWorkbook()
    for (name in names(tables)) {
      openxlsx::addWorksheet(wb, name)
      openxlsx::writeData(wb, name, tables[[name]])
    }
    openxlsx::saveWorkbook(wb, file.path(dir, filename), overwrite = TRUE)
  } else {
    for (name in names(tables)) {
      write.table(tables[[name]], file.path(dir, paste0(name, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  invisible(tables)
}
