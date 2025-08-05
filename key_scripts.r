# =============================================================================
# scripts/01_define_genomic_region.R
# =============================================================================

#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  library(GenomicRanges)
  library(biomaRt)
  library(rtracklayer)
})

# Get parameters from Snakemake
gene_name <- snakemake@params[["gene"]]
window_size <- as.numeric(snakemake@params[["window"]])
genome_build <- snakemake@params[["genome"]]

# Connect to Ensembl
if (genome_build == "hg19") {
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                     dataset = "hsapiens_gene_ensembl",
                     host = "grch37.ensembl.org")
} else {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
}

# Get gene coordinates
gene_info <- getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                                  "start_position", "end_position", "strand"),
                   filters = "hgnc_symbol",
                   values = gene_name,
                   mart = ensembl)

if (nrow(gene_info) == 0) {
  stop(paste("Gene", gene_name, "not found in Ensembl"))
}

# Take the first entry if multiple
gene_info <- gene_info[1, ]

# Define analysis region
chr <- paste0("chr", gene_info$chromosome_name)
tss <- ifelse(gene_info$strand == 1, gene_info$start_position, gene_info$end_position)
region_start <- max(1, tss - window_size)
region_end <- tss + window_size

# Create GRanges object
analysis_region <- GRanges(
  seqnames = chr,
  ranges = IRanges(start = region_start, end = region_end),
  strand = "*",
  gene = gene_name,
  tss = tss
)

# Save as RDS for R scripts
saveRDS(analysis_region, snakemake@output[["granges"]])

# Save as BED for other tools
bed_df <- data.frame(
  chr = as.character(seqnames(analysis_region)),
  start = start(analysis_region) - 1,  # Convert to 0-based
  end = end(analysis_region),
  name = paste0(gene_name, "_region"),
  score = 1000,
  strand = "."
)

write.table(bed_df, snakemake@output[["bed"]], 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Analysis region defined:", chr, ":", region_start, "-", region_end, "\n")

# =============================================================================
# scripts/03_conservation_analysis.R  
# =============================================================================

#!/usr/bin/env Rscript

suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
})

# Load conservation scores
conservation_bed <- fread(snakemake@input[["scores"]], 
                         col.names = c("chr", "start", "end", "score"))

# Load genomic region
region <- readRDS(snakemake@input[["region"]])  # Assuming we also save as RDS

# Parameters
threshold <- as.numeric(snakemake@params[["threshold"]])
tile_size <- as.numeric(snakemake@params[["tile_size"]])

# Create tiles across the region
region_chr <- as.character(seqnames(region))
region_start <- start(region)
region_end <- end(region)

# Generate tiles
tile_starts <- seq(region_start, region_end - tile_size, by = tile_size)
tile_ends <- tile_starts + tile_size - 1
tile_ends[length(tile_ends)] <- min(tile_ends[length(tile_ends)], region_end)

tiles <- GRanges(
  seqnames = region_chr,
  ranges = IRanges(start = tile_starts, end = tile_ends)
)

# Calculate average conservation score per tile
conservation_gr <- GRanges(
  seqnames = conservation_bed$chr,
  ranges = IRanges(start = conservation_bed$start + 1, end = conservation_bed$end),
  score = conservation_bed$score
)

# Find overlaps and calculate average scores
overlaps <- findOverlaps(tiles, conservation_gr)
tile_scores <- aggregate(
  conservation_gr$score[subjectHits(overlaps)], 
  by = list(queryHits(overlaps)), 
  FUN = mean
)

# Add scores to tiles
tiles$avg_score <- 0
tiles$avg_score[tile_scores$Group.1] <- tile_scores$x

# Identify high conservation tiles
high_conservation_tiles <- tiles[tiles$avg_score > threshold]

# Save results
tiles_df <- data.frame(
  chr = as.character(seqnames(tiles)),
  start = start(tiles) - 1,
  end = end(tiles),
  score = tiles$avg_score
)

high_cons_df <- data.frame(
  chr = as.character(seqnames(high_conservation_tiles)),
  start = start(high_conservation_tiles) - 1,
  end = end(high_conservation_tiles),
  name = paste0("high_cons_", seq_along(high_conservation_tiles)),
  score = high_conservation_tiles$avg_score,
  strand = "."
)

write.table(tiles_df, snakemake@output[["tiles"]], 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(high_cons_df, snakemake@output[["high_cons"]], 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Identified", nrow(high_cons_df), "high conservation regions\n")

# =============================================================================
# scripts/04_motif_scanning.R
# =============================================================================

#!/usr/bin/env Rscript

suppressMessages({
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)  # or hg38
  library(motifmatchr)
  library(TFBSTools)
  library(