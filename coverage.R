suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

cat("[INFO] Script started\n")

# ---- Inputs ----
args <- commandArgs(trailingOnly = TRUE)

bed_file <- args[1]
bam <- args[2]
cnv <- args[3]
out <- args[4]

# ---- sanity check ----
if (!file.exists(bed_file)) stop("[ERROR] BED file not found: ", bed_file)
if (!file.exists(bam))      stop("[ERROR] BAM file not found: ", bam)
if (!file.exists(cnv))      stop("[ERROR] CNV file not found: ", cnv)
if (is.null(out))           stop("[ERROR] Output path (args[4]) is missing")

cat("[INFO] Input files OK\n")

# ---- read BED ----
cat("[INFO] Reading BED file...\n")
bed <- fread(bed_file, header = FALSE)
if (nrow(bed) == 0) stop("[ERROR] BED file is empty")

# ---- merge overlapping BED ----
cat("[INFO] Merging overlapping BED intervals...\n")

bed2 <- bed %>%
  arrange(V1, V2) %>%
  mutate(
    overlap = if_else(
      V1 == lag(V1) & V2 < lag(V3),
      TRUE,
      FALSE,
      missing = FALSE
    )
  )

bed3 <- bed2 %>%
  mutate(new_block = if_else(overlap, 0L, 1L)) %>%
  mutate(block_id = cumsum(new_block))

bed_merged <- bed3 %>%
  group_by(V1, block_id) %>%
  summarise(
    start = min(V2),
    end   = max(V3),
    .groups = "drop"
  ) %>%
  select(1, start, end)

colnames(bed_merged) <- c("chr", "start", "end")

if (nrow(bed_merged) == 0)
  stop("[ERROR] Merged BED is empty")

cat("[INFO] Merged BED intervals:", nrow(bed_merged), "\n")

# ---- samtools depth ----
temp_bed <- tempfile()
temp_cov <- tempfile()

write.table(bed_merged, file = temp_bed, sep="\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

cmd <- sprintf(
  "samtools depth -b %s -Q 20 -G 3332 %s > %s",
  temp_bed, bam, temp_cov
)

cat("[INFO] Running samtools depth...\n")
ret <- system(cmd)

if (ret != 0) stop("[ERROR] samtools depth failed")

if (!file.exists(temp_cov) || file.size(temp_cov) == 0)
  stop("[ERROR] samtools depth produced empty output")

average_bed <- fread(temp_cov)
mean <- mean(average_bed$V3)

cat("[INFO] Mean depth over capture regions:", mean, "\n")

# ---- samtools bedcov ----
cnv_cov_temp <- tempfile()
cmd_cnv <- sprintf(
  "samtools bedcov -Q 20 -G 3332 %s %s > %s",
  cnv, bam, cnv_cov_temp
)

cat("[INFO] Running samtools bedcov...\n")
ret <- system(cmd_cnv)

if (ret != 0) stop("[ERROR] samtools bedcov failed")

cnv_cov <- fread(cnv_cov_temp, header = FALSE) %>%
  select(1:3, last_col())

if (nrow(cnv_cov) == 0)
  stop("[ERROR] CNV coverage table is empty")

colnames(cnv_cov) <- c("chr", "start", "end", "coverage")

cat("[INFO] CNV events:", nrow(cnv_cov), "\n")

# ---- overlap calculation ----
cat("[INFO] Calculating overlap between CNV and capture regions...\n")

overlap <- c()

for (i in seq_len(nrow(cnv_cov))) {
  
  if (i %% 100 == 0)
    cat("[INFO] Processing CNV", i, "/", nrow(cnv_cov), "\n")
  
  event <- list(
    chr   = cnv_cov$chr[i],
    start = cnv_cov$start[i],
    end   = cnv_cov$end[i],
    depth = cnv_cov$coverage[i]
  )
  
  overlap_inside <- bed_merged %>%
    filter(chr == event$chr) %>%
    mutate(
      overlap = pmax(0L, pmin(end, event$end) - pmax(start, event$start))
    ) %>%
    pull(overlap)
  
  overlap <- c(overlap, sum(overlap_inside))
}

if (any(overlap == 0))
  warning("[WARN] Some CNVs have zero overlap with capture regions")

cnv_cov <- cnv_cov %>%
  mutate(overlap = overlap,
         K = coverage / overlap)

# ---- write output ----
write.table(
  cnv_cov,
  file = out,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

cat("[INFO] Finished successfully\n")
cat("[INFO] Output written to:", out, "\n")
