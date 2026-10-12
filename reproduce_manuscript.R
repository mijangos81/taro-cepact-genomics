# ============================================================================
# Reproduction script for
#   Mijangos, Restall, Kilian & Campbell. Genome-wide SNPs reveal clonal
#   redundancy and connectivity in the world's largest taro collection.
#   Horticulture Research.
#
# Runs from the repository root and regenerates every figure, table and
# reported number of the main text and the Supplementary Information, in the
# order the analyses depend on each other:
#   1  Setup: packages, seed, directories, shared helpers
#   2  Dataset assembly (Methods) and Fig. 1
#   3  Ploidy callers (SI S1; Figs S1, S2)
#   4  Diploid working set (SI S2; Figs S3, S4; Supplementary Data 4)
#   5  Qploidy (SI S1.4)
#   6  Clonal lineages, cultivar synonymy and the candidate core
#      (Figs 2, 3, 5, S6; Supplementary Table 3; Supplementary Data 1, 2, 5, 6)
#   7  Ploidy threshold and redundancy (SI S2.1; Table S4 rows)
#   8  Population structure (Figs 6, 7, 8, S7, S8, S9; Supplementary Data 3)
#   9  Diversity by country and cluster, manuscript tables
#      (Table 1, Fig. 4; Tables S1, S2, S4, S5; Supplementary Data 7, 8)
#  10  Selection scan and gene context (Fig. S5; SI S3)
#
# Inputs: the two DArTseq SNP reports and the three count reports in
# input_files/, the sample metadata and country coordinates in info/, and the
# flow cytometry cytotypes in outputs/ploidy.csv. Blocks that need files not
# deposited with the paper (read alignments, BAM files) are marked.
# Binaries on PATH: plink, qpas, nemeco, selscan and convert (OHANA; these two
# names clash with the iHS selscan and ImageMagick, so the OHANA folder must
# come first on PATH), fitch and consense (PHYLIP), blastn, blastp,
# makeblastdb, samtools, bcftools, PopClusterMac / PopClusterLnx /
# PopClusterWin, curl. gbs2ploidy is installed from the CRAN archive
# (https://cran.r-project.org/src/contrib/Archive/gbs2ploidy/) and
# gl.find.genes.for.loci needs dartR.popgen from the dev branch (PR #87).
# ============================================================================

library(ape)
library(colorspace)
library(dartR.base)
library(dartR.popgen)
library(dartRverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(ggtree)
library(gridExtra)
library(LEA)
library(nQuack)
library(patchwork)
library(poppr)
library(Qploidy2)
library(R.utils)
library(rjags)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(sf)
library(tidyr)
library(vcfR)

set.seed(42)

# Output directories, including those the external binaries write into
for (d in c("input_files", "outputs", "figures", "data/annotation",
            "outputs/bam_ploidy", "outputs/runnQuack_output/02_prepared/temp",
            "outputs/runnQuack_output/03_processed", "outputs/runnQuack_output/04_output",
            "outputs/ohana", "outputs/ohana_declone", "outputs/ohana_multistart/runs",
            "outputs/ohana_multistart/best", "outputs/snmf_kselect", "outputs/snmf_declone",
            "outputs/popcluster", "outputs/ohana_consensus",
            paste0("outputs/ohana_selection_", rep(c("genets", "consensus", "accessions"), each = 3),
                   "/", c("multistart", "best", "component_scans")),
            "outputs/ohana_candidate_genes/blastdb", "outputs/ohana_candidate_genes_ahi"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- Shared helpers ---

# Display names for the run-together country labels stored in pop()
.country_full <- c(
  CookIslands = "Cook Islands", FederatedStatesofMicronesia = "Micronesia (FSM)",
  IndiaNorthen = "India (northern)", NewCaledonia = "New Caledonia",
  PapuaNewGuinea = "Papua New Guinea", SolomonIslands = "Solomon Islands",
  "Tahiti&FrenchPolynesia" = "Tahiti", UnitedStatesofAmerica = "United States",
  WSamoa = "Western Samoa")
.country_abbr <- c(
  CookIslands = "Cook Is.", FederatedStatesofMicronesia = "FSM",
  IndiaNorthen = "N. India", NewCaledonia = "New Caledonia",
  PapuaNewGuinea = "PNG", SolomonIslands = "Solomon Is.",
  "Tahiti&FrenchPolynesia" = "Tahiti", UnitedStatesofAmerica = "USA",
  WSamoa = "W. Samoa")
pretty_country <- function(x, short = FALSE) {
  m   <- if (short) .country_abbr else .country_full
  x   <- as.character(x)
  out <- m[x]
  ifelse(is.na(out), x, unname(out))
}

# Polychrome palette (darkened 20%), recycled to n colours
polychrome <- function(n) {
  pal <- c("#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF",
           "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356",
           "#85660D", "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF",
           "#F7E1A0", "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C",
           "#B5EFB5", "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF",
           "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3B00FB")
  names(pal) <- c(
    "Dark_Purplish_Gray", "Purplish_White", "Vivid_Red", "Vivid_Purple",
    "Vivid_Yellowish_Green", "Strong_Purplish_Blue", "Vivid_Orange_Yellow",
    "Vivid_Purplish_Red", "Brilliant_Green", "Vivid_Yellow_Green", "Vivid_Blue",
    "Brilliant_Purple", "Vivid_Violet", "Strong_Pink", "Strong_Blue",
    "Strong_Reddish_Orange", "Vivid_Green", "Light_Olive_Brown",
    "Vivid_Reddish_Purple", "Vivid_Greenish_Yellow", "Vivid_Yellowish_Green",
    "Vivid_Red", "Vivid_Purplish_Red", "Pale_Yellow", "Strong_Reddish_Purple",
    "Vivid_Violet", "Vivid_Yellow_Green", "Very_Light_Blue",
    "Strong_Reddish_Brown", "Very_Light_Yellowish_Green",
    "Very_Light_Bluish_Green", "Deep_Greenish_Blue", "Vivid_Purple",
    "Deep_Purple", "Brilliant_Blue", "Vivid_Violet")
  pal2 <- darken(pal, amount = 0.2)
  rep(pal2, ceiling(n / length(pal2)))[1:n]
}

# Colour-blind-friendly population palette; HCL fallback beyond 16 colours
pop_palette <- function(n) {
  base_cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E41A1C",
                 "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  if (n <= length(base_cols)) {
    base_cols[seq_len(n)]
  } else {
    h <- seq(0, 360, length.out = n + 1)[-1]
    grDevices::hcl(h = h, c = 70, l = 55)
  }
}

# ============================================================================
# 2. Dataset assembly (Methods) and Fig. 1
# ============================================================================

## Methods: read the two DArTseq SNP reports (same samples, disjoint marker panels)
g4 <- gl.read.dart("input_files/Report_DTa25-11061_3_moreOrders_SNP_mapping_4.csv", verbose = 0)
saveRDS(g4, "input_files/taro_report4_raw.rds")
g5 <- gl.read.dart("input_files/Report_DTa25-11061_3_moreOrders_SNP_mapping_5.csv", verbose = 0)
saveRDS(g5, "input_files/taro_report5_raw.rds")

# Per-individual call rate from SNPbin NA positions
callrate_ind <- function(g) 1 - vapply(g@gen, function(x) length(x@NA.posi), 1) / nLoc(g)

# Collapse technical replicates: keep the highest-call-rate instance per canonical name
collapse_reps <- function(g, canonize) {
  cn <- canonize(indNames(g))
  cr <- callrate_ind(g)
  keep <- tapply(seq_along(cn), cn, function(ix) ix[which.max(cr[ix])])
  gk <- g[as.integer(keep), ]
  indNames(gk) <- names(keep)
  gk
}
canon4 <- function(x) sub("_$", "", sub("_[a-dA-D]$", "", x))   # report 4: _a/_b/_c/_d suffix
canon5 <- function(x) sub("\\.[0-9]+$", "", x)                   # report 5: gl.read.dart .N uniquifier
g4c <- collapse_reps(g4, canon4)
g5c <- collapse_reps(g5, canon5)

# Common individuals in identical order, then join the two marker panels side by side
common <- intersect(indNames(g4c), indNames(g5c))
g4c <- g4c[match(common, indNames(g4c)), ]
g5c <- g5c[match(common, indNames(g5c)), ]
gj <- gl.join(g4c, g5c, verbose = 0)
gj <- gl.compliance.check(gj, verbose = 0)
saveRDS(gj, "input_files/taro_joined_raw.rds")

## Methods: attach metadata and keep Colocasia esculenta (1,282 accessions)
t1 <- readRDS("input_files/taro_joined_raw.rds")
meta <- read.csv("info/metadata_DTa23-8697.csv", stringsAsFactors = FALSE, check.names = FALSE)
meta$canon <- sub("_$", "", sub("_[a-dA-D]$", "", trimws(meta$id)))
meta1 <- meta[!duplicated(meta$canon), ]
met <- meta1[match(indNames(t1), meta1$canon), ]
met$id <- indNames(t1)
t1$other$ind.metrics <- met
pop(t1) <- t1$other$ind.metrics$Species
t2 <- gl.keep.pop(t1, pop.list = "Colocasia esculenta", verbose = 0)
pop(t2) <- t2$other$ind.metrics$pop

# Country coordinates (left join keeps every sample)
coord <- read.csv("info/coor_countries.csv")
coord$pop <- trimws(coord$pop, which = "both")
t2$other$ind.metrics$order <- 1:nInd(t2)
met_tmp <- t2$other$ind.metrics
met_tmp$pop <- as.character(met_tmp$pop)
met_tmp2 <- merge(met_tmp, coord, by = "pop", all.x = TRUE)
met_tmp2 <- met_tmp2[order(met_tmp2$order, decreasing = FALSE), ]
t2$other$ind.metrics <- met_tmp2
t2@other$latlon <- data.frame(lat = met_tmp2[, "lat.y"], lon = met_tmp2[, "lon.y"])

## Methods: QC filters (read depth, reproducibility >= 0.98, one SNP per tag, call rate ind >= 0.75, loc >= 0.80)
t3 <- gl.filter.rdepth(t2, verbose = 0)
t3 <- gl.filter.reproducibility(t3, threshold = 0.98, verbose = 0)
t3 <- gl.filter.secondaries(t3, method = "best", verbose = 0)
t3 <- gl.filter.callrate(t3, threshold = 0.75, method = "ind", verbose = 0)
t3 <- gl.filter.callrate(t3, method = "loc", threshold = 0.8, verbose = 0)
t3$other$loc.metrics$TrimmedSequence <- t3$other$loc.metrics$TrimmedSequenceSnp

# Clean population and individual names
popNames(t3) <- gsub(" ", "", popNames(t3))
popNames(t3) <- gsub("\\(", "_", popNames(t3))
popNames(t3) <- gsub("\\)", "", popNames(t3))
popNames(t3) <- gsub("\\/", "&", popNames(t3))
popNames(t3) <- gsub(",", "_", popNames(t3))
popNames(t3) <- gsub("\\.", "_", popNames(t3))
popNames(t3) <- gsub("_", "", popNames(t3))
indNames(t3) <- gsub(" ", "_", indNames(t3))
indNames(t3) <- gsub("_", "", indNames(t3))
indNames(t3) <- sub("a$", "", indNames(t3))

## Methods: keep loci on the 14 assembled chromosomes, drop unplaced, duplicate-coordinate and monomorphic loci -> 1,212 x 16,910
idx <- grep("chromosome", t3$other$loc.metrics$Chrom_Colocasia_esculenta_JAAS_v1.0)
t4 <- gl.keep.loc(t3, loc.list = locNames(t3)[idx], verbose = 0)
chr_tmp <- sub("^[^_]*_[^_]*_", "", t4$other$loc.metrics$Chrom_Colocasia_esculenta_JAAS_v1.0)
t4$chromosome <- as.factor(chr_tmp)
t4$position <- t4$other$loc.metrics$ChromPosSnp_Colocasia_esculenta_JAAS_v1.0
pos0 <- which(t4$position == 0)
if (length(pos0)) t4 <- gl.drop.loc(t4, loc.list = locNames(t4)[pos0], verbose = 0)
t4$other$loc.metrics$loc <- paste0(as.character(t4$chromosome), "_", t4$position)
locNames(t4) <- t4$other$loc.metrics$loc
dup <- which(duplicated(locNames(t4)) == TRUE)
if (length(dup)) t4 <- gl.drop.loc(t4, loc.list = locNames(t4)[dup], verbose = 0)
t4 <- gl.filter.monomorphs(t4, verbose = 0)
t4 <- t4[, order(as.numeric(t4$chromosome), t4$position, decreasing = FALSE)]
saveRDS(t4, "input_files/taro.rds")

## Fig. 1: geographic sampling of the taro collection (figures/Fig_taro_map.png)
gl  <- readRDS("input_files/taro.rds")
pts <- data.frame(lon = gl@other$latlon$lon, lat = gl@other$latlon$lat,
                  location = as.character(pop(gl)))
pts <- pts[!is.na(pts$lon) & !is.na(pts$lat), ]

# Pacific-centred longitude; one point per country centroid sized by sample count
wrap360 <- function(x) ifelse(x < 0, x + 360, x)
pts$lon360 <- wrap360(pts$lon)
locs <- pts %>%
  group_by(location) %>%
  summarise(lon360 = mean(lon360), lat = mean(lat), n = n(), .groups = "drop") %>%
  arrange(desc(n))

# Indo-Pacific extent; countries outside it are not drawn
xlim <- c(78, 213); ylim <- c(-30, 40)
off <- locs %>% filter(lon360 < xlim[1] | lon360 > xlim[2] | lat < ylim[1] | lat > ylim[2])
locs_in <- locs %>% filter(!(location %in% off$location))
loc_cols <- setNames(polychrome(nrow(locs_in)), locs_in$location)

# Natural Earth land, cropped into two dateline-safe boxes and the eastern box shifted to [180, 213]
sf::sf_use_s2(FALSE)
land <- st_make_valid(ne_countries(scale = "medium", returnclass = "sf"))
west <- st_crop(land, xmin = xlim[1], ymin = ylim[1], xmax = 180, ymax = ylim[2])
east <- st_crop(land, xmin = -180, ymin = ylim[1], xmax = xlim[2] - 360, ymax = ylim[2])
east <- st_shift_longitude(east)
land <- rbind(west, east)

p <- ggplot() +
  geom_sf(data = land, fill = "#efeee9", colour = "#b9b6ad", linewidth = 0.2) +
  geom_point(data = locs_in, aes(x = lon360, y = lat, colour = location, size = n),
             alpha = 0.85) +
  geom_text_repel(data = locs_in,
                  aes(x = lon360, y = lat, label = location, colour = location),
                  fontface = "bold", size = 3.3, max.overlaps = Inf,
                  box.padding = 0.5, min.segment.length = 0,
                  segment.colour = "grey50", segment.size = 0.3, show.legend = FALSE) +
  scale_colour_manual(values = loc_cols) +
  scale_size_area(max_size = 11, name = "Samples (n)", breaks = c(10, 50, 100, 200)) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   pad_x = unit(0.4, "cm"), pad_y = unit(0.4, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1.0, "cm"), width = unit(1.0, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  guides(colour = "none") +
  theme_bw(base_size = 12) +
  theme(panel.background  = element_rect(fill = "#dbe9f2"),
        panel.grid.major  = element_line(colour = "#ffffff", linewidth = 0.2),
        axis.title        = element_blank(),
        legend.position   = c(0.08, 0.22),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA))
ggsave("figures/Fig_taro_map.png", p, width = 9.0, height = 6.0, units = "in",
       dpi = 300, bg = "white")

# ============================================================================
# 3. Ploidy callers (SI S1; Figs S1, S2)
# ============================================================================

# read a DArT two-row counts report (odd row = reference allele, even row = alternate allele) for the requested individuals
read.dart.counts <- function(counts.file, lastmetric = NULL, ind.vec, id = 4) {
  tdummy <- read.csv(counts.file, na.strings = "-", check.names = FALSE, nrows = 20,
                     header = FALSE, stringsAsFactors = TRUE)
  colskip <- sum(tdummy[1, ] == "*") + 1
  topskip <- sum(tdummy[, 1] == "*")
  snpraw <- read.csv(counts.file, na.strings = "-", skip = topskip, check.names = FALSE,
                     stringsAsFactors = TRUE)
  counts_names <- as.character(unname(unlist(tdummy[id, colskip:ncol(tdummy)])))
  if (is.null(lastmetric)) {
    top <- read.csv(counts.file, header = FALSE, nrows = 20, stringsAsFactors = FALSE)
    lastmetric <- top[max(which(top[, 1] == "*")) + 1, max(which(top[1, ] == "*"))]
  }
  lmet <- which(colnames(snpraw) == lastmetric)
  covmetrics <- snpraw[, 1:lmet]
  datas <- snpraw[, (lmet + 1):ncol(snpraw)]
  dup.names <- which(!duplicated(counts_names))
  datas <- datas[, dup.names]
  colnames(datas) <- counts_names[dup.names]
  datas <- datas[, which(colnames(datas) %in% make.unique(ind.vec))]
  list(counts = datas, loc_metrics = covmetrics)
}

# split the allele rows into two individuals x loci count matrices (allele_1 = even rows, allele_2 = odd rows)
CountsSplit <- function(count.data) {
  cdata <- count.data$counts
  esl <- seq(2, nrow(cdata), 2)
  list(allele_1 = as.matrix(t(cdata[esl, ])), allele_2 = as.matrix(t(cdata[-esl, ])))
}

# expected fraction of true 0.5 heterozygotes with allele ratio in [0.42, 0.58], binomial at the observed depths
exp_frac_half <- function(depths) {
  depths <- depths[is.finite(depths) & depths > 0]
  if (!length(depths)) return(NA_real_)
  mean(sapply(depths, function(d) {
    d <- round(d); lo <- ceiling(0.42 * d); hi <- floor(0.58 * d)
    if (hi < lo) 0 else sum(dbinom(lo:hi, d, 0.5))
  }))
}

# flow cytometry truth set (65 individuals) and the three DArT count reports
truth <- read.csv("outputs/ploidy.csv", stringsAsFactors = FALSE)
colnames(truth) <- c("id", "ploidy_fc")
truth$ploidy_fc <- trimws(truth$ploidy_fc)
COUNT_FILES <- c("input_files/SEQ_SNPs_counts_0_Target.csv",
                 "input_files/SEQ_SNPs_counts_1_Target.csv",
                 "input_files/SEQ_SNPs_counts_2_Target.csv")

## ---- gbs2ploidy (Gompert & Mock 2017) on the DArT count reports ----

# Colocasia esculenta samples (canonical IDs) from the sample metadata
meta <- read.csv("info/metadata_DTa23-8697.csv", stringsAsFactors = FALSE)
meta$canon <- sub("_$", "", sub("_[a-dA-D]$", "", trimws(meta$id)))
ind_vec <- unique(meta$canon[meta$Species == "Colocasia esculenta"])

# allele-1 / allele-2 count matrices over the three reports
parts <- lapply(COUNT_FILES, function(f) read.dart.counts(f, ind.vec = ind_vec))
count.data_all <- list(counts = do.call(rbind, lapply(parts, `[[`, "counts")),
                       loc_metrics = do.call(rbind, lapply(parts, `[[`, "loc_metrics")))
dat <- CountsSplit(count.data_all)
cov1_mat <- dat$allele_1
cov2_mat <- dat$allele_2

# locus QC: mean depth 10-60, call rate >= 0.9, MAF >= 0.05
tot_depth <- cov1_mat + cov2_mat
mean_depth <- colMeans(tot_depth, na.rm = TRUE)
call_rate <- colMeans(!is.na(tot_depth) & tot_depth > 0)
ref_depth <- colSums(cov1_mat, na.rm = TRUE)
alt_depth <- colSums(cov2_mat, na.rm = TRUE)
allele2_freq <- alt_depth / (ref_depth + alt_depth)
maf <- pmin(allele2_freq, 1 - allele2_freq)
keep_loci <- which(mean_depth >= 10 & mean_depth <= 60 & call_rate >= 0.9 & maf >= 0.05)
cov1_use <- cov1_mat[, keep_loci, drop = FALSE]
cov2_use <- cov2_mat[, keep_loci, drop = FALSE]

# keep heterozygous genotypes only: depth >= 10 and allele proportion in (0.1, 0.9)
tot_use <- cov1_use + cov2_use
p_alt <- cov2_use / tot_use
drop_gt <- tot_use < 10 | is.na(tot_use) | p_alt <= 0.1 | p_alt >= 0.9
cov1_use[drop_gt] <- NA
cov2_use[drop_gt] <- NA

# posterior weights on five dosage categories per individual (rjags MCMC)
props <- gbs2ploidy::estprops(cov1 = t(cov1_use), cov2 = t(cov2_use),
                              props = c(0.25, 0.33, 0.5, 0.66, 0.75),
                              mcmc.nchain = 2, mcmc.steps = 10000, mcmc.burnin = 1000, mcmc.thin = 2)

# heterozygosity and mean depth per individual; drop individuals without informative loci
hx <- apply(!is.na(cov1_use + cov2_use), 1, mean)
dx <- apply(cov1_use + cov2_use, 1, mean, na.rm = TRUE)
keep_ind <- which(!is.na(dx))

## Results (S1.1): gbs2ploidy cytotype classes, unsupervised, three classes on PC1-PC2
pl <- gbs2ploidy::estploidy(alphas = props[keep_ind], het = hx[keep_ind], depth = dx[keep_ind],
                            train = FALSE, pl = NA, set = NA, nclasses = 3, pcs = 1:2)
m <- pl$pp
m[, 2:4] <- round(m[, 2:4], 2)
m <- as.data.frame(m)
m$ids <- rownames(m)
m$ploidy <- apply(m[, 2:4], 1,
                  function(x) if (max(x, na.rm = TRUE) < 0.7) NA_integer_ else which.max(x) + 1)
write.csv(m, "outputs/m_ploidy.csv", row.names = FALSE)

## Results (S1.1): gbs2ploidy class (2 = diploid, 3 = triploid, 4 = tetraploid, NA = uncertain) vs flow cytometry
gbs <- merge(truth, m[, c("ids", "ploidy")], by.x = "id", by.y = "ids")
gbs2ploidy_tab <- table(flow_cytometry = gbs$ploidy_fc, gbs2ploidy = gbs$ploidy, useNA = "ifany")
write.csv(as.data.frame.matrix(gbs2ploidy_tab), "outputs/gbs2ploidy_vs_flow_cytometry.csv")

## ---- nQuack (Gaynor et al. 2024) on the repeat-masked alignments; BAM files are not deposited ----
# alignments: bam_files/repeats_removed/<TargetID>.FASTQ.gz_sort_removed.bam (bowtie2 to GCA_014218235.1, RepeatMasker regions removed)
bam_dir <- "bam_files/repeats_removed"
qdir <- "outputs/runnQuack_output"
bam_list <- sub("[.]bam$", "", list.files(bam_dir, pattern = "[.]bam$"))
for (s in bam_list)
  prepare_data(name = s, inpath = paste0(bam_dir, "/"), outpath = paste0(qdir, "/02_prepared/"),
               tempfolder = paste0(qdir, "/02_prepared/temp"))

# per-site allele counts: total depth >= 2, depth <= 90th percentile, error 0.01, no allele-frequency truncation
for (f in list.files(paste0(qdir, "/02_prepared"), pattern = "[.]txt$")) {
  xm <- process_data(paste0(qdir, "/02_prepared/", f), min.depth = 2, max.depth.quantile.prob = 0.9,
                     error = 0.01, trunc = c(0, 0))
  write.csv(xm, paste0(qdir, "/03_processed/", sub("[.]txt$", ".csv", f)), row.names = FALSE)
}

# normal, beta and beta-binomial mixtures per sample
for (f in list.files(paste0(qdir, "/03_processed"), pattern = "[.]csv$")) {
  xm <- as.matrix(read.csv(paste0(qdir, "/03_processed/", f)))
  allout <- rbind(quackNormal(xm = xm, samplename = f, cores = 10, parallel = FALSE),
                  quackBeta(xm = xm, samplename = f, cores = 10, parallel = FALSE),
                  quackBetaBinom(xm = xm, samplename = f, cores = 10, parallel = FALSE))
  write.csv(allout, paste0(qdir, "/04_output/", f), row.names = FALSE)
}

# lowest-BIC mixture per sample, joined to the flow cytometry set through the TargetID crosswalk
top_rows <- lapply(list.files(paste0(qdir, "/04_output"), pattern = "[.]csv$"), function(f) {
  df <- read.csv(paste0(qdir, "/04_output/", f), stringsAsFactors = FALSE)
  df <- df[order(df$BIC, decreasing = FALSE), ][1, , drop = FALSE]
  df$file_name <- sub("\\..*$", "", f)
  df
})
final_df <- do.call(rbind, top_rows)
write.csv(final_df, "outputs/final_combined_top_rows.csv", row.names = FALSE)
r1 <- read.csv("info/ID_targets.csv")
r2 <- truth
colnames(r2) <- c("Genotype", "ploidy")
r3 <- merge(r2, r1, by = "Genotype")
r4 <- final_df[, c("mixture", "file_name")]
colnames(r4) <- c("mixture", "TargetID")
r5 <- merge(r3, r4, by = "TargetID")
write.csv(r5, "outputs/ploidy_taro.csv")

## Results (S1.1): nQuack cytotype per verified record vs flow cytometry
nq <- r5[, c("Genotype", "ploidy", "mixture")]
colnames(nq) <- c("Genotype", "Reported ploidy", "nQuack")
write.csv(nq, "outputs/ploidy_taro_nQuack.csv", row.names = FALSE)
nquack_tab <- table(flow_cytometry = nq[["Reported ploidy"]], nQuack = nq$nQuack)
write.csv(as.data.frame.matrix(nquack_tab), "outputs/nquack_vs_flow_cytometry.csv")

## ---- Allele balance at heterozygous sites of the flow cytometry set (Fig. S1) ----

# allele-1 / allele-2 read matrices (loci x individuals) for the 65 verified individuals
ind_fc <- unique(truth$id)
parts <- lapply(COUNT_FILES, function(f) read.dart.counts(f, ind.vec = ind_fc))
counts_all <- do.call(rbind, lapply(parts, `[[`, "counts"))
nr <- nrow(counts_all)
a1 <- as.matrix(counts_all[seq(1, nr, 2), , drop = FALSE])
a2 <- as.matrix(counts_all[seq(2, nr, 2), , drop = FALSE])
storage.mode(a1) <- "numeric"; storage.mode(a2) <- "numeric"
ids <- colnames(counts_all)
tot <- a1 + a2
ratio <- a2 / tot

## Results (S1.2): median read depth per genotype
median_depth <- median(tot[!is.na(tot) & tot > 0])

## Results (S1.2): fraction of heterozygous sites in [0.42, 0.58] per depth threshold vs the binomial expectation
THRESHOLDS <- c(20, 40, 60, 100)
summ <- list()
for (TH in THRESHOLDS) {
  het <- !is.na(tot) & tot >= TH & a1 > 0 & a2 > 0
  idx <- which(het, arr.ind = TRUE)
  long <- data.frame(id = ids[idx[, 2]], ratio = ratio[het], depth = tot[het])
  long <- merge(long, truth, by = "id")
  long <- long[long$ploidy_fc %in% c("Diploid", "Triploid"), ]
  band <- long %>%
    filter(ratio > 0.15, ratio < 0.85) %>%
    group_by(id, ploidy_fc) %>%
    summarise(n_het = n(), frac_half = mean(ratio >= 0.42 & ratio <= 0.58), .groups = "drop")
  grp <- band %>% group_by(ploidy_fc) %>%
    summarise(n_ind = n(), med_n_het = median(n_het), med_frac_half = median(frac_half))
  dip_depths <- long$depth[long$ploidy_fc == "Diploid" & long$ratio > 0.15 & long$ratio < 0.85]
  grp$binom_expected_half <- ifelse(grp$ploidy_fc == "Diploid", exp_frac_half(dip_depths), NA)
  grp$threshold <- TH
  summ[[as.character(TH)]] <- grp
}
summ_df <- do.call(rbind, summ)
write.csv(summ_df, "outputs/ploidy_signal_depth_summary.csv", row.names = FALSE)

## Fig. S1: half-band fraction of verified diploids and triploids against the depth threshold
summ_df$ploidy_fc <- factor(summ_df$ploidy_fc, levels = c("Diploid", "Triploid"))
exp_df <- summ_df[summ_df$ploidy_fc == "Diploid", c("threshold", "binom_expected_half")]
p <- ggplot(summ_df, aes(threshold, med_frac_half, colour = ploidy_fc)) +
  geom_line(data = exp_df, aes(threshold, binom_expected_half),
            inherit.aes = FALSE, linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n_ind)) +
  annotate("text", x = max(exp_df$threshold), y = exp_df$binom_expected_half[nrow(exp_df)],
           label = "binomial expectation\n(true 0.5 heterozygote)",
           hjust = 1, vjust = 1.4, size = 3, colour = "grey40") +
  scale_colour_manual(values = c(Diploid = "#0072B2", Triploid = "#D55E00"), name = NULL) +
  scale_size_area(max_size = 6, name = "Individuals (N)") +
  scale_x_continuous(breaks = THRESHOLDS) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Minimum read-depth threshold", y = "Fraction of heterozygous sites in [0.42, 0.58]") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave("figures/ploidy_signal_depth.pdf", p, width = 7.5, height = 5)
ggsave("figures/ploidy_signal_depth.png", p, width = 7.5, height = 5, dpi = 300)

## ---- HDplot paralog screen (McKinney et al. 2017) on the same read matrices (Fig. S2) ----
GENO_MINDEPTH <- 8
HET_LO <- 0.10; HET_HI <- 0.90
H_MAX <- 0.50
D_MAX <- 7
nloc <- nrow(a1)

# per-locus heterozygosity H and read-ratio deviation D over the heterozygous genotypes
geno_ok <- !is.na(tot) & tot >= GENO_MINDEPTH
het <- geno_ok & ratio > HET_LO & ratio < HET_HI
hetN <- het; hetN[is.na(hetN)] <- FALSE
a1z <- a1; a1z[!hetN] <- 0; a1z[is.na(a1z)] <- 0
a2z <- a2; a2z[!hetN] <- 0; a2z[is.na(a2z)] <- 0
n_geno <- rowSums(geno_ok, na.rm = TRUE)
n_het <- rowSums(hetN)
H <- ifelse(n_geno > 0, n_het / n_geno, NA)
A <- rowSums(a1z); B <- rowSums(a2z)
D <- ifelse((A + B) > 0, (A - B) / sqrt(A + B), NA)
mean_ratio_het <- ifelse(n_het > 0, rowSums(ratio * hetN, na.rm = TRUE) / n_het, NA)
hd <- data.frame(locus = seq_len(nloc), n_geno, n_het, H, A, B, D, mean_ratio_het)

## Results (S1.3): evaluable loci (n_het >= 5, n_geno >= 20) classed as singleton (H <= 0.5, |D| <= 7) or paralog
evaluable <- hd$n_het >= 5 & hd$n_geno >= 20
hd$class <- "unevaluable"
hd$class[evaluable] <- ifelse(hd$H[evaluable] <= H_MAX & abs(hd$D[evaluable]) <= D_MAX,
                              "singleton", "paralog")
write.csv(hd, "outputs/ploidy_hdplot_loci.csv", row.names = FALSE)
singleton_loci <- which(hd$class == "singleton")

## Fig. S2: HDplot, D against H for the evaluable loci
hdplot_fig <- function() {
  ev <- which(evaluable)
  plot(hd$H[ev], hd$D[ev], pch = 16, cex = 0.5,
       col = ifelse(hd$class[ev] == "singleton", "#2c7fb888", "#d9544288"),
       xlab = "H (proportion heterozygous)", ylab = "D (read-ratio z-score)")
  abline(h = c(-D_MAX, D_MAX), v = H_MAX, lty = 2, col = "grey40")
}
pdf("figures/ploidy_hdplot.pdf", width = 10, height = 6); hdplot_fig(); dev.off()
png("figures/ploidy_hdplot.png", width = 10, height = 6, units = "in", res = 300); hdplot_fig(); dev.off()

## Results (S1.3): diploid convergence on all loci vs HDplot-singleton loci
conv <- function(loc_idx, label) {
  out <- list()
  for (TH in c(20, 40, 60)) {
    sub_tot <- tot[loc_idx, , drop = FALSE]
    sub_a1 <- a1[loc_idx, , drop = FALSE]
    sub_a2 <- a2[loc_idx, , drop = FALSE]
    sub_rat <- sub_a2 / sub_tot
    keep <- !is.na(sub_tot) & sub_tot >= TH & sub_a1 > 0 & sub_a2 > 0
    idx <- which(keep, arr.ind = TRUE)
    long <- data.frame(id = ids[idx[, 2]], ratio = sub_rat[keep], depth = sub_tot[keep])
    long <- merge(long, truth, by = "id")
    dip <- long[long$ploidy_fc == "Diploid" & long$ratio > 0.15 & long$ratio < 0.85, ]
    band <- dip %>% group_by(id) %>%
      summarise(n = n(), fh = mean(ratio >= 0.42 & ratio <= 0.58), .groups = "drop")
    out[[as.character(TH)]] <- data.frame(
      set = label, threshold = TH, n_ind = nrow(band),
      med_n_het = median(band$n), med_frac_half = median(band$fh),
      binom_expected = exp_frac_half(dip$depth))
  }
  do.call(rbind, out)
}
res <- rbind(conv(seq_len(nloc), "all_loci"), conv(singleton_loci, "singleton_loci"))
res <- res[order(res$threshold, res$set), ]
write.csv(res, "outputs/ploidy_hdplot_convergence.csv", row.names = FALSE)

# pooled verified-diploid allele ratios at depth >= 40, all loci vs singleton loci
pooled <- function(loc_idx, TH = 40) {
  st <- tot[loc_idx, , drop = FALSE]; sa1 <- a1[loc_idx, , drop = FALSE]
  sa2 <- a2[loc_idx, , drop = FALSE]; sr <- sa2 / st
  keep <- !is.na(st) & st >= TH & sa1 > 0 & sa2 > 0
  idx <- which(keep, arr.ind = TRUE)
  long <- merge(data.frame(id = ids[idx[, 2]], ratio = sr[keep]), truth, by = "id")
  long$ratio[long$ploidy_fc == "Diploid"]
}
pdf("figures/ploidy_hdplot_convergence.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
for (cfg in list(c("all loci", "all"), c("singleton loci", "sing"))) {
  r <- if (cfg[2] == "all") pooled(seq_len(nloc)) else pooled(singleton_loci)
  hist(r, breaks = 50, freq = FALSE,
       main = sprintf("Verified diploids, depth>=40\n%s (n=%d)", cfg[1], length(r)),
       xlab = "alt-allele read ratio", col = "grey85", border = "grey60")
  abline(v = c(1/3, 1/2, 2/3), col = c("red", "blue", "red"), lty = 2, lwd = 1.3)
}
dev.off()

## ---- Allele depths re-called from the alignments; BAM files are not deposited ----
# one bowtie2 alignment (bam_files/<TargetID>.FASTQ.gz_sort.bam, @RG SM = sample ID) per verified individual sequenced once
hdr <- read.csv(COUNT_FILES[1], header = FALSE, nrows = 6, check.names = FALSE, stringsAsFactors = FALSE)
colskip <- sum(hdr[1, ] == "*") + 1
crosswalk <- data.frame(id = as.character(unlist(hdr[4, colskip:ncol(hdr)])),
                        TargetID = as.character(unlist(hdr[6, colskip:ncol(hdr)])),
                        stringsAsFactors = FALSE)
crosswalk <- crosswalk[crosswalk$id %in% truth$id, ]
crosswalk <- crosswalk[!crosswalk$id %in% crosswalk$id[duplicated(crosswalk$id)], ]
writeLines(file.path("bam_files", paste0(crosswalk$TargetID, ".FASTQ.gz_sort.bam")),
           "outputs/bam_ploidy/bamlist.txt")

# reference/ref_bamnamed.fasta: reference FASTA with headers renamed to the alignment contig names (CM024753.1 -> CM024753.1_chromosome_1), indexed with samtools faidx
# joint genotyping over the 14 chromosomes with allele depths, MAPQ >= 20, base quality >= 20; biallelic SNPs kept
chroms <- paste0("CM0247", 53:66, ".1_chromosome_", 1:14)
system2("bcftools", c("mpileup", "-f", "reference/ref_bamnamed.fasta", "-a", "FORMAT/AD,FORMAT/DP",
                      "-b", "outputs/bam_ploidy/bamlist.txt", "-r", paste(chroms, collapse = ","),
                      "-q", "20", "-Q", "20", "--threads", "4", "-Ou", "|",
                      "bcftools", "call", "-m", "-v", "--threads", "4", "-Oz",
                      "-o", "outputs/bam_ploidy/chroms.vcf.gz"))
system2("bcftools", c("view", "-m2", "-M2", "-v", "snps", "-Oz",
                      "-o", "outputs/bam_ploidy/chroms.biallelic.vcf.gz", "outputs/bam_ploidy/chroms.vcf.gz"))
system2("bcftools", c("index", "outputs/bam_ploidy/chroms.biallelic.vcf.gz"))

## Results (S1.3): biallelic SNP count and AD-based allele ratios of the verified alignments
v <- read.vcfR("outputs/bam_ploidy/chroms.biallelic.vcf.gz", verbose = FALSE)
ad <- extract.gt(v, element = "AD")
samp <- colnames(ad)
n_bam_snps <- nrow(ad)
vals <- vapply(strsplit(ad, ",", fixed = TRUE),
               function(p) if (length(p) < 2) c(NA, NA) else suppressWarnings(as.numeric(p[1:2])),
               numeric(2))
refd <- matrix(vals[1, ], nrow(ad), dimnames = dimnames(ad))
altd <- matrix(vals[2, ], nrow(ad), dimnames = dimnames(ad))
tot_ad <- refd + altd
ratio_ad <- altd / tot_ad

## Results (S1.3): median BAM allele depth per genotype
median_depth_bam <- median(tot_ad[!is.na(tot_ad) & tot_ad > 0])

# heterozygous genotypes at a depth threshold, joined to flow cytometry
make_long <- function(TH) {
  keep <- !is.na(tot_ad) & tot_ad >= TH & refd > 0 & altd > 0
  idx <- which(keep, arr.ind = TRUE)
  long <- data.frame(id = samp[idx[, 2]], ratio = ratio_ad[keep], depth = tot_ad[keep])
  merge(long, truth, by = "id")
}

## Results (S1.3): half-band fraction vs depth threshold on the BAM allele depths
summ <- list()
for (TH in c(20, 40, 60)) {
  long <- make_long(TH)
  long <- long[long$ploidy_fc %in% c("Diploid", "Triploid"), ]
  band <- long %>% filter(ratio > 0.15, ratio < 0.85) %>%
    group_by(id, ploidy_fc) %>%
    summarise(n_het = n(), frac_half = mean(ratio >= 0.42 & ratio <= 0.58), .groups = "drop")
  grp <- band %>% group_by(ploidy_fc) %>%
    summarise(n_ind = n(), med_n_het = median(n_het), med_frac_half = median(frac_half), .groups = "drop")
  dipd <- long$depth[long$ploidy_fc == "Diploid" & long$ratio > 0.15 & long$ratio < 0.85]
  grp$binom_expected <- ifelse(grp$ploidy_fc == "Diploid", exp_frac_half(dipd), NA)
  grp$threshold <- TH
  summ[[as.character(TH)]] <- grp
}
res <- do.call(rbind, summ)
write.csv(res, "outputs/bam_ploidy/bam_ad_convergence.csv", row.names = FALSE)

# pooled allele-ratio histograms per cytotype and depth threshold
pdf("figures/ploidy_bam_ad_check.pdf", width = 11, height = 7)
par(mfrow = c(2, 3))
for (g in c("Diploid", "Triploid")) {
  for (TH in c(20, 40, 60)) {
    long <- make_long(TH); long <- long[long$ploidy_fc == g, ]
    r <- long$ratio
    hist(r, breaks = 50, freq = FALSE,
         main = paste0(g, "  BAM-AD depth>=", TH, "\nn=", length(r)),
         xlab = "alt-allele read ratio", col = "grey85", border = "grey60", cex.main = 0.9)
    abline(v = c(1/3, 1/2, 2/3), col = c("red", "blue", "red"), lty = 2, lwd = 1.3)
  }
}
dev.off()

# per-individual histograms at depth >= 20 (up to 16 diploids and 16 triploids)
pdf("figures/ploidy_bam_ad_individuals.pdf", width = 11, height = 8)
long20 <- make_long(20)
plot_inds <- function(group, max_n) {
  these <- head(unique(long20$id[long20$ploidy_fc == group]), max_n)
  if (!length(these)) return(invisible())
  par(mfrow = c(4, 4), mar = c(2.3, 2.3, 1.8, 0.8))
  for (s in these) {
    r <- long20$ratio[long20$id == s]
    hist(r, breaks = 40, freq = FALSE, main = paste0(s, " (", group, ", n=", length(r), ")"),
         xlab = "", ylab = "", col = "grey85", border = "grey60", cex.main = 0.8)
    abline(v = c(1/3, 1/2, 2/3), col = c("red", "blue", "red"), lty = 2, lwd = 1.1)
  }
}
plot_inds("Diploid", 16)
plot_inds("Triploid", 16)
dev.off()

# ============================================================================
# 4. Diploid working set (SI S2; Figs S3, S4; Supplementary Data 4)
# ============================================================================

## Fig. S3 / Results: individual heterozygosity separates diploids from polyploids
# flow-cytometry cytotypes (outputs/ploidy.csv: id, "Flow cytometry verified")
g <- readRDS("input_files/taro.rds")
truth <- read.csv("outputs/ploidy.csv", stringsAsFactors = FALSE)
colnames(truth) <- c("id", "ploidy_fc"); truth$ploidy_fc <- trimws(truth$ploidy_fc)
truth$class <- ifelse(truth$ploidy_fc == "Diploid", "Diploid", "Polyploid")

# per-individual observed heterozygosity and call rate from the 0/1/2 calls
M <- as.matrix(g)
nonmiss <- rowSums(!is.na(M))
Ho <- rowSums(M == 1, na.rm = TRUE) / nonmiss
ind <- data.frame(id = indNames(g), Ho = Ho, n_called = nonmiss,
                  call_rate = nonmiss / ncol(M), stringsAsFactors = FALSE)
ind <- merge(ind, truth[, c("id", "ploidy_fc", "class")], by = "id", all.x = TRUE)
ind$labelled <- !is.na(ind$class)
write.csv(ind, "outputs/ploidy_individual_Ho.csv", row.names = FALSE)

## Results: median Ho by cytotype (0.168 polyploid vs 0.033 diploid), AUC 0.997
lab <- ind[ind$labelled, ]
dip <- lab$Ho[lab$class == "Diploid"]
pol <- lab$Ho[lab$class == "Polyploid"]
n_dip <- length(dip); n_poly <- length(pol)
ho_by_class <- aggregate(Ho ~ class, lab, function(x)
  round(c(n = length(x), median = median(x), min = min(x), max = max(x), mean = mean(x)), 4))
auc <- mean(outer(pol, dip, ">") + 0.5 * outer(pol, dip, "=="))

# gate A (max diploid Ho), gate B (min polyploid Ho) and the Youden-optimal threshold (0.0923)
thrA <- max(dip)
thrB <- min(pol)
cand <- sort(unique(lab$Ho))
youden <- sapply(cand, function(t) mean(pol > t) + mean(dip <= t) - 1)
tstar <- cand[which.max(youden)]

# Fig. S3: labelled Ho by class (upper), all individuals with labelled rugs (lower)
draw_ho_gate <- function() {
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  br <- seq(0, max(lab$Ho) * 1.05, length.out = 40)
  hd <- hist(dip, breaks = br, plot = FALSE)
  hp <- hist(pol, breaks = br, plot = FALSE)
  plot(hd, col = rgb(0.2, 0.5, 0.8, 0.6), border = "white", xlim = range(br),
       ylim = c(0, max(hd$counts, hp$counts)), main = "",
       xlab = "Individual observed heterozygosity (Ho)")
  plot(hp, col = rgb(0.85, 0.3, 0.2, 0.6), border = "white", add = TRUE)
  abline(v = c(thrA, thrB), lty = c(2, 3), col = c("blue", "red"))
  legend("topright", c("Diploid", "Polyploid", "max-dip threshold", "min-poly threshold"),
         fill = c(rgb(0.2, 0.5, 0.8, 0.6), rgb(0.85, 0.3, 0.2, 0.6), NA, NA),
         border = NA, lty = c(NA, NA, 2, 3), col = c(NA, NA, "blue", "red"), bty = "n")
  hist(ind$Ho, breaks = 60, col = "grey85", border = "grey60", main = "",
       xlab = "Individual observed heterozygosity (Ho)")
  rug(dip, col = rgb(0.2, 0.5, 0.8, 0.9), lwd = 2)
  rug(pol, col = rgb(0.85, 0.3, 0.2, 0.9), lwd = 2, ticksize = 0.06)
  abline(v = c(thrA, thrB), lty = c(2, 3), col = c("blue", "red"))
}
pdf("figures/ploidy_diploid_filter_Ho.pdf", width = 10, height = 7); draw_ho_gate(); dev.off()
png("figures/ploidy_diploid_filter_Ho.png", width = 10, height = 7, units = "in", res = 200)
draw_ho_gate(); dev.off()

## SI S2.1: AUC bootstrap 95% CI (0.984-1.000); sensitivity 11/11 and specificity 51/52 at Ho 0.0923 with exact CIs
THRESH <- 0.0923
mixo <- lab$Ho[lab$ploidy_fc == "Mixoploid"]
trip <- lab$Ho[lab$ploidy_fc == "Triploid"]
auc_fn <- function(a, b) {
  gt <- outer(b, a, ">"); eq <- outer(b, a, "==")
  (sum(gt) + 0.5 * sum(eq)) / (length(a) * length(b))
}
set.seed(42)
B <- 5000
boot <- replicate(B, auc_fn(sample(dip, replace = TRUE), sample(pol, replace = TRUE)))
auc_ci <- quantile(boot, c(0.025, 0.975))
sens <- mean(pol > THRESH); spec <- mean(dip <= THRESH)
k_sens <- sum(pol > THRESH); k_spec <- sum(dip <= THRESH)
sens_ci <- binom.test(k_sens, n_poly)$conf.int
spec_ci <- binom.test(k_spec, n_dip)$conf.int
write.csv(data.frame(
  metric = c("n_diploid", "n_polyploid", "n_triploid", "n_mixoploid",
             "AUC", "AUC_CI_low", "AUC_CI_high",
             "sensitivity", "sens_CI_low", "sens_CI_high",
             "specificity", "spec_CI_low", "spec_CI_high", "youden_threshold"),
  value = round(c(n_dip, n_poly, length(trip), length(mixo),
                  auc, auc_ci[1], auc_ci[2],
                  sens, sens_ci[1], sens_ci[2],
                  spec, spec_ci[1], spec_ci[2], tstar), 4)),
  "outputs/ploidy_gate_ci.csv", row.names = FALSE)

## SI S2.1: contamination and hybridisation cross-check of the Ho screen
ind <- ind[match(indNames(g), ind$id), ]
ind$cat <- with(ind, ifelse(!is.na(class) & class == "Diploid", "verif-diploid",
                     ifelse(!is.na(class) & class == "Polyploid", "verif-polyploid",
                     ifelse(Ho >= THRESH, "unlab-highHo", "unlab-lowHo"))))
# (A) contamination: Ho against call rate (median call rate 0.92 in the high-Ho group)
cc <- cor(ind$Ho, ind$call_rate, use = "complete.obs")
callrate_by_cat <- round(tapply(ind$call_rate, ind$cat, median), 3)
# (B) hybridisation: PCA on mean-imputed calls; high-Ho individuals sit with the verified polyploids
Mi <- M
cm <- colMeans(Mi, na.rm = TRUE)
ix <- which(is.na(Mi), arr.ind = TRUE)
Mi[ix] <- cm[ix[, 2]]
Mi <- Mi[, apply(Mi, 2, function(x) var(x) > 0)]
pc <- prcomp(Mi, center = TRUE, scale. = FALSE)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
ind$PC1 <- pc$x[, 1]; ind$PC2 <- pc$x[, 2]; ind$PC3 <- pc$x[, 3]
ho_pc_cor <- sapply(1:3, function(k) cor(ind$Ho, pc$x[, k]))
grp <- ind$cat %in% c("verif-polyploid", "unlab-highHo")
pc1_median <- c(highHo_poly = median(ind$PC1[grp]), diploid_lowHo = median(ind$PC1[!grp]))

cols <- c("verif-diploid" = "#1f78b4", "verif-polyploid" = "#e31a1c",
          "unlab-lowHo" = "#a6cee3", "unlab-highHo" = "#fb9a99")
pchs <- c("verif-diploid" = 17, "verif-polyploid" = 17, "unlab-lowHo" = 16, "unlab-highHo" = 16)
pdf("figures/ploidy_diploid_crosscheck.pdf", width = 11, height = 5.5)
par(mfrow = c(1, 2))
ord <- order(ind$cat %in% c("verif-diploid", "verif-polyploid"))
plot(ind$PC1[ord], ind$PC2[ord], col = cols[ind$cat[ord]], pch = pchs[ind$cat[ord]],
     cex = ifelse(grepl("verif", ind$cat[ord]), 1.2, 0.6),
     xlab = sprintf("PC1 (%.1f%%)", ve[1]), ylab = sprintf("PC2 (%.1f%%)", ve[2]),
     main = paste0("PCA — all ", nrow(ind), " individuals"))
legend("topright", names(cols), col = cols, pch = pchs[names(cols)], pt.cex = 1, bty = "n", cex = 0.8)
plot(ind$call_rate, ind$Ho, col = cols[ind$cat], pch = pchs[ind$cat],
     cex = ifelse(grepl("verif", ind$cat), 1.2, 0.6),
     xlab = "Call rate (1 - missingness)", ylab = "Individual Ho",
     main = sprintf("Contamination check (r=%.2f)", cc))
abline(h = THRESH, lty = 2, col = "grey40")
dev.off()
write.csv(ind[, c("id", "Ho", "call_rate", "n_called", "ploidy_fc", "class", "cat", "PC1", "PC2", "PC3")],
          "outputs/ploidy_crosscheck.csv", row.names = FALSE)

## Supplementary Data 4: per-accession diploid decision (flow cytometry where available, Ho < 0.0923 otherwise)
ind$decision <- with(ind, ifelse(
  !is.na(class) & class == "Diploid", "keep-verified-diploid",
  ifelse(!is.na(class) & class == "Polyploid", "drop-verified-polyploid",
  ifelse(Ho < THRESH, "keep-Ho-diploid", "drop-Ho-polyploid"))))
ind$is_diploid <- ind$decision %in% c("keep-verified-diploid", "keep-Ho-diploid")
write.csv(ind[, c("id", "Ho", "call_rate", "ploidy_fc", "class", "decision", "is_diploid")],
          "outputs/diploid_set_decisions.csv", row.names = FALSE)

## Results: 1,064 of 1,212 accessions retained as the diploid working set
gd <- gl.keep.ind(g, ind.list = ind$id[ind$is_diploid], verbose = 0)
gd <- gl.recalc.metrics(gd, verbose = 0)
saveRDS(gd, "input_files/taro_diploid.rds")

## Fig. S4: PCA of all accessions before the diploid filter, coloured by ploidy
set.seed(42)
dec <- read.csv("outputs/diploid_set_decisions.csv", stringsAsFactors = FALSE)
fc <- trimws(as.character(dec$ploidy_fc))
lab_ploidy <- ifelse(!is.na(fc) & fc != "" & fc != "NA", fc,
                     ifelse(dec$is_diploid, "Diploid (Ho)", "Polyploid (Ho)"))
names(lab_ploidy) <- dec$id
glp <- g
pop(glp) <- factor(lab_ploidy[indNames(glp)])
# impute on the original country populations, then PCA
t1 <- gl.gen2fbm(g, verbose = 0)
t2 <- gl.impute(t1, verbose = 0)
pca <- gl.pcoa(t2, verbose = 0)
saveRDS(pca, "outputs/pca_preploidy.rds")
p <- gl.pcoa.plot(pca, glp, pt.size = 3, verbose = 0)
ggplot2::ggsave("figures/pca_preploidy_by_ploidy.png", p, width = 8, height = 6, dpi = 200)
ggplot2::ggsave("figures/pca_preploidy_by_ploidy.pdf", p, width = 8, height = 6)

## SI S2.1: leave-one-out cross-validation of the Youden threshold (11/11 polyploids, 51/52 diploids)
youden_thr <- function(a, b) {
  cands <- sort(unique(c(a, b)))
  grid <- c(cands[1] - 1e-6, (head(cands, -1) + tail(cands, -1)) / 2, tail(cands, 1) + 1e-6)
  J <- sapply(grid, function(t) mean(b > t) + mean(a <= t) - 1)
  grid[which.max(J)]
}
labeled <- rbind(data.frame(Ho = dip, y = "Diploid"), data.frame(Ho = pol, y = "Polyploid"))
pred <- character(nrow(labeled))
for (i in seq_len(nrow(labeled))) {
  a <- labeled$Ho[labeled$y == "Diploid" & seq_len(nrow(labeled)) != i]
  b <- labeled$Ho[labeled$y == "Polyploid" & seq_len(nrow(labeled)) != i]
  pred[i] <- ifelse(labeled$Ho[i] > youden_thr(a, b), "Polyploid", "Diploid")
}
cm_cv <- table(true = labeled$y, pred = factor(pred, levels = c("Diploid", "Polyploid")))
cv_sens <- cm_cv["Polyploid", "Polyploid"] / n_poly
cv_spec <- cm_cv["Diploid", "Diploid"] / n_dip
# separation gap: most heterozygous verified diploid (0.133) vs least heterozygous verified polyploid (0.120)
gap_lo <- max(dip); gap_hi <- min(pol)

## Results: 131 (12%) working-set accessions within 0.02 Ho below the cutoff, mostly Indonesia, Vietnam, Malaysia
ind$called_dip <- ifelse(!is.na(ind$class), ind$class == "Diploid", ind$Ho <= THRESH)
ws <- ind[ind$called_dip, ]
ctry <- setNames(as.character(pop(gd)), indNames(gd))
near_cut_country <- sort(table(ctry[intersect(ws$id[ws$Ho > THRESH - 0.02], names(ctry))]), decreasing = TRUE)

# ============================================================================
# 5. Qploidy (SI S1.4)
# ============================================================================

## ---- Qploidy (Taniguti et al. 2025): standardised B-allele frequency from the DArT counts ----

# write the DArT read counts of the filtered genlight's loci into the AD field of a GT:AD VCF (gl2vcf writes GT only)
gl <- readRDS("input_files/taro.rds")
cd <- read.dart.counts(COUNT_FILES[1], ind.vec = indNames(gl), id = 4)
cm <- as.matrix(cd$counts)
cm <- matrix(as.integer(as.character(cm)), nrow = nrow(cm), dimnames = dimnames(cm))
odd <- seq(1, nrow(cm), 2); even <- seq(2, nrow(cm), 2)
ref_mat <- cm[odd, , drop = FALSE]
alt_mat <- cm[even, , drop = FALSE]
key <- as.character(cd$loc_metrics$AlleleID[even])
rownames(ref_mat) <- rownames(alt_mat) <- key

# match loci by AlleleID; CHROM/POS from the genlight slots, REF/ALT from the SNP field ("pos:REF>ALT")
glm <- gl$other$loc.metrics
gkey <- as.character(glm$AlleleID)
common_loci <- intersect(gkey, key)
samps <- intersect(indNames(gl), colnames(cm))
gi <- match(common_loci, gkey)
chr <- as.character(gl$chromosome[gi])
pos <- as.integer(gl$position[gi])
ra <- sub(".*:", "", sub("/", ">", as.character(glm$SNP[gi])))
REF <- sub(">.*", "", ra); ALT <- sub(".*>", "", ra)
nt <- c("A", "C", "G", "T")
ok <- !is.na(chr) & chr != "" & !is.na(pos) & REF %in% nt & ALT %in% nt
ok[is.na(ok)] <- FALSE
ord <- order(as.integer(chr[ok]), pos[ok])
common_loci <- common_loci[ok][ord]; gi <- gi[ok][ord]
chr <- chr[ok][ord]; pos <- pos[ok][ord]; REF <- REF[ok][ord]; ALT <- ALT[ok][ord]

# GT from the genlight calls, AD from the counts
GT <- t(as.matrix(gl)[, gi, drop = FALSE])[, samps, drop = FALSE]
gt <- matrix("./.", nrow(GT), ncol(GT))
gt[GT == 0] <- "0/0"; gt[GT == 1] <- "0/1"; gt[GT == 2] <- "1/1"
rs <- ref_mat[common_loci, samps, drop = FALSE]
as_ <- alt_mat[common_loci, samps, drop = FALSE]
rs[is.na(rs)] <- 0L; as_[is.na(as_)] <- 0L
cell <- matrix(paste0(gt, ":", rs, ",", as_), nrow(gt))
body <- paste(chr, pos, common_loci, REF, ALT, ".", "PASS", ".", "GT:AD",
              apply(cell, 1, paste, collapse = "\t"), sep = "\t")
header <- c("##fileformat=VCFv4.2", "##source=qploidy_from_dart_counts.R",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">',
            paste0("##contig=<ID=", unique(chr), ">"),
            paste0("#", paste(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samps),
                              collapse = "\t")))
writeLines(c(header, body), "outputs/taro_qploidy_all.vcf")

# standardise BAF per marker against a diploid reference (three dosage clusters, count mode)
qp_data <- qploidy_read_vcf("outputs/taro_qploidy_all.vcf")
qp_genos <- qploidy_read_vcf("outputs/taro_qploidy_all.vcf", geno = TRUE)
qp_pos <- qploidy_read_vcf("outputs/taro_qploidy_all.vcf", geno.pos = TRUE)
std <- standardize(data = qp_data, genos = qp_genos, geno.pos = qp_pos,
                   ploidy.standardization = 2, threshold.n.clusters = 3,
                   type = "counts", n.cores = 4, verbose = FALSE)
saveRDS(std, "outputs/qploidy_std_all.rds")

## Results (S1.4): Qploidy copy number (grid 2-4) for the flow cytometry set vs flow cytometry
dec <- read.csv("outputs/diploid_set_decisions.csv", stringsAsFactors = FALSE)
ver <- dec[!is.na(dec$ploidy_fc) & dec$ploidy_fc != "" & dec$ploidy_fc != "NA" &
           dec$id %in% unique(std$data$SampleName), ]
calls <- do.call(rbind, lapply(ver$id, function(s) {
  b <- select_best_baf_model(qploidy_standardization = std, sample = s, cn_grid = 2:4, plot = FALSE)$best
  data.frame(id = s, qploidy_cn = b$best_cn, n_obs = b$n_obs, prob = round(b$prob, 3),
             stringsAsFactors = FALSE)
}))
ver <- merge(ver[, c("id", "ploidy_fc", "class", "Ho", "call_rate")], calls, by = "id")
write.csv(ver, "outputs/qploidy_verified_calls.csv", row.names = FALSE)
ver$qp_class <- ifelse(ver$qploidy_cn == 2, "Diploid", "Polyploid")
qploidy_tab <- table(flow_cytometry = ver$ploidy_fc, qploidy_cn = ver$qploidy_cn)
qploidy_class_tab <- table(flow_cytometry = ver$class, qploidy = ver$qp_class)
write.csv(as.data.frame.matrix(qploidy_tab), "outputs/qploidy_vs_flow_cytometry.csv")

## Figure: standardised BAF of two verified diploids, triploids and mixoploids (figures/qploidy_baf_examples.png)
pick <- function(cls, n = 2) head(ver$id[ver$ploidy_fc == cls], n)
ex <- data.frame(id = c(pick("Diploid"), pick("Triploid"), pick("Mixoploid")))
ex$ploidy_fc <- ver$ploidy_fc[match(ex$id, ver$id)]
ex$lab <- paste0(ex$id, " (", ex$ploidy_fc, ")")
bd <- std$data[std$data$SampleName %in% ex$id, c("SampleName", "baf")]
names(bd) <- c("id", "baf")
bd <- bd[!is.na(bd$baf), ]
bd$lab <- ex$lab[match(bd$id, ex$id)]
p <- ggplot(bd, aes(baf)) +
  geom_histogram(bins = 40, fill = "#2E7D32", colour = "white", linewidth = .1) +
  geom_vline(xintercept = c(1/3, 1/2, 2/3), linetype = 2, colour = "red", linewidth = .3) +
  facet_wrap(~ lab, scales = "free_y", ncol = 2) +
  labs(x = "Standardized B-allele frequency (dashed: 1/3, 1/2, 2/3)", y = "Markers",
       title = "Qploidy standardized BAF: verified diploids vs polyploids (taro)") +
  theme_bw(base_size = 10)
ggsave("figures/qploidy_baf_examples.png", p, width = 9, height = 7, dpi = 200, bg = "white")

# ============================================================================
# 6. Clonal lineages, cultivar synonymy and the candidate core
# ============================================================================

# ---------------------------------------------------------------------------
## Results: 1,064 diploid accessions collapse to 392 clonal lineages (Arnaud-Haond et al. 2007)
# ---------------------------------------------------------------------------
gl      <- readRDS("input_files/taro_diploid.rds")
country <- as.character(pop(gl))
sc      <- as.snpclone(gl)

# pairwise distance = fraction of differing calls; complete linkage = farthest neighbour
D  <- bitwise.dist(sc, percent = TRUE, missing_match = TRUE, threads = 1L)
saveRDS(D, "outputs/mlg_bitwise_dist.rds")
Dv <- as.numeric(D)
hc <- hclust(D, method = "complete")

## Methods: cutoff_predictor (largest gap in the filter_stats merge thresholds) supports 0.008
fstats   <- filter_stats(sc, distance = D, plot = FALSE, threads = 1L)
auto_thr <- vapply(c("farthest", "average", "nearest"),
                   function(a) cutoff_predictor(fstats[[a]]$THRESHOLDS), numeric(1))
THRESH <- 0.008

# pairwise-distance distribution (bimodal; lineages cut at the valley)
zoom_hi <- max(0.05, THRESH * 4)
p_full <- ggplot(data.frame(d = Dv), aes(d)) +
  geom_histogram(bins = 100, fill = "#2E7D32", colour = "white", linewidth = .1) +
  geom_vline(xintercept = THRESH, colour = "red", linetype = 2) +
  labs(x = "Pairwise distance (fraction of differing loci)", y = "Pairs",
       title = "Taro diploid set: pairwise genetic distance",
       subtitle = sprintf("red = threshold %.4f", THRESH)) +
  theme_bw(base_size = 11)
p_zoom <- ggplot(data.frame(d = Dv[Dv <= zoom_hi]), aes(d)) +
  geom_histogram(bins = 80, fill = "#1565C0", colour = "white", linewidth = .1) +
  geom_vline(xintercept = THRESH, colour = "red", linetype = 2) +
  labs(x = "Pairwise distance (zoom on low-distance / clone peak)", y = "Pairs",
       title = sprintf("Zoom 0 - %.3f", zoom_hi)) +
  theme_bw(base_size = 11)
ggsave("figures/mlg_distance_distribution.png",
       gridExtra::arrangeGrob(p_full, p_zoom, ncol = 1),
       width = 8, height = 7, dpi = 200, bg = "white")

# exact multilocus genotypes (distance 0) and lineages (complete linkage cut at 0.008)
mlg_exact <- as.integer(cutree(hc, h = 0))
n_exact   <- length(unique(mlg_exact))
mll_id    <- as.integer(cutree(hc, h = THRESH))
n_mll     <- length(unique(mll_id))

# one representative per lineage: the ramet with the highest call rate
M        <- as.matrix(gl)
callrate <- 1 - rowMeans(is.na(M))
ord    <- order(mll_id, -callrate)
reps   <- rownames(M)[ord][!duplicated(mll_id[ord])]
is_rep <- indNames(gl) %in% reps

# clone-censored allele frequencies (one ramet per lineage) and F_IS per locus
Mr  <- M[reps, , drop = FALSE]
nrp <- length(reps)
q   <- colMeans(Mr, na.rm = TRUE) / 2
fl  <- 1 / (2 * nrp)
q   <- pmin(pmax(q, fl), 1 - fl)
Ho  <- colMeans(Mr == 1, na.rm = TRUE)
He  <- 2 * q * (1 - q)
FIS <- ifelse(He > 0, 1 - Ho / He, 0)
FIS[!is.finite(FIS)] <- 0
FIS <- pmin(pmax(FIS, -0.999), 0.999)

lq <- log10(q); l1mq <- log10(1 - q); l2 <- log10(2)
N  <- nrow(M); L <- ncol(M)
idx1 <- which(M == 1L); idx2 <- which(M == 2L); idxNA <- which(is.na(M))
col1 <- ((idx1 - 1L) %/% N) + 1L
col2 <- ((idx2 - 1L) %/% N) + 1L

# p_gen under HWE: product over loci (2008 erratum), computed as a sum of log10
G <- matrix(2 * l1mq, nrow = N, ncol = L, byrow = TRUE)
G[idx1] <- (l2 + lq + l1mq)[col1]
G[idx2] <- (2 * lq)[col2]
G[idxNA] <- 0
log10_pgen_hwe <- rowSums(G)
l_scored <- rowSums(!is.na(M))

# p_gen with the F_IS correction: homozygotes x (1 + FIS), heterozygotes x (1 - FIS)
adj_hom <- log10(1 + FIS); adj_het <- log10(1 - FIS)
Gf <- matrix(2 * l1mq + adj_hom, nrow = N, ncol = L, byrow = TRUE)
Gf[idx1] <- (l2 + lq + l1mq + adj_het)[col1]
Gf[idx2] <- (2 * lq + adj_hom)[col2]
Gf[idxNA] <- 0
log10_pgen_fis <- rowSums(Gf)

## Results: p_sex, the probability that a genotype recurs n times by sexual reproduction
mll_n <- as.integer(table(mll_id)[as.character(mll_id)])
psex_log10 <- function(log10p, n, Ntot) {
  p <- 10^log10p
  log10_1mp <- if (p <= 0 || !is.finite(p)) 0 else log10(1 - p)
  is <- n:min(Ntot, n + 50L)
  lt <- lchoose(Ntot, is) / log(10) + is * log10p + (Ntot - is) * log10_1mp
  m  <- max(lt)
  m + log10(sum(10^(lt - m)))
}
psex_hwe <- mapply(psex_log10, log10_pgen_hwe, mll_n, MoreArgs = list(Ntot = N))
psex_fis <- mapply(psex_log10, log10_pgen_fis, mll_n, MoreArgs = list(Ntot = N))

## Supplementary Data 5: lineage membership of every diploid accession
im <- gl$other$ind.metrics
membership <- data.frame(
  id              = indNames(gl),
  country         = country,
  variety         = as.character(im$Variety),
  common_name     = as.character(im$Common_name),
  mlg_exact       = mlg_exact,
  mll             = mll_id,
  mll_size        = mll_n,
  loci_scored     = l_scored,
  call_rate       = round(callrate, 4),
  log10_pgen_hwe  = round(log10_pgen_hwe, 2),
  log10_psex_hwe  = round(psex_hwe, 2),
  log10_pgen_fis  = round(log10_pgen_fis, 2),
  log10_psex_fis  = round(psex_fis, 2),
  is_representative = is_rep,
  stringsAsFactors = FALSE)
write.csv(membership, "outputs/mlg_membership.csv", row.names = FALSE)

## Results: clonal richness R = 0.37, Simpson D*, evenness V and Pareto beta = 0.44, overall and per country
pareto_beta <- function(sizes) {
  tx  <- sort(unique(sizes[sizes > 0]))
  Nge <- vapply(tx, function(x) sum(sizes[sizes >= x]), numeric(1))
  d   <- data.frame(X = tx, Nge = Nge)
  d   <- d[d$X > 0 & d$Nge > 0, , drop = FALSE]
  if (nrow(d) < 3 || length(unique(d$X)) < 3)
    return(list(beta = NA, r2 = NA, p = NA, npts = nrow(d), data = d))
  fit <- lm(log10(Nge) ~ log10(X), data = d); s <- summary(fit)
  list(beta = unname(-coef(fit)[2]), r2 = s$r.squared,
       p = coef(s)[2, 4], npts = nrow(d), data = d)
}
clonal_indices <- function(mll_sub) {
  Nn <- length(mll_sub)
  ni <- as.integer(table(mll_sub))
  G  <- length(ni)
  pi <- ni / Nn
  Pd <- G / Nn
  R  <- if (Nn > 1) (G - 1) / (Nn - 1) else NA
  Go <- 1 / sum(pi^2)
  Dstar <- if (Nn > 1) 1 - sum(ni * (ni - 1)) / (Nn * (Nn - 1)) else NA
  Hp <- -sum(pi * log(pi))
  Jp <- if (G > 1) Hp / log(G) else NA
  if (Nn > 1 && G > 1) {
    Dmin <- ((2 * Nn - G) * (G - 1) / Nn^2) * (Nn / (Nn - 1))
    Dmax <- ((G - 1) / G) * (Nn / (Nn - 1))
    V <- if (Dmax > Dmin) (Dstar - Dmin) / (Dmax - Dmin) else NA
  } else V <- NA
  pb <- pareto_beta(ni)
  data.frame(N = Nn, G_MLL = G, Pd = round(Pd, 4), R = round(R, 4),
             Go = round(Go, 3), Dstar = round(Dstar, 4),
             Hp = round(Hp, 3), J = round(Jp, 4), V = round(V, 4),
             pareto_beta = round(pb$beta, 3), pareto_r2 = round(pb$r2, 3),
             pareto_p = signif(pb$p, 3))
}
idx_overall <- cbind(level = "OVERALL", clonal_indices(mll_id))
idx_by_pop  <- do.call(rbind, lapply(sort(unique(country)), function(cc)
  cbind(level = cc, clonal_indices(mll_id[country == cc]))))
g_exact_pop <- tapply(mlg_exact, country, function(x) length(unique(x)))
idx_by_pop$G_MLG_exact  <- as.integer(g_exact_pop[idx_by_pop$level])
idx_overall$G_MLG_exact <- n_exact
diversity <- rbind(idx_overall, idx_by_pop)
write.csv(diversity, "outputs/mlg_diversity_indices.csv", row.names = FALSE)

## Results: 53 lineages occur in more than one country (518 accessions; one lineage of 135 across 11 countries)
share <- do.call(rbind, lapply(sort(unique(mll_id)), function(m) {
  rows <- which(mll_id == m)
  cty  <- sort(unique(country[rows]))
  if (length(cty) < 2) return(NULL)
  vr <- membership$variety[rows]; vr <- unique(vr[!is.na(vr) & vr != ""])
  data.frame(mll = m, n_ramets = length(rows), n_countries = length(cty),
             countries = paste(cty, collapse = "; "),
             varieties = paste(utils::head(vr, 8), collapse = "; "),
             stringsAsFactors = FALSE)
}))
share <- share[order(-share$n_countries, -share$n_ramets), ]
write.csv(share, "outputs/mlg_crosscountry_sharing.csv", row.names = FALSE)

# clone-censored genlight: 392 genets, one highest-call-rate ramet per lineage
gl_mll <- gl.keep.ind(gl, ind.list = reps, verbose = 0)
gl_mll <- gl.recalc.metrics(gl_mll, verbose = 0)
saveRDS(gl_mll, "input_files/taro_diploid_mll.rds")

# ---------------------------------------------------------------------------
## Results, SI S4.1: single, average and complete linkage give 358, 376 and 392 lineages at 0.008
# ---------------------------------------------------------------------------
CUT <- 0.008
n <- attr(D, "Size")
link <- sapply(c(complete = "complete", average = "average", single = "single"),
               function(mth) length(unique(cutree(hclust(D, method = mth), h = CUT))))
red <- round(100 * (1 - link / n), 1)
valley <- sapply(c(0.006, 0.008, 0.010), function(h)
  sapply(c("complete", "average", "single"),
         function(m) length(unique(cutree(hclust(D, method = m), h = h)))))
ids <- indNames(gl)
linkage <- data.frame(
  metric = c("genets_complete", "genets_average", "genets_single",
             "redundancy_complete", "redundancy_average", "redundancy_single",
             "n_accessions", "n_unique_ids"),
  value  = c(link["complete"], link["average"], link["single"],
             red["complete"], red["average"], red["single"],
             length(ids), length(unique(ids))))
write.csv(linkage, "outputs/linkage_sensitivity.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
## SI S4.1, Methods: genotyping error from 56 technical replicate pairs (23 samples run 2-4 times)
# ---------------------------------------------------------------------------
g4     <- readRDS("input_files/taro_report4_raw.rds")
cn     <- canon4(indNames(g4))

# replicate discordance on the report-4 loci that survive into the analysis set
keep <- which(g4@other$loc.metrics$AlleleID %in% gl@other$loc.metrics$AlleleID)
g4f  <- g4[, keep]
M4   <- as.matrix(g4f)
rownames(M4) <- indNames(g4)
sets  <- names(which(table(cn) >= 2))
pairs <- do.call(rbind, lapply(sets, function(s) t(combn(which(cn == s), 2))))
disc <- apply(pairs, 1, function(p) {
  a <- M4[p[1], ]; b <- M4[p[2], ]
  both <- !is.na(a) & !is.na(b)
  n <- sum(both); d <- sum(a[both] != b[both])
  c(n_both = n, n_diff = d, discord = d / n)
})
err <- data.frame(canon = cn[pairs[, 1]],
                  inst_a = rownames(M4)[pairs[, 1]], inst_b = rownames(M4)[pairs[, 2]],
                  t(disc))
write.csv(err, "outputs/somatic_error_floor.csv", row.names = FALSE)
disc_all <- median(err$discord)

## SI S4.1: discordance 2.49% at all analysis loci (RepAvg >= 0.98) and 0.31% at RepAvg = 1.00
repL <- g4f@other$loc.metrics$RepAvg
disc_band <- function(mask) median(apply(pairs, 1, function(p) {
  a <- M4[p[1], mask]; b <- M4[p[2], mask]; both <- !is.na(a) & !is.na(b)
  sum(a[both] != b[both]) / sum(both)
}))
disc_repavg <- sapply(c(0.98, 0.99, 0.995, 1.0), function(th) disc_band(repL >= th))
names(disc_repavg) <- c("0.98", "0.99", "0.995", "1.00")

## SI S4.1: replicate samples have a lower call rate (0.76) than the collection (0.97)
cr4    <- 1 - vapply(g4@gen, function(x) length(x@NA.posi), 1) / nLoc(g4)
cr_rep <- median(cr4[cn %in% sets])
cr_all <- median(callrate)

## Methods: median pairwise distance within lineages (0.0035)
Dm <- as.matrix(D)
within <- unlist(lapply(unique(mll_id[duplicated(mll_id)]), function(m) {
  idx <- which(mll_id == m); sub <- Dm[idx, idx]; sub[lower.tri(sub)]
}))
within_med <- median(within)

# ---------------------------------------------------------------------------
# allele-presence matrix: REF carried if genotype 0/1, ALT if 1/2 (columns 1..L REF, L+1..2L ALT)
# ---------------------------------------------------------------------------
ref <- M <= 1L; ref[is.na(ref)] <- FALSE
alt <- M >= 1L; alt[is.na(alt)] <- FALSE
P   <- cbind(ref, alt); rm(ref, alt)
present_full <- colSums(P) > 0
A_full <- sum(present_full)
Pp <- P[, present_full, drop = FALSE]

# ---------------------------------------------------------------------------
## Fig. S6: lineage count, redundancy and core allele capture across the lineage threshold
# ---------------------------------------------------------------------------
CHOSEN <- 0.008
grid <- sort(unique(c(0.002, 0.004, 0.005, 0.006, 0.007, CHOSEN, 0.009,
                      0.010, 0.012, 0.015, 0.020)))
sens <- do.call(rbind, lapply(grid, function(t) {
  mll <- cutree(hc, h = t)
  g   <- length(unique(mll))
  ord <- order(mll, -callrate)
  reps_idx <- ord[!duplicated(mll[ord])]
  cap <- sum(colSums(Pp[reps_idx, , drop = FALSE]) > 0)
  data.frame(threshold = t, n_genets = g,
             pct_redundancy = round(100 * (N - g) / N, 1),
             pct_allele_capture = round(100 * cap / A_full, 2))
}))
write.csv(sens, "outputs/mlg_threshold_sensitivity.csv", row.names = FALSE)

sc_y <- max(sens$n_genets) / 100
p_s6 <- ggplot(sens, aes(threshold)) +
  geom_vline(xintercept = CHOSEN, linetype = 2, colour = "grey40") +
  geom_line(aes(y = n_genets, colour = "Distinct genets"), linewidth = 0.8) +
  geom_point(aes(y = n_genets, colour = "Distinct genets")) +
  geom_line(aes(y = pct_allele_capture * sc_y, colour = "Allele capture (%)"),
            linewidth = 0.8) +
  geom_point(aes(y = pct_allele_capture * sc_y, colour = "Allele capture (%)")) +
  scale_y_continuous(name = "Distinct genets (MLLs)",
                     sec.axis = sec_axis(~ . / sc_y, name = "Core allele capture (%)")) +
  scale_colour_manual(values = c("Distinct genets" = "#2E7D32",
                                 "Allele capture (%)" = "#D55E00"), name = NULL) +
  annotate("text", x = CHOSEN, y = max(sens$n_genets),
           label = sprintf("chosen = %.3f", CHOSEN), hjust = -0.05, size = 3,
           colour = "grey40") +
  labs(x = "MLL distance threshold (fraction of differing loci)") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.98, 0.5), legend.justification = c(1, 0.5),
        panel.grid.minor = element_blank())
ggsave("figures/mlg_threshold_sensitivity.png", p_s6, width = 7.5, height = 5,
       dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
## Fig. 3: sharing of clonal lineages among countries (Pacific-centred map)
# ---------------------------------------------------------------------------
mm <- read.csv("outputs/mlg_membership.csv", stringsAsFactors = FALSE)
d <- data.frame(
  id      = indNames(gl),
  country = as.character(pop(gl)),
  lon360  = wrap360(gl@other$latlon$lon),
  lat     = gl@other$latlon$lat,
  stringsAsFactors = FALSE)
d <- merge(d, mm[, c("id", "mll")], by = "id")
d <- d[!is.na(d$lon360) & !is.na(d$lat), ]

## Results: private (single-country) lineages per country, e.g. 63 of 72 in Indonesia
lin_countries <- tapply(d$country, d$mll, function(x) unique(x))
countries <- sort(unique(d$country))
csum <- do.call(rbind, lapply(countries, function(k) {
  rows <- d[d$country == k, ]
  lin  <- unique(rows$mll)
  priv <- sum(vapply(lin, function(L) length(lin_countries[[as.character(L)]]) == 1, logical(1)))
  data.frame(country = k, N = nrow(rows), G = length(lin),
             private = priv, private_frac = priv / length(lin),
             lon360 = mean(rows$lon360), lat = mean(rows$lat),
             stringsAsFactors = FALSE)
}))

# map extent (Indo-Pacific core); the United States falls outside and is dropped
xlim <- c(78, 213); ylim <- c(-30, 40)
in_ext <- with(csum, lon360 >= xlim[1] & lon360 <= xlim[2] &
                     lat >= ylim[1] & lat <= ylim[2])
csum <- csum[in_ext, ]
write.csv(csum, "outputs/mlg_country_summary.csv", row.names = FALSE)

# edges: country pairs sharing at least one lineage, weighted by the number shared
pair_acc <- list()
for (L in names(lin_countries)) {
  cs <- intersect(lin_countries[[L]], csum$country)
  if (length(cs) < 2) next
  cmb <- combn(sort(cs), 2)
  for (j in seq_len(ncol(cmb))) {
    key <- paste(cmb[1, j], cmb[2, j], sep = "||")
    pair_acc[[key]] <- c(pair_acc[[key]], as.integer(L))
  }
}
edges <- do.call(rbind, lapply(names(pair_acc), function(key) {
  ab <- strsplit(key, "\\|\\|")[[1]]
  data.frame(a = ab[1], b = ab[2], n_shared = length(pair_acc[[key]]),
             stringsAsFactors = FALSE)
}))
edges$ax <- csum$lon360[match(edges$a, csum$country)]
edges$ay <- csum$lat[match(edges$a, csum$country)]
edges$bx <- csum$lon360[match(edges$b, csum$country)]
edges$by <- csum$lat[match(edges$b, csum$country)]
edges <- edges[order(edges$n_shared), ]
write.csv(edges[, c("a", "b", "n_shared")],
          "outputs/mlg_country_sharing_edges.csv", row.names = FALSE)

# basemap: crop the world into two dateline-safe boxes, then shift the eastern box to [180, 213]
sf::sf_use_s2(FALSE)
land <- st_make_valid(ne_countries(scale = "medium", returnclass = "sf"))
west <- suppressWarnings(st_crop(land, xmin = xlim[1], ymin = ylim[1],
                                 xmax = 180, ymax = ylim[2]))
east <- suppressWarnings(st_crop(land, xmin = -180, ymin = ylim[1],
                                 xmax = xlim[2] - 360, ymax = ylim[2]))
east <- st_shift_longitude(east)
land <- rbind(west, east)

p_map <- ggplot() +
  geom_sf(data = land, fill = "#efeee9", colour = "#b9b6ad", linewidth = 0.2) +
  geom_segment(data = edges,
               aes(x = ax, y = ay, xend = bx, yend = by,
                   linewidth = n_shared, alpha = n_shared),
               colour = "#7d3c98", lineend = "round") +
  scale_linewidth_continuous(range = c(0.2, 3.0), name = "Shared lineages",
                             breaks = c(1, 3, 6, 9)) +
  scale_alpha_continuous(range = c(0.15, 0.85), guide = "none") +
  geom_point(data = csum,
             aes(x = lon360, y = lat, size = N, fill = private_frac),
             shape = 21, colour = "grey20", stroke = 0.4) +
  scale_size_area(max_size = 12, name = "Samples (N)",
                  breaks = c(10, 50, 100, 185)) +
  scale_fill_gradient(low = "#e67300", high = "#1a6b3b",
                      name = "Endemic fraction\n(private lineages / G)",
                      limits = c(0, 1)) +
  geom_text_repel(data = csum,
                  aes(x = lon360, y = lat, label = pretty_country(country)),
                  size = 3.1, fontface = "bold", colour = "grey15",
                  max.overlaps = Inf, box.padding = 0.5,
                  min.segment.length = 0, segment.colour = "grey55",
                  segment.size = 0.3) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#dbe9f2"),
    panel.grid.major = element_line(colour = "#ffffff", linewidth = 0.2),
    axis.title       = element_blank(),
    legend.position  = "right",
    legend.box       = "vertical",
    plot.subtitle    = element_text(size = 8.5, colour = "grey30"))
ggsave("figures/mlg_clone_sharing_map.png", p_map, width = 10.5, height = 6.5,
       units = "in", dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
## Fig. 5: clonal richness R by country (N >= 6) and richness vs Simpson evenness V (N >= 20)
# ---------------------------------------------------------------------------
dv <- read.csv("outputs/mlg_diversity_indices.csv", stringsAsFactors = FALSE)
overall_R <- dv$R[dv$level == "OVERALL"]
pp <- dv[dv$level != "OVERALL" & !is.na(dv$R) & dv$N >= 6, ]
pp$country <- pp$level
region <- c(
  Indonesia = "SE/E Asia", Malaysia = "SE/E Asia", Thailand = "SE/E Asia",
  Vietnam = "SE/E Asia", Philippines = "SE/E Asia", Japan = "SE/E Asia",
  IndiaNorthen = "SE/E Asia",
  PapuaNewGuinea = "Melanesia", SolomonIslands = "Melanesia",
  Vanuatu = "Melanesia", NewCaledonia = "Melanesia", Fiji = "Melanesia",
  Samoa = "Polynesia", WSamoa = "Polynesia", Tonga = "Polynesia",
  Niue = "Polynesia", CookIslands = "Polynesia",
  "Tahiti&FrenchPolynesia" = "Polynesia", Hawaii = "Polynesia",
  Palau = "Micronesia", FederatedStatesofMicronesia = "Micronesia",
  Australia = "Australia", UnitedStatesofAmerica = "Other")
pp$region <- factor(region[pp$country],
                    levels = c("SE/E Asia", "Melanesia", "Polynesia",
                               "Micronesia", "Australia"))
pp$label <- pretty_country(pp$country, short = TRUE)
reg_cols <- c("SE/E Asia" = "#E69F00", "Melanesia" = "#009E73",
              "Polynesia" = "#0072B2", "Micronesia" = "#CC79A7",
              "Australia" = "#9E4A06")

pA <- ggplot(pp, aes(x = R, y = reorder(label, R))) +
  geom_segment(aes(x = 0, xend = R, yend = reorder(label, R)),
               colour = "grey80", linewidth = 0.6) +
  geom_vline(xintercept = overall_R, linetype = 2, colour = "grey40") +
  annotate("text", x = overall_R, y = 0.7,
           label = sprintf("overall R = %.2f", overall_R),
           hjust = -0.05, vjust = 0, size = 2.7, colour = "grey40") +
  geom_point(aes(colour = region, size = N)) +
  scale_colour_manual(values = reg_cols, name = "Region") +
  scale_size_area(max_size = 9, name = "Samples (N)",
                  breaks = c(10, 50, 100, 185)) +
  scale_x_continuous(limits = c(0, 1.0), expand = expansion(c(0, 0.05))) +
  labs(x = "Clonal richness  R = (G − 1)/(N − 1)", y = NULL, tag = "A") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        panel.grid.minor = element_blank())

pB_dat <- pp[pp$N >= 20, ]
pB <- ggplot(pB_dat, aes(x = R, y = V)) +
  geom_point(aes(size = N, colour = region), alpha = 0.85) +
  scale_colour_manual(values = reg_cols, name = "Region", guide = "none") +
  geom_text_repel(aes(label = label), size = 2.9, colour = "grey15",
                  max.overlaps = Inf, box.padding = 0.45,
                  min.segment.length = 0, segment.colour = "grey60",
                  segment.size = 0.3) +
  scale_size_area(max_size = 9, breaks = c(20, 50, 100, 185),
                  name = "Samples (N)", guide = "none") +
  scale_x_continuous(limits = c(0, 1.0)) +
  scale_y_continuous(limits = c(0, 1.0)) +
  annotate("text", x = 0.02, y = 0.03, hjust = 0, vjust = 0, size = 2.9,
           fontface = "italic", colour = "#b0392b",
           label = "few dominant clones") +
  annotate("text", x = 0.98, y = 0.99, hjust = 1, vjust = 1, size = 2.9,
           fontface = "italic", colour = "#1f6f3f",
           label = "many even genets") +
  labs(x = "Clonal richness (R)", y = "Clonal evenness (V)", tag = "B") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        panel.grid.minor = element_blank())
fig5 <- arrangeGrob(pA, pB, ncol = 2, widths = c(1.08, 1))
ggsave("figures/mlg_diversity_by_population.png", fig5,
       width = 11, height = 5.6, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
## Supplementary Data 1 and 2: cultivar synonymy (one genet, many names) and homonymy (one name, many genets)
# ---------------------------------------------------------------------------
memb <- read.csv("outputs/mlg_membership.csv", stringsAsFactors = FALSE,
                 fileEncoding = "latin1")
MISSING <- c("", "n/a", "na", "unknown", "?", "-", "--", "none", "variety")
norm_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", " ", x)
  ifelse(tolower(x) %in% MISSING, NA_character_, x)
}
memb$variety_clean <- norm_name(memb$variety)
memb$name_key      <- tolower(memb$variety_clean)
named <- memb[!is.na(memb$variety_clean), ]

## Results: 98 of 258 named lineages carry two or more names; one lineage carries 99 names in 9 countries
uniq_display <- function(x) paste(sort(unique(x)), collapse = "; ")
syn <- do.call(rbind, lapply(sort(unique(named$mll)), function(m) {
  s <- named[named$mll == m, ]
  keys <- unique(s$name_key)
  if (length(keys) < 2) return(NULL)
  data.frame(
    mll              = m,
    n_ramets         = nrow(s),
    n_distinct_names = length(keys),
    n_countries      = length(unique(s$country)),
    names            = uniq_display(s$variety_clean),
    countries        = uniq_display(s$country),
    common_names     = uniq_display(norm_name(s$common_name)[!is.na(norm_name(s$common_name))]),
    stringsAsFactors = FALSE)
}))
syn <- syn[order(-syn$n_distinct_names, -syn$n_ramets), ]
write.csv(syn, "outputs/cultivar_synonymy.csv", row.names = FALSE)

## Results: 41 names are each attached to two or more lineages
hom <- do.call(rbind, lapply(sort(unique(named$name_key)), function(k) {
  s <- named[named$name_key == k, ]
  mlls <- unique(s$mll)
  if (length(mlls) < 2) return(NULL)
  data.frame(
    name        = s$variety_clean[which.max(nchar(s$variety_clean))][1],
    n_lineages  = length(mlls),
    n_ramets    = nrow(s),
    n_countries = length(unique(s$country)),
    mlls        = paste(sort(mlls), collapse = "; "),
    countries   = uniq_display(s$country),
    stringsAsFactors = FALSE)
}))
hom <- hom[order(-hom$n_lineages, -hom$n_ramets), ]
write.csv(hom, "outputs/cultivar_homonymy.csv", row.names = FALSE)

summ <- data.frame(
  metric = c("accessions_named", "distinct_names", "genets_with_a_name",
             "synonymous_genets", "accessions_in_synonymous_genets",
             "max_names_one_genet", "homonymous_names",
             "max_genets_one_name"),
  value  = c(nrow(named), length(unique(named$name_key)),
             length(unique(named$mll)),
             nrow(syn), sum(syn$n_ramets), max(syn$n_distinct_names),
             nrow(hom), max(hom$n_lineages)))
write.csv(summ, "outputs/cultivar_synonymy_summary.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
## Fig. 2, Supplementary Data 6: candidate core of one ramet per lineage
# ---------------------------------------------------------------------------
memb <- read.csv("outputs/mlg_membership.csv", stringsAsFactors = FALSE)
memb <- memb[match(indNames(gl), memb$id), ]
N_full   <- nInd(gl)
core_ids <- memb$id[memb$is_representative]
N_core   <- length(core_ids)
is_core  <- indNames(gl) %in% core_ids

## Results: the 392-genet core holds 28,809 of 29,718 alleles (96.9%)
Pc     <- Pp[is_core, , drop = FALSE]
A_core <- sum(colSums(Pc) > 0)

## Results: mean He 0.053 in the core against 0.059 in the full set
he_mean <- function(mat) {
  q <- colMeans(mat, na.rm = TRUE) / 2
  mean(2 * q * (1 - q), na.rm = TRUE)
}
poly_frac <- function(mat) {
  q <- colMeans(mat, na.rm = TRUE) / 2
  mean(q > 0 & q < 1, na.rm = TRUE)
}
He_full   <- he_mean(M);   He_core   <- he_mean(M[is_core, , drop = FALSE])
poly_full <- poly_frac(M); poly_core <- poly_frac(M[is_core, , drop = FALSE])

## Results: greedy selection reaches 95% of the alleles with 41 genets and 99% with 167
greedy_order <- function(Pmat) {
  n <- nrow(Pmat)
  captured <- logical(ncol(Pmat))
  order_idx <- integer(0); gained <- integer(0)
  remaining <- seq_len(n)
  repeat {
    newc <- as.integer(Pmat[remaining, !captured, drop = FALSE] %*%
                         rep(1L, sum(!captured)))
    best <- which.max(newc)
    if (length(best) == 0L || newc[best] == 0L) break
    pick <- remaining[best]
    order_idx <- c(order_idx, pick); gained <- c(gained, newc[best])
    captured <- captured | Pmat[pick, ]
    remaining <- remaining[-best]
    if (all(captured) || length(remaining) == 0L) break
  }
  list(order = order_idx, gained = gained, cum = cumsum(gained))
}
gr      <- greedy_order(Pc)
A_reach <- max(gr$cum)
n_for   <- function(frac) which(gr$cum >= frac * A_reach)[1]

# accumulation curves: greedy genets vs random genets vs random accessions (50 permutations)
rand_curve <- function(Pmat, nperm = 50L, npts = 60L) {
  n <- nrow(Pmat)
  xs <- unique(round(seq(1, n, length.out = min(npts, n))))
  mat <- matrix(NA_real_, nperm, length(xs))
  for (b in seq_len(nperm)) {
    ord <- sample.int(n)
    cap <- logical(ncol(Pmat)); got <- integer(n)
    for (i in seq_len(n)) {
      cap <- cap | Pmat[ord[i], ]; got[i] <- sum(cap)
    }
    mat[b, ] <- got[xs]
  }
  data.frame(x = xs, frac = colMeans(mat) / A_full,
             lo = apply(mat, 2, quantile, 0.025) / A_full,
             hi = apply(mat, 2, quantile, 0.975) / A_full)
}
set.seed(42)
cur_greedy <- data.frame(x = seq_along(gr$cum), frac = gr$cum / A_full,
                         lo = gr$cum / A_full, hi = gr$cum / A_full,
                         set = "Greedy (one ramet per lineage)")
rg <- rand_curve(Pc, nperm = 50L); rg$set <- "Random genets"
ra <- rand_curve(Pp, nperm = 50L); ra$set <- "Random accessions (with clones)"
curves <- rbind(cur_greedy, rg[, c("x", "frac", "lo", "hi", "set")],
                ra[, c("x", "frac", "lo", "hi", "set")])
curves$set <- factor(curves$set,
  levels = c("Greedy (one ramet per lineage)", "Random genets",
             "Random accessions (with clones)"))

p_core <- ggplot(curves, aes(x, frac, colour = set)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = set), colour = NA, alpha = 0.18,
              show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = N_core, linetype = 2, colour = "grey40") +
  annotate("text", x = N_core, y = 0.15, hjust = -0.05, size = 3, colour = "grey40",
           label = sprintf("core = %d genets", N_core)) +
  geom_hline(yintercept = A_core / A_full, linetype = 3, colour = "grey60") +
  scale_colour_manual(values = c("#2E7D32", "#0072B2", "#D55E00"), name = NULL) +
  scale_fill_manual(values = c("#2E7D32", "#0072B2", "#D55E00")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Accessions retained", y = "SNP alleles captured (% of full set)") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.7)),
        panel.grid.minor = element_blank())
ggsave("figures/core_collection_accumulation.png", p_core,
       width = 7.5, height = 5.2, dpi = 300, bg = "white")

core <- memb[memb$is_representative,
             c("id", "country", "variety", "common_name", "mll", "mll_size", "call_rate")]
core <- core[order(-core$mll_size, core$country), ]
names(core)[names(core) == "mll_size"] <- "ramets_represented"
write.csv(core, "outputs/core_collection.csv", row.names = FALSE)

summary_tab <- data.frame(
  metric = c("accessions_full", "genets_core", "duplicates_removed",
             "maintenance_saving_pct", "genets_retained_pct",
             "alleles_full", "alleles_core", "allele_retention_pct",
             "alleles_lost_with_duplicates", "He_full", "He_core",
             "He_retention_pct", "polymorphic_pct_full", "polymorphic_pct_core",
             "greedy_core_95pct", "greedy_core_99pct", "greedy_core_100pct"),
  value = c(N_full, N_core, N_full - N_core,
            round(100 * (N_full - N_core) / N_full, 1), 100,
            A_full, A_core, round(100 * A_core / A_full, 2),
            A_full - A_core, round(He_full, 4), round(He_core, 4),
            round(100 * He_core / He_full, 1), round(100 * poly_full, 1),
            round(100 * poly_core, 1), n_for(0.95), n_for(0.99), length(gr$order)))
write.csv(summary_tab, "outputs/core_collection_summary.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
## Results: the 909 alleles the core does not capture are rare (median frequency 0.0009); 684 are singletons
# ---------------------------------------------------------------------------
present_core <- colSums(P[is_core, , drop = FALSE]) > 0
lost <- which(present_full & !present_core)
mll  <- memb$mll
carrier_n <- integer(length(lost))
lineage_n <- integer(length(lost))
for (j in seq_along(lost)) {
  who <- which(P[, lost[j]])
  carrier_n[j] <- length(who)
  lineage_n[j] <- length(unique(mll[who]))
}
singleton  <- carrier_n == 1
multi_1lin <- carrier_n >= 2 & lineage_n == 1
multi_xlin <- lineage_n >= 2
maf_lost   <- colSums(P[, lost, drop = FALSE]) / nrow(M)
lost_tab <- data.frame(
  metric = c("alleles_full", "alleles_core", "alleles_lost",
             "lost_singletons", "lost_singleton_pct",
             "lost_multi_one_lineage", "lost_multi_one_lineage_pct",
             "lost_multi_across_lineages", "lost_multi_across_pct",
             "lost_median_freq"),
  value  = round(c(A_full, sum(present_core), length(lost),
                   sum(singleton), 100 * mean(singleton),
                   sum(multi_1lin), 100 * mean(multi_1lin),
                   sum(multi_xlin), 100 * mean(multi_xlin),
                   median(maf_lost)), 4))
write.csv(lost_tab, "outputs/lost_alleles.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
## Supplementary Table 3: core allele capture as a function of the minimum allele count
# ---------------------------------------------------------------------------
# allele counts as carriers (accessions carrying the allele) and as copies
ref_cop <- 2L - M; ref_cop[is.na(ref_cop)] <- 0L
alt_cop <- M;      alt_cop[is.na(alt_cop)] <- 0L
Ccopies <- c(colSums(ref_cop), colSums(alt_cop)); rm(ref_cop, alt_cop)
carriers <- colSums(P)
keep_any <- carriers > 0
core_carriers <- colSums(P[is_core, , drop = FALSE])

mac_row <- function(mac, basis = c("carriers", "copies")) {
  basis <- match.arg(basis)
  cnt <- if (basis == "carriers") carriers else Ccopies
  sel <- keep_any & cnt >= mac
  n_all <- sum(sel)
  n_cap <- sum(sel & core_carriers > 0)
  data.frame(basis = basis, min_allele_count = mac,
             alleles_in_full_set = n_all,
             alleles_in_core = n_cap,
             alleles_missed = n_all - n_cap,
             capture_pct = round(100 * n_cap / n_all, 2))
}
mac_tab <- rbind(
  do.call(rbind, lapply(1:5, mac_row, basis = "carriers")),
  do.call(rbind, lapply(c(1, 2, 3, 4, 5, 10), mac_row, basis = "copies")))
write.csv(mac_tab, "outputs/brad03_core_mac.csv", row.names = FALSE)

## Supplementary Table 3: genets needed for 95% and 99% of the alleles at MAC 1, 2 and 3
greedy_curve <- function(sel_alleles, ids_idx) {
  Ps <- P[ids_idx, sel_alleles, drop = FALSE]
  reachable <- colSums(Ps) > 0
  Ps <- Ps[, reachable, drop = FALSE]
  total <- ncol(Ps)
  remaining <- rep(TRUE, total)
  avail <- rep(TRUE, nrow(Ps))
  captured <- integer(0)
  repeat {
    idx_rem <- which(remaining)
    gains <- rowSums(Ps[, idx_rem, drop = FALSE])
    gains[!avail] <- -1
    best <- which.max(gains)
    if (gains[best] <= 0) break
    avail[best] <- FALSE
    remaining[idx_rem[Ps[best, idx_rem]]] <- FALSE
    captured <- c(captured, total - sum(remaining))
    if (!any(remaining)) break
  }
  data.frame(n_selected = seq_along(captured),
             alleles = captured, pct = 100 * captured / total,
             reachable_alleles = total)
}
core_idx <- which(is_core)
mac_curves <- do.call(rbind, lapply(c(1, 2, 3), function(mac) {
  sel <- which(keep_any & carriers >= mac)
  cv <- greedy_curve(sel, core_idx)
  cv$min_allele_count <- mac
  cv$alleles_in_full_set <- length(sel)
  cv
}))
write.csv(mac_curves, "outputs/brad03_core_accumulation.csv", row.names = FALSE)

# ============================================================================
# 7. Ploidy threshold and redundancy (SI S2.1; Table S4 rows)
# ============================================================================

## Results: stricter thresholds; 0.0923 -> 0.06 gives 810 accessions, redundancy 63.2% -> 63.6%, core allele capture 96.9% -> 96.5%
THRESH <- 0.0923
Md <- as.matrix(gd); rownames(Md) <- indNames(gd)
Dmat <- as.matrix(readRDS("outputs/mlg_bitwise_dist.rds"))
CUT <- 0.008
alleles_present <- function(ids) {
  m <- Md[ids, , drop = FALSE]
  sum(colSums(m <= 1L, na.rm = TRUE) > 0) + sum(colSums(m >= 1L, na.rm = TRUE) > 0)
}
redundancy_at <- function(ids) {
  hc <- hclust(as.dist(Dmat[ids, ids, drop = FALSE]), method = "complete")
  mll <- cutree(hc, h = CUT)
  cr <- rowSums(!is.na(Md[ids, , drop = FALSE]))
  reps <- tapply(seq_along(ids), mll, function(ix) ids[ix][which.max(cr[ix])])
  ng <- length(unique(mll)); n <- length(ids)
  data.frame(threshold = NA, N = n, genets = ng, redundancy = round(100 * (1 - ng / n), 1),
             allele_capture = round(100 * alleles_present(unlist(reps)) / alleles_present(ids), 2))
}
sweep <- c(0.0923, 0.08, 0.07, 0.06, 0.05)
res <- do.call(rbind, lapply(sweep, function(thr) {
  ids <- ind$id[(!is.na(ind$class) & ind$class == "Diploid") | (is.na(ind$class) & ind$Ho <= thr)]
  r <- redundancy_at(intersect(ids, indNames(gd))); r$threshold <- thr; r
}))
# countries of the accessions dropped when tightening 0.0923 -> 0.06
dropped_ids <- setdiff(ws$id, ind$id[(!is.na(ind$class) & ind$class == "Diploid") | (is.na(ind$class) & ind$Ho <= 0.06)])
dropped_country <- sort(table(ctry[intersect(dropped_ids, names(ctry))]), decreasing = TRUE)
write.csv(data.frame(
  metric = c("cv_sensitivity", "cv_specificity", "loo_misclassified",
             "gap_max_diploid", "gap_min_polyploid", "classes_overlap",
             "working_set_n", "near_cut_0.02_n", "near_cut_0.02_pct",
             "redundancy_at_0.0923", "redundancy_at_0.06",
             "capture_at_0.0923", "capture_at_0.06"),
  value = round(c(cv_sens, cv_spec, sum(labeled$y != pred),
                  gap_lo, gap_hi, as.numeric(gap_lo > gap_hi),
                  nrow(ws), sum(ws$Ho > THRESH - 0.02), 100 * mean(ws$Ho > THRESH - 0.02),
                  res$redundancy[res$threshold == 0.0923], res$redundancy[res$threshold == 0.06],
                  res$allele_capture[res$threshold == 0.0923], res$allele_capture[res$threshold == 0.06]), 4)),
  "outputs/ploidy_gate_cv.csv", row.names = FALSE)

# ============================================================================
# 8. Population structure (Figs 6, 7, 8, S7, S8, S9; Supplementary Data 3)
# ============================================================================

# Clone-censored genets (392), loci monomorphic after decloning removed
gl_mll <- gl.filter.monomorphs(readRDS("input_files/taro_diploid_mll.rds"), verbose = 0)

# Country order west -> east by mean longitude (dateline-unwrapped)
lon360 <- gl_mll@other$latlon$lon
lon360 <- ifelse(is.finite(lon360) & lon360 < 0, lon360 + 360, lon360)
pop_order <- names(sort(tapply(lon360, as.character(pop(gl_mll)), mean, na.rm = TRUE)))

# Country -> region
region <- c(
  Indonesia = "SE/E Asia", Malaysia = "SE/E Asia", Thailand = "SE/E Asia",
  Vietnam = "SE/E Asia", Philippines = "SE/E Asia", Japan = "SE/E Asia",
  PapuaNewGuinea = "Melanesia", SolomonIslands = "Melanesia",
  Vanuatu = "Melanesia", NewCaledonia = "Melanesia", Fiji = "Melanesia",
  Samoa = "Polynesia", WSamoa = "Polynesia", Tonga = "Polynesia",
  Niue = "Polynesia", CookIslands = "Polynesia",
  "Tahiti&FrenchPolynesia" = "Polynesia", Hawaii = "Polynesia",
  Palau = "Micronesia", FederatedStatesofMicronesia = "Micronesia")
region_levels <- c("SE/E Asia", "Melanesia", "Polynesia", "Micronesia")

# genlight -> PLINK bed -> 1/2-coded ped -> OHANA .dgm
gl2dgm <- function(gl, out) {
  gl2plink(gl, outpath = out, bed.files = TRUE,
           plink.bin.path = dirname(Sys.which("plink")), verbose = 0)
  system2("plink", c("--bfile", file.path(out, "gl_plink"), "--recode", "12", "tab",
                     "--allow-extra-chr", "--out", file.path(out, "plink")))
  system2("convert", c("ped2dgm", file.path(out, "plink.ped"), file.path(out, "g.dgm")))
}

# Final iteration, log likelihood and delta from a qpas log
parse_qpas <- function(path) {
  rows <- grep("^[[:space:]]*[0-9]+[[:space:]]", readLines(path, warn = FALSE), value = TRUE)
  f <- strsplit(trimws(tail(rows, 1)), "[[:space:]]+")[[1]]
  list(iteration = as.integer(f[1]), loglik = as.numeric(f[3]), delta = as.numeric(f[4]))
}

# --- OHANA on the 392 genets -------------------------------------------------

## Methods: OHANA qpas K = 2-8 on the genets (outputs/ohana_declone/q.kK.matrix, f.kK.matrix)
ohana_dir <- "outputs/ohana_declone"
gl2dgm(gl_mll, ohana_dir)
loglik <- data.frame(K = 2:8, loglik = NA_real_)
for (k in 2:8) {
  log_k <- file.path(ohana_dir, sprintf("qpas.k%d.log", k))
  system2("qpas", c(file.path(ohana_dir, "g.dgm"), "-k", k,
                    "-qo", file.path(ohana_dir, sprintf("q.k%d.matrix", k)),
                    "-fo", file.path(ohana_dir, sprintf("f.k%d.matrix", k)),
                    "-mi", 200, "-e", "1e-4"), stdout = log_k)
  loglik$loglik[loglik$K == k] <- parse_qpas(log_k)$loglik
}
write.csv(loglik, file.path(ohana_dir, "qpas_loglik_by_K.csv"), row.names = FALSE)

# Component covariance and component tree at K = 6
system2("nemeco", c(file.path(ohana_dir, "g.dgm"), file.path(ohana_dir, "f.k6.matrix"),
                    "-co", file.path(ohana_dir, "c.matrix"), "-mi", 200))
system2("convert", c("cov2nwk", file.path(ohana_dir, "c.matrix"), file.path(ohana_dir, "tree.nwk")))

## Fig. 6: OHANA K = 6 admixture of the 392 genets, countries west to east (figures/ohana_declone_admixture_k6.png)
fam <- read.table(file.path(ohana_dir, "gl_plink.fam"), stringsAsFactors = FALSE)
Q <- as.matrix(read.table(file.path(ohana_dir, "q.k6.matrix"), skip = 1))
colnames(Q) <- paste0("K", seq_len(ncol(Q)))
ctry <- setNames(as.character(pop(gl_mll)), indNames(gl_mll))
df <- data.frame(ind = fam$V2, pop = factor(ctry[fam$V2], levels = pop_order),
                 row = seq_len(nrow(Q)), stringsAsFactors = FALSE)
# genets within a country ordered by that country's dominant component
pop_dom <- tapply(seq_len(nrow(df)), df$pop,
                  function(i) which.max(colMeans(Q[i, , drop = FALSE])))
df$sortkey <- Q[cbind(seq_len(nrow(df)), pop_dom[as.integer(df$pop)])]
df <- df[order(df$pop, -df$sortkey), ]
df$x <- seq_len(nrow(df))
levels(df$pop) <- pretty_country(levels(df$pop), short = TRUE)
long <- pivot_longer(cbind(df[, c("x", "pop")], as.data.frame(Q)[df$row, ]),
                     starts_with("K"), names_to = "comp", values_to = "prop")
long$comp <- factor(long$comp, levels = colnames(Q))
p6 <- ggplot(long, aes(x = x, y = prop, fill = comp)) +
  geom_col(width = 1, position = "stack") +
  facet_grid(. ~ pop, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = setNames(unname(polychrome(ncol(Q))), colnames(Q)),
                    name = "Component") +
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Ancestry proportion") +
  theme_minimal(base_size = 11) +
  theme(panel.spacing = unit(0.28, "lines"), panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
        strip.placement = "outside", legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, title.position = "left", label.position = "right"))
ggsave("figures/ohana_declone_admixture_k6.png", p6, width = 13, height = 4.6,
       dpi = 200, bg = "white")

# --- Number of ancestry components (all 1,064 accessions) --------------------

## Methods: OHANA input for all 1,064 accessions (outputs/ohana/g.dgm)
gl_dip <- readRDS("input_files/taro_diploid.rds")
gl2dgm(gl_dip, "outputs/ohana")

## Fig. S7A: 30 qpas starts per K (seeds 0-28 plus the seed of the original run), best kept (outputs/ohana_multistart/)
ms_dir <- "outputs/ohana_multistart"
run_dir <- file.path(ms_dir, "runs"); best_dir <- file.path(ms_dir, "best")
legacy_seeds <- data.frame(K = 2:8,
                           seed = c(77527064, 3834932940, 2644120585, 1497834555,
                                    3915235949, 1285561651, 1282352437))
jobs <- rbind(expand.grid(seed = 0:28, K = 8:2),
              legacy_seeds[order(legacy_seeds$K, decreasing = TRUE), ])
seed_stem <- function(k, seed)
  sprintf("k%d.seed%s", k, if (seed <= 31) sprintf("%02.0f", seed) else sprintf("%.0f", seed))
run_one <- function(k, seed) {
  stem <- seed_stem(k, seed)
  log_file <- file.path(run_dir, paste0(stem, ".log"))
  status <- system2("qpas", c("outputs/ohana/g.dgm", "-k", k, "-s", sprintf("%.0f", seed),
                              "-qo", file.path(run_dir, paste0(stem, ".q.matrix")),
                              "-fo", file.path(run_dir, paste0(stem, ".f.matrix")),
                              "-mi", 200, "-e", "1e-4"),
                    stdout = log_file, stderr = log_file)
  p <- parse_qpas(log_file)
  data.frame(K = k, seed = seed, loglik = p$loglik, final_iteration = p$iteration,
             final_delta = p$delta, reached_limit = p$iteration >= 200, status = status)
}
results <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)),
                                             function(i) run_one(jobs$K[i], jobs$seed[i]),
                                             mc.cores = 6, mc.preschedule = FALSE))
results <- results[order(results$K, results$seed), ]
write.csv(results, file.path(ms_dir, "all_runs.csv"), row.names = FALSE)
best <- do.call(rbind, lapply(split(results, results$K), function(x) x[which.max(x$loglik), ]))
best <- best[order(best$K), ]
best$marginal_gain <- c(NA, diff(best$loglik))
write.csv(best, file.path(ms_dir, "best_runs_by_K.csv"), row.names = FALSE)
for (i in seq_len(nrow(best))) {
  stem <- seed_stem(best$K[i], best$seed[i])
  file.copy(file.path(run_dir, paste0(stem, c(".log", ".q.matrix", ".f.matrix"))),
            file.path(best_dir, sprintf(c("qpas.k%d.log", "q.k%d.matrix", "f.k%d.matrix"), best$K[i])),
            overwrite = TRUE)
}
write.csv(best[c("K", "loglik")], file.path(ms_dir, "qpas_loglik_by_K.csv"), row.names = FALSE)

## Fig. S7B: sNMF cross-entropy K = 1-10, five runs, on the same accessions (outputs/snmf_kselect/)
gl_ce <- gl.filter.monomorphs(gl_dip, verbose = 0)
gl2geno(gl_ce, outfile = "taro_dip", outpath = "outputs/snmf_kselect", verbose = 0)
proj <- snmf("outputs/snmf_kselect/taro_dip.geno", K = 1:10, entropy = TRUE,
             repetitions = 5, project = "new", CPU = 6, seed = 42)
ce_df <- data.frame(K = 1:10, cross_entropy = sapply(1:10, function(k)
  mean(cross.entropy(proj, K = k), na.rm = TRUE)))
ce_df$dCE <- c(NA, diff(ce_df$cross_entropy))

## Fig. S7: OHANA log likelihood and sNMF cross-entropy by K (figures/k_selection.png)
CHOSEN <- 6
pA <- ggplot() +
  geom_vline(xintercept = CHOSEN, linetype = 2, colour = "grey40") +
  geom_point(data = results, aes(K, loglik), colour = "grey65", alpha = 0.3) +
  geom_line(data = best, aes(K, loglik), colour = "#2E7D32", linewidth = 0.8) +
  geom_point(data = best, aes(K, loglik), colour = "#2E7D32") +
  scale_x_continuous(breaks = best$K) +
  labs(x = "K (ancestry components)", y = "OHANA log likelihood", tag = "A") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
pB <- ggplot(ce_df, aes(K, cross_entropy)) +
  geom_vline(xintercept = CHOSEN, linetype = 2, colour = "grey40") +
  geom_line(colour = "#0072B2", linewidth = 0.8) + geom_point(colour = "#0072B2") +
  scale_x_continuous(breaks = ce_df$K) +
  labs(x = "K (ancestry components)", y = "sNMF cross entropy", tag = "B") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("figures/k_selection.png", arrangeGrob(pA, pB, ncol = 2), width = 10,
       height = 4.5, dpi = 300, bg = "white")

## SI S4.2: log-likelihood gain and cross-entropy change per added K (outputs/k_selection.csv)
ksel <- data.frame(K = best$K, best_seed = best$seed, ohana_loglik = round(best$loglik, 1),
                   ohana_dLL = round(best$marginal_gain, 1))
ksel <- merge(ksel, ce_df, by = "K", all = TRUE)
write.csv(ksel[order(ksel$K), ], "outputs/k_selection.csv", row.names = FALSE)

# --- PCA of the genets -------------------------------------------------------

## Fig. 7: PCA of the 392 genets, PC1 vs PC2 and PC1 vs PC3 (figures/pca_declone_combined.png, outputs/pca_declone.rds)
set.seed(42)
region_pca <- c(region, IndiaNorthen = "SE/E Asia", Australia = "Australia",
                UnitedStatesofAmerica = "Other")
glr <- gl_mll
pop(glr) <- factor(region_pca[as.character(pop(gl_mll))],
                   levels = c(region_levels, "Australia", "Other"))
# impute globally (one dummy population) so imputation does not inflate the structure shown
imp <- gl_mll
pop(imp) <- factor(rep("all", nInd(imp)))
imp <- gl.gen2fbm(imp, verbose = 0)
imp <- gl.impute(imp, verbose = 0)
pca <- gl.pcoa(imp, verbose = 0)
saveRDS(pca, "outputs/pca_declone.rds")
p12 <- gl.pcoa.plot(pca, glr, xaxis = 1, yaxis = 2, pt.size = 3, verbose = 0)
p13 <- gl.pcoa.plot(pca, glr, xaxis = 1, yaxis = 3, pt.size = 3, verbose = 0)
ggsave("figures/pca_declone_combined.png",
       p12 + p13 + patchwork::plot_annotation(tag_levels = "A"),
       width = 12, height = 5.2, dpi = 200)

# --- Country tree ------------------------------------------------------------

## Fig. 8: F81 distances among countries from the tag sequences, Fitch-Margoliash tree (outputs/dist_phylo_declone.rds, outputs/tree_fitch_declone.rds)
set.seed(42)
gl_tree <- gl_mll
popNames(gl_tree) <- gsub("[^A-Za-z0-9]", "", popNames(gl_tree))
gl_tree@other$loc.metrics$TrimmedSequence <- gl_tree@other$loc.metrics$TrimmedSequenceSnp
Dphylo <- gl.dist.phylo(gl_tree, verbose = 0)
saveRDS(Dphylo, "outputs/dist_phylo_declone.rds")
pdf("figures/taro_tree_declone.pdf", width = 8, height = 8)
tr <- gl.tree.fitch(Dphylo, x = gl_tree, phylip.path = dirname(Sys.which("fitch")),
                    out.path = tempdir(), tree.method = "FM", plot.type = "phylogram",
                    verbose = 0)
dev.off()
saveRDS(tr, "outputs/tree_fitch_declone.rds")

## Fig. 8: country tree with full names and region colours (figures/taro_tree_declone_ggtree.png)
# PHYLIP truncates tip labels to 10 alphanumeric characters
orig <- names(region)
map <- data.frame(trunc = substr(gsub("[^A-Za-z0-9]", "", orig), 1, 10),
                  pretty = unname(pretty_country(orig)), region = unname(region),
                  stringsAsFactors = FALSE)
tr_plot <- tr
tr_plot$tip.label <- map$pretty[match(tr_plot$tip.label, map$trunc)]
tipdat <- data.frame(label = map$pretty, region = factor(map$region, levels = region_levels))
pal <- c("SE/E Asia" = "#F8766D", "Melanesia" = "#7CAE00",
         "Polynesia" = "#00BFC4", "Micronesia" = "#C77CFF")
p8 <- ggtree(tr_plot, linewidth = 0.6) %<+% tipdat +
  geom_tippoint(aes(colour = region), size = 2.6) +
  geom_tiplab(aes(colour = region), size = 3.2,
              offset = max(tr_plot$edge.length) * 0.03, show.legend = FALSE) +
  scale_colour_manual(values = pal, name = "Region", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 3.5))) +
  ggtree::hexpand(0.35) +
  theme_tree2() +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.background = element_blank(), legend.key = element_blank(),
        axis.title.x = element_text(size = 9)) +
  labs(x = "Branch length (F81 distance)")
ggsave("figures/taro_tree_declone_ggtree.png", p8, width = 8, height = 7, dpi = 200, bg = "white")

## Supplementary Data 3: bootstrap support of the country tree, 1,000 locus resamples (outputs/tree_support.csv)
# Vectorised F81 distance among individuals from the tags, averaged between countries (reproduces gl.dist.phylo)
gl_bs <- gl.filter.overshoot(gl_tree, verbose = 0)
M <- as.matrix(gl_bs)
called <- !is.na(M)
ref_hom <- M == 0 & called
alt_hom <- M == 2 & called
homozygous <- ref_hom | alt_hom
tag <- as.character(gl_bs@other$loc.metrics$TrimmedSequence)
tag_length <- nchar(tag)
base_levels <- c("A", "C", "G", "T")
# invariant bases per tag after removing the SNP position
flank_counts <- t(vapply(tag, function(s) {
  tabulate(match(strsplit(toupper(s), "", fixed = TRUE)[[1]], base_levels), nbins = 4)
}, numeric(4)))
snp_index <- gl_bs@other$loc.metrics$SnpPosition + 1L
snp_base_index <- match(substring(tag, snp_index, snp_index), base_levels)
flank_cell <- cbind(seq_len(nLoc(gl_bs)), snp_base_index)
flank_counts[flank_cell] <- flank_counts[flank_cell] - 1
# heterozygous SNPs are ambiguity codes that ape excludes from base frequencies
alleles <- do.call(rbind, strsplit(as.character(gl_bs@loc.all), "/", fixed = TRUE))
n_called <- colSums(called); n_ref <- colSums(ref_hom); n_alt <- colSums(alt_hom)
locus_base_counts <- flank_counts * n_called
for (b in seq_along(base_levels)) {
  locus_base_counts[, b] <- locus_base_counts[, b] +
    n_ref * (alleles[, 1] == base_levels[b]) + n_alt * (alleles[, 2] == base_levels[b])
}
country <- as.character(pop(gl_bs))
country_levels <- levels(pop(gl_bs))
country_index <- lapply(country_levels, function(z) which(country == z))
calculate_country_f81 <- function(weights) {
  ref_alt <- tcrossprod(sweep(ref_hom, 2, weights, `*`), alt_hom)
  mismatch <- ref_alt + t(ref_alt)
  comparable <- tcrossprod(sweep(called, 2, weights * (tag_length - 1L), `*`), called) +
    tcrossprod(sweep(homozygous, 2, weights, `*`), homozygous)
  base_count <- colSums(locus_base_counts * weights)
  base_frequency <- base_count / sum(base_count)
  e_f81 <- 1 - sum(base_frequency^2)
  individual_f81 <- -e_f81 * log1p(-(mismatch / comparable) / e_f81)
  diag(individual_f81) <- 0
  country_f81 <- matrix(0, length(country_levels), length(country_levels),
                        dimnames = list(country_levels, country_levels))
  for (i in seq_len(length(country_levels) - 1L)) {
    for (j in (i + 1L):length(country_levels)) {
      country_f81[i, j] <- country_f81[j, i] <-
        mean(individual_f81[country_index[[i]], country_index[[j]], drop = FALSE], na.rm = TRUE)
    }
  }
  as.dist(country_f81)
}
# one fixed seed per replicate, drawn from the master seed
B <- 1000
set.seed(42)
replicate_seeds <- sample.int(.Machine$integer.max, B)
bootstrap_distances <- lapply(seq_len(B), function(i) {
  set.seed(replicate_seeds[i])
  weights <- tabulate(sample.int(nLoc(gl_bs), nLoc(gl_bs), replace = TRUE), nbins = nLoc(gl_bs))
  calculate_country_f81(weights)
})
# all replicates as one PHYLIP multiple-dataset infile, fitted with the same settings as the main tree
fitch_dir <- tempfile("tree_fitch_bootstrap_")
dir.create(fitch_dir)
phylip_labels <- sprintf("%-10s", substr(country_levels, 1L, 10L))
con <- file(file.path(fitch_dir, "infile"), open = "wt")
for (distance in bootstrap_distances) {
  md <- as.matrix(distance)[country_levels, country_levels]
  writeLines(as.character(length(country_levels)), con)
  writeLines(vapply(seq_len(nrow(md)), function(i) {
    paste(phylip_labels[i], paste(format(md[i, ], digits = 16, scientific = TRUE, trim = TRUE),
                                  collapse = " "))
  }, character(1)), con)
}
close(con)
# multiple-dataset mode turns jumbling on; toggle it off to keep the input order of the main tree
writeLines(c("M", as.character(B), "331", "J", "2", "3", "Y"), file.path(fitch_dir, "fitch.cmd"))
old_wd <- setwd(fitch_dir)
system2("fitch", stdin = "fitch.cmd", stdout = "fitch.log", stderr = "fitch.log")
setwd(old_wd)
bootstrap_trees <- ape::read.tree(file.path(fitch_dir, "outtree"))
saveRDS(bootstrap_trees, "outputs/tree_fitch_bootstrap_trees.rds")
support <- ape::prop.clades(tr, bootstrap_trees, rooted = FALSE) / B * 100
supported_tree <- tr
supported_tree$node.label <- ifelse(is.na(support), NA_character_, sprintf("%.1f", support))
saveRDS(supported_tree, "outputs/tree_fitch_declone_supported.rds")
number_tips <- length(tr$tip.label)
support_table <- do.call(rbind, lapply(seq_len(tr$Nnode), function(i) {
  tips <- ape::extract.clade(tr, number_tips + i)$tip.label
  data.frame(support = round(support[i], 1), n_tips = length(tips),
             clade = paste(sort(tips), collapse = ","), stringsAsFactors = FALSE)
}))
write.csv(support_table[order(-support_table$support), ], "outputs/tree_support.csv",
          row.names = FALSE)

# --- Admixture validation with sNMF and PopCluster ---------------------------

## Fig. S8: sNMF K = 6 on the 392 genets, 10 runs, regularisation 1, dendrogram-ordered (figures/snmf_declone_K6_den.png)
set.seed(42)
gl_s <- gl_mll
levels(gl_s@pop) <- unname(pretty_country(levels(gl_s@pop), short = TRUE))
sr <- gl.run.snmf(gl_s, filename = "taro_declone", minK = 6, maxK = 6, rep = 10,
                  regularization = 1, ncores = 4, plot.out = FALSE, verbose = 0)
saveRDS(sr, "outputs/snmf_declone/snmf_declone_res.rds")
region_s <- setNames(region, pretty_country(names(region), short = TRUE))
qc <- paste0("Pop_", 1:6)

## SI S4.3: mean sNMF ancestry by region (outputs/snmf_declone/mean_ancestry_by_region.csv)
reg <- factor(region_s[as.character(pop(gl_s))], levels = region_levels)
write.csv(aggregate(sr$matrix[["K6"]][qc], list(region = reg), function(x) round(mean(x), 2)),
          "outputs/snmf_declone/mean_ancestry_by_region.csv", row.names = FALSE)

clust_cols <- setNames(unname(polychrome(6)), qc)
reg_cols <- c("SE/E Asia" = "#E69F00", "Melanesia" = "#009E73",
              "Polynesia" = "#0072B2", "Micronesia" = "#CC79A7")
pop_cols <- unname(reg_cols[region_s[levels(pop(gl_s))]])
png("figures/snmf_declone_K6_den.png", width = 16, height = 7.5, units = "in", res = 300)
invisible(gl.plot.snmf(snmf.result = sr, plot.K = 6, den = TRUE, inverse.den = TRUE,
                       x = gl_s, color.clusters = clust_cols, plot.colors.pop = pop_cols,
                       ind.name = FALSE, plot.out = TRUE, plot.file = NULL, verbose = 0))
dev.off()

## Fig. S9: PopCluster K = 1-6, admixture model 2, five runs, on the 392 genets (figures/popcluster_taro_K6.png)
gl_pc <- readRDS("input_files/taro_diploid_mll.rds")
ploidy(gl_pc) <- rep(2L, nInd(gl_pc))
gl_pc <- gl.filter.callrate(gl_pc, method = "ind", threshold = 0.80, verbose = 0)
gl_pc <- gl.filter.allna(gl_pc, verbose = 0)
popcluster_bin <- Sys.which(c("PopClusterMac", "PopClusterLnx", "PopClusterWin"))
res <- gl.run.popcluster(gl_pc, popcluster.path = dirname(popcluster_bin[popcluster_bin != ""][1]),
                         output.path = normalizePath("outputs/popcluster"),
                         filename = "taro_declone", minK = 1, maxK = 6, rep = 5,
                         model = 2, parallel = FALSE, plot.out = FALSE, verbose = 0)
saveRDS(res, "outputs/popcluster/popcluster_taro_result.rds")
Qp <- res$matrix[["taro_declone.popcluster_K_6_R_1"]]
Qp$Pop <- pretty_country(Qp$Pop, short = TRUE)
kc <- paste0("Pop_", 1:6)
Mp <- as.matrix(Qp[, kc])
dom <- max.col(Mp, ties.method = "first")
dfp <- data.frame(pop = factor(Qp$Pop, levels = intersect(pretty_country(pop_order, short = TRUE),
                                                          unique(Qp$Pop))),
                  dom = dom, domval = Mp[cbind(seq_len(nrow(Mp)), dom)], Mp, check.names = FALSE)
dfp <- dfp[order(dfp$pop, dfp$dom, -dfp$domval), ]
dfp$x <- seq_len(nrow(dfp))
longp <- pivot_longer(dfp[, c("x", "pop", kc)], cols = all_of(kc),
                      names_to = "comp", values_to = "prop")
longp$comp <- factor(longp$comp, levels = kc)
p9 <- ggplot(longp, aes(x = x, y = prop, fill = comp)) +
  geom_col(width = 1, position = "stack") +
  facet_grid(. ~ pop, scales = "free_x", space = "free_x", switch = "x") +
  scale_fill_manual(values = setNames(unname(polychrome(6)), kc), name = "Cluster") +
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  labs(y = "Ancestry proportion") + theme_minimal(base_size = 11) +
  theme(panel.spacing = unit(0.28, "lines"), panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
        strip.placement = "outside", legend.position = "right")
ggsave("figures/popcluster_taro_K6.png", p9, width = 13, height = 4.5, dpi = 200, bg = "white")

# ============================================================================
# 9. Diversity by country and cluster, manuscript tables
# ============================================================================

MIN_N <- 5L   # minimum genets for a group to be reported

# diversity estimators: unbiased He (Nei 1978), FIS = 1 - Ho/He, rarefied Ar (El Mousadik & Petit 1996)
ar_locus_vec <- function(n0, n1, twoN, k) {
  out <- rep(NA_real_, length(twoN))
  ok <- twoN >= k & twoN > 0
  if (!any(ok)) return(out)
  p0 <- exp(lchoose(twoN[ok] - n0[ok], k) - lchoose(twoN[ok], k))
  p1 <- exp(lchoose(twoN[ok] - n1[ok], k) - lchoose(twoN[ok], k))
  out[ok] <- (1 - p0) + (1 - p1)
  out
}
group_div <- function(M, idx, k_rare) {
  m  <- M[idx, , drop = FALSE]
  nn <- colSums(!is.na(m))
  het <- colSums(m == 1L, na.rm = TRUE)
  alt <- colSums(m, na.rm = TRUE)
  twoN <- 2 * nn
  q <- ifelse(nn > 0, alt / twoN, NA)
  Ho_l <- ifelse(nn > 0, het / nn, NA)
  He_l <- ifelse(nn > 1, (twoN / (twoN - 1)) * 2 * q * (1 - q), NA)
  Ar_l <- ar_locus_vec(twoN - alt, alt, twoN, k_rare)
  data.frame(N = length(idx),
             poly_pct = round(100 * mean(q > 0 & q < 1, na.rm = TRUE), 1),
             Ho = round(mean(Ho_l, na.rm = TRUE), 4),
             He = round(mean(He_l, na.rm = TRUE), 4),
             FIS = round(1 - mean(Ho_l, na.rm = TRUE) / mean(He_l, na.rm = TRUE), 4),
             Ar = round(mean(Ar_l, na.rm = TRUE), 4))
}
# allele-presence matrix: REF (genotype 0/1) in columns 1..L, ALT (genotype 1/2) in L+1..2L
pres_matrix <- function(M) {
  ref <- M <= 1L; ref[is.na(ref)] <- FALSE
  alt <- M >= 1L; alt[is.na(alt)] <- FALSE
  cbind(ref, alt)
}
# per-group table for groups with >= MIN_N genets; Ar rarefied to 2 x smallest retained group
build_table <- function(M, groups) {
  Pres <- pres_matrix(M)
  tab  <- sort(table(groups), decreasing = TRUE)
  keep <- names(tab)[tab >= MIN_N]
  k_rare <- 2L * min(tab[keep])
  levs <- names(table(groups))
  GP <- t(vapply(levs, function(g) colSums(Pres[groups == g, , drop = FALSE]) > 0,
                 logical(ncol(Pres))))
  priv_cols <- which(colSums(GP) == 1L)
  owner <- levs[max.col(t(GP[, priv_cols, drop = FALSE]), ties.method = "first")]
  priv <- setNames(as.integer(table(factor(owner, levels = levs))), levs)
  res <- do.call(rbind, lapply(keep, function(g) {
    cbind(group = g, group_div(M, which(groups == g), k_rare),
          private_alleles = priv[[g]])
  }))
  res[order(-res$He), ]
}

## Table 1 / Supplementary Data 7: diversity by country on one ramet of every lineage recorded in the country
gl   <- readRDS("input_files/taro_diploid.rds")
memb <- read.csv("outputs/mlg_membership.csv", stringsAsFactors = FALSE)
memb <- memb[match(indNames(gl), memb$id), ]
Mfull <- as.matrix(gl); rownames(Mfull) <- indNames(gl)
country_full <- memb$country
idx_current <- which(memb$is_representative)            # one ramet per lineage overall (392 genets)
key <- paste(memb$country, memb$mll, sep = "|")
ord <- order(key, -memb$call_rate)
idx_percountry <- sort(ord[!duplicated(key[ord])])       # highest-call-rate ramet per lineage per country
by_country_pc <- build_table(Mfull[idx_percountry, , drop = FALSE], country_full[idx_percountry])
write.csv(by_country_pc, "outputs/brad01_diversity_percountry.csv", row.names = FALSE)

## Results / SI S4.5: 53 lineages shared among countries hold 518 accessions; 339 single-country lineages
n_cty_per_mll <- tapply(memb$country, memb$mll, function(x) length(unique(x)))
shared     <- names(n_cty_per_mll)[n_cty_per_mll > 1]
single_mll <- as.integer(names(n_cty_per_mll)[n_cty_per_mll == 1L])
idx_single <- idx_current[memb$mll[idx_current] %in% single_mll]
write.csv(data.frame(
  analysis = c("n_single_country_lineages", "n_genets_single_country",
               "n_shared_lineages", "n_accessions_in_shared_lineages",
               "n_genets_current", "n_genets_percountry"),
  value = c(length(single_mll), length(idx_single),
            length(shared), sum(memb$mll %in% as.integer(shared)),
            length(idx_current), length(idx_percountry))),
  "outputs/brad01_ibd_singlecountry.csv", row.names = FALSE)

## Table S1 / Supplementary Data 8: diversity by K = 6 ancestry cluster on the 392 genets
gl_mll <- readRDS("input_files/taro_diploid_mll.rds")
M_mll  <- as.matrix(gl_mll)
country_mll <- as.character(pop(gl_mll))
fam <- read.table("outputs/ohana_declone/gl_plink.fam", stringsAsFactors = FALSE)
Q   <- as.matrix(read.table("outputs/ohana_declone/q.k6.matrix", skip = 1))
clust_by_id <- setNames(paste0("C", max.col(Q, ties.method = "first")), fam$V2)
cluster <- unname(clust_by_id[indNames(gl_mll)])
by_cluster <- build_table(M_mll, cluster)
write.csv(by_cluster, "outputs/diversity_by_cluster.csv", row.names = FALSE)

## Table S2: bootstrap 95% CIs for He, Ho, FIS and Ar by country (1,000 resamples of genets)
set.seed(42)
gl_pc <- gl.keep.ind(gl, ind.list = indNames(gl)[idx_percountry], verbose = 0)
cmap  <- setNames(memb$country, memb$id)
pop(gl_pc) <- factor(unname(cmap[indNames(gl_pc)]))   # gl.keep.ind keeps genlight order: look up by name
M_pc <- as.matrix(gl_pc)
country_pc <- as.character(pop(gl_pc))
B <- 1000L
keep_pc   <- names(which(table(country_pc) >= MIN_N))
k_rare_pc <- 2L * min(table(country_pc)[keep_pc])
stats_of <- function(M, idx, k_rare) {
  m <- M[idx, , drop = FALSE]; nn <- colSums(!is.na(m))
  het <- colSums(m == 1L, na.rm = TRUE); alt <- colSums(m, na.rm = TRUE); twoN <- 2 * nn
  q  <- ifelse(nn > 0, alt / twoN, NA)
  Ho <- mean(ifelse(nn > 0, het / nn, NA), na.rm = TRUE)
  He <- mean(ifelse(nn > 1, (twoN / (twoN - 1)) * 2 * q * (1 - q), NA), na.rm = TRUE)
  Ar <- mean(ar_locus_vec(twoN - alt, alt, twoN, k_rare), na.rm = TRUE)
  c(Ho = Ho, He = He, FIS = 1 - Ho / He, Ar = Ar)
}
rows <- lapply(keep_pc, function(g) {
  idx <- which(country_pc == g); n <- length(idx)
  pt <- stats_of(M_pc, idx, k_rare_pc)
  bt <- replicate(B, stats_of(M_pc, sample(idx, n, replace = TRUE), k_rare_pc))
  ci <- apply(bt, 1, quantile, c(.025, .975), na.rm = TRUE)
  data.frame(country = g, N = n,
             He = round(pt["He"], 4),   He_lo = round(ci[1, "He"], 4),   He_hi = round(ci[2, "He"], 4),
             Ho = round(pt["Ho"], 4),   Ho_lo = round(ci[1, "Ho"], 4),   Ho_hi = round(ci[2, "Ho"], 4),
             FIS = round(pt["FIS"], 4), FIS_lo = round(ci[1, "FIS"], 4), FIS_hi = round(ci[2, "FIS"], 4),
             Ar = round(pt["Ar"], 4),   Ar_lo = round(ci[1, "Ar"], 4),   Ar_hi = round(ci[2, "Ar"], 4))
})
ci_tab <- do.call(rbind, rows)
ci_tab <- ci_tab[order(-ci_tab$He), ]
write.csv(ci_tab, "outputs/brad01_diversity_ci.csv", row.names = FALSE)

## SI S4.4: raw private-allele count vs N (Spearman rho = 0.34, P = 0.29) on the 392 genets
Pres_mll <- pres_matrix(M_mll)
keep_mll <- names(which(table(country_mll) >= MIN_N))
priv_raw <- function(groups, levs) {
  GP <- vapply(levs, function(g) colSums(Pres_mll[groups == g, , drop = FALSE]) > 0,
               logical(ncol(Pres_mll)))
  owner_ok <- rowSums(GP) == 1L
  sapply(levs, function(g) sum(owner_ok & GP[, g]))
}
raw <- priv_raw(country_mll, keep_mll)
Ns  <- table(country_mll)[keep_mll]
priv_vs_N <- cor.test(raw, as.numeric(Ns), method = "spearman")

## SI S4.4: rarefied private allelic richness, balanced panel of K = 10 genets, 500 resamples
set.seed(42)
K <- 10L; R <- 500L
pool <- names(which(table(country_mll) >= K))
acc <- matrix(0, nrow = R, ncol = length(pool), dimnames = list(NULL, pool))
for (r in seq_len(R)) {
  idx <- unlist(lapply(pool, function(g) sample(which(country_mll == g), K)))
  grp <- rep(pool, each = K)
  GP  <- vapply(pool, function(g) colSums(Pres_mll[idx[grp == g], , drop = FALSE]) > 0,
                logical(ncol(Pres_mll)))
  only1 <- rowSums(GP) == 1L
  acc[r, ] <- sapply(pool, function(g) sum(only1 & GP[, g]))
}
rare <- data.frame(country = pool,
                   private_raw = raw[pool],
                   N = as.integer(Ns[pool]),
                   private_rarefied_K10 = round(colMeans(acc), 1),
                   rare_lo = round(apply(acc, 2, quantile, .025), 1),
                   rare_hi = round(apply(acc, 2, quantile, .975), 1))
rare <- rare[order(-rare$private_rarefied_K10), ]
write.csv(rare, "outputs/private_rarefied.csv", row.names = FALSE)

## Fig. 4: expected heterozygosity and FIS by country (figures/diversity_by_country.png)
region <- c(
  Indonesia = "SE/E Asia", Malaysia = "SE/E Asia", Thailand = "SE/E Asia",
  Vietnam = "SE/E Asia", Philippines = "SE/E Asia", Japan = "SE/E Asia",
  IndiaNorthen = "SE/E Asia",
  PapuaNewGuinea = "Melanesia", SolomonIslands = "Melanesia",
  Vanuatu = "Melanesia", NewCaledonia = "Melanesia", Fiji = "Melanesia",
  Samoa = "Polynesia", WSamoa = "Polynesia", Tonga = "Polynesia",
  Niue = "Polynesia", CookIslands = "Polynesia",
  "Tahiti&FrenchPolynesia" = "Polynesia", Hawaii = "Polynesia",
  Palau = "Micronesia", FederatedStatesofMicronesia = "Micronesia",
  Australia = "Australia", UnitedStatesofAmerica = "Other")
reg_cols <- c("SE/E Asia" = "#E69F00", "Melanesia" = "#009E73",
              "Polynesia" = "#0072B2", "Micronesia" = "#CC79A7",
              "Australia" = "#9E4A06")
pc <- by_country_pc[order(-by_country_pc$He), ]
pc$region <- factor(region[pc$group],
                    levels = c("SE/E Asia", "Melanesia", "Polynesia", "Micronesia", "Australia"))
pc$label <- pretty_country(pc$group, short = TRUE)
pA <- ggplot(pc, aes(x = He, y = reorder(label, He))) +
  geom_segment(aes(x = 0, xend = He, yend = reorder(label, He)),
               colour = "grey80", linewidth = 0.6) +
  geom_point(aes(colour = region, size = N)) +
  scale_colour_manual(values = reg_cols, name = "Region") +
  scale_size_area(max_size = 8, name = "Genets (N)") +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  labs(x = "Unbiased expected heterozygosity (He)", y = NULL, tag = "A") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        panel.grid.minor = element_blank())
pB <- ggplot(pc, aes(x = FIS, y = reorder(label, He))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_segment(aes(x = 0, xend = FIS, yend = reorder(label, He)),
               colour = "grey80", linewidth = 0.6) +
  geom_point(aes(colour = region, size = N)) +
  scale_colour_manual(values = reg_cols, name = "Region", guide = "none") +
  scale_size_area(max_size = 8, guide = "none") +
  labs(x = expression("Inbreeding coefficient  " * F[IS]), y = NULL, tag = "B") +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        axis.text.y = element_blank(), panel.grid.minor = element_blank())
g <- arrangeGrob(pA, pB, ncol = 2, widths = c(1.25, 1))
ggsave("figures/diversity_by_country.png", g, width = 11, height = 5.6, dpi = 300, bg = "white")

## Table 1 and Table S1 as Markdown text (outputs/ms_tables.md)
clon <- read.csv("outputs/mlg_diversity_indices.csv", stringsAsFactors = FALSE)
m <- merge(by_country_pc, clon[, c("level", "N", "G_MLL", "R")],
           by.x = "group", by.y = "level", suffixes = c("_genet", "_ram"))
m <- m[order(-m$He), ]
con <- file("outputs/ms_tables.md", "w")
w <- function(...) cat(..., "\n", file = con, sep = "")
w("## Table 1. Per-country diversity of the diploid taro collection\n")
w("Accessions (ramets), distinct genets (G), clonal richness R, and gene-diversity")
w("statistics on one ramet of every lineage recorded in the country: percent")
w("polymorphic loci, observed and unbiased expected heterozygosity (Ho, He), FIS,")
w("rarefied allelic richness (Ar, to 10 gene copies) and private alleles.")
w("Countries with >= 5 genets.\n")
w("| Country | Accessions | Genets | R | Poly % | Ho | He | FIS | Ar | Private |")
w("|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
for (i in seq_len(nrow(m))) {
  r <- m[i, ]
  w(sprintf("| %s | %d | %d | %.2f | %.0f | %.3f | %.3f | %.3f | %.3f | %d |",
            pretty_country(r$group), r$N_ram, r$G_MLL, r$R, r$poly_pct,
            r$Ho, r$He, r$FIS, r$Ar, r$private_alleles))
}
w("\n## Table S. Per-K=6-cluster diversity (clone-censored genets)\n")
w("| Cluster | Genets | Poly % | Ho | He | FIS | Ar | Private |")
w("|---|--:|--:|--:|--:|--:|--:|--:|")
for (i in seq_len(nrow(by_cluster))) {
  r <- by_cluster[i, ]
  w(sprintf("| %s | %d | %.0f | %.3f | %.3f | %.3f | %.3f | %d |",
            r$group, r$N, r$poly_pct, r$Ho, r$He, r$FIS, r$Ar, r$private_alleles))
}
close(con)

## Results: He declines with distance from Indonesia (rho = -0.61, P = 0.010) and with longitude (rho = -0.59, P = 0.012)
# per-country centroids (one lon/lat per country in ind.metrics), great-circle distance in km
im <- gl@other$ind.metrics
cc_lon <- tapply(im$lon.y, country_full, function(x) x[1])
cc_lat <- tapply(im$lat.y, country_full, function(x) x[1])
haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; d2r <- pi / 180
  dlat <- (lat2 - lat1) * d2r; dlon <- (lon2 - lon1) * d2r
  a <- sin(dlat / 2)^2 + cos(lat1 * d2r) * cos(lat2 * d2r) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}
div <- by_country_pc
lon_east <- ifelse(cc_lon < 0, cc_lon + 360, cc_lon)   # dateline-unwrapped longitude (degrees east)
div$lon <- lon_east[div$group]
div$lat <- cc_lat[div$group]
div <- div[order(div$lon), ]
sp_he <- cor.test(div$He, div$lon, method = "spearman")
sp_ar <- cor.test(div$Ar, div$lon, method = "spearman")
lm_he <- lm(He ~ lon, data = div)
src <- "Indonesia"
div$dist_src <- haversine(cc_lat[src], cc_lon[src], div$lat,
                          ifelse(div$lon > 180, div$lon - 360, div$lon))
sp_ds <- cor.test(div$He, div$dist_src, method = "spearman")
write.csv(data.frame(
  metric = c("n_countries_reg", "He_lon_spearman_rho", "He_lon_spearman_p",
             "Ar_lon_spearman_rho", "Ar_lon_spearman_p",
             "He_lon_ols_slope", "He_lon_ols_R2", "He_lon_ols_p",
             "He_distsrc_rho", "He_distsrc_p"),
  value = signif(c(nrow(div), sp_he$estimate, sp_he$p.value,
                   sp_ar$estimate, sp_ar$p.value,
                   coef(lm_he)[2], summary(lm_he)$r.squared, summary(lm_he)$coefficients[2, 4],
                   sp_ds$estimate, sp_ds$p.value), 4)),
  "outputs/cline_test.csv", row.names = FALSE)

## Results: pairwise regional FST on the 392 genets (SE/E Asia vs Melanesia 0.17; Pacific regions 0.04 to 0.10)
set.seed(42)
region_fst <- c(
  Indonesia = "SEEAsia", Malaysia = "SEEAsia", Thailand = "SEEAsia", Vietnam = "SEEAsia",
  Philippines = "SEEAsia", Japan = "SEEAsia", IndiaNorthen = "SEEAsia",
  PapuaNewGuinea = "Melanesia", SolomonIslands = "Melanesia", Vanuatu = "Melanesia",
  NewCaledonia = "Melanesia", Fiji = "Melanesia",
  Samoa = "Polynesia", WSamoa = "Polynesia", Tonga = "Polynesia", Niue = "Polynesia",
  CookIslands = "Polynesia", "Tahiti&FrenchPolynesia" = "Polynesia", Hawaii = "Polynesia",
  Palau = "Micronesia", FederatedStatesofMicronesia = "Micronesia",
  Australia = "Other", UnitedStatesofAmerica = "Other")
reg <- region_fst[country_mll]
keep_reg <- !is.na(reg) & reg != "Other"
gl_reg <- gl_mll[keep_reg, ]
pop(gl_reg) <- factor(reg[keep_reg], levels = c("SEEAsia", "Melanesia", "Polynesia", "Micronesia"))
fst <- gl.fst.pop(gl_reg, nboots = 1000, percent = 95, nclusters = 1, verbose = 0)
write.csv(fst$Bootstraps, "outputs/fst_regions.csv", row.names = FALSE)

## Table S5: isolation by distance under four designs (design A is the main-text Mantel r = 0.37, P = 0.0002)
set.seed(42)
keymap <- c(TahitiFrenchPolynesia = "Tahiti&FrenchPolynesia")   # popNames stripped of "&" back to centroid names
# per-country F81 distance matrix of a genlight whose pop() is country
f81_by_country <- function(g) {
  g <- gl.filter.monomorphs(g, verbose = 0)
  popNames(g) <- gsub("[^A-Za-z0-9]", "", popNames(g))
  g@other$loc.metrics$TrimmedSequence <- g@other$loc.metrics$TrimmedSequenceSnp
  as.matrix(gl.dist.phylo(g))
}
# Mantel test of F81 vs great-circle distance between country centroids, 9,999 permutations
mantel_geo <- function(Dg, B = 9999) {
  rn <- rownames(Dg); rn2 <- ifelse(rn %in% names(keymap), keymap[rn], rn)
  have <- rn2 %in% names(cc_lon)
  Dg <- Dg[have, have]; ctys <- rn2[have]
  lat <- cc_lat[ctys]; lon <- cc_lon[ctys]
  Geo <- outer(seq_along(ctys), seq_along(ctys),
               Vectorize(function(i, j) haversine(lat[i], lon[i], lat[j], lon[j])))
  gd <- Dg[lower.tri(Dg)]; ge <- Geo[lower.tri(Geo)]
  r <- cor(gd, ge); n <- nrow(Dg)
  perm <- replicate(B, { p <- sample(n); Dp <- Dg[p, p]
                         cor(Dp[lower.tri(Dp)], ge) })
  list(r = r, p = (1 + sum(perm >= r)) / (B + 1), n = n)
}
gl_of <- function(idx) {   # idx ascending, so pop can be assigned positionally
  g <- gl.keep.ind(gl, ind.list = indNames(gl)[idx], verbose = 0)
  pop(g) <- factor(memb$country[idx])
  g
}
idx_all <- idx_current
tab_s <- table(memb$country[idx_single])
cty_D <- names(tab_s)[tab_s >= MIN_N]                    # countries with >= 5 single-country genets
ibd_A <- mantel_geo(f81_by_country(gl_of(idx_all)))
ibd_B <- mantel_geo(f81_by_country(gl_of(idx_single)))
ibd_C <- mantel_geo(f81_by_country(gl.keep.pop(gl_of(idx_all), pop.list = cty_D, verbose = 0)))
ibd_D <- mantel_geo(f81_by_country(gl.keep.pop(gl_of(idx_single), pop.list = cty_D, verbose = 0)))
write.csv(data.frame(
  design = c("A all lineages / all countries",
             "B single-country lineages / all countries",
             "C all lineages / countries of D",
             "D single-country lineages / countries with >= 5 genets"),
  genets = c(length(idx_all), length(idx_single), length(idx_all), length(idx_single)),
  countries = c(ibd_A$n, ibd_B$n, ibd_C$n, ibd_D$n),
  mantel_r = signif(c(ibd_A$r, ibd_B$r, ibd_C$r, ibd_D$r), 4),
  mantel_p = signif(c(ibd_A$p, ibd_B$p, ibd_C$p, ibd_D$p), 4)),
  "outputs/brad04_ibd_matched.csv", row.names = FALSE)

## Table S4: sensitivity of diversity ranking, PCA axes, regional FST and Mantel test to the ploidy threshold
set.seed(42)
CUT      <- 0.008     # main-text lineage cut height (complete linkage)
THR_MAIN <- 0.0923
THR_CONS <- 0.06
gl_all   <- readRDS("input_files/taro.rds")
ho_calls <- read.csv("outputs/ploidy_individual_Ho.csv", stringsAsFactors = FALSE)
# diploid working set: flow-cytometry cytotype where available, individual Ho otherwise
diploid_ids <- function(thr) {
  ids <- ho_calls$id[(!is.na(ho_calls$class) & ho_calls$class == "Diploid") |
                     (is.na(ho_calls$class) & ho_calls$Ho < thr)]
  intersect(ids, indNames(gl_all))
}
ids_main <- diploid_ids(THR_MAIN)
ids_cons <- diploid_ids(THR_CONS)
Dbit <- as.matrix(readRDS("outputs/mlg_bitwise_dist.rds"))
ctry_all <- setNames(as.character(pop(gl)), indNames(gl))
# lineages recut at the main-text height inside each working set; highest-call-rate ramet per lineage
reps_for <- function(ids) {
  hc  <- hclust(as.dist(Dbit[ids, ids, drop = FALSE]), method = "complete")
  mll <- cutree(hc, h = CUT)
  cr  <- rowSums(!is.na(Mfull[ids, , drop = FALSE]))
  unname(unlist(tapply(seq_along(ids), mll, function(ix) ids[ix][which.max(cr[ix])])))
}
reps_main <- reps_for(ids_main)
reps_cons <- reps_for(ids_cons)
div_main <- build_table(Mfull[reps_main, , drop = FALSE], ctry_all[reps_main])
div_cons <- build_table(Mfull[reps_cons, , drop = FALSE], ctry_all[reps_cons])
cmp <- merge(div_main, div_cons, by = "group", suffixes = c("_main", "_cons"))
cmp <- cmp[order(-cmp$He_main), ]
cmp$dHe <- round(cmp$He_cons - cmp$He_main, 4)
cmp$dAr <- round(cmp$Ar_cons - cmp$Ar_main, 4)
write.csv(cmp, "outputs/brad02_threshold_diversity.csv", row.names = FALSE)
common <- intersect(div_main$group, div_cons$group)
r_He <- cor(div_main$He[match(common, div_main$group)],
            div_cons$He[match(common, div_cons$group)], method = "spearman")
r_Ar <- cor(div_main$Ar[match(common, div_main$group)],
            div_cons$Ar[match(common, div_cons$group)], method = "spearman")
# PCA: correlation of PC scores over the genets shared by the two working sets
mk_gl <- function(ids) {
  g <- gl.keep.ind(gl, ind.list = ids, verbose = 0)
  pop(g) <- factor(unname(ctry_all[indNames(g)]))       # gl.keep.ind keeps genlight order: look up by name
  g
}
pcoa_scores <- function(g, nf = 3) {
  p <- gl.pcoa(g, nfactors = nf, verbose = 0)
  s <- p$scores[, seq_len(nf), drop = FALSE]
  rownames(s) <- indNames(g); s
}
sc_main <- pcoa_scores(mk_gl(reps_main))
sc_cons <- pcoa_scores(mk_gl(reps_cons))
shared_ind <- intersect(rownames(sc_main), rownames(sc_cons))
pc_r <- vapply(1:3, function(k) abs(cor(sc_main[shared_ind, k], sc_cons[shared_ind, k])), numeric(1))
# pairwise regional FST under each threshold
fst_regional <- function(ids) {
  g <- mk_gl(ids)
  rg <- region_fst[as.character(pop(g))]
  keep <- !is.na(rg) & rg != "Other"
  g2 <- g[keep, ]
  pop(g2) <- factor(rg[keep], levels = c("SEEAsia", "Melanesia", "Polynesia", "Micronesia"))
  gl.fst.pop(g2, nboots = 100, percent = 95, nclusters = 1, verbose = 0)$Fst
}
F_main <- fst_regional(reps_main)
F_cons <- fst_regional(reps_cons)
fm <- as.numeric(F_main[lower.tri(F_main)]); fc <- as.numeric(F_cons[lower.tri(F_cons)])
ok <- is.finite(fm) & is.finite(fc)
# isolation by distance under each threshold, countries with >= MIN_N genets
mantel_for <- function(ids) {
  g  <- mk_gl(ids)
  tb <- table(pop(g))
  g  <- gl.keep.pop(g, pop.list = names(tb)[tb >= MIN_N], verbose = 0)
  mantel_geo(f81_by_country(g))
}
mt_main <- mantel_for(reps_main)
mt_cons <- mantel_for(reps_cons)
summ <- data.frame(
  result = c("Individuals retained", "Genets", "Countries with N >= 5",
             "Per-country He, rank correlation with main",
             "Per-country Ar, rank correlation with main",
             "PCA PC1 |r| with main", "PCA PC2 |r| with main", "PCA PC3 |r| with main",
             "Regional FST, correlation with main",
             "Regional FST, largest absolute difference",
             "Mantel r", "Mantel P", "Mantel countries"),
  main = c(length(ids_main), length(reps_main), nrow(div_main),
           NA, NA, NA, NA, NA, NA, NA,
           signif(mt_main$r, 3), signif(mt_main$p, 4), mt_main$n),
  conservative = c(length(ids_cons), length(reps_cons), nrow(div_cons),
                   round(r_He, 3), round(r_Ar, 3),
                   round(pc_r[1], 3), round(pc_r[2], 3), round(pc_r[3], 3),
                   round(cor(fm[ok], fc[ok]), 3), round(max(abs(fm[ok] - fc[ok])), 4),
                   signif(mt_cons$r, 3), signif(mt_cons$p, 4), mt_cons$n))
write.csv(summ, "outputs/brad02_threshold_summary.csv", row.names = FALSE)

# ============================================================================
# 10. Selection scan and gene context (Fig. S5; SI S3)
# ============================================================================

# OHANA matrix helpers (first line holds the dimensions)
read_matrix <- function(path, integer_values = FALSE) {
  con <- file(path, "r")
  on.exit(close(con))
  dims <- scan(con, what = integer(), n = 2L, quiet = TRUE)
  values <- scan(con, what = if (integer_values) integer() else numeric(), quiet = TRUE)
  matrix(values, nrow = dims[1L], ncol = dims[2L], byrow = TRUE)
}
write_matrix <- function(x, path) {
  writeLines(sprintf("%d %d", nrow(x), ncol(x)), path)
  write.table(x, path, append = TRUE, row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = "\t")
}

# ---- Consensus genotype per clonal lineage (sensitivity panel) ----
## Methods: consensus panel, modal called genotype per lineage, ties set to missing
g <- read_matrix("outputs/ohana/g.dgm", integer_values = TRUE)
bim <- read.table("outputs/ohana/gl_plink.bim", stringsAsFactors = FALSE)
fam <- read.table("outputs/ohana/gl_plink.fam", stringsAsFactors = FALSE)
membership <- read.csv("outputs/mlg_membership.csv", stringsAsFactors = FALSE)
membership <- membership[match(fam[[2L]], membership$id), ]
groups <- split(seq_len(nrow(g)), membership$mll)

consensus_one <- function(index) {
  if (length(index) == 1L) return(g[index, ])
  x <- g[index, , drop = FALSE]
  counts <- rbind(colSums(x == 0L), colSums(x == 1L), colSums(x == 2L))
  maxima <- apply(counts, 2L, max)
  tied <- colSums(counts == rep(maxima, each = 3L)) > 1L
  value <- max.col(t(counts), ties.method = "first") - 1L
  value[maxima == 0L | tied] <- 3L
  value
}
consensus <- do.call(rbind, lapply(groups, consensus_one))
storage.mode(consensus) <- "integer"
lineages <- do.call(rbind, lapply(seq_along(groups), function(i) {
  index <- groups[[i]]
  data.frame(mll = as.integer(names(groups)[i]), n_ramets = length(index),
             representative_id = membership$id[index][which.max(membership$call_rate[index])],
             stringsAsFactors = FALSE)
}))
write_matrix(consensus, "outputs/ohana_consensus/g.dgm")
write.table(bim, "outputs/ohana_consensus/markers.bim", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(lineages, "outputs/ohana_consensus/lineages.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
rm(g, bim, fam, membership, groups, consensus, lineages)

# ---- Component-specific selection scan: one function, three panels ----
k <- 6L
n_starts <- 30L
max_iterations <- 500L
tolerance <- 1e-4
min_call_rate <- 0.80
min_maf <- 0.05
scale_increment <- 10
alpha <- 0.05
workers <- 6L

seed_name <- function(seed) {
  paste0("seed", if (seed <= 99) sprintf("%02.0f", seed) else sprintf("%.0f", seed))
}
parse_qpas_log <- function(path) {
  rows <- grep("^[[:space:]]*[0-9]+[[:space:]]", readLines(path, warn = FALSE), value = TRUE)
  fields <- strsplit(trimws(tail(rows, 1L)), "[[:space:]]+")[[1L]]
  list(iteration = as.integer(fields[1L]), loglik = as.numeric(fields[3L]),
       delta = as.numeric(fields[4L]))
}

run_component_scan <- function(panel, dgm, bim, out, seeds) {
  run_dir <- file.path(out, "multistart")
  best_dir <- file.path(out, "best")
  scan_dir <- file.path(out, "component_scans")

  # Selection panel: call rate >= 0.80 and MAF >= 0.05 recomputed within the panel
  g_all <- read_matrix(dgm, integer_values = TRUE)
  bim_all <- read.table(bim, stringsAsFactors = FALSE)
  called <- g_all != 3L
  n_called <- colSums(called)
  allele_frequency <- colSums(ifelse(called, g_all, 0L)) / (2 * n_called)
  maf <- pmin(allele_frequency, 1 - allele_frequency)
  call_rate <- n_called / nrow(g_all)
  keep <- is.finite(maf) & maf >= min_maf & call_rate >= min_call_rate
  g <- g_all[, keep, drop = FALSE]
  markers <- data.frame(locus = bim_all[keep, 2L], chr = bim_all[keep, 1L],
                        pos = bim_all[keep, 4L], call_rate = call_rate[keep],
                        maf = maf[keep], stringsAsFactors = FALSE)
  selection_dgm <- file.path(out, "selection.g.dgm")
  write_matrix(g, selection_dgm)
  write.table(markers, file.path(out, "selection_markers.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
  rm(g_all, called)

  # 30 qpas starts at K = 6; keep the run with the highest log likelihood
  run_one <- function(seed) {
    stem <- file.path(run_dir, seed_name(seed))
    status <- system2("qpas", c(selection_dgm, "-k", k, "-s", format(seed, scientific = FALSE),
                                "-qo", paste0(stem, ".q.matrix"), "-fo", paste0(stem, ".f.matrix"),
                                "-mi", max_iterations, "-e", format(tolerance, scientific = TRUE)),
                      stdout = paste0(stem, ".log"), stderr = paste0(stem, ".log"))
    parsed <- parse_qpas_log(paste0(stem, ".log"))
    data.frame(seed = seed, loglik = parsed$loglik, final_iteration = parsed$iteration,
               final_delta = parsed$delta, reached_limit = parsed$iteration >= max_iterations,
               status = status)
  }
  runs <- do.call(rbind, parallel::mclapply(seeds, run_one, mc.cores = workers, mc.preschedule = FALSE))
  runs <- runs[order(runs$seed), ]
  write.csv(runs, file.path(out, "qpas_all_runs.csv"), row.names = FALSE)
  valid <- runs[runs$status == 0L & is.finite(runs$loglik), ]
  best <- valid[which.max(valid$loglik), , drop = FALSE]
  best$K <- k
  best$n_starts <- n_starts
  best$n_loci <- ncol(g)
  write.csv(best, file.path(out, "qpas_best_run.csv"), row.names = FALSE)
  best_stem <- file.path(run_dir, seed_name(best$seed))
  best_f <- file.path(best_dir, "f.k6.matrix")
  file.copy(paste0(best_stem, ".q.matrix"), file.path(best_dir, "q.k6.matrix"), overwrite = TRUE)
  file.copy(paste0(best_stem, ".f.matrix"), best_f, overwrite = TRUE)
  file.copy(paste0(best_stem, ".log"), file.path(best_dir, "qpas.k6.log"), overwrite = TRUE)

  # Genome-wide component covariance (nemeco) and component tree
  c_path <- file.path(out, "c.matrix")
  system2("nemeco", c(selection_dgm, best_f, "-co", c_path, "-mi", max_iterations),
          stdout = file.path(out, "nemeco.log"), stderr = file.path(out, "nemeco.log"))
  system2("convert", c("cov2nwk", c_path, file.path(out, "tree.nwk")))

  # Six component-specific selscan tests: +h on all cells targets the root component 0,
  # +h on one diagonal cell targets components 1..K-1
  c_global <- read_matrix(c_path)
  scan_one <- function(component) {
    c_scale <- c_global
    if (component == 0L) {
      c_scale <- c_scale + scale_increment
    } else {
      c_scale[component, component] <- c_scale[component, component] + scale_increment
    }
    scale_path <- file.path(scan_dir, sprintf("c_scale.component%d.matrix", component))
    scan_path <- file.path(scan_dir, sprintf("selscan.component%d.txt", component))
    write_matrix(c_scale, scale_path)
    system2("selscan", c(selection_dgm, best_f, c_path, "-cs", scale_path),
            stdout = scan_path, stderr = scan_path)
    z <- read.table(scan_path, header = TRUE, check.names = FALSE)
    data.frame(panel = panel, component = component, markers, lrt = z[["lle-ratio"]],
               global_loglik = z[["global-lle"]], local_loglik = z[["local-lle"]],
               stringsAsFactors = FALSE)
  }
  selection <- do.call(rbind, lapply(0:(k - 1L), scan_one))

  # 50:50 point-mass/chi-square(1) null; one Bonferroni correction across loci x components
  selection$p_value <- ifelse(selection$lrt <= 0, 1,
                              0.5 * pchisq(selection$lrt, df = 1, lower.tail = FALSE))
  total_tests <- nrow(selection)
  selection$p_bonferroni <- pmin(1, selection$p_value * total_tests)
  selection$minus_log10_p <- -log10(pmax(selection$p_value, .Machine$double.xmin))
  selection$significant <- selection$p_bonferroni <= alpha
  selection <- selection[order(selection$component, selection$chr, selection$pos), ]
  write.csv(selection, file.path(out, "component_selection.csv"), row.names = FALSE)

  theory_q90 <- qchisq(0.8, df = 1)
  summary_rows <- do.call(rbind, lapply(split(selection, selection$component), function(z) {
    top <- z[which.max(z$lrt), ]
    data.frame(panel = panel, component = top$component, n_loci = nrow(z),
               lambda_90 = as.numeric(quantile(z$lrt, 0.90, na.rm = TRUE)) / theory_q90,
               n_significant = sum(z$significant, na.rm = TRUE),
               top_locus = top$locus, top_chr = top$chr, top_pos = top$pos,
               top_lrt = top$lrt, top_p = top$p_value, top_p_bonferroni = top$p_bonferroni,
               stringsAsFactors = FALSE)
  }))
  write.csv(summary_rows, file.path(out, "component_selection_summary.csv"), row.names = FALSE)

  # Six-panel Manhattan plot, marker order taken from the BIM file
  plot_dat <- selection
  plot_dat$component_label <- factor(plot_dat$component, levels = 0:(k - 1L),
                                     labels = paste("Ancestry component", seq_len(k)))
  plot_dat$chr <- factor(plot_dat$chr, levels = unique(markers$chr))
  chr_info <- aggregate(pos ~ chr, markers, max)
  chr_info$offset <- c(0, head(cumsum(chr_info$pos), -1L))
  offsets <- setNames(chr_info$offset, chr_info$chr)
  plot_dat$x <- plot_dat$pos + offsets[as.character(plot_dat$chr)]
  axis_pos <- chr_info$offset + chr_info$pos / 2
  bonf_line <- -log10(alpha / total_tests)
  p <- ggplot(plot_dat, aes(x = x, y = minus_log10_p, colour = chr)) +
    geom_point(size = 0.45, alpha = 0.75) +
    geom_hline(yintercept = bonf_line, linetype = 2, colour = "#B2182B") +
    facet_wrap(~ component_label, ncol = 2) +
    scale_x_continuous(breaks = axis_pos, labels = chr_info$chr) +
    scale_colour_manual(values = rep(c("#2C7FB8", "#7FCDBB"), length.out = nrow(chr_info))) +
    labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(sprintf("figures/ohana_selection_components_%s.png", panel), p,
         width = 11, height = 9, dpi = 300, bg = "white")
  invisible(selection)
}

## Fig. S5: OHANA selection scans for the six ancestry components (392 genets)
## Results: three candidates after Bonferroni across 19,044 tests (3,174 loci x 6 components)
run_component_scan("genets", "outputs/ohana_declone/g.dgm", "outputs/ohana_declone/gl_plink.bim",
                   "outputs/ohana_selection_genets", c(0:28, 3011845641))
## Methods: consensus-genotype panel (3,115 loci)
run_component_scan("consensus", "outputs/ohana_consensus/g.dgm", "outputs/ohana_consensus/markers.bim",
                   "outputs/ohana_selection_consensus", 0:29)
## Methods / SI S3: all 1,064 accessions (3,690 loci)
run_component_scan("accessions", "outputs/ohana/g.dgm", "outputs/ohana/gl_plink.bim",
                   "outputs/ohana_selection_accessions", c(0:28, 3915235949))

# ---- Align components across panels and compare the candidates ----
## SI S3: matched-component allele frequency correlations (r = 0.967-0.986; accession best match 0.249)
## SI S3 table: adjusted P of the candidates in the genet and consensus panels
panels <- c("genets", "consensus", "accessions")
reference_panel <- "genets"
read_panel <- function(panel) {
  path <- file.path("outputs", paste0("ohana_selection_", panel))
  markers <- read.delim(file.path(path, "selection_markers.tsv"), stringsAsFactors = FALSE)
  frequencies <- t(read_matrix(file.path(path, "best", "f.k6.matrix")))
  rownames(frequencies) <- markers$locus
  selection <- read.csv(file.path(path, "component_selection.csv"), stringsAsFactors = FALSE)
  list(markers = markers, frequencies = frequencies, selection = selection)
}
all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], all_permutations(x[-i]))))
}
data <- setNames(lapply(panels, read_panel), panels)
permutations <- all_permutations(seq_len(k))

# One-to-one component matching that maximises the summed allele frequency correlation
align_one <- function(comparison_panel) {
  reference <- data[[reference_panel]]$frequencies
  comparison <- data[[comparison_panel]]$frequencies
  shared <- intersect(rownames(reference), rownames(comparison))
  correlation <- cor(reference[shared, , drop = FALSE], comparison[shared, , drop = FALSE])
  scores <- apply(permutations, 1L, function(p) sum(correlation[cbind(seq_len(k), p)]))
  best <- permutations[which.max(scores), ]
  do.call(rbind, lapply(seq_len(k), function(component) {
    reference_top <- shared[order(reference[shared, component], decreasing = TRUE)[1:50]]
    comparison_top <- shared[order(comparison[shared, best[component]], decreasing = TRUE)[1:50]]
    data.frame(reference_panel = reference_panel, reference_component = component - 1L,
               comparison_panel = comparison_panel, comparison_component = best[component] - 1L,
               allele_frequency_correlation = correlation[component, best[component]],
               n_shared_loci = length(shared),
               top_50_overlap = length(intersect(reference_top, comparison_top)),
               stringsAsFactors = FALSE)
  }))
}
alignment <- do.call(rbind, lapply(setdiff(panels, reference_panel), align_one))
write.csv(alignment, "outputs/ohana_selection_component_alignment.csv", row.names = FALSE)

# Every locus significant in any panel, tested in the matched component of the other panels
significant <- do.call(rbind, lapply(panels, function(panel) {
  z <- data[[panel]]$selection
  z <- z[z$significant, c("panel", "component", "locus", "chr", "pos")]
  if (panel == reference_panel) {
    z$reference_component <- z$component
  } else {
    map <- alignment[alignment$comparison_panel == panel, ]
    z$reference_component <- map$reference_component[match(z$component, map$comparison_component)]
  }
  z
}))
candidate_keys <- unique(significant[c("locus", "chr", "pos", "reference_component")])
candidate_comparison <- do.call(rbind, lapply(seq_len(nrow(candidate_keys)), function(i) {
  key <- candidate_keys[i, ]
  do.call(rbind, lapply(panels, function(panel) {
    component <- if (panel == reference_panel) {
      key$reference_component
    } else {
      map <- alignment[alignment$comparison_panel == panel &
                         alignment$reference_component == key$reference_component, ]
      map$comparison_component
    }
    z <- data[[panel]]$selection
    z <- z[z$locus == key$locus & z$component == component, ]
    if (!nrow(z)) {
      return(data.frame(locus = key$locus, chr = key$chr, pos = key$pos,
                        reference_component = key$reference_component, panel = panel,
                        panel_component = component, passed_panel_filters = FALSE,
                        lrt = NA_real_, p_value = NA_real_, p_bonferroni = NA_real_,
                        significant = FALSE, stringsAsFactors = FALSE))
    }
    data.frame(locus = key$locus, chr = key$chr, pos = key$pos,
               reference_component = key$reference_component, panel = panel,
               panel_component = component, passed_panel_filters = TRUE,
               lrt = z$lrt, p_value = z$p_value, p_bonferroni = z$p_bonferroni,
               significant = z$significant, stringsAsFactors = FALSE)
  }))
}))
write.csv(candidate_comparison, "outputs/ohana_selection_candidate_comparison.csv", row.names = FALSE)

## SI S3: mean ancestry of the candidate component by country (Philippines 0.64, Palau 0.73)
fam_g <- read.table("outputs/ohana_declone/gl_plink.fam", stringsAsFactors = FALSE)
q_g <- read_matrix("outputs/ohana_selection_genets/best/q.k6.matrix")
sel_g <- data$genets$selection
comp_g <- unique(sel_g$component[sel_g$significant]) + 1L
ancestry_g <- aggregate(q_g[, comp_g, drop = FALSE], list(country = fam_g$V1), mean)
write.csv(ancestry_g, "outputs/ohana_selection_genets/candidate_component_ancestry.csv", row.names = FALSE)

## SI S3: the accession-panel component holding the five candidates is one clonal lineage of 23 ramets
fam_a <- read.table("outputs/ohana/gl_plink.fam", stringsAsFactors = FALSE)
q_a <- read_matrix("outputs/ohana_selection_accessions/best/q.k6.matrix")
sel_a <- data$accessions$selection
comp_a <- unique(sel_a$component[sel_a$significant]) + 1L
memb <- read.csv("outputs/mlg_membership.csv", stringsAsFactors = FALSE)
mll_a <- memb$mll[match(fam_a$V2, memb$id)]
ancestry_by_mll <- tapply(q_a[, comp_a[1]], mll_a, mean)
clone <- as.integer(names(which.max(ancestry_by_mll)))
in_clone <- mll_a == clone
# alternate allele frequency of the five loci inside and outside the lineage; presence in the genet panel
g_a <- read_matrix("outputs/ohana/g.dgm", integer_values = TRUE)
bim_a <- read.table("outputs/ohana/gl_plink.bim", stringsAsFactors = FALSE)
loci_a <- sel_a$locus[sel_a$significant]
cols_a <- match(loci_a, bim_a$V2)
alt_freq <- function(rows) {
  x <- g_a[rows, cols_a, drop = FALSE]; x[x == 3L] <- NA
  colMeans(x, na.rm = TRUE) / 2
}
clone_check <- data.frame(locus = loci_a, clone_lineage = clone, clone_ramets = sum(in_clone),
                          clone_mean_ancestry = round(ancestry_by_mll[as.character(clone)], 3),
                          alt_freq_clone = round(alt_freq(in_clone), 3),
                          alt_freq_others = round(alt_freq(!in_clone), 3),
                          in_genet_panel = loci_a %in% data$genets$markers$locus)
write.csv(clone_check, "outputs/ohana_selection_accessions/clone_component_check.csv", row.names = FALSE)

# ---- Gene context of the candidates on the annotated Ahi assembly (GCA_009445465.1) ----
## Results: Taro_020035 (intronic, S/T kinase), Taro_012402 (0.9 kb upstream, PP2C), Taro_020717 (intronic, no hit)
threads <- 8L
annot_dir <- "data/annotation"
out_dir <- "outputs/ohana_candidate_genes_ahi"
genome_gz <- "data/GCA_009445465.1_ASM944546v1_genomic.fna.gz"
gff_gz <- "data/GCA_009445465.1_ASM944546v1_genomic.gff.gz"
genome_fa <- file.path(annot_dir, "GCA_009445465.1_ASM944546v1_genomic.fna")
ahi_zip <- file.path(annot_dir, "GCA_009445465.1_ahi.zip")
protein_fa <- file.path(annot_dir, "GCA_009445465.1_protein.faa")
sprot_gz <- file.path(annot_dir, "uniprot_sprot.fasta.gz")
sprot_fa <- file.path(annot_dir, "uniprot_sprot.fasta")
sprot_db <- "outputs/ohana_candidate_genes/blastdb/sprot"

# Downloads: Ahi assembly and GFF (NCBI FTP), Ahi proteins (NCBI datasets), Swiss-Prot (UniProt)
ncbi_ftp <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/445/465/GCA_009445465.1_ASM944546v1/"
system2("curl", c("-sSL", "--max-time", "900", "-o", genome_gz, paste0(ncbi_ftp, basename(genome_gz))))
system2("curl", c("-sSL", "--max-time", "900", "-o", gff_gz, paste0(ncbi_ftp, basename(gff_gz))))
system2("curl", c("-sSL", "--max-time", "900", "-o", ahi_zip,
                  shQuote("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_009445465.1/download?include_annotation_type=PROT_FASTA")))
unzip(ahi_zip, files = "ncbi_dataset/data/GCA_009445465.1/protein.faa",
      exdir = file.path(annot_dir, "ahi_tmp"))
file.copy(file.path(annot_dir, "ahi_tmp/ncbi_dataset/data/GCA_009445465.1/protein.faa"),
          protein_fa, overwrite = TRUE)
unlink(file.path(annot_dir, "ahi_tmp"), recursive = TRUE)
system2("curl", c("-sSL", "--max-time", "900", "-o", sprot_gz,
                  "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"))
R.utils::gunzip(genome_gz, destname = genome_fa, remove = FALSE)
R.utils::gunzip(sprot_gz, destname = sprot_fa, remove = FALSE)
system2("samtools", c("faidx", genome_fa))

# Candidate loci: every locus significant in any panel
x <- readRDS("input_files/taro_diploid_mll.rds")
comparison <- read.csv("outputs/ohana_selection_candidate_comparison.csv", stringsAsFactors = FALSE)
sig <- comparison[comparison$significant %in% TRUE, ]
candidates <- unique(comparison[c("locus", "chr", "pos")])
candidates$panels_significant <- vapply(candidates$locus, function(l) {
  paste(sort(unique(sig$panel[sig$locus == l])), collapse = ";")
}, character(1))
candidates <- candidates[order(candidates$chr, candidates$pos), ]

# Map every DArTseq tag to the Ahi assembly (megablast, one best hit per tag)
xb <- gl.blast(x, ref_genome = genome_fa, task = "megablast", Percentage_identity = 70,
               Percentage_overlap = 0.8, bitscore = 50, number_of_threads = threads, verbose = 0)
lmb <- xb$other$loc.metrics
hit_cols <- c("sacc", "pident", "length", "evalue", "bitscore", "qstart", "qend",
              "sstart", "send", "qlen", "slen", "PercentageOverlap")
all_hits <- data.frame(locus = locNames(x), lmb[hit_cols], stringsAsFactors = FALSE)
write.csv(all_hits, file.path(out_dir, "tag_blast_all_loci.csv"), row.names = FALSE)
# All hits passing the filters, from the unfiltered BLAST table gl.blast() leaves in tempdir()
raw <- read.table(file.path(tempdir(), "output_blast.txt"), sep = "\t", quote = "\"",
                  comment.char = "", stringsAsFactors = FALSE)
names(raw) <- c("qseqid", "sacc", "stitle", "qseq", "sseq", "nident", "mismatch", "pident",
                "length", "evalue", "bitscore", "qstart", "qend", "sstart", "send",
                "gapopen", "gaps", "qlen", "slen")
raw$overlap <- raw$length / pmin(raw$qlen, raw$slen)
raw <- raw[raw$overlap >= 0.8 & raw$bitscore >= 50, ]
raw$locus <- locNames(x)[raw$qseqid]
write.csv(raw[c("locus", "sacc", "pident", "length", "evalue", "bitscore", "qstart", "qend",
                "sstart", "send", "qlen", "slen")],
          file.path(out_dir, "tag_blast_filtered_hits.csv"), row.names = FALSE)
n_hits <- table(raw$locus)
lm <- x$other$loc.metrics

# SNP position on the assembly from the hit coordinates (DArT SnpPosition is 0-based)
snp_in_tag <- as.integer(lm$SnpPosition) + 1L
strand <- ifelse(is.na(all_hits$sstart), NA, ifelse(all_hits$send >= all_hits$sstart, "+", "-"))
in_alignment <- !is.na(all_hits$qstart) & snp_in_tag >= all_hits$qstart & snp_in_tag <= all_hits$qend
offset <- snp_in_tag - all_hits$qstart
snp_pos <- ifelse(in_alignment,
                  ifelse(strand == "+", all_hits$sstart + offset, all_hits$sstart - offset),
                  NA_integer_)
xb$chromosome <- factor(ifelse(is.na(snp_pos), NA, as.character(all_hits$sacc)))
xb$position <- as.integer(snp_pos)

ci <- match(candidates$locus, locNames(xb))
candidate_hits <- data.frame(
  candidates, scaffold = as.character(all_hits$sacc[ci]), strand = strand[ci],
  snp_position = snp_pos[ci], pident = all_hits$pident[ci], aln_length = all_hits$length[ci],
  tag_length = all_hits$qlen[ci], evalue = all_hits$evalue[ci], bitscore = all_hits$bitscore[ci],
  n_filtered_hits = as.integer(n_hits[candidates$locus]), snp_within_alignment = in_alignment[ci],
  jaas_strand = as.character(lm[[grep("^Strand_", names(lm))[1]]][ci]), stringsAsFactors = FALSE)
candidate_hits$n_filtered_hits[is.na(candidate_hits$n_filtered_hits)] <- 0L
write.csv(candidate_hits, file.path(out_dir, "candidate_hits.csv"), row.names = FALSE)

# Overlapping or nearest gene from the NCBI GFF
genes <- as.data.frame(gl.find.genes.for.loci(xb, gff.file = gff_gz, loci = candidates$locus,
                                              verbose = 0))

# Exon/CDS status and codon change for SNPs inside a gene
gff <- as.data.table(ape::read.gff(gzfile(gff_gz)))
attr_of <- function(a, key) {
  r <- regexpr(paste0("(^|;)", key, "=[^;]+"), a)
  out <- rep(NA_character_, length(a))
  out[r > 0] <- sub(paste0("^;?", key, "="), "", regmatches(a, r))
  out
}
gff[, ID := attr_of(attributes, "ID")]
gff[, Parent := attr_of(attributes, "Parent")]
gff[, Name := attr_of(attributes, "Name")]
faidx <- function(region) {
  toupper(paste(system2("samtools", c("faidx", genome_fa, region), stdout = TRUE)[-1], collapse = ""))
}
codon_table <- c(
  TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
  ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
  TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
  ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
  TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
  AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
  TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
  AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G")
comp <- function(s) chartr("ACGT", "TGCA", s)
# DArT SNP field is "pos:ref>alt" in tag orientation; complemented for minus-strand hits
tag_alleles <- do.call(rbind, strsplit(sub("^[0-9]+:", "", as.character(lm$SNP)), ">"))

gene_context <- function(j) {
  g <- genes[j, ]
  row <- list(exon = NA, cds = NA, codon_ref = NA, codon_alt = NA, aa_ref = NA, aa_alt = NA, effect = NA)
  if (is.na(g$nearest_side) || g$nearest_side != "inside") return(row)
  i <- match(g$locus, locNames(xb))
  pos <- xb$position[i]; scaffold <- as.character(xb$chromosome[i])
  mrna <- gff[type == "mRNA" & Parent == g$gene_id]
  if (!nrow(mrna)) return(row)
  parts <- gff[Parent %in% mrna$ID & seqid == scaffold]
  exons <- parts[type == "exon"]; cds <- parts[type == "CDS"]
  row$exon <- any(exons$start <= pos & exons$end >= pos)
  row$cds <- any(cds$start <= pos & cds$end >= pos)
  if (!row$cds) return(row)
  tid <- cds$Parent[cds$start <= pos & cds$end >= pos][1]
  cds <- cds[Parent == tid]
  minus <- cds$strand[1] == "-"
  cds <- cds[order(start, decreasing = minus)]
  seg <- which(cds$start <= pos & cds$end >= pos)
  before <- if (seg > 1) sum(cds$end[seq_len(seg - 1)] - cds$start[seq_len(seg - 1)] + 1L) else 0L
  within <- if (minus) cds$end[seg] - pos else pos - cds$start[seg]
  phase <- as.integer(cds$phase[1]); if (is.na(phase)) phase <- 0L
  cpos <- before + within - phase
  if (cpos < 0) return(row)
  codon_start <- cpos - cpos %% 3
  coding_to_genome <- function(c0) {
    c0 <- c0 + phase; acc <- 0L
    for (s in seq_len(nrow(cds))) {
      len <- cds$end[s] - cds$start[s] + 1L
      if (c0 < acc + len) return(if (minus) cds$end[s] - (c0 - acc) else cds$start[s] + (c0 - acc))
      acc <- acc + len
    }
    NA_integer_
  }
  gpos <- vapply(codon_start + 0:2, coding_to_genome, numeric(1))
  bases <- vapply(gpos, function(p) faidx(sprintf("%s:%d-%d", scaffold, p, p)), character(1))
  if (minus) bases <- comp(bases)
  ref_codon <- paste(bases, collapse = "")
  alt_base <- tag_alleles[i, 2]; ref_base <- tag_alleles[i, 1]
  if (strand[i] == "-") { alt_base <- comp(alt_base); ref_base <- comp(ref_base) }
  if (minus) { alt_base <- comp(alt_base); ref_base <- comp(ref_base) }
  idx <- which(gpos == pos)
  alt_codon <- ref_codon; substr(alt_codon, idx, idx) <- alt_base
  row$codon_ref <- ref_codon; row$codon_alt <- alt_codon
  row$aa_ref <- unname(codon_table[ref_codon]); row$aa_alt <- unname(codon_table[alt_codon])
  row$effect <- if (substr(ref_codon, idx, idx) != ref_base) "reference allele mismatch" else
    if (row$aa_ref == row$aa_alt) "synonymous" else
      if (row$aa_alt == "*") "stop gained" else "nonsynonymous"
  row
}
ctx <- do.call(rbind, lapply(seq_len(nrow(genes)), function(j) as.data.frame(gene_context(j))))
genes <- cbind(genes, ctx)

# Putative function of each reported gene model: BLASTP of its protein against Swiss-Prot
system2("makeblastdb", c("-in", sprot_fa, "-dbtype", "prot", "-out", sprot_db), stdout = FALSE)
headers <- grep("^>", readLines(protein_fa), value = TRUE)
prot_acc <- sub("^>([^ ]+) .*", "\\1", headers)
prot_gene <- sub(".*hypothetical protein ([^ ,]+).*", "\\1", headers)
genes$protein_id <- prot_acc[match(genes$gene_name, prot_gene)]
ids <- unique(genes$protein_id[!is.na(genes$protein_id)])
system2("samtools", c("faidx", protein_fa))
q <- file.path(out_dir, "gene_proteins.query.fa")
system2("samtools", c("faidx", protein_fa, ids), stdout = q)
o <- file.path(out_dir, "gene_proteins.query.sprot.tsv")
fields <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle")
system2("blastp", c("-query", q, "-db", sprot_db, "-evalue", "1e-10", "-max_target_seqs", "5",
                    "-num_threads", threads,
                    "-outfmt", shQuote(paste("6", paste(fields, collapse = " "))), "-out", o))
hits <- read.delim(o, header = FALSE, col.names = fields, quote = "", stringsAsFactors = FALSE)
missing <- setdiff(ids, hits$qseqid)
hits <- rbind(hits, data.frame(qseqid = missing, sseqid = "no_hit", pident = NA, length = NA,
                               evalue = NA, bitscore = 0, stitle = "no Swiss-Prot hit",
                               stringsAsFactors = FALSE))
write.table(hits, file.path(out_dir, "gene_proteins.sprot.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
best <- hits[order(hits$qseqid, -hits$bitscore), ]
best <- best[!duplicated(best$qseqid), ]
fn <- ifelse(best$sseqid == "no_hit", "no Swiss-Prot hit",
             sprintf("%s [%.0f%% id; e=%s]", sub(" OS=.*", "", best$stitle), best$pident, best$evalue))
genes$sprot_function <- unname(setNames(fn, best$qseqid)[genes$protein_id])

genes <- merge(candidate_hits[c("locus", "panels_significant", "strand", "pident", "n_filtered_hits")],
               genes, by = "locus", all.x = TRUE)
genes <- genes[order(match(genes$locus, candidates$locus)), ]
write.csv(genes, file.path(out_dir, "candidate_genes.csv"), row.names = FALSE)
