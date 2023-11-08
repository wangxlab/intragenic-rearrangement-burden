# ----- 1.1 Prepare Data (DO NOT RUN) -----
# this intragenic rearrangement (IGR) is calculated from WGS data
icgc_igr <- read.table("data_for_manuscript/ICGC_PANCAN_20180104_StSM-allfiles.GeneRearrangements3kbPromoter.UCSCAndGenCodeMergeGRCh37.PassOnly.V9.svcountv4.bysvType-HanUpdate.tsv", 
                       sep = "\t", head = TRUE, check.names = FALSE)
icgc_mut <- read.table("data_for_manuscript/ICGC_Variations_All.tsv", 
                       sep = "\t", head = TRUE, check.names = FALSE)
# combine IGR and TMB
icgc_igr_mut <- inner_join(icgc_igr, icgc_mut, by = c("SampleID" = "icgc_specimen_id"))
icgc_igr_mut$sqrt_igr <- sqrt(icgc_igr_mut$Intragenic.DEL.n. + icgc_igr_mut$Intragenic.DUP.n.)
icgc_igr_mut$sqrt_tmb <- sqrt(icgc_igr_mut$missense_variant)
icgc_igr_mut$sqrt_indel <- sqrt(icgc_igr_mut$frameshift_variant)
# ----- Somatic CNV -----
setwd("~/Documents/Research/WangLab/Manuscript/rebuttal_files/")
a = list.files("consensus.20170119.somatic.cna.icgc.public/")
a2 = list.files("consensus.20170119.somatic.cna.tcga.public/")
purity = read.table("consensus.20170217.purity.ploidy.txt", sep = "\t", head = T)
aa = c(a, a2)
b = unlist(sapply(aa, function(x) unlist(strsplit(x, "\\."))[1]))

purity = purity[match(as.character(b), purity$samplename), ]

project = read.table("~/Documents/Research/WangLab/Projects/IGR/data_for_manuscript/pcawg_sample_sheet.v1.4.2016-09-14.tsv", 
                     head = T, sep = "\t")
a = paste0("consensus.20170119.somatic.cna.icgc.public/", a)
a2 = paste0("consensus.20170119.somatic.cna.tcga.public/", a2)
a = c(a, a2)
a = lapply(a, read.table, head = T)
cnv = sapply(a, function(x) {
    nrow(x[x$total_cn != 2 | x$major_cn != 1 | x$minor_cn != 1, ])
})

adj_cnv = c()
for (i in 1:2778) {
    tmp = a[[i]]
    pp = purity$purity[i]
    nominator = tmp$total_cn*pp + 2*(1-pp)
    #demoninator = pp*3 + 2*(1-pp)
    demoninator = pp*purity$ploidy[i] + 2*(1-pp)
    tmp = nominator/demoninator
    tmp = log2(abs(na.omit(tmp)))
    tmp = tmp[tmp >= 0.35]
    tmp = ifelse(tmp > 1, 2, 1)
    adj_cnv[i] = sum(tmp)
}
adj_cnv2 = adj_cnv

adj_cnv = c()
for (i in 1:2778) {
    tmp = a[[i]]
    pp = purity$purity[i]
    nominator = tmp$total_cn*pp + 2*(1-pp)
    demoninator = pp*3 + 2*(1-pp)
    # demoninator = pp*purity$ploidy[i] + 2*(1-pp)
    tmp = nominator/demoninator
    tmp = abs(na.omit(tmp))
    tmp = tmp[tmp >= 0.35]
    tmp = ifelse(tmp > 1, 2, 1)
    adj_cnv[i] = sum(tmp)
}

cnv = data.frame(id = b, cnv = cnv, scna = adj_cnv, scna2 = adj_cnv2)
cnv_db= inner_join(cnv, project, by = c("id" = "aliquot_id"))
cnv_db= inner_join(cnv_db, purity, by = c("id" = "samplename"))

icgc = left_join(icgc, cnv_db, by = c("SampleID" = "icgc_specimen_id"))
icgc$sqrt_scna = sqrt(icgc$scna)
icgc_igr_mut$catype <- icgc_igr_mut$project_code.x
write.table(icgc_igr_mut, "data_for_manuscript/icgc_igr_mut.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

icgc_expr <- ReadExpr("data_for_manuscript/tophat_star_fpkm.v2_aliquot_gl.tsv.gz")
icgc_smp_info <- read.table("data_for_manuscript/pcawg_sample_sheet.v1.4.2016-09-14.tsv", 
                            sep = "\t", head = TRUE, check.names = FALSE)
icgc_msi <- readxl::read_xlsx("data_for_manuscript/MS_analysis.PCAWG_release_v1.RIKEN.xlsx")
icgc_msi_smp <- inner_join(icgc_msi, icgc_smp_info, by = c("ID" = "aliquot_id"))
icgc_msi_smp <- icgc_msi_smp[, c("icgc_specimen_id", "mutation_rate")]

# the raw column names of expr matrix are "aliquot_id". We need to map them back to 
# "icgc_specimen_id"
colnames(icgc_expr) <- colnames(icgc_expr) %>%
    match(icgc_smp_info$aliquot_id) %>%
    icgc_smp_info$icgc_specimen_id[.]
icgc_expr <- icgc_expr[, !is.na(colnames(icgc_expr))]
icgc_expr_tpm <- fpkm2tpm(icgc_expr) %>%
    as.data.frame()
# "gene symbol" = gsub("\\..*", "", rownames(icgc_expr_tpm)
# look up "gene symbol" in [https://www.biotools.fr/human/ensembl_symbol_converter]
# the returned mapping list is saved in ENSG_to_Gene.txt
ensg_gene_list <- read.table("data_for_manuscript/ENSG_to_Gene.txt", fill = TRUE)
valid_idx <- which(ensg_gene_list$V2 != "")
icgc_expr_tpm <- icgc_expr_tpm[valid_idx, ]
rownames(icgc_expr_tpm) <- make.names(ensg_gene_list$V2[valid_idx], unique = TRUE)
write.table(icgc_expr_tpm, "data_for_manuscript/icgc_expr_tpm.tsv", 
            quote = FALSE, sep = "\t")

# ----- 1.2 Calculate Signatures (DO NOT RUN) -----
icgc_expr_tpm <- ReadExpr("data_for_manuscript/icgc_expr_tpm.tsv")
# cibersort
icgc_cibersort <- CalCiber(icgc_expr_tpm, "cibersort")
# inflame signature
icgc_inflame <- CalImmuSig(icgc_expr_tpm)
# cell cycle signature
icgc_cellcycle <- CalCellCycleSig(icgc_expr_tpm)
icgc_sig <- inner_join(icgc_cibersort, icgc_inflame, by = c("samples" = "SampleName"))
icgc_sig <- inner_join(icgc_sig, icgc_cellcycle, by = c("samples" = "icgc_specimen_id"))
write.table(icgc_sig, "data_for_manuscript/icgc_signatures.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# merge IGR, TMB, Indel, Signatures and Microsatellite mutation rate
icgc_all <- inner_join(icgc_igr_mut, icgc_sig, by = c("SampleID" = "samples"))
icgc_all <- left_join(icgc_all, icgc_msi_smp, by = c("SampleID" = "icgc_specimen_id"))
write.table(icgc_all, "data_for_manuscript/icgc_all.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
cancer_types <- unique(icgc_all$catype)
cancer_igr <- array()
cancer_tmb <- array()
cancer_indel <- array()
cancer_scnv <- array()
for (i in 1:length(cancer_types)) {
    tmp = icgc_all[icgc_all$catype == cancer_types[i], ]
    cancer_igr[i] <- median(tmp[, "sqrt_igr"])
    cancer_tmb[i] <- median(tmp[, "sqrt_tmb"])
    cancer_indel[i] <- median(tmp[, "sqrt_indel"])
    cancer_scnv[i] <- median(tmp[, "sqrt_scna2"]) 
}
cancer_info = data.frame(sqrt_igr   = cancer_igr, 
                         sqrt_tmb   = cancer_tmb, 
                         sqrt_indel = cancer_indel, 
                         sqrt_scna2 = cancer_scnv,
                         catype     = cancer_types)

for (i in 1:length(cancer_types)) {
    tmp = icgc_all[icgc_all$catype == cancer_types[i], ]
    cancer_igr[i] <- mean(tmp[, "sqrt_igr"])
    cancer_tmb[i] <- mean(tmp[, "sqrt_tmb"])
    cancer_indel[i] <- mean(tmp[, "sqrt_indel"])
    cancer_scnv[i] <- mean(tmp[, "sqrt_scna2"], na.rm = T) 
}
cancer_info_mean = data.frame(sqrt_igr   = cancer_igr, 
                         sqrt_tmb   = cancer_tmb, 
                         sqrt_indel = cancer_indel, 
                         sqrt_scna2 = cancer_scnv,
                         catype     = cancer_types)
write.table(cancer_info_mean, "data_for_manuscript/icgc_catype_mean.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)

# ----- 1.3 Label TMB and IGR  -----
setwd("~/Documents/Research/WangLab/Projects/IGR/")
library(dplyr)
icgc <- read.table("data_for_manuscript/icgc_all.tsv", 
                   sep = "\t", head = TRUE)
icgc$inflame <- icgc$T_effector_signature_TotalScore
icgc$msi <- log10(icgc$mutation_rate)
n <- array()
for (i in 1:nrow(icgc)) {
    val <- icgc$catype[i]
    n[i] <- which(icgc$catype == val) %>% length()
}
icgc$catype_n <- paste0(icgc$catype, " (n=", n, ")")
cancer_info <- read.table("data_for_manuscript/icgc_catype.tsv", 
                          sep = "\t", head = TRUE)
cancer_info_norm <- apply(cancer_info[, 1:3], 2, function (x) { 
    (x - mean(x)) / sd(x)})
skcm_idx <- which(cancer_info$catype == "SKCM")
cancer_info_norm <- data.frame(cancer_info_norm, 
                               catype = cancer_info$catype)
# SKCM is an outlier, keeping SKCM will enlarger the sd, thus remove it
cancer_info_norm$sqrt_tmb <- 
    (cancer_info$sqrt_tmb - mean(cancer_info$sqrt_tmb[-skcm_idx])) / sd(cancer_info$sqrt_tmb[-skcm_idx])

# ----- Figure 1a -----
# Landscape violion plot of IGR Burden across all Cancer Types
# see PRISM result ["~/Projects/IGR/data_for_manuscript/IGR.pzfx"]
# Fig1. Pancancer IGR violion plot

# ----- Figure 1b -----
# Scatter plot Showing IGR burden and mutation burden are not associated
fig1b <- ggplot(icgc, aes(x = sqrt_igr, y = sqrt_tmb, col = inflame)) + 
    geom_point(size = 4.5, alpha = 0.4) + 
    labs(x = "IGR Burden", y = "Tumor Mutation Burden", col = "T-Inflamed Sig.") + 
    scale_color_gradient2(low = "#3300FF", high = "#FF0000", mid = "white", 
                          midpoint = median(icgc$inflame), 
                          oob = scales::squish, 
                          limits = quantile(icgc$inflame, c(0.2, 0.8)),
                          breaks = quantile(icgc$inflame, c(0.2, 0.5, 0.8)), 
                          labels = c("Low (<20%)", "Median", "High (>80%)")) + 
    theme_classic() + 
    theme(legend.position = c(0.95, 0.6), 
          legend.title = element_text(size = 18, face = "bold"), 
          legend.text = element_text(size = 15),
          legend.justification = c(1, 0), 
          axis.text = element_text(size = 18), 
          axis.title = element_text(size = 18))
#    annotate(geom = "curve", x = 30, y = 90, xend = 200, yend = 20, 
#             curvature = 0.3, lwd = 1.5) + 
#    annotate("text", x = 70, y = 90, label = "TMB-driven", size = 8) + 
#    annotate("text", x = 200, y = 30, label = "IGR-driven", size = 8)
#ggsave("figures/fig1b.pdf", fig1b, units = "in", width = 8, height = 6)
ggsave("~/Documents/fig1b.pdf", fig1b, units = "in", width = 7, height = 5)
# ----- Figure S2a -----
# Scatter plot showing frameshift and mutation burden are positively associated
figS2a <- ggplot(icgc, aes(x = sqrt_indel, y = sqrt_tmb, col = inflame)) + 
    geom_point(size = 4.5, alpha = 0.4) + 
    labs(x = "Frameshift Burden", y = "Tumor Mutation Burden", col = "Inflame Sig.") + 
    scale_color_gradient2(low = "#3300FF", high = "#FF0000", mid = "white", 
                          midpoint = median(icgc$inflame), 
                          oob = scales::squish, 
                          limits = quantile(icgc$inflame, c(0.2, 0.8)),
                          breaks = quantile(icgc$inflame, c(0.2, 0.5, 0.8)), 
                          labels = c("Low (<20%)", "Median", "High (>80%)")) + 
    theme_classic() + 
    theme(legend.position = c(0.6, 0.6), 
          legend.title = element_text(size = 15, face = "bold"), 
          legend.text = element_text(size = 12),
          legend.justification = c(1, 0), 
          axis.text = element_text(size = 18), 
          axis.title = element_text(size = 18))
ggsave("figures/figS2a.pdf", figS2a, units = "in", width = 8, height = 6)

# ----- Figure S2b -----
# Scatter plot showing IGR burden and frameshift are not associated
figS2b <- ggplot(icgc, aes(x = sqrt_indel, y = sqrt_igr, col = inflame)) + 
    geom_point(size = 4.5, alpha = 0.4) + 
    labs(x = "Frameshift Burden", y = "IGR Burden", col = "Inflame Sig.") + 
    scale_color_gradient2(low = "#3300FF", high = "#FF0000", mid = "white", 
                          midpoint = median(icgc$inflame), 
                          oob = scales::squish, 
                          limits = quantile(icgc$inflame, c(0.2, 0.8)),
                          breaks = quantile(icgc$inflame, c(0.2, 0.5, 0.8)), 
                          labels = c("Low (<20%)", "Median", "High (>80%)")) + 
    theme_classic() + 
    theme(legend.position = c(0.95, 0.6), 
          legend.title = element_text(size = 18, face = "bold"), 
          legend.text = element_text(size = 15),
          legend.justification = c(1, 0), 
          axis.text = element_text(size = 18), 
          axis.title = element_text(size = 18))
ggsave("figures/figS2b.pdf", figS2b, units = "in", width = 8, height = 6)

# ----- Figure S3 -----
# Figure 1b by cancer type
figS3 <- ggplot(icgc, aes(x = sqrt_igr, y = sqrt_scna, col = inflame)) + 
    geom_point(size = 2.5, alpha = 0.4) + 
    labs(x = "IGR Burden", y = "Tumor Mutation Burden", col = "Inflame Sig.") + 
    scale_color_gradient2(low = "#3300FF", high = "#FF0000", mid = "white", 
                          midpoint = median(icgc$inflame), 
                          oob = scales::squish, 
                          limits = quantile(icgc$inflame, c(0.2, 0.8)),
                          breaks = quantile(icgc$inflame, c(0.2, 0.5, 0.8)), 
                          labels = c("Low (<20%)", "Median", "High (>80%)")) + 
    theme_classic() + 
    theme(legend.title = element_text(size = 15, face = "bold"), 
          legend.text = element_text(size = 12),
          strip.text.x = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 18)) + 
    facet_wrap(~catype_n, nc = 7)
ggsave("figures/figS3.pdf", figS3, units = "in", width = 12, height = 5)

# ----- Figure 1c -----
# see PRISM result ["~/Projects/IGR/data_for_manuscript/IGR.pzfx"]
# Fig 1c. Cancer Types (IGR Burden vs Mutation Burden)
write.table(cancer_info_norm, "~/Documents/cancer_info_norm.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

# ----- Figure 1d -----
igr_catype <- cancer_info_norm[cancer_info_norm$sqrt_igr >= 1 & 
                                   cancer_info_norm$sqrt_tmb < 1, "catype"]
tmb_catype <- cancer_info_norm[cancer_info_norm$sqrt_igr < 1 & 
                                   cancer_info_norm$sqrt_tmb >= 1, "catype"]
interested <- list()
for (i in igr_catype) {
    tmp <- icgc[icgc$catype == i, ]
    tmp$igr_level_in_cancer <- ifelse(tmp$sqrt_igr >= quantile(tmp$sqrt_igr, .5), 
                                      "IGR_High", "IGR_Low")
    tmp$tmb_level_in_cancer <- ifelse(tmp$sqrt_tmb >= quantile(tmp$sqrt_tmb, .5), 
                                      "TMB_High", "TMB_Low")
    interested[[i]] = tmp
}
interested <- bind_rows(interested)
interested$total_level_same_cancer <- 
    paste0(interested$igr_level_in_cancer, "/", interested$tmb_level_in_cancer)

ggplot(interested, aes(x    = total_level_same_cancer, 
                       y    = inflame, 
                       fill = total_level_same_cancer)) + 
    geom_boxplot(lwd = .75) + labs(x = "", y = "Inflammatory Signature", fill = "") + 
    stat_compare_means(comparisons = list(c("IGR_High/TMB_High", "IGR_High/TMB_Low"),
                                          c("IGR_High/TMB_Low", "IGR_Low/TMB_Low"),
                                          c("IGR_High/TMB_High", "IGR_Low/TMB_High"),
                                          c("IGR_High/TMB_High", "IGR_Low/TMB_Low")), 
                       method.args = list(alternative = "greater")) + 
    scale_fill_manual(values = c("#d11141", "#ffc425", "#00aedb", "#00b159")) + 
    theme_pubclean() + 
    theme(legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 14),
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 18),
          axis.ticks.x = element_blank()) + 
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
