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
icgc_igr_mut$sqrt_cnv <- sqrt(icgc_igr_mut$gain + icgc_igr_mut$loss)
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
# ----- Figure S1a -----
# Scatter plot showing frameshift and mutation burden are positively associated
figS1a <- ggplot(icgc, aes(x = sqrt_indel, y = sqrt_tmb, col = inflame)) + 
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
ggsave("figures/figS1a.pdf", figS1a, units = "in", width = 8, height = 6)

# ----- Figure S1b -----
# Scatter plot showing IGR burden and frameshift are not associated
figS1b <- ggplot(icgc, aes(x = sqrt_indel, y = sqrt_igr, col = inflame)) + 
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
ggsave("figures/figS1b.pdf", figS1b, units = "in", width = 8, height = 6)

# ----- Figure S2 -----
# Figure 1b by cancer type
figS2 <- ggplot(icgc, aes(x = sqrt_igr, y = sqrt_tmb, col = inflame)) + 
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
ggsave("figures/figS2.pdf", figS2, units = "in", width = 12, height = 5)

# ----- Figure 1c -----
# see PRISM result ["~/Projects/IGR/data_for_manuscript/IGR.pzfx"]
# Fig 1c. Cancer Types (IGR Burden vs Mutation Burden)
write.table(cancer_info_norm, "~/Documents/cancer_info_norm.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

# ----- Fig S3 -----
# see PRISM result ["~/Projects/IGR/data_for_manuscript/IGR.pzfx"]
# Fig S3. Cancer Types (Frameshift Burden vs Mutation Burden)

# ----- Figure 1d -----
igr_catype <- cancer_info_norm[cancer_info_norm$sqrt_igr >= 1 & 
                                   cancer_info_norm$sqrt_tmb < 1, "catype"]
tmb_catype <- cancer_info_norm[cancer_info_norm$sqrt_igr < 1 & 
                                   cancer_info_norm$sqrt_tmb >= 1, "catype"]
catype3 = c("KIRC", "KIRP", "KICH")
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

# measure R^2 (variance explained in predicting inflame sig)
# igr significant, p = 0.024
m1 <- lm(inflame ~ sqrt_tmb*catype, data = interested)
m2 <- lm(inflame ~ sqrt_igr*catype, data = interested)
m3 <- lm(inflame ~ sqrt_igr*catype + sqrt_tmb*catype, data = interested)
rs <- c(summary(m1)$r.squared, summary(m2)$r.squared, summary(m3)$r.squared)
anova(m1, m3)

m1 <- lm(mutation_rate ~ sqrt_tmb*catype, data = interested)
m2 <- lm(mutation_rate ~ sqrt_igr*catype, data = interested)
m3 <- lm(mutation_rate ~ sqrt_igr*catype + sqrt_tmb*catype, data = interested)
rs <- c(summary(m1)$r.squared, summary(m2)$r.squared, summary(m3)$r.squared)
anova(m1, m3)
#anova(m2, m3)
#rs

#write.table(interested$cell_cycle_sig[interested$tmb_level_in_cancer == "TMB_Low"],
#            "~/Documents/cc_tmb.txt", 
#            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

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
for (i in c("IGR_High/TMB_High", "IGR_High/TMB_Low",
            "IGR_Low/TMB_High", "IGR_Low/TMB_Low")) {
    write.table(interested[interested$total_level_same_cancer == i, "inflame"],
                paste0("~/Documents/", gsub("/", "_", i), ".txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# IGR vs TMB vs Indel vs CNV
na <- which(! is.na(icgc$sqrt_cnv))
icgc_na <- droplevels(icgc[na, ])
figS1_igr <- ggplot(icgc_na, aes(x = sqrt_igr, y = sqrt_cnv)) + 
    geom_point(size = 4.5, alpha = 0.4) + 
    labs(x = "IGR Burden", y = "CNV Burden", 
         title = paste0("Pearson=", 
                        round(cor(icgc_na$sqrt_igr, icgc_na$sqrt_cnv), digit = 3), 
                        ", pvalue=", 
                        round(cor.test(icgc_na$sqrt_igr, icgc_na$sqrt_cnv)$p.value, digit = 3))) + 
    theme_classic() +
    geom_smooth(method = "lm", se = FALSE) + 
    theme(axis.text = element_text(size = 18), 
          axis.title = element_text(size = 18), 
          title = element_text(size = 12))
# ggsave("figures/figS1_igr.pdf", figS1_igr)
figS1_tmb <- ggplot(icgc_na, aes(x = sqrt_tmb, y = sqrt_cnv)) + 
    geom_point(size = 4.5, alpha = 0.4) + 
    labs(x = "Mutation Burden", y = "CNV Burden", 
         title = paste0("Pearson=", 
                        round(cor(icgc_na$sqrt_tmb, icgc_na$sqrt_cnv), digit = 3), 
                        ", pvalue=", 
                        round(cor.test(icgc_na$sqrt_tmb, icgc_na$sqrt_cnv)$p.value, digit = 3))) + 
    theme_classic() +
    geom_smooth(method = "lm", se = FALSE) + 
    theme(axis.text = element_text(size = 18), 
          axis.title = element_text(size = 18), 
          title = element_text(size = 12))
figS1_indel <- ggplot(icgc_na, aes(x = sqrt_indel, y = sqrt_cnv)) + 
    geom_point(size = 4.5, alpha = 0.4) + 
    labs(x = "Indel Burden", y = "CNV Burden", 
         title = paste0("Pearson=", 
                        round(cor(icgc_na$sqrt_indel, icgc_na$sqrt_cnv), digit = 3), 
                        ", pvalue=", 
                        round(cor.test(icgc_na$sqrt_indel, icgc_na$sqrt_cnv)$p.value, digit = 3))) + 
    theme_classic() +
    geom_smooth(method = "lm", se = FALSE) + 
    theme(axis.text = element_text(size = 18), 
          axis.title = element_text(size = 18), 
          title = element_text(size = 12))
figS1_f = ggarrange(figS1_igr, figS1_tmb, figS1_indel, ncol = 3)
ggsave("figures/figS1_f.pdf", figS1_f, units = "in", width = 12, height = 4)

# lung squamous
lusc = icgc[icgc$catype == "LUSC", ] %>% droplevels()
