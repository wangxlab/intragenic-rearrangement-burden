# ----- 1.1 Prepare Data (DO NOT RUN) -----
# get log2 FPKM data from below
# data_for_manuscript/nature17676-s3/Supplementary Table 7.Transcriptomic.342.txt
wgs560_log2_fpkm <- ReadExpr("data_for_manuscript/wgs560-expr-log-fpkm.tsv")
wgs560_expr_fpkm <- 2 ^ wgs560_log2_fpkm
wgs560_expr_fpkm[wgs560_expr_fpkm == 1] <- 0
wgs560_expr_tpm <- fpkm2tpm(wgs560_expr_fpkm)
write.table(wgs560_expr_fpkm, "data_for_manuscript/wgs560_expr_fpkm.tsv", 
            quote = FALSE, sep = "\t")
write.table(wgs560_expr_tpm, "data_for_manuscript/wgs560_expr_tpm.tsv", 
            quote = FALSE, sep = "\t")

# identify triple negative breast cancer
sample_info <- read.csv("data_for_manuscript/nature17676-s3/Supplementary Table 1 CLINICAL.PATHOLOGY.DATA.FREEZE.ANALYSIS.v4.032015.csv", 
                        check.names = FALSE)
names(sample_info) <- sample_info[1, ]
sample_info <- sample_info[-1, ]
sample_info$TNBC <- ifelse(sample_info$final.ER == "negative" & 
                               sample_info$final.HER2 == "negative" & 
                               sample_info$final.PR == "negative", 
                           "Yes", "No")
# collect mutation information
mut <- readxl::read_xlsx("data_for_manuscript/nature17676-s3/Supplementary Table 3.Summary.Somatic.Catalogue.v1.xlsx")
mut <- as.data.frame(mut)
mut <- mut[-561, ] # grand total
mut$id <- c(substr(mut$sample[1:296], 1, 7), substr(mut$sample[297:560], 1, 6))

# collect igr
igr <- read.table("data_for_manuscript/wgs560_igr.tsv", sep = "\t", head = TRUE)
igr <- igr[grep("^PD", igr$SampleID), ]

# hrd score
hrd <- readxl::read_xlsx("data_for_manuscript/nature17676-s3/Supplementary Table 17.HRD.intermediary.09042015.v1.xlsx")
hrd <- as.data.frame(hrd)
hrd <- hrd[, c("Sample", "HRD")]

# subtype
subtype <- readxl::read_xlsx("data_for_manuscript/nature17676-s3/Supplementary Table 18.Expression.Subtyping.xlsx")
subtype <- as.data.frame(subtype)
subtype <- subtype[, c("sample_name", "PAM50.subtype")]

# merge datasets
wgs560_db <- left_join(sample_info, mut, by = c("sample_name" = "id"))
wgs560_db <- left_join(wgs560_db, igr, by = c("sample" = "SampleID"))
wgs560_db$sqrt_igr <- sqrt(wgs560_db$Intragenic.DEL.n. + 
                               wgs560_db$Intragenic.DUP.n. + 
                               wgs560_db$Intragenic.INV3.n. + 
                               wgs560_db$Intragenic.INV5.n.)
# before remove all other samples
write.table(wgs560_db[wgs560_db$TNBC == "Yes", "sqrt_igr"], "~/Documents/tnbc.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(wgs560_db[wgs560_db$TNBC == "No", "sqrt_igr"], "~/Documents/nontnbc.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
wgs560_db <- left_join(wgs560_db, hrd, by = c("sample_name" = "Sample"))
wgs560_db <- left_join(wgs560_db, subtype, by = c("sample_name" = "sample_name"))
wgs560_tnbc_db <- wgs560_db[wgs560_db$TNBC == "Yes", ]

# ----- 1.2 Calculate Signatures (DO NOT RUN) -----
sample_name_in_expr <- c(substr(colnames(wgs560_expr_tpm)[1:149], 1, 7), 
                         substr(colnames(wgs560_expr_tpm)[150:342], 1, 6))
sample_name_in_expr <- gsub("R", "D", sample_name_in_expr)
idx <- which(sample_name_in_expr %in% wgs560_tnbc_db$sample_name)
wgs560_tnbc_expr_tpm <- wgs560_expr_tpm[, idx]
write.table(wgs560_tnbc_expr_tpm, "data_for_manuscript/wgs560_tnbc_expr_tpm.tsv", 
            quote = FALSE, sep = "\t")
colnames(wgs560_tnbc_expr_tpm) <- sample_name_in_expr[idx]

wgs560_tnbc_log2_fpkm <- wgs560_log2_fpkm[, idx]
colnames(wgs560_tnbc_log2_fpkm) <- sample_name_in_expr[idx]
write.table(wgs560_tnbc_log2_fpkm, "data_for_manuscript/wgs560_tnbc_expr_log2_fpkm.tsv", 
            quote = FALSE, sep = "\t")

# cibersort
wgs560_tnbc_cibersort <- CalCiber(wgs560_tnbc_expr_tpm, "cibersort")
# inflame signature
wgs560_tnbc_inflame <- CalImmuSig(wgs560_tnbc_expr_tpm)
# cell cycle signature
wgs560_tnbc_cellcycle <- CalCellCycleSig(wgs560_tnbc_expr_tpm)

wgs560_tnbc_sig <- inner_join(wgs560_tnbc_cibersort, wgs560_tnbc_inflame, by = c("samples" = "SampleName"))
wgs560_tnbc_sig <- inner_join(wgs560_tnbc_sig, wgs560_tnbc_cellcycle, by = c("samples" = "icgc_specimen_id"))

wgs560_tnbc_all <- left_join(wgs560_tnbc_db, wgs560_tnbc_sig, by = c("sample_name" = "samples"))
wgs560_tnbc_all$inflame <- wgs560_tnbc_all$T_effector_signature_TotalScore
write.table(wgs560_tnbc_all, "data_for_manuscript/wgs560_tnbc_all.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)

# ----- 1.3 Differential Expressed Genes -----
setwd("~/Documents/Research/WangLab/Projects/IGR/")
tnbc <- read.table("data_for_manuscript/wgs560_tnbc_all.tsv", sep = "\t", head = TRUE)
tnbc <- tnbc[! is.na(tnbc$sqrt_igr), ] %>% droplevels()
tnbc$igr_level <- ifelse(tnbc$sqrt_igr >= quantile(tnbc$sqrt_igr, 0.5), "High", "Low")
tnbc$tmb_level <- ifelse(tnbc$substitutions >= quantile(tnbc$substitutions, 0.5), "TMB_High", "TMB_Low")
tnbc$levels <- paste0(tnbc$igr_level, "/", tnbc$tmb_level)
tnbc$sqrt_scna = sqrt(tnbc$Intergenic.Intro.chr.n. + tnbc$Intergenic.AGT.n.)
tnbc$sqrt_scna2 = sqrt(tnbc$rearrangements)
tnbc$sqrt_tmb = sqrt(tnbc$substitutions)
tnbc$sqrt_indel = sqrt(tnbc$insertion.deletions)
m1 = lm(inflame ~ sqrt_tmb, data = tnbc)
m2 = lm(inflame ~ sqrt_tmb+sqrt_indel, data = tnbc)
m3 = lm(inflame ~ sqrt_tmb+sqrt_indel+sqrt_igr, data = tnbc)
m4 = lm(inflame ~ sqrt_tmb+sqrt_indel+sqrt_igr+sqrt_scna2, data = tnbc)
m5 = lm(inflame ~ sqrt_tmb+sqrt_indel+sqrt_igr+sqrt_scna2+HRD, data = tnbc)
anova(m1, m2, m3, m4, m5)

tnbc_subtype = read.csv("~/Downloads/TNBCType.c36afc01-11a3-4f59-9459-df47e4975af7_result.csv")

basal = tnbc$sqrt_igr[tnbc$PAM50.subtype == "Basal"]
nonbasal = tnbc$sqrt_igr[tnbc$PAM50.subtype != "Basal"]
write.table(na.omit(basal), "~/Documents/basal.txt", quote = F, row.names = F, col.names = F)
write.table(nonbasal, "~/Documents/nonbasal.txt", quote = F, row.names = F, col.names = F)

library(bnlearn)
tmp = tnbc[, c("HRD", "inflame", "sqrt_igr", "sqrt_scna", "sqrt_tmb", "sqrt_indel")]
for (i in 1:ncol(tmp)) tmp[, i] = as.numeric(tmp[, i])


# DE analysis
# since no count data available, we will use log2 fpkm
tnbc_log2_fpkm <- ReadExpr("data_for_manuscript/wgs560_tnbc_expr_log2_fpkm.tsv")
tnbc_de <- tnbc[match(names(tnbc_log2_fpkm), tnbc$sample_name), ] %>% droplevels()
all(tnbc_de$sample_name == names(tnbc_log2_fpkm))
num_na <- apply(tnbc_log2_fpkm2, 1, function(x) length(which(x == 0)))
# 10% of 73 = 8
tnbc_log2_fpkm <- tnbc_log2_fpkm2[num_na <= 8, ]

d0 <- tnbc_log2_fpkm
tnbc_tpm <- fpkm2tpm(2^tnbc_log2_fpkm)
group <- tnbc_de$igr_level
mm <- model.matrix(~ 0 + group)
#y <- voom(2^d0, mm, plot = F)
fit <- lmFit(log2(tnbc_tpm+0.01), mm)
#fit <- lmFit(tnbc_log2_fpkm, mm)
grp_contrast <- paste0("groupHigh-groupLow")
contr <- makeContrasts(contrasts = grp_contrast, levels = colnames(coef(fit)))
fit2 <- contrasts.fit(fit, contr)
tmp <- eBayes(fit2)
top_table <- topTable(tmp, sort.by = "none", n = Inf, coef = 1)
names(top_table) <- c("logFC", "AveExpr", "t", "p.value", "q.value", "B")
top_table <- top_table[order(top_table$logFC), ]

# Pathway analysis
library(fgsea)
library(msigdbr)
atab <- top_table

gene_pool <- atab[atab$p.value < 0.05, ]
geneRanks <- gene_pool$logFC
names(geneRanks) <- rownames(gene_pool)

hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list_hallmark <- split(x = hallmark_gene_sets$gene_symbol, 
                              f = hallmark_gene_sets$gs_name)
hallmark <- fgsea(pathways = msigdbr_list_hallmark, geneRanks, minSize = 15, 
                 maxSize = 500, nperm = 10000)

h <- hallmark[hallmark$padj < 0.05, ]

# ----- Figure 2a -----
# mean barplots with 95% CI shows why choose TNBC over non-TNBC
# see PRISM result ["~/Projects/IGR/data_for_manuscript/IGR.pzfx"]

# ----- Figure 2b -----
# stacked barplots showing IGR is associated with cisplatin and PARP inhibitors?
# see PRISM result ["~/Projects/IGR/data_for_manuscript/IGR.pzfx"]
tnbc$Lymphocyte_infiltration2 <- ifelse(tnbc$Lymphocyte_infiltration == "nil" | tnbc$Lymphocyte_infiltration == "mild", 
                                        "nil+mild", tnbc$Lymphocyte_infiltration)
fig2b <- table(tnbc$igr_level, tnbc$Lymphocyte_infiltration)
apply(fig2b, 1, function(x) x/sum(x))

tnbc$mitotic_score2 <- ifelse(tnbc$mitotic_score == "1" | tnbc$mitotic_score == "2", 
                                        "<3", tnbc$mitotic_score)
fig2b <- table(tnbc$igr_level, tnbc$mitotic_score2)
apply(fig2b, 1, function(x) x/sum(x))

# ----- Figure 2c -----
# boxplots
for (i in c("HRD", "T_cell_CD8.", "inflame", "Macrophage_M1", "total_mitoses", "T_cell_CD4._memory_activated", "Macrophage_M2")) {
    for (j in unique(tnbc$levels)) {
        tmp <- tnbc[tnbc$levels == j, i]
        write.table(tmp, paste0("~/Documents/", i, "_", gsub("/", "-", j), ".txt"),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
}
idx2 <- which(tnbc$total_mitoses != "no_data_supplied")
wilcox.test(as.numeric(tnbc$total_mitoses[idx2]) ~ tnbc$igr_level[idx2])

# ----- Figure S5c -----
# igr vs HRD + inflame vs HRD
# see PRISM result ["~/Projects/IGR/data_for_manuscript/IGR.pzfx"]
write.table(tnbc$HRD, "~/Documents/HRD.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
write.table(tnbc$inflame, "~/Documents/inflame.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
write.table(tnbc$sqrt_igr, "~/Documents/sqrt_igr.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
cor(tnbc$HRD, tnbc$inflame)
cor(tnbc$HRD, tnbc$sqrt_igr)
figS5C <- data.frame(hrd = tnbc$HRD,
                    inflame = tnbc$inflame,
                    igr = tnbc$sqrt_igr)
figS5C <- ggplot(figS5C, aes(x = hrd, y = inflame)) + 
    geom_point(size = 3, col = "lightgrey") + 
    geom_smooth(method = "lm", se = FALSE, lty = 2, col = "black") + 
    labs(x = "HRD Score", y = "Inflame Signature") + 
    theme_classic() + 
    theme(legend.title = element_text(size = 15, face = "bold"), 
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 15), 
          axis.title = element_text(size = 18)) + 
    annotate("text", 0, 0, label = "cor=0.108\np=0.361", size = 7, hjust = 0)
ggsave("figures/figS5C.pdf", fig2d1, units = "in", width = 7, height = 7)

# ----- Figure S5f -----
# heatmap
maria = list()
maria$CD8 = c("CD8A", "GZMA", "GZMB", "PRF1", "CXCL9", "CXCL10", "TBX21")
maria$APM = c("TAP1", "TAP2", "B2M", "HLA-A", "HLA-B", "HLA-C")
maria$ckpt = c("CD274", "PDCD1LG2", "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT")
maria$ddr = c("BRCA2", "ERCC2", "ERCC4", "FANCA", "FANCB", "FANCD2", "PALB2", "POLE", "RAD51C")
maria$cc = c("MKI67", "CCNE1", "BUB1", "BUB1B", "CCNB2", "CDC25C", "CDK2", "MCM4", "MCM6", "MCM2")
maria$angio = c("TEK", "CDH5", "SOX17", "SOX18")
maria$emt = c("CLDN3", "CLDN7", "CLDN4", "CDH1", "VIM", "TWIST1", "ZEB1", "ZEB2")
maria$TGFb = c("TGFB1", "TGFBR2")
maria = as.character(unlist(maria))
lm22 <- read.table("script/LM22.txt", head = TRUE, sep = "\t")
lm22 <- c(lm22$Gene.symbol, maria)
tnbc_tpm_lm22 <- tnbc_tpm[rownames(tnbc_tpm) %in% lm22, ]
plot.matrix = t(scale(t(tnbc_tpm_lm22), center = TRUE, scale = TRUE))
all(colnames(plot.matrix) == tnbc_de$sample_name)
quantile.range <- quantile(as.matrix(plot.matrix), probs = seq(0, 1, 0.01))
myBreaks <- seq(quantile.range["10%"], quantile.range["90%"], 0.1)
myColor  <- colorRampPalette(c("skyblue", "white", "red"))(length(myBreaks) - 1)
annotation = data.frame(IGR = tnbc_de$igr_level,
                        TIL = ifelse(tnbc_de$Lymphocyte_infiltration == "no_data_supplied", "NA", as.character(tnbc_de$Lymphocyte_infiltration)),
                        Grade = ifelse(tnbc_de$tumour_grade == "III", "III", "II or NA"),
                        Mitotic = ifelse(tnbc_de$mitotic_score == "no_data_supplied", "NA", tnbc_de$mitotic_score))
rownames(annotation) = tnbc_de$sample_name
var1 =  c("#FB0106", "#0F7FFE")
names(var1) = c("High", "Low")
var2 = c("#00798c", "#66a182", "#edae49", "#d1495b", "grey")
names(var2) = c("nil", "mild", "moderate", "severe", "NA")
var3 = c("#2e4057", "#8d96a3")
names(var3) = c("III", "II or NA")
var4 = c("#CCFFCC", "#66FF33", "#336600", "grey")
names(var4) = c("1", "2", "3", "NA")
anno_colors = list(IGR = var1, TIL = var2, Grade = var3, Mitotic = var4)

my_rank = array()
for (i in 1:nrow(tnbc_de)) {
    if (tnbc_de$igr_level[i] == "High" & tnbc_de$Lymphocyte_infiltration[i] == "severe") {
        my_rank[i] = 1
    } else if (tnbc_de$igr_level[i] == "High" & tnbc_de$Lymphocyte_infiltration[i] == "moderate") {
        my_rank[i] = 2
    } else if (tnbc_de$igr_level[i] == "High" & tnbc_de$Lymphocyte_infiltration[i] == "mild") {
        my_rank[i] = 3
    } else if (tnbc_de$igr_level[i] == "High" & tnbc_de$Lymphocyte_infiltration[i] == "nil") {
        my_rank[i] = 4
    } else if (tnbc_de$igr_level[i] == "High" & tnbc_de$Lymphocyte_infiltration[i] == "no_data_supplied") {
        my_rank[i] = 5
    } else if (tnbc_de$igr_level[i] == "Low" & tnbc_de$Lymphocyte_infiltration[i] == "severe") {
        my_rank[i] = 6
    } else if (tnbc_de$igr_level[i] == "Low" & tnbc_de$Lymphocyte_infiltration[i] == "moderate") {
        my_rank[i] = 7
    } else if (tnbc_de$igr_level[i] == "Low" & tnbc_de$Lymphocyte_infiltration[i] == "mild") {
        my_rank[i] = 8
    } else if (tnbc_de$igr_level[i] == "Low" & tnbc_de$Lymphocyte_infiltration[i] == "nil") {
        my_rank[i] = 9
    } else if (tnbc_de$igr_level[i] == "Low" & tnbc_de$Lymphocyte_infiltration[i] == "no_data_supplied") {
        my_rank[i] = 10
    }
}
tnbc_de$rank = my_rank
tnbc2 <- tnbc_de[sort.list(tnbc_de$rank), ]
#tnbc2 <- tnbc_de[sort.list(tnbc_de$sqrt_igr, decreasing = TRUE), ]
plot.matrix2 <- plot.matrix[, match(tnbc2$sample_name, colnames(plot.matrix))]
all(colnames(plot.matrix2) == tnbc2$sample_name)
plot.matrix.p <- array()
for (i in 1:nrow(plot.matrix2))
    plot.matrix.p[i] <- wilcox.test(as.numeric(tnbc_tpm_lm22[i, ]) ~ tnbc_de$igr_level, 
                                    alternative = "greater")$p.value
names(plot.matrix.p) <- rownames(plot.matrix)
#interested_genes <- c(names(plot.matrix.p[p.adjust(plot.matrix.p) < 1]), 
#                      "CD8A", "CXCL10", "PDCD1LG2", "CD274")
interested_genes <- names(plot.matrix.p[maria2][plot.matrix.p[maria2] < 0.05])
interested_genes <- c("CD274","CTLA4","LAG3","HAVCR2","PDCD1","PDCD1LG2","CD160","CD244","TIGIT","ENTPD1","BTLA","IFNG","FOXP3","CD8A","GZMB")
interested_genes = plot.matrix.p[plot.matrix.p < 0.1 & names(plot.matrix.p) %in% interested_genes] %>% names()
interested_genes <- interested_genes[interested_genes != "CTLA4"]
plot.matrix2 <- plot.matrix2[rownames(plot.matrix2) %in% interested_genes, ]
my_rownames <- ifelse(rownames(plot.matrix2) %in% interested_genes, 
                      rownames(plot.matrix2), "")

# ----- Figure 2e -----
# see PRISM result ["~/Projects/IGR/data_for_manuscript/IGR.pzfx"]
# pathway barplot
pathway_p <- -log10(h$padj) * sign(h$NES)
names(pathway_p) <- h$pathway


# ----- Figure S5g -----
# gsea results of significant pathways
for (i in h$pathway) {
    tmp <- paste0("p.adj=", round(h[h$pathway == i, "padj"], digits = 4), 
                  "; NES=", round(h[h$pathway == i, "NES"], digits = 4))
    a <- plotEnrichment(msigdbr_list_hallmark[[i]], geneRanks) + 
        labs(title = gsub("_", " ", gsub("HALLMARK_", "", i)), x = "Gene Rank", 
             subtitle = tmp) + 
        theme(plot.title = element_text(size = 16), 
              axis.title = element_text(size = 16))
    ggsave(paste0("figures/gsea/", i, ".pdf"), a, units = "in", width = 5, height = 3)
}


# 
fig2e <- ggplot(tnbc, aes(x = HRD, y = inflame, col = sqrt_igr)) + 
    geom_point(size = 4.5, alpha = 0.5) + 
    labs(x = "HRD Score", y = "T-Inflamed Signature", col = "IGR Burden") + 
    scale_color_gradient2(low = "#0F7FFE", high = "#FB0106", mid = "white", 
                          midpoint = median(tnbc$sqrt_igr),
                          oob = scales::squish,
                          limits = quantile(tnbc$sqrt_igr, c(0.2, 0.8)),
                          breaks = quantile(tnbc$sqrt_igr, c(0.2, 0.5, 0.8)),
                          labels = c("Lowest", "Median", "Highest")) + 
    theme_classic() + 
    theme(legend.position = c(1, 0.05),
          legend.title = element_text(size = 18, face = "bold"), 
          legend.text = element_text(size = 15),
          legend.justification = c(1, 0), 
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18))
ggsave("figures/fig2e.pdf", fig2e, units = "in", width = 8, height = 4.5)
