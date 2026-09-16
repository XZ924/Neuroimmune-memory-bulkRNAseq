library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(ggrepel)
library(biomaRt)
library(fgsea)
library(writexl)
library(readxl)
library(glue)
library(patchwork)


# read in and format data -------------------------------------------------

## read in and format count data
cts <- data.frame(read_csv("data/all_expected_counts.csv", show_col_types = FALSE))
cat("Number of duplicated genes: ", sum(duplicated(cts$ensembl_id))) # check id's are unique
gene_names <- cts[1:2] # store gene name mappings
cts <- cts[3:ncol(cts)]
cts <- round(cts) # round counts to integers for future normalization
rownames(cts) <- gene_names$ensembl_id

## read in and format sample information 
sample_info <- data.frame(read_csv("data/sokol_sample_names_updated.csv", show_col_types = FALSE))
treatment <- str_sub(sample_info$Sample_Plate, end=-3)
coldata <- data.frame(treatment, row.names= sample_info$Sample_Name)
coldata$treatment <- factor(coldata$treatment)

# check that columns in count matrix match rows in metadata
cat('\nMetadata aligns with counts: ')
all(rownames(coldata) %in% colnames(cts)) & all(rownames(coldata) == colnames(cts))


# run PCA -----------------------------------------------------------------

# create deseq2 object
dds <- DESeqDataSetFromMatrix(countData=cts, 
                              colData=coldata,
                              design = ~treatment)

# normalize counts
dds <- estimateSizeFactors(dds)
norm_cts <- counts(dds, normalized=TRUE)

# moderate variance across mean of norm counts for improved clustering distances
vsd <- vst(dds, blind=TRUE) 
vsd_mat <- assay(vsd)

# run PCA
pca <- prcomp(t(vsd_mat))
pca_df <- cbind(coldata, pca$x) 
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pca_df <- pca_df %>%
  rownames_to_column(var = "Sample_Name") %>%
  as_tibble() %>%
  left_join(sample_info)
  
# plot by treatment type
ggplot(pca_df, aes(x=PC1, y=PC2, color=treatment)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample_Plate), show.legend = FALSE) +
  theme_bw() +
  scale_color_discrete(name="Treatment") + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  ggtitle("PC1 vs PC2 by Treatment Type")
ggsave('pca/treatment_pc1v2.png', width = 7, height = 5)

ggplot(pca_df, aes(x=PC3, y=PC4, color=treatment)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample_Plate), show.legend = FALSE) +
  theme_bw() +
  scale_color_discrete(name="Treatment") + 
  xlab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) +
  ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) +
  ggtitle("PC3 vs PC4 by Treatment Type")
ggsave('pca/treatment_pc3v4.png', width = 7, height = 5)


# run DE analysis ---------------------------------------------------------

# For NC reference level: OVA, Papain, Papain/QX314
dds$treatment <- relevel(dds$treatment, ref='NC') 
dds <- DESeq(dds)
plotDispEsts(dds)
resultsNames(dds)

ova_res_nc_ref <- results(dds, name = "treatment_OVA_vs_NC")
papain_res_nc_ref <- results(dds, name = "treatment_Papain_vs_NC")
papain_qx314_res_nc_ref <- results(dds, name = "treatment_PapainQX314_vs_NC")

# for Papain reference level: Papain/QX314
dds$treatment <- relevel(dds$treatment, ref='Papain') 
dds <- DESeq(dds)
resultsNames(dds)

papain_qx314_res_papain_ref <- results(dds, name = "treatment_PapainQX314_vs_Papain")
ova_res_papain_ref <- results(dds, name = "treatment_OVA_vs_Papain")

# For OVA reference level: Papain, Papain/QX314
dds$treatment <- relevel(dds$treatment, ref='OVA')
dds <- DESeq(dds)
resultsNames(dds)

papain_qx314_res_ova_ref <- results(dds, name = "treatment_PapainQX314_vs_OVA")

# save results to excel
format_sig_deres <- function(de_res, sig_level = 0.1){
  de_res %>%
    as.data.frame() %>%
    rownames_to_column(var = 'ensembl_id') %>%
    as_tibble() %>%
    left_join(gene_names) %>%
    relocate(ensembl_id, gene_symbol) %>%
    arrange(padj) %>%
    filter(padj <= sig_level)
}

write_xlsx(list("OVA_vs_NC" = format_sig_deres(ova_res_nc_ref),
                "Papain_vs_NC" = format_sig_deres(papain_res_nc_ref),
                "PapainQX314_vs_NC" = format_sig_deres(papain_qx314_res_nc_ref),
                "PapainQX314_vs_Papain" = format_sig_deres(papain_qx314_res_papain_ref),
                "OVA_vs_Papain" = format_sig_deres(ova_res_papain_ref),
                "PapainQX314_vs_OVA" = format_sig_deres(papain_qx314_res_ova_ref)), 
           "de_res/aggregated_de_results_padj_0.1.xlsx")

write_xlsx(list("OVA_vs_NC" = format_sig_deres(ova_res_nc_ref, 1),
                "Papain_vs_NC" = format_sig_deres(papain_res_nc_ref, 1),
                "PapainQX314_vs_NC" = format_sig_deres(papain_qx314_res_nc_ref, 1),
                "PapainQX314_vs_Papain" = format_sig_deres(papain_qx314_res_papain_ref, 1),
                "OVA_vs_Papain" = format_sig_deres(ova_res_papain_ref, 1),
                "PapainQX314_vs_OVA" = format_sig_deres(papain_qx314_res_ova_ref, 1)), 
           "de_res/aggregated_de_results_all.xlsx")


# make volcano plots ------------------------------------------------------

plot_volcano <- function(res, comparison, title_name, padj_cutoff = 0.1, top_n=20){
  # label significant genes and get gene symbols
  res <- res %>%
    data.frame() %>%
    dplyr::select(pvalue, padj, log2FoldChange) %>%
    rownames_to_column(var="ensembl_id") %>%
    as_tibble() %>%
    drop_na(padj) %>%
    mutate(case = factor(case_when(
      (padj < padj_cutoff) & (log2FoldChange < 0) ~ 'ref_up',
      (padj < padj_cutoff) & (log2FoldChange > 0) ~ 'comp_up',
      TRUE ~ "non_sig"), 
      levels= c('ref_up','comp_up','non_sig'))) %>%
    arrange(padj) %>%
    left_join(gene_names) %>% 
    mutate(genelabels ="")
  
  # label the top genes by padj
  sig_gene_names <- res$gene_symbol[1:top_n]
  res$genelabels[1:top_n] <- sig_gene_names 
  
  # additional annotation
  colors <- brewer.pal(5, "Set1")
  group_colors <- c('Papain' = colors[3], 
                    'OVA' = colors[2], 
                    'PapainQX314' = colors[4])
  comp_group <- comparison[1]
  ref_group <- comparison[2]

  # plot
  ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(colour = case)) +
    geom_text_repel(aes(label = genelabels), max.overlaps = Inf) +
    scale_color_manual(values = c(unname(group_colors[ref_group]), 
                                  unname(group_colors[comp_group]), 
                                  "darkgrey"),
                       labels = c(paste0("Up in ", ref_group), 
                                  paste0("Up in ", comp_group), 
                                  paste0("padj > ", padj_cutoff))) +
    theme_bw() +
    theme(legend.title = element_blank())  +
    labs(title = title_name, 
         xlab = "log2 fold change", 
         ylab = "-log10 p-value")
}

v1_papain_vs_ova <- plot_volcano(ova_res_papain_ref, 
             c('OVA', 'Papain'), 
             "Papain vs OVA", 
             padj_cutoff = 0.1, 
             top_n = sum(ova_res_papain_ref$padj <= 0.1, na.rm=T))
pdf("de_res/volcanoplot_papain_vs_ova.pdf", width = 7, height = 6)
print(v1_papain_vs_ova)
dev.off()

v2_papain_vs_qx314 <- plot_volcano(papain_qx314_res_papain_ref, 
             c('PapainQX314', 'Papain'), 
             "Papain vs PapainQX314", 
             padj_cutoff = 0.1, 
             top_n = sum(papain_qx314_res_papain_ref$padj <= 0.1, na.rm=T))
pdf("de_res/volcanoplot_papain_vs_qx314.pdf", width = 7, height = 6)
print(v2_papain_vs_qx314)
dev.off()


# get mouse to human gene map for gsea ------------------------------------

# select ensembl database and get links to human and mouse datasets (need to use old host to connect)
mouse_mart <- useMart(biomart='ensembl', 
                      dataset='mmusculus_gene_ensembl', 
                      host="https://dec2021.archive.ensembl.org")
human_mart <- useMart(biomart="ensembl", 
                      dataset="hsapiens_gene_ensembl",
                      host="https://dec2021.archive.ensembl.org")

# get mouse gene symbols from mouse ensembl dataset
mouse_genes <- getBM('mgi_symbol', mart=mouse_mart)  

# link the mouse and human ensembl datasets (translating to homology mapping) and get the hngc_symbol associated with mgi_symbol
mouse_to_human <- getLDS(attributes=c('ensembl_gene_id', 'mgi_symbol'),
                         filters = 'mgi_symbol',
                         values = mouse_genes$mgi_symbol,
                         mart = mouse_mart,
                         attributesL = c('ensembl_gene_id', 'hgnc_symbol'),
                         martL= human_mart,
                         uniqueRows=TRUE)
colnames(mouse_to_human) <- c('mouse_ensembl_id', 'mgi_symbol', 'human_ensembl_id', 'hgnc_symbol')


# run gsea ----------------------------------------------------------------

# load GSEA pathways into a named list
pathways <- gmtPathways("data/msigdb_symbols.gmt")

# get hallmark pathways
hallmark_pathways <- pathways[grep('HALLMARK', names(pathways))]
hallmark_df <- tibble(pathway_name = names(hallmark_pathways), 
                      pathway_genes = unlist(lapply(hallmark_pathways, str_c, collapse=',')))

run_gsea <- function(res){
  # get mgi symbols and stat statistics to rank by
  ranked_df <-  res %>%
    as.data.frame() %>%
    rownames_to_column(var = 'ensembl_id') %>%
    inner_join(gene_names) %>%
    left_join(mouse_to_human, 
              by=c('gene_symbol'='mgi_symbol'),
              relationship = "many-to-many") %>%
    dplyr::select(hgnc_symbol, stat) %>%
    na.omit() %>%
    filter(hgnc_symbol != "") %>%
    distinct() %>%
    group_by(hgnc_symbol) %>%
    summarize(stat=mean(stat))
  
  # format ranked dataset for GSEA
  ranks <- deframe(ranked_df)
  
  # look through hallmakr gene sets
  fgseaRes <- fgsea(pathways=hallmark_pathways, 
                      stats=ranks,
                      minSize=10,
                      nPermSimple=10000) %>%
      as_tibble() %>%
      arrange(desc(NES))
    
  # format pathway names
  fgseaRes$pathway_name <- sapply(str_split(fgseaRes$pathway, "_", n=2), tail, 1)
  return(fgseaRes)
}

format_gsea_result <- function(gsea_res){
  gsea_res %>% 
    left_join(hallmark_df, by = c('pathway' = 'pathway_name')) %>%
    dplyr::select(pathway_name, pval, padj, log2err, ES, NES, leadingEdge, size, pathway_genes) %>%
    rowwise() %>% 
    mutate(leadingEdge = str_c(leadingEdge, collapse = ',')) %>%
    arrange(padj)
}

# run for papain-vs-ova and papain-vs-qx314 comparisons
gsea_res_papain_vs_ova <- run_gsea(ova_res_papain_ref)
gsea_res_papain_vs_ova_formatted  <- format_gsea_result(gsea_res_papain_vs_ova)

gsea_res_papain_vs_qx314 <- run_gsea(papain_qx314_res_papain_ref)
gsea_res_papain_vs_qx314_formatted  <- format_gsea_result(gsea_res_papain_vs_qx314)

write_xlsx(list(papain_vs_ova_gsea = gsea_res_papain_vs_ova_formatted, 
                papain_vs_papainQX314 = gsea_res_papain_vs_qx314_formatted), 
           path = "gsea/all_hallmark_gsea_results.xlsx")


# make gsea dotplot -------------------------------------------------------

# read in gsea results
read_sheet <- function(file_name, sheet_name){
  read_xlsx(file_name, sheet = sheet_name) %>%
    mutate(sheet = sheet_name) %>%
    dplyr::select(pathway_name, padj, NES, sheet) %>%
    mutate(
      case = case_when(
        NES < 0 ~ 'down',
        NES > 0 ~ 'up'
      ),
      case = factor(case, levels = c('down', 'up')),
      NES_abs = abs(NES)) %>%
    arrange(sheet, padj, NES_abs) %>%
    mutate(pathway_name = factor(pathway_name, levels = .$pathway_name))
}

gsea_res_papain_vs_ova_formatted <- read_sheet("gsea/all_hallmark_gsea_results.xlsx", 'papain_vs_ova_gsea')
gsea_res_papain_vs_qx314_formatted <- read_sheet("gsea/all_hallmark_gsea_results.xlsx", 'papain_vs_papainQX314')

# set p-value range for plots
p_min <- min(c(gsea_res_papain_vs_ova_formatted$padj, 
               gsea_res_papain_vs_qx314_formatted$padj))
p_max <- 0.1

# make plots 
g1_pap_vs_ova <- ggplot(gsea_res_papain_vs_ova_formatted %>% 
         filter(padj <= 0.1) %>%
         mutate(pathway_name = as.character(pathway_name)), 
       aes(x = NES,
           y = pathway_name, 
           color = padj)) +
  geom_point(size=3) +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_discrete(limits = rev) +
  xlim(c(-3, 3)) +
  scale_size(range = c(1, 6)) +
  scale_color_viridis_c(limits = c(p_min, p_max), 
                        direction = -1) +
  theme_classic() +
  labs(title = 'Pap vs. Sham',
       y = 'Hallmark Pathway',
       x = "Normalized Enrichment Score")

g2_pap_vs_qx314 <- ggplot(gsea_res_papain_vs_qx314_formatted %>% 
                          filter(padj <= 0.1) %>%
                          mutate(pathway_name = as.character(pathway_name)), 
                        aes(x = NES,
                            y = pathway_name, 
                            color = padj)) +
  geom_point(size=3) +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_discrete(limits = rev) +
  xlim(c(-3, 3)) +
  scale_size(range = c(1, 6)) +
  scale_color_viridis_c(limits = c(p_min, p_max), 
                        direction = -1) +
  theme_classic() +
  labs(title = 'Pap vs. Pap/QX314',
       y = 'Hallmark Pathway',
       x = "Normalized Enrichment Score")

# combine plots
pdf('gsea/gsea_hallmark_dotplot.pdf', width = 10, height = 8)
g1_pap_vs_ova + 
  g2_pap_vs_qx314 + 
  plot_layout(guides = "collect", axis_titles = "collect_y")
dev.off()


# combine bulk-seq spreadsheets -------------------------------------------

# for gsea
gsea_res_papain_vs_ova_formatted <- read_xlsx("gsea/all_hallmark_gsea_results.xlsx", sheet='papain_vs_ova_gsea')
gsea_res_papain_vs_qx314_formatted <- read_xlsx("gsea/all_hallmark_gsea_results.xlsx", sheet='papain_vs_papainQX314')
all_gsea_res <- gsea_res_papain_vs_ova_formatted %>% 
  mutate(comparison = 'papain_vs_ova') %>%
  relocate(comparison) %>%
  bind_rows(gsea_res_papain_vs_qx314_formatted %>%
              mutate(comparison = 'papain_vs_papainQX314') %>%
              relocate(comparison))

# for de
de_res_papain_vs_ova <- read_xlsx("de_res/aggregated_de_results_all.xlsx", sheet = "OVA_vs_Papain")
de_res_papain_vs_qx314 <- read_xlsx("de_res/aggregated_de_results_all.xlsx", sheet = "PapainQX314_vs_Papain")
all_de_res <- de_res_papain_vs_ova %>%
  mutate(comparison = 'papain_vs_ova') %>%
  relocate(comparison) %>%
  bind_rows(de_res_papain_vs_qx314 %>%
              mutate(comparison = 'papain_vs_papainQX314') %>%
              relocate(comparison))

write_xlsx(list(de_results = all_de_res, 
                gsea_results = all_gsea_res), 
           path = "bulkseq_combined_results.xlsx")

