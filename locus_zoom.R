args <- commandArgs(trailingOnly = TRUE)
locus.frame <- read.delim(args[1], stringsAsFactors = F)
start.pos <- args[2]
end.pos <- args[3]
out.file.prefix <- args[4]

#Annotate the gene region using biomaRt
library(biomaRt) 
library(patchwork)
library(forcats)
library(dplyr)
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) # we will need an additional mart for genes
out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'), 
  filters = c('chromosome_name','start','end'), 
  values = list(1, start.pos, end.pos), 
  mart = gene.ensembl)
out.bm.genes.region <- out.bm.genes.region %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype), 
                                                                                     "protein_coding"), external_gene_name = fct_reorder2(external_gene_name, 
                                                                                                                                          start_position, gene_biotype_fac, .desc = TRUE))
plot.range <- c(start.pos, end.pos)

#Create the plots
library(ggplot2)
pdf(paste0(out.file.prefix, "_allele_dose.pdf"), width = 10, height = 10)
ggplot(data = locus.frame) + 
  geom_point(aes(position, -log10(allele_dose_p)), shape = 1) +
  ggplot(data = out.bm.genes.region) + 
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
  coord_flip() + ylab("") +
  ylim(plot.range) + 
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
  labs(title = "Allele", subtitle = paste0("Genes"), caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom", 
        panel.grid.major.y = element_blank()) + 
  expand_limits(y=c(-1, 1)) +
  scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) +
  plot_layout(ncol = 1, heights = c(6, 2))
dev.off()

pdf(paste0(out.file.prefix, "_afr_dose.pdf"), width = 10, height = 10)
ggplot(data = locus.frame) + 
  geom_point(aes(position, -log10(afr_dose_p)), shape = 1) +
  ggplot(data = out.bm.genes.region) + 
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
  coord_flip() + ylab("") +
  ylim(plot.range) + 
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
  labs(title = "African", subtitle = paste0("Genes"), caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom", 
        panel.grid.major.y = element_blank()) + 
  expand_limits(y=c(-1, 1)) +
  scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) +
  plot_layout(ncol = 1, heights = c(6, 2))
dev.off()

pdf(paste0(out.file.prefix, "_allele_afr_dose.pdf"), width = 10, height = 10)
ggplot(data = locus.frame) + 
  geom_point(aes(position, -log10(allele_afr_dose_p)), shape = 1) +
  ggplot(data = out.bm.genes.region) + 
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
  coord_flip() + ylab("") +
  ylim(plot.range) + 
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
  labs(title = "Allele-African", subtitle = paste0("Genes"), caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        strip.text.y = element_text(angle = 0),
        legend.position="bottom", 
        panel.grid.major.y = element_blank()) + 
  expand_limits(y=c(-1, 1)) +
  scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) +
  plot_layout(ncol = 1, heights = c(6, 2))
dev.off()