args <- commandArgs(trailingOnly = TRUE)
locus.frame <- read.delim(args[1], stringsAsFactors = F)
out.file.prefix <- args[2]
chr <- args[3]
start.pos <- args[4]
end.pos <- args[5]
width <- args[6]
height <- args[7]
base.font.size <- args[8]
gene.text.size <- args[9]
min.p <- args[10]

#Set optional parameters
if (is.na(start.pos) | is.na(end.pos)) {
  start.pos <- locus.frame$position[1]
  end.pos <- locus.frame$position[nrow(locus.frame)]
}
if (is.na(width)) {
  width <- 14
}
if (is.na(height)) {
  height <- 14
}
if (is.na(base.font.size)) {
  base.font.size <- 20
}
if (is.na(gene.text.size)) {
  gene.text.size <- 5
}
if (is.na(min.p)) {
  min.p < 0.001
}

#Annotate the gene region using biomaRt
library(biomaRt)
library(patchwork)
library(forcats)
library(dplyr)
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) # we will need an additional mart for genes
out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'),
  filters = c('chromosome_name','start','end'),
  values = list(chr, start.pos, end.pos),
  mart = gene.ensembl)
out.bm.genes.region <- out.bm.genes.region %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype),
                                                                                     "protein_coding"), external_gene_name = fct_reorder2(external_gene_name,
                                                                                                                                          start_position, gene_biotype_fac, .desc = TRUE))
plot.range <- c(start.pos, end.pos)

#Create the plots
library(ggplot2)
pdf(paste0(out.file.prefix, "_allele_dose.pdf"), width = width, height = height)
theme_set(theme_classic(base_size=base.font.size))
significant <- rep("grey", nrow(locus.frame))
significant[locus.frame$allele_dose_p <= min.p] <- "red"
ggplot(data = locus.frame) +
  geom_point(aes(position, -log10(allele_dose_p)), shape = 19, colour=significant) +
  ggplot(data = out.bm.genes.region) +
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
  coord_flip() + ylab("") +
  ylim(plot.range) +
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size=gene.text.size) +
  labs(title = "Allele", color = "Gene Biotype") +
  expand_limits(y=c(-1, 1)) +
  scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position="bottom",
        panel.grid.major.y = element_blank()) +
  plot_layout(ncol = 1, heights = c(6, 3))
dev.off()

pdf(paste0(out.file.prefix, "_afr_dose.pdf"), width = width, height = height)
theme_set(theme_classic(base_size=base.font.size))
significant <- rep("grey", nrow(locus.frame))
significant[locus.frame$afr_dose_p <= min.p] <- "red"
ggplot(data = locus.frame) +
  geom_point(aes(position, -log10(afr_dose_p)), shape = 19, colour=significant) +
  ggplot(data = out.bm.genes.region) +
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
  coord_flip() + ylab("") +
  ylim(plot.range) +
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size=gene.text.size) +
  labs(title = "African", color = "Gene Biotype") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position="bottom",
        panel.grid.major.y = element_blank()) +
  expand_limits(y=c(-1, 1)) +
  scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) +
  plot_layout(ncol = 1, heights = c(6, 3))
dev.off()

pdf(paste0(out.file.prefix, "_allele_afr_dose.pdf"), width = width, height = height)
theme_set(theme_classic(base_size=base.font.size))
significant <- rep("grey", nrow(locus.frame))
significant[locus.frame$allele_afr_dose_p <= min.p] <- "red"
ggplot(data = locus.frame) +
  geom_point(aes(position, -log10(allele_afr_dose_p)), shape = 19, colour=significant) +
  ggplot(data = out.bm.genes.region) +
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
  coord_flip() + ylab("") +
  ylim(plot.range) +
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size=gene.text.size) +
  labs(title = "Allele-African", color = "Gene Biotype") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position="bottom",
        panel.grid.major.y = element_blank()) +
  expand_limits(y=c(-1, 1)) +
  scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) +
  plot_layout(ncol = 1, heights = c(6, 3))
dev.off()


allele.frame <- locus.frame[,c("position", "allele_dose_p")]
names(allele.frame)[2] <- "p"
allele.frame$dose_type <- "Allele Dose"
afr.frame <- locus.frame[,c("position", "afr_dose_p")]
names(afr.frame)[2] <- "p"
afr.frame$dose_type <- "African Dose"
allele.afr.frame <- locus.frame[,c("position", "allele_afr_dose_p")]
names(allele.afr.frame)[2] <- "p"
allele.afr.frame$dose_type <- "Allele-African Dose"
merged.locus.frame <- merge(allele.frame, afr.frame, all.x=T, all.y=T)
merged.locus.frame <- merge(merged.locus.frame, allele.afr.frame, all.x=T, all.y=T)
merged.locus.frame$dose_type <- factor(merged.locus.frame$dose_type)
pdf(paste0(out.file.prefix, "_combined_dose.pdf"), width = width, height = height)
theme_set(theme_classic(base_size=base.font.size))
ggplot(data = merged.locus.frame) +
  geom_point(aes(position, -log10(p), colour=dose_type), shape = 19) +
  ggplot(data = out.bm.genes.region) +
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
  coord_flip() + ylab("") +
  ylim(plot.range) +
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size=gene.text.size) +
  labs(title = out.file.prefix, subtitle = paste0("Genes"), caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position="bottom",
        panel.grid.major.y = element_blank()) +
  expand_limits(y=c(-1, 1)) +
  scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) +
  plot_layout(ncol = 1, heights = c(6, 3))
dev.off()
