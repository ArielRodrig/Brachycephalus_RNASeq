library("devtools")
install_github("nkurniansyah/Olivia")
library(Olivia)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tableone)
library(ggrepel)

setwd() #Set this to your working forlder

expression<-read.csv("Gene TPM.csv",stringsAsFactors=FALSE) # Take normalized expression data,as outputed by 3DRNASeq package
rownames(expression)<-expression$gene
# Import coordinates in Peafowl perceptual space.
data<-read.csv("dorsal_coordinates_in_perceptually-corrected_colorspace.csv")
rownames(data)<-data$indiv


# Transcriptome-wide association analysis using Olivia (https://github.com/nkurniansyah/Olivia)


expression<-read.csv("Gene TPM.csv",stringsAsFactors=FALSE)
rownames(expression)<-expression$gene
# Import coordinates in Peafowl perceptual space.
data<-read.csv("dorsal_coordinates_in_perceptually-corrected_colorspace.csv")
rownames(data)<-data$indiv
rnaseq_count_matrix<-expression[-1]
metadata<-read.csv("samples_metadata.csv")
# Merge metdata and average predator space coordinates in one data frame
phenotype<-merge(metadata,data[1:4], by.x="field.nr", by.y="indiv")
rownames(phenotype)<-phenotype$sample
summary_phen<- summarize_phenotypes(pheno = phenotype,
                                    numeric_variables = c("x","y","lum","admix.prop"),
                                    strata = "color")
summary_phen

#phenotype<-phenotype[c(2:9)]									
IDs_both <- intersect(rownames(phenotype), colnames(rnaseq_count_matrix))
rnaseq_matrix <- rnaseq_count_matrix[, IDs_both] 
phenotypes <- phenotype[match(IDs_both,rownames(phenotype)),]

covariates_string <- "admix.prop"
resid_plot<- residual_plot(pheno = phenotype, 
                           traits = "x",
                           covariates_string = covariates_string)
resid_plot
median_norm <- median_normalization(rnaseq_matrix)
raw.range <- range(colSums(rnaseq_matrix))
norm.range <- range(colSums(median_norm))

par(mfrow = c(1,2))
hist(colSums(rnaseq_matrix), 
     xlim = raw.range, 
     main = "",
     xlab = "Before median normalization",
     ylab = "Library size distribution")
hist(colSums(median_norm), 
     xlim = norm.range, 
     main = "",
     xlab = "After median normalization",
     ylab = "Library size distribution")

# As we already have TPM normalized expression values, the previous step can be skipped.

# Apply a filter to remove genes with zero counts in <50% of samples

clean_count_matrix <- apply_filters(count_matrix = rnaseq_matrix, 
                                    median_min = 0, 
                                    expression_sum_min = 10, 
                                    max_min = 0, 
                                    range_min = 0, 
                                    prop_zero_max = 0.75)
									
pdf("densi-plots_before-and-after_normalization.pdf")
par(mfrow=c(1,2))
plot(density(log_replace_half_min(rnaseq_matrix)[,"ARO01_RNA1037"]),
     main="Before filtering", xlab="log trascript")
plot(density(log_replace_half_min(clean_count_matrix)[,"ARO01_RNA1037"]), 
     main="After filtering", xlab="log trascript")
dev.off()	 

log_counts<- log_replace_half_min(clean_count_matrix)
log_counts<-melt(log_counts)

pdf("box-plot_sample_expression_loghalf-min_replacements.pdf")
box_plot<- ggplot(log_counts, aes(x = Var2, y = value)) + 
                  stat_boxplot(aes(Var2, value), 
                  geom='errorbar', linetype=1, width=0.5)+ 
                  xlab("Sample")+ 
                  ylab("log(Transcripts)")+ 
                  geom_boxplot( aes(Var2, value),outlier.shape=1)+
                  stat_summary(fun = mean, geom = "errorbar", 
                                aes(ymax = ..y.., ymin = ..y..),
                                width = .75, linetype = "dotted") +
                 theme_bw()+
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
dev.off()

# Perform association analysis
set.seed(12)
covariates_string="admix.prop"
quantile_emp_trascript<-lm_count_mat_emp_pval(clean_count_matrix, 
                                              pheno=phenotypes, 
                                              trait="x",
                                              covariates_string=covariates_string,
                                              n_permute=100,
                                              log_transform = "log_replace_half_min",
                                              outcome_type ="continuous",
                                              gene_IDs=rownames(clean_count_matrix))
tophits <- quantile_emp_trascript[which(quantile_emp_trascript$bh_emp_pvals< 0.05),]
head(tophits)
tophits[order(tophits$bh_emp_pvals),]
write.csv(tophits,"DE_genes_Olivia.x-with-admix-prop_as_covariate.csv")
volcano <- volcano_plot(deg_res = quantile_emp_trascript,
                        significant_threshold = 0.05)
pdf("volcano_plot_olivia.x-scores.with-admix-prop_as_covariate.pdf")
volcano
dev.off()
ggsave("volcano_plot_olivia.x-scores.with-admix-prop_as_covariate.png")

# Similar analysis with the luminance values
covariates_string="admix.prop" 
quantile_emp_lum<-lm_count_mat_emp_pval(clean_count_matrix, 
                                              pheno=phenotypes, 
                                              trait="lum",
                                              covariates_string=covariates_string,
                                              n_permute=100,
                                              log_transform = "log_replace_half_min",
                                              outcome_type ="continuous",
                                              gene_IDs=rownames(clean_count_matrix))
tophits <- quantile_emp_lum[which(quantile_emp_lum$bh_emp_pvals< 0.05),]
head(tophits)
tophits[order(tophits$bh_emp_pvals),]
write.csv(tophits,"DE_genes_Olivia.x-with-admix-prop_as_covariate.csv")
volcano2 <- volcano_plot(deg_res = quantile_emp_multi,
                        significant_threshold = 0.05)
pdf("volcano_plot_olivia.multivariate.with-admix-prop_as_covariate.pdf")
volcano2
dev.off()

# Perform association analysis by multiple regression
covariates_string="admix.prop" 
quantile_emp_multi<-lm_mult_count_mat_emp_pval(clean_count_matrix, 
       pheno=phenotypes, traits=c("x","lum"), covariates_string=covariates_string,
       n_permute=100, log_transform = "log_replace_half_min", outcome_type ="continuous",
       gene_IDs=rownames(clean_count_matrix))

top_emp_multi<-quantile_emp_multi[quantile_emp_multi$fdr_bh< 0.05,]
top_emp_multi[order(top_emp_multi$adjLogFC_x),]
write.csv(top_emp_multi,"DE_genes_multiple-regression.csv")

