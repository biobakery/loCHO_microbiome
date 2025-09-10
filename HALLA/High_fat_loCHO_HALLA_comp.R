#################################
# Difference in correlations 
# from Halla
#
# Adapted code from Hanseul Kim
#
#################################

rm(list = ls(all.names = TRUE))

library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggcorrplot)
library(pheatmap)
library(gplots)
library(plotrix)
library(ggridges)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggcorrplot")
#devtools::install_github("nicolash2/ggdendroplot")

#Read in the two files from halla
shift_colors <- c("+,+"="#7F58Af",
                  "-,+"="#64C5EB",
                  "+,-"="#E84D8A",
                  "-,-"="#FEB326")

High_carb_assocations <-read.table('HALLA/RES_FILT/Carb_res/all_associations.txt', header=T, sep = "\t", quote="",
                                   comment.char="")


High_fat_assocations <-read.table('HALLA/RES_FILT/Fat_res/all_associations.txt', header=T, sep = "\t",
                                  quote="", comment.char = "")


#join the datasets and get absolute difference and signs
merged_data <- full_join(High_carb_assocations,High_fat_assocations,by=c("X_features","Y_features"))

#set na to 0 
merged_data$association.x[which(is.na(merged_data$association.x))] <- 0
merged_data$association.y[which(is.na(merged_data$association.y))] <- 0

merged_data <- merged_data %>%
  mutate(absdiff=abs(association.x-association.y)) %>% mutate(sign=association.x*association.y) %>%
  filter(!is.na(absdiff)) %>% 
  mutate(group=case_when(association.x>=0 & association.y>=0 ~"+,+", association.x<=0 & association.y<=0 ~"-,-",
                         association.x>=0 & association.y<=0 ~"+,-", association.x<=0 & association.y>=0 ~"-,+")) %>%
  mutate(groupinnum=case_when(association.x>=0 & association.y>=0 ~1, association.x<=0 & association.y<=0 ~2,
                              association.x>=0 & association.y<=0 ~3, association.x<=0 & association.y>=0 ~4)) 

#significance calculation
#https://janhove.github.io/analysis/2014/10/28/assessing-differences-of-significance
library(psych)
for(i in 1:nrow(merged_data)){ 
  summary_test<-r.test(n = 35, r12 = merged_data$association.x[i], 
                       n2 = 35, r34 = merged_data$association.y[i])
  merged_data$difference_p[i]<-as.numeric(summary_test$p)
}

merged_data$fdr_qval<-p.adjust(merged_data$difference_p,method="fdr")

to_spread<-merged_data %>% select(X_features,Y_features,groupinnum)

to_cluster <- spread(to_spread, Y_features, groupinnum) %>% column_to_rownames("X_features") %>% as.matrix()

rowclust1 = hclust(dist(to_cluster))
colclust1 = hclust(dist(t(to_cluster)))
clustered_df = to_cluster[rowclust1$order,colclust1$order]

merged_data_sig <- merged_data %>% filter(difference_p<0.05) 


#20*12
#difference heatmap with all
ggplot(merged_data, aes(x = X_features, y = Y_features, color = group)) + 
  geom_count(aes(size = absdiff)) +
  scale_size_area(max_size = 3) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12)) + 
  scale_y_discrete(limits = colnames(clustered_df)) + 
  scale_x_discrete(limits = rownames(clustered_df))

## Save the data
saveRDS(merged_data, "HALLA/RES_FILT/Fat_res/Fat_merged_res.RDS")

#difference heatmap with only significant values
ggplot(merged_data_sig, aes(x = X_features,y = Y_features,color=group)) + 
  geom_count(aes(size = absdiff)) +
  scale_size_area(max_size = 3) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_discrete(limits = colnames(clustered_df)) + scale_x_discrete(limits = rownames(clustered_df))

### we can filter out metabolites that are not significant in 


merged_data_sig_effect_filt <- merged_data_sig %>% filter(absdiff > 0.7)

keep_metas <- merged_data_sig_effect_filt %>% group_by(Y_features) %>% summarize(number=n()) %>% filter(number > 5) %>% pull(Y_features)


merged_data_sig_effect_filt <- merged_data_sig_effect_filt %>% filter(Y_features %in% keep_metas)


keep_taxa <- merged_data_sig_effect_filt %>% group_by(X_features) %>% summarize(number=n()) %>% filter(number > 5) %>%
  pull(X_features)

merged_data_sig_effect_filt <- merged_data_sig_effect_filt %>% filter(X_features %in% keep_taxa)

## need to cluster these
to_spread<-merged_data_sig_effect_filt %>% select(X_features,Y_features,groupinnum)

to_cluster <- spread(to_spread, Y_features, groupinnum) %>% column_to_rownames("X_features") %>% as.matrix()
to_cluster[is.na(to_cluster)] <- 0

rowclust1 = hclust(dist(to_cluster))
colclust1 = hclust(dist(t(to_cluster)))
clustered_df = to_cluster[rowclust1$order,colclust1$order]


final_plot <- ggplot(merged_data_sig_effect_filt, aes(x = X_features,y = Y_features,color=group)) + 
  geom_count(aes(size = absdiff)) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_size_continuous(name="Absolute difference in correlation") +
  scale_color_manual(values=shift_colors, name="Sign difference") +
  ylab("Metabolites") +
  xlab("Microbes") +
  scale_y_discrete(limits = colnames(clustered_df)) + scale_x_discrete(limits = rownames(clustered_df))


final_plot
ggsave(filename = "HALLA/RES_FILT/multi_omics_overview_High_fat.pdf",
       height=6, width=12)

