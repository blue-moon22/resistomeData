# Load packages
library(resistomeAnalysis)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(purrr)
library(grid)
library(stringr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(metafor)
library(openxlsx)

#### Read mapping data ####
# Read non-subsampled mapping data
df_map <- readMappingData("db/MAPPING_DATA/nonsubsampled_merged_bedtools.csv", without_US_duplicates = TRUE)
df_map_dup <- readMappingData("db/MAPPING_DATA/nonsubsampled_merged_bedtools.csv", without_US_duplicates = FALSE)
# Read subsampled mapping data
df_map_subsampled_all <- readMappingData("db/MAPPING_DATA/subsampled_all_merged_bedtools.csv", without_US_duplicates = TRUE)
# Read subsampled mapping data for ARG richness
df_map_subsampled_argrich <- readMappingData("db/MAPPING_DATA/subsampled_argrich_merged_bedtools.csv", without_US_duplicates = TRUE)
# Read subsampled saliva mapping data
df_map_sub_saliva <- readMappingData("db/MAPPING_DATA/subsampled_saliva_merged_bedtools.csv", without_US_duplicates = TRUE)
# Read subsampled dental mapping data
df_map_sub_dental <- readMappingData("db/MAPPING_DATA/subsampled_dental_merged_bedtools.csv", without_US_duplicates = TRUE)
# Read subsampled stool mapping data
df_map_sub_stool <- readMappingData("db/MAPPING_DATA/subsampled_stool_merged_bedtools.csv", without_US_duplicates = TRUE)

#### Figure 1a ####
# Select colours for graph
cols <- c("white", "grey20", "grey40", "grey60", "grey80")
names(cols) <- c("China", "Fiji", "Philippines", "Western Europe", "USA")

# Percentage saliva samples containing ARG class
df_map_pb_saliva_class <- joinProportionAndBootstrap(df_map_sub_saliva, "Drug.Class", B = 20)

# Generate figure
tiff(filename = "figures/Figure1a.tiff", width = 2200, height = 1000, res = 180)
plotPercentages(df_map_pb_saliva_class, cols) + xlab("ARG class") + ylab("% saliva samples")
dev.off()

#### Figure 1b ####
# Percentage saliva samples containing ARG mechanism
df_map_pb_saliva_mech <- joinProportionAndBootstrap(df_map_sub_saliva, "Resistance.Mechanism", B = 20)

# Generate figure
tiff(filename = "figures/Figure1b.tiff", width = 1000, height = 1000, res = 180)
plotPercentages(df_map_pb_saliva_mech, cols) + xlab("ARG mechanism") + ylab("% saliva samples")
dev.off()

#### Figure 1c ####
# Percentage dental samples contraining ARG class
df_map_pb_dental_class <- joinProportionAndBootstrap(df_map_sub_dental, "Drug.Class", B = 20)

# Generate figure
tiff(filename = "figures/Figure1c.tiff", width = 1250, height = 1000, res = 170)
plotPercentages(df_map_pb_dental_class, cols) + xlab("ARG class") + ylab("% dental samples")
dev.off()

#### Figure 1d ####
# Percentage dental samples containing ARG mechanism
df_map_pb_dental_mech <- joinProportionAndBootstrap(df_map_sub_dental, "Resistance.Mechanism", B = 20)

# Generate figure
tiff(filename = "figures/Figure1d.tiff", width = 700, height = 1000, res = 180)
plotPercentages(df_map_pb_dental_mech, cols) + xlab("ARG mechanism") + ylab("% dental samples")
dev.off()

#### Figure 1e ####
# Percentage stool samples contraining ARG class
df_map_pb_stool_class <- joinProportionAndBootstrap(df_map_sub_stool, "Drug.Class", B = 20)

# Generate figure
tiff(filename = "figures/Figure1e.tiff", width = 2500, height = 1000, res = 170)
plotPercentages(df_map_pb_stool_class, cols) + xlab("ARG class") + ylab("% stool samples")
dev.off()

#### Figure 1f ####
# Percentage stool samples containing ARG mechanism
df_map_pb_stool_mech <- joinProportionAndBootstrap(df_map_sub_stool, "Resistance.Mechanism", B = 20)

# Generate figure
tiff(filename = "figures/Figure1f.tiff", width = 1000, height = 1000, res = 180)
plotPercentages(df_map_pb_stool_mech, cols) + xlab("ARG mechanism") + ylab("% stool samples")
dev.off()

#### Supplementary Figure 1 ####
# Draw core ARGs for saliva samples
tiff("figures/Supplementary_Figure1a.tiff", width=2800, height=2400, res=150)
drawCoreARGs(df_map_sub_saliva, c(0.5, rep(1, 5), 0.5, rep(0, 4)), bar_label_size = 3, group_label_size = 3.8, B = 20)
dev.off()

# Draw core ARGs for dental samples
tiff("figures/Supplementary_Figure1b.tiff", width=1600, height=1200, res=150)
drawCoreARGs(df_map_sub_dental, c(0.5, rep(1, 4), 0.5, rep(0, 3)), bar_label_size = 2.5, group_label_size = 2.5, B = 20)
dev.off()

# Draw core ARGs for stool samples
tiff("figures/Supplementary_Figure1c.tiff", width=2500, height=2000, res=150)
drawCoreARGs(df_map_sub_stool, c(0.5, rep(1, 5), 0.5, rep(0, 4)), bar_label_size = 3, group_label_size = 3, B = 20)
dev.off()

#### Supplementary Figure 2 ####
# Create dendrogram of US duplicates
top_col <- rev(brewer.pal(8, "Spectral"))
coloured_labels <- c("stool" = top_col[1], "dorsum of tongue" = top_col[2], "buccal mucosa" = top_col[3], "dental" = top_col[4])

# Select US longitudinal samples
df_map_hmp <- df_map_dup[df_map_dup$Location == "USA",]
body_sites <- unique(df_map_hmp$sample_type)
for(i in 1:length(body_sites)){
  df_map_hmp_body_site <- df_map_hmp[df_map_hmp$sample_type == body_sites[i],]
  df_samples_one <- unique(df_map_hmp_body_site$Sample.name[df_map_hmp_body_site$Visit_Number == 1])
  df_samples_two <- unique(df_map_hmp_body_site$Sample.name[df_map_hmp_body_site$Visit_Number >= 2])
  df_samples_rm <- df_samples_one[!(df_samples_one %in% df_samples_two)]
  df_ids_rm <- df_map_hmp_body_site$ID[df_map_hmp_body_site$Sample.name %in% df_samples_rm &
                                         df_map_hmp_body_site$sample_type == body_sites[i]]
  df_map_hmp <- df_map_hmp[!(df_map_hmp$ID %in% df_ids_rm),]
}

tiff(filename = "figures/Supplementary_Figure2.tiff", width = 3000, height = 3000, res = 200)
createUSDendrogram(df_map_hmp, coloured_labels)
dev.off()

#### Figure 2a ####
# Run principle coordinate analysis
mds <- runPrincipleCoordinateAnalysis(df_map_subsampled_all)

# Plot resisostypes
tiff("figures/Figure2a.tiff", width = 1600, height = 1000, res = 220)
plotResistotypes(mds) + theme(axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=10))
dev.off()


#### Figure 2b ####
# Calculate proportion of samples
cluster_res <- mds %>% group_by(clusters, Location, sample_type) %>% summarise(n = n()) %>%
  group_by(clusters) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop_cluster = n/total_n*100, Location_sampletype = paste(sample_type, "-", Location))

# Label by sample types
cols <- c("grey", brewer.pal(9, "Blues")[c(4,8)], "yellowgreen", brewer.pal(9, "YlOrRd")[c(3,5,7,9)], brewer.pal(9, "RdPu")[c(3,5,7,9)])
tiff("figures/Figure2b.tiff", width = 1400, height = 1000, res = 220)
ggplot(cluster_res, aes(clusters, prop_cluster, fill=Location_sampletype)) +
  geom_bar(stat="identity") +
  scale_fill_manual("Body Site - Geographical Location", values = cols) +
  xlab("Resistotypes") + ylab("Percentage") +
  theme_classic()
dev.off()

#### Figure 2c ####
# Read resistance map use data
rm_use_card <- readResistanceMapUseData("db/ResistanceMap/resistanceMap_use_190221.csv")

df_map_unnest <- df_map %>%
  mutate(Drug.Class = strsplit(Drug.Class, ";")) %>%
  unnest(Drug.Class)

df_map_use <- df_map_unnest %>% group_by(Location, sample_type, Drug.Class) %>%
  summarise(average_rpkm = sum(rpkm)/length(rpkm), std_error = sqrt(var(rpkm)/length(rpkm))) %>%
  left_join(rm_use_card, by = c("Location", "Drug.Class" = "class")) %>%
  filter(!is.na(sum.DDD.Per.1000.Pop))

cols <- brewer.pal(length(unique(df_map_use$Drug.Class)), "Paired")
names(cols) <- unique(df_map_use$Drug.Class)
shape_values <- c(0, 1, 2, 4, 6)
names(shape_values) <- unique(df_map_use$sample_type)

plotUseGraph <- function(df_map_use, Location, cols, shape_values) {
  ggplot(df_map_use, aes(sum.DDD.Per.1000.Pop, average_rpkm, colour = Drug.Class, shape = sample_type)) +
    geom_point() +
    geom_errorbar(aes(ymin=average_rpkm-std_error, ymax=average_rpkm+std_error), width=.1, position=position_dodge(0.1)) +
    theme_classic() +
    scale_colour_manual("ARG Class", values = cols) +
    scale_shape_manual("Sample type", values = shape_values) +
    xlab("Defined Daily Dose Per 1000 Individuals") + ylab("mean RPKM") + ggtitle(Location)
}

locations <- unique(df_map_use$Location)
g <- list()
linearMods <- data.frame()
count <- 1
for (i in 1:length(locations)){
  g[[i]] <- plotUseGraph(df_map_use[df_map_use$Location == locations[i],], locations[i], cols, shape_values) +
    theme(legend.position = "none", axis.title.y = element_text(size=14), plot.title = element_text(size=18))

  df_map_location <- df_map_use[df_map_use$Location == locations[i],]
  sample_types <- unique(df_map_location$sample_type)
  for (j in 1:length(sample_types)){
    # Linear model
    linRes <- lm(log10(average_rpkm) ~ sum.DDD.Per.1000.Pop, data = df_map_location[df_map_location$sample_type == sample_types[j],])
    linearMods[count,1] <- summary(linRes)$r.squared
    linearMods[count,2] <- summary(linRes)$coefficients[1,4]
    linearMods[count,3] <- summary(linRes)$coefficients[2,4]
    linearMods[count,4] <- locations[i]
    linearMods[count,5] <- sample_types[j]
    count <- count + 1
  }
}
names(linearMods) <- c("r2 value", "intercept", "p-value", "Country", "Body Site")
write.csv(linearMods, "db/ResistanceMap/linear_regression_results.csv", row.names = FALSE)

g[[length(g)+1]] <- g_legend(plotUseGraph(df_map_use, "", cols, shape_values) +
                               theme(legend.title = element_text(size=16), legend.text = element_text(size=12)))

# Generate figure
tiff("figures/Figure2c.tiff", width = 1500, height = 750, res = 140)
grid.arrange(grobs = g, layout_matrix = rbind(c(1,2,5),c(3,4,5)))
dev.off()

#### Figure 3a ####
# Create dataset of paired samples
pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"),
                  c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"),
                  c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
df_map_pairs <- createPairedData(df_map, pair_list) # Remove efflux pumps and regulatory elements

# Plot boxplots which includes comparison with stool
plotMultipleRPKM <- function(df_map_pairs){
  graphs <- list()
  count <- 0
  unique_Location <- unique(df_map_pairs$Location)
  for(i in 1:length(unique_Location)){
    tmp <- df_map_pairs[df_map_pairs$Location == unique_Location[i],]
    unique_groups <- unique(tmp$group)
    for(j in 1:length(unique_groups)){
      count <- count + 1
      df_map_pairs_groups <- tmp[tmp$group == unique_groups[j],]
      graphs[[count]] <- plotRPKM(df_map_pairs_groups) + ylim(c(0, max(df_map_pairs$rpkm_log)))
    }
  }
  return(graphs)
}

stool_comp_graphs <- plotMultipleRPKM(df_map_pairs[df_map_pairs$group_mod %in% c("stool \nvs. dental", "stool \nvs. saliva", "stool \nvs. buccal mucosa", "stool \nvs. dorsum of tongue"),])
lay <- rbind(c(1, 2, 3), c(4, 5, 6), c(7, NA, NA))
tiff("figures/Figure3a.tiff", width = 1800, height = 1500, res = 230)
grid.arrange(grobs = stool_comp_graphs, layout_matrix = lay)
dev.off()

#### Supplementary Figure 3 ####
oral_comp_graphs <- plotMultipleRPKM(df_map_pairs[df_map_pairs$group_mod %in% c("dorsum of tongue \nvs. buccal mucosa", "dorsum of tongue \nvs. dental", "buccal mucosa \nvs. dental"),])
tiff("figures/Supplementary_Figure3.tiff", width = 1800, height = 500, res = 200)
grid.arrange(grobs = oral_comp_graphs, nrow = 1)
dev.off()

#### Figure 3b ####
# Get relative abundance for each country and sample type with more than one sample type
df_map_rel <- getRelativeAbundance(df_map[df_map$Location %in% c("China", "USA", "Fiji", "Western Europe"),])

# Create muliple graphs of relative abundance
uniq_location <- unique(df_map_rel$Location)
cols <- brewer.pal(length(unique(df_map_rel$Drug.Class.Alt)), "Set3")
names(cols) <- unique(df_map_rel$Drug.Class.Alt)
g <- list()
for(i in 1:length(uniq_location)){
  g[[i]] <- plotARGClassAbundance(df_map_rel[df_map_rel$Location == uniq_location[i],], cols) + ggtitle(gsub(" - ", "\n", uniq_location[i])) +
    theme(legend.position = "none", plot.title = element_text(size=18), axis.title = element_text(size=16), axis.text.x = element_text(size=14))
}
g[[length(g)+1]] <- g_legend(plotARGClassAbundance(df_map_rel, cols) + theme(legend.title = element_text(size=18), legend.text = element_text(size=12)))

# Plot relative abundance
lay <- rbind(c(1, 2, 4), c(3, 5, 5))
tiff("figures/Figure3b.tiff", width = 2800, height = 2000, res = 250)
grid.arrange(grobs = g, layout_matrix = lay)
dev.off()

#### Supplementary Figure 4 ####
# Get individual relative abundance of ARGs by ARG class
df_map_rel_ind <- getRelativeAbundanceIndividuals(df_map[df_map$Location %in% c("China", "USA", "Fiji", "Western Europe"),])

# Create multiple plots of relative abundance
uniq_location <- unique(df_map_rel_ind$Location)
cols <- brewer.pal(length(unique(df_map_rel_ind$Drug.Class.Alt)), "Set3")
names(cols) <- unique(df_map_rel_ind$Drug.Class.Alt)

# Plot muliple graphs
g <- list()
for(i in 1:length(uniq_location)){
  g[[i]] <- plotIndividualAbundance(df_map_rel_ind[df_map_rel_ind$Location == uniq_location[i],], cols) + ggtitle(gsub(" - ", "\n", uniq_location[i])) +
    theme(legend.position = "none", plot.title = element_text(size=18), axis.title = element_text(size=16))
}
g[[length(g)+1]] <- g_legend(plotIndividualAbundance(df_map_rel_ind, cols) + theme(legend.title = element_text(size=18), legend.text = element_text(size=12)))

# Plot relative abundance
lay <- rbind(c(1,2), c(3,4), c(5,5))
tiff("figures/Supplementary_Figure4.tiff", width = 2000, height = 1500, res = 140)
grid.arrange(grobs = g, layout_matrix = lay)
dev.off()

#### Supplementary Figure 5a ####
# Create RPKM heatmap for China
top_col <- brewer.pal(6, "Paired")
col_vector <- c("stool" = top_col[1], "saliva" = top_col[3], "dental" = top_col[5])

tiff(filename = "figures/Supplementary_Figure5a.tiff", width = 3000, height = 1500, res = 300)
createRPKMHeatmap(df_map, "China", col_vector)
dev.off()

#### Supplementary Figure 5b ####
# Create RPKM heatmap for US (including duplicates)
top_col <- brewer.pal(4, "Paired")
col_vector <- c("stool" = top_col[1], "saliva" = top_col[3])

tiff(filename = "figures/Supplementary_Figure5b.tiff", width = 3000, height = 1500, res = 300)
createRPKMHeatmap(df_map_dup, "Fiji", col_vector)
dev.off()

#### Supplementary Figure 5c ####
# Create RPKM heatmap for US (including duplicates)
top_col <- rev(brewer.pal(8, "Spectral"))
col_vector <- c("stool" = top_col[1], "dorsum of tongue" = top_col[2], "buccal mucosa" = top_col[3], "dental" = top_col[4])

tiff(filename = "figures/Supplementary_Figure5c.tiff", width = 3000, height = 1500, res = 300)
createRPKMHeatmap(df_map_dup, "USA", col_vector)
dev.off()

#### Supplementary Figure 5d ####
# Create RPKM heatmap for Western Europe
top_col <- brewer.pal(4, "Paired")
col_vector <- c("stool" = top_col[1], "saliva" = top_col[3])

tiff(filename = "figures/Supplementary_Figure5d.tiff", width = 3000, height = 1500, res = 300)
createRPKMHeatmap(df_map, "Western Europe", col_vector)
dev.off()

#### Supplementary Figure 6 ####
list_deseq <- list(China = list(sample_type = c("saliva", "stool", "dental")),
                   USA = list(sample_type = c("stool", "buccal mucosa", "dorsum of tongue", "dental")),
                   "Western Europe" = list(sample_type = c("saliva", "stool")),
                   Fiji = list(sample_type = c("saliva", "stool")))

# DeSeq2 analysis through each comparison
results_list <- list()
name <- list()
sample_comps <- list()
dds_list <- list()
count <- 0
for(i in 1:length(list_deseq)){
  for(k in 1:(length(list_deseq[[i]]$sample_type)-1)){
    for(l in (k+1):length(list_deseq[[i]]$sample_type)){
      count = count + 1
      deseq_output <- runDESeq2(df_map[df_map$Location %in% c("China", "USA", "Western Europe", "Fiji"),], Location = names(list_deseq[i]),
                                compare_samples = c(list_deseq[[i]]$sample_type[k], list_deseq[[i]]$sample_type[l]))
      results_list[[count]] <- deseq_output$results
      sample_comps[[count]] <- deseq_output$sample_comparison
      name[[count]] <- paste(names(list_deseq[i]), list_deseq[[i]]$sample_type[k], list_deseq[[i]]$sample_type[l], sep= "_")
      dds_list[[count]] <- deseq_output$dds
    }
  }
}

# Create volcano plots
class_names <- unique(df_map$Drug.Class_mod)
n <- length(class_names)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
class_colours = rev(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
names(class_colours) <- class_names
class_names <- class_names[!is.na(names(class_colours))]
locations <- sapply(unlist(name), function(x) strsplit(x, "_")[[1]][1])
sample_comps <- gsub("_", " ", gsub("\\.", " ", sample_comps))
titles <- map2_chr(.x=locations, .y=sample_comps, .f = function(.x, .y) paste0(.x, "\n", .y))
p <- list() # Print plots out individually and save in grob list
for(i in 1:length(results_list)){
  p[[i]] <- plotVolcano(results_list[[i]], df_map[df_map$Location %in% c("China", "USA", "Western Europe", "Fiji"),], titles[i], class_colours, label_size = 3) +
    theme(legend.position = "none") + ylim(c(0,max(sapply(results_list, function(x) max(-log10(x$padj[x$padj < 0.05 & x$padj > 0 & abs(x$log2FoldChange) > 2]), na.rm = TRUE)))))
}

# Create legend
g_plot <- ggplot(df_map, aes(Drug.Class_mod, rpkm, fill = Drug.Class_mod)) +
  geom_bar(stat="identity") +
  scale_fill_manual("ARG Class", breaks = names(class_colours), values = class_colours) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=20))

p[[length(p)+1]] <- g_legend(g_plot)

# Print volcano plots
lay <- rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,NA), c(12, 12, 12,12))
tiff("figures/Supplementary_Figure6.tiff", width = 4000, height = 5000, res = 100)
grid.arrange(grobs = p, layout_matrix = lay)
dev.off()

#### Figure 3c ####

# Meta-analysis
randomEffectModel <- function(results_list, sample_comps, name, dds_list, group){
  # Get shared differentially abundant genes between stool and oral
  shrd_gn <- Reduce(intersect, sapply(results_list[grep(group, sample_comps)], function(x) row.names(x)[x$padj < 0.05 & !is.na(x$padj)]))

  # Get effect sizes and variances to oral and stool comparisons
  meta_val <- list()
  for (i in 1:length(shrd_gn)){

    meta_val[[i]] <- map_df(.x=grep(group, sample_comps), function(.x) {
      results <- results_list[[.x]]
      dds <- dds_list[[.x]]
      norm_counts <- counts(dds, normalized=TRUE)

      variance <- 1/mean(norm_counts[rownames(norm_counts) == shrd_gn[i],]) + dispersions(dds)[row.names(results) == shrd_gn[i]]
      effct_sz <- results$log2FoldChange[row.names(results) == shrd_gn[i]]
      comparison <- name[.x]
      valence <- ifelse(effct_sz > 0, "pos", "neg")
      return(c(variance = variance, effct_sz = effct_sz, comparison = comparison, valence = valence))
    })
  }

  # Run random-effect model
  res <- lapply(meta_val, function(x) rma(yi=effct_sz, vi=variance, data=x, method = "REML"))
  names(res) <- shrd_gn
  res <- res[sapply(res, function(x) x$pval < 0.05)]
  res_df <- as.data.frame(t(map_df(res, function(x) c(estimate = as.numeric(x$b), ci.lb = x$ci.lb, ci.ub = x$ci.ub))))
  names(res_df) <- c("estimate", "ci.lb", "ci.ub")
  res_df$ARG <- row.names(res_df)

  return(res_df)
}

stool_saliva_res <- randomEffectModel(results_list, sample_comps, name, dds_list, "stool vs saliva")
stool_dental_res <- randomEffectModel(results_list, sample_comps, name, dds_list, "stool vs dental")

## Plot results
stool_saliva_res$sample_type <- ifelse(stool_saliva_res$estimate < 0, "saliva", "stool")
stool_saliva_res$ARG <- gsub("part_of.*", "part of efflux pump complex]", stool_saliva_res$ARG)
stool_saliva_res$ARG <- gsub("regulates.*", "regulates efflux pump]", stool_saliva_res$ARG)
stool_saliva_res$comparison <- "saliva vs. stool"

stool_dental_res$sample_type <- ifelse(stool_dental_res$estimate < 0, "dental", "stool")
stool_dental_res$ARG <- gsub("part_of.*", "part of efflux pump complex]", stool_dental_res$ARG)
stool_dental_res$ARG <- gsub("regulates.*", "regulates efflux pump]", stool_dental_res$ARG)
stool_dental_res$comparison <- "dental vs. stool"

all_res <- rbind(stool_saliva_res, stool_dental_res)
all_res$ARG <- gsub(" \\[.*", "", all_res$ARG)

tiff("figures/Figure3c.tiff", width = 2800, height = 1500, res = 250)
ggplot(all_res, aes(ARG, abs(estimate), colour = sample_type)) +
  geom_point() +
  geom_errorbar(aes(ARG, ymin=abs(ci.lb), ymax = abs(ci.ub))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, face="italic")) +
  ylab("Estimated average log2 fold change") +
  scale_colour_manual("Sample Type", values = brewer.pal(3,"Dark2")) +
  facet_grid(~comparison, scale = "free", space = "free")
dev.off()

#### Supplementary Figure 7 ####
# Get ARGs found exclusively in countries
stool_saliva_exclusive_args <- map_df(.x = grep("stool vs saliva", sample_comps), function(.x) data.frame(results_list[[.x]][results_list[[.x]]$padj < 0.05 & !is.na(results_list[[.x]]$padj),],
                                                                                                          ARG = row.names(results_list[[.x]])[results_list[[.x]]$padj < 0.05 & !is.na(results_list[[.x]]$padj)],
                                                                                                          Location = gsub("_.*", "", name[[.x]]),
                                                                                                          comparison = "saliva vs. stool"))
stool_saliva_exclusive_args <- stool_saliva_exclusive_args[!duplicated(stool_saliva_exclusive_args$ARG),]
stool_saliva_exclusive_args$sample_type <- ifelse(stool_saliva_exclusive_args$log2FoldChange < 0, "saliva", "stool")

stool_dental_exclusive_args <- map_df(.x = grep("stool vs dental", sample_comps), function(.x) data.frame(results_list[[.x]][results_list[[.x]]$padj < 0.05 & !is.na(results_list[[.x]]$padj),],
                                                                                                          ARG = row.names(results_list[[.x]])[results_list[[.x]]$padj < 0.05 & !is.na(results_list[[.x]]$padj)],
                                                                                                          Location = gsub("_.*", "", name[[.x]]),
                                                                                                          comparison = "dental vs. stool"))
stool_dental_exclusive_args <- stool_dental_exclusive_args[!duplicated(stool_dental_exclusive_args$ARG),]
stool_dental_exclusive_args$sample_type <- ifelse(stool_dental_exclusive_args$log2FoldChange < 0, "dental", "stool")

stool_saliva_dental_excl <- rbind(stool_saliva_exclusive_args, stool_dental_exclusive_args)
stool_saliva_dental_excl$ARG <- gsub("part_of.*", "part of efflux pump complex]", stool_saliva_dental_excl$ARG)
stool_saliva_dental_excl$ARG <- gsub("regulates.*", "regulates efflux pump]", stool_saliva_dental_excl$ARG)

# Plot results
stool_saliva_dental_excl$ARG <- gsub(" \\[.*", "", stool_saliva_dental_excl$ARG)
tiff("figures/Supplementary_Figure7.tiff", width = 5000, height = 1500, res = 250)
ggplot(stool_saliva_dental_excl, aes(ARG, abs(log2FoldChange), colour = sample_type, shape = Location)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, face = "italic")) +
  ylab("Log2 Fold Change") +
  scale_colour_manual("Sample Type", values = brewer.pal(3,"Dark2")) +
  scale_shape_manual("Only found in:", values = c(1,2,3,4)) +
  facet_grid(~comparison, scale = "free", space = "free")
dev.off()

#### Figure 4 ####
metadata_file <- "db/SAMPLES/metadata/supplementary_metadata.csv"
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Remove longitudinal metadata
metadata <- metadata[metadata$Visit_Number== 1,]

metadata_pairs <- data.frame()
for(i in 1:length(pair_list)){
  samples_one <- unique(metadata$Sample.name[metadata$sample_type == pair_list[[i]][1]])
  samples_two <- unique(metadata$Sample.name[metadata$sample_type == pair_list[[i]][2]])
  metadata_pair <- metadata[(metadata$Sample.name %in% Reduce(intersect,list(samples_one, samples_two))) & (metadata$sample_type %in% pair_list[[i]]),]

  metadata_pair <- cbind(metadata_pair, group = paste(pair_list[[i]][1], "vs.", pair_list[[i]][2]))
  metadata_pairs <- rbind(metadata_pairs, metadata_pair)
}

# Find lowest sequence number per group for subsampling (for Arg richness)
options(scipen=999)
subsample_summary <- metadata_pairs %>%
  filter(seq_num > 1) %>%
  group_by(Location, group) %>%
  summarise(num = n_distinct(ID), seq_num.x = signif(floor(min(seq_num)/0.1)*0.1, 2)*1e6) %>%
  mutate(sample_type = strsplit(as.character(group), " vs. ")) %>%
  unnest(sample_type) %>%
  mutate_if(is.factor, as.character)


# Create dataset of paired samples with efflux pumps
df_map_subsampled_argrich_pairs <- left_join(df_map_subsampled_argrich, subsample_summary, by = c("Location", "seq_num.x", "sample_type"))
sample_pairs_summary <- df_map_subsampled_argrich_pairs %>%
  group_by(Location, group, sample_type, Sample.name) %>%
  summarise(n = n_distinct(Sample.name)) %>%
  group_by(Location, group, Sample.name) %>%
  summarise(total_n = sum(n))
df_map_subsampled_argrich_pairs <- left_join(df_map_subsampled_argrich_pairs, sample_pairs_summary, by = c("Location", "group", "Sample.name")) %>%
  filter(total_n == 2)

summ <- df_map_subsampled_argrich_pairs %>%
  group_by(Location, group, sample_type, seq_num.x) %>%
  summarise(n = n_distinct(Sample.name))

# Run t-test
ttest_groups <- runTtest(df_map_subsampled_argrich_pairs)

# Plot boxplots which includes comparison with stool
graphs_comp <- plotMultipleARGRichnessGraphs(ttest_groups, df_map_subsampled_argrich_pairs)
lay <- rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,NA))
tiff("figures/Figure4.tiff", width = 2100, height = 1500, res = 200)
grid.arrange(grobs = graphs_comp, layout_matrix = lay)
dev.off()

#### Supplementary Figure 8 ####
# Create dataset of paired samples without efflux pumps
df_map_subsampled_argrich_pairs_noefflux <- df_map_subsampled_argrich_pairs[is.na(df_map_subsampled_argrich_pairs$Multi.Efflux.Role),]

sample_pairs_summary_noefflux <- df_map_subsampled_argrich_pairs_noefflux %>%
  group_by(Location, group, sample_type, Sample.name) %>%
  summarise(n = n_distinct(Sample.name)) %>%
  group_by(Location, group, Sample.name) %>%
  summarise(total_n_efflux = sum(n))
df_map_subsampled_argrich_pairs_noefflux <- left_join(df_map_subsampled_argrich_pairs_noefflux, sample_pairs_summary_noefflux, by = c("Location", "group", "Sample.name")) %>%
  filter(total_n_efflux == 2)

# Run t-test
ttest_groups_noefflux <- runTtest(df_map_subsampled_argrich_pairs_noefflux)

# Plot boxplots which includes comparison with stool
graphs_comp_noefflux <- plotMultipleARGRichnessGraphs(ttest_groups_noefflux, df_map_subsampled_argrich_pairs_noefflux)
lay <- rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,NA))
tiff("figures/Supplementary_Figure8.tiff", width = 2100, height = 1500, res = 200)
grid.arrange(grobs = graphs_comp_noefflux, layout_matrix = lay)
dev.off()

#### Figure 5a ####
# Combine metaphlan and rpkm for each sample
metaphlan <- read.csv("db/METAPHLAN/metaphlan.csv", stringsAsFactors = FALSE, row.names = 1)
metaphlan_rpkm <- combineMetaphlanandARG(df_map, metaphlan)

# Chinese healthy saliva with paired stool
china_sample_ids <- Reduce(intersect, list(unique(df_map$Sample.name[df_map$Location == "China" & df_map$sample_type == "saliva"]),
                                     unique(df_map$Sample.name[df_map$Location == "China" & df_map$sample_type == "stool"])))


# Get spearman's correlation
high_cor_china_saliva <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Sample.name %in% china_sample_ids & df_map$sample_type == "saliva"]), taxon_level = "t", taxon_ignore = "@@@")

# Plot heatmap
phyla <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Candidatus_Saccharibacteria", "Fusobacteria", "Spirochaetes",
            "Verrucomicrobia", "Ascomycota", "Synergistetes")
tiff("figures/Figure5a.tiff", width = 4500, height = 3500, res = 400)
drawCorrelationHeatmap(high_cor_china_saliva[high_cor_china_saliva$phylum %in% phyla,], 150, phyla = phyla, left_margin = 60)
dev.off()

#### Figure 5b ####
# Spearman's correlation of Philippines saliva
high_cor_philippines_saliva <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Location == "Philippines" & df_map$sample_type == "saliva"]), taxon_level = "t", taxon_ignore = "@@@")

# Plot heatmap
tiff("figures/Figure5b.tiff", width = 2500, height = 1800, res = 300)
drawCorrelationHeatmap(high_cor_philippines_saliva[high_cor_philippines_saliva$phylum %in% phyla,], 150, phyla = phyla, left_margin = 60)
dev.off()

#### Supplementary Figure 9 ####
# Get spearman's correlation
high_cor_china_stool <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Sample.name %in% china_sample_ids & df_map$sample_type == "stool"]), taxon_level = "t", taxon_ignore = "@@@")

# Plot heatmap
tiff("figures/Supplementary_Figure9.tiff", width = 3000, height = 5000, res = 400)
drawCorrelationHeatmap(high_cor_china_stool[high_cor_china_stool$phylum %in% phyla,], 150, 200, phyla = phyla)
dev.off()

#### Save source data ####
# Use openxlsx instead?
source_data = "Source Data.xlsx"
source_datasets <- list(df_map_pb_saliva_class, df_map_pb_saliva_mech, df_map_pb_dental_class, 
                        df_map_pb_dental_mech, df_map_pb_stool_class, df_map_pb_stool_mech, 
                        df_map_sub_saliva, df_map_sub_dental, df_map_sub_stool, df_map_hmp, mds, cluster_res,
                        df_map_use, 
                        df_map_pairs[df_map_pairs$group_mod %in% c("stool \nvs. dental", "stool \nvs. saliva", "stool \nvs. buccal mucosa", "stool \nvs. dorsum of tongue"),],
                        df_map_pairs[df_map_pairs$group_mod %in% c("dorsum of tongue \nvs. buccal mucosa", "dorsum of tongue \nvs. dental", "buccal mucosa \nvs. dental"),],
                        df_map_rel, df_map_rel_ind,
                        df_map[df_map$Location == "China",], df_map[df_map$Location == "Fiji",],
                        df_map[df_map$Location == "USA",], df_map[df_map$Location == "Western Europe",],
                        do.call("rbind", lapply(1:length(name), function(x) cbind(as.data.frame(results_list[[x]]), cohort = name[[x]]))),
                        all_res, stool_saliva_dental_excl, df_map_subsampled_argrich_pairs,
                        df_map_subsampled_argrich_pairs_noefflux,
                        high_cor_china_saliva[high_cor_china_saliva$phylum %in% phyla,],
                        high_cor_philippines_saliva[high_cor_philippines_saliva$phylum %in% phyla,],
                        high_cor_china_stool[high_cor_china_stool$phylum %in% phyla,])
#saveRDS(source_datasets, "source_data.RDS")

sheet_names <- c(paste0("Figure1", c("a", "b", "c", "d", "e", "f")), 
                 paste0("SupplementaryFigure1", c("a", "b", "c")),
                 "SupplementaryFigure2", paste0("Figure2", c("a", "b", "c")),
                 "Figure3a", "SupplementaryFigure3", "Figure3b", "SupplementaryFigure4",
                 paste0("SupplementaryFigure5", c("a", "b", "c", "d")), "SupplementaryFigure6",
                 "Figure3c", "SupplementaryFigure7", "Figure4", "SupplementaryFigure8", 
                 paste0("Figure5", c("a", "b")), "SupplementaryFigure9")

wb <- createWorkbook()
for (i in 1:length(sheet_names)){
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet = sheet_names[i], source_datasets[[i]])
}
saveWorkbook(wb, source_data, overwrite = TRUE)
