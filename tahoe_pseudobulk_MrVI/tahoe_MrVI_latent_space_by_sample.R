library(tidyverse)

col_data <- read_csv("tahoe_pseudobulk_colData.csv")

u_space <- read_csv("tahoe_pseudobulk_MrVI_u_factor.csv",
                    col_names = F)

u_space_for_print <- col_data %>%
  select(cell_drug) %>%
  bind_cols(u_space)
write_csv(u_space_for_print, "tahoe_pseudobulk_u-space_by_sample.csv")

col_data <- read_csv("tahoe_pseudobulk_colData.csv")

z_space <- read_csv("tahoe_pseudobulk_MrVI_z_factor.csv",
                    col_names = F)

z_space_for_print <- col_data %>%
  select(cell_drug) %>%
  bind_cols(z_space)
write_csv(z_space_for_print, "tahoe_pseudobulk_z-space_by_sample.csv")

# ================================================================================
# UMAPs by tissue
# ================================================================================
library(monocle3)
col_data <- read_csv("tahoe_pseudobulk_colData.csv")
cell_line_meta <- read_csv("cell_line_metadata.csv") %>%
  select(cell_name, Organ, Cell_ID_Cellosaur) %>%
  distinct()

col_data <- left_join(col_data, cell_line_meta, by = c("Cell_ID_Cellosaur"))

u_space <- read_csv("tahoe_pseudobulk_MrVI_u_factor.csv",
                    col_names = F)

cds_temp <- new_cell_data_set(expression_data = t(as.matrix(u_space)),
                              cell_metadata = col_data)

row.names(u_space) <- col_data %>% select(cell_drug) %>% pull()

reducedDims(cds_temp)[["PCA"]] <- u_space

cds_temp <- reduce_dimension(cds_temp, preprocess_method = "PCA",
                             reduction_method = "UMAP",
                             verbose = T)

row.names(colData(cds_temp)) <- paste0(colData(cds_temp)$cell_drug,"_",colData(cds_temp)$`_indices`)

plot_cells(cds_temp,
           reduction_method = "UMAP",
           color_cells_by = "Organ",
          label_cell_groups = F) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave("tahoe_pseudobulk_MrVI_UMAP_u-space.png",
       dpi = 600,
       height = 4, width = 6)


# ==============================================================================
# 
# ==============================================================================

col_data <- read_csv("tahoe_pseudobulk_colData.csv")
cell_line_meta <- read_csv("cell_line_metadata.csv") %>%
  select(cell_name, Organ, Cell_ID_Cellosaur) %>%
  distinct()
trial_success_meta <- read_csv("tahoe_drugs_clinical_trials.csv")

col_data <- left_join(col_data, cell_line_meta, by = c("Cell_ID_Cellosaur"))

z_space <- read_csv("tahoe_pseudobulk_MrVI_z_factor.csv",
                    col_names = F)

gene_meta <- data.frame("id" = colnames(z_space),
                        "gene_short_name" = colnames(z_space))
row.names(gene_meta) <- gene_meta$id

cds_temp <- new_cell_data_set(expression_data = t(as.matrix(z_space)),
                              cell_metadata = col_data,
                              gene_metadata = gene_meta)

row.names(z_space) <- col_data %>% select(cell_drug) %>% pull()

reducedDims(cds_temp)[["PCA"]] <- z_space

cds_temp <- reduce_dimension(cds_temp, preprocess_method = "PCA",
                             reduction_method = "UMAP",
                             verbose = T)

row.names(colData(cds_temp)) <- paste0(colData(cds_temp)$cell_drug,"_",colData(cds_temp)$`_indices`)

plot_cells(cds_temp,
           reduction_method = "UMAP",
           color_cells_by = "Organ",
           label_cell_groups = F) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave("tahoe_pseudobulk_MrVI_UMAP_z-space.png",
       dpi = 600,
       height = 4, width = 6)


colData(cds_temp)$drug_of_interest <- case_when(
  grepl(pattern = "stat", colData(cds_temp)$drug) ~ colData(cds_temp)$drug,
  TRUE ~ NA
)


plot_cells(cds_temp,
           reduction_method = "UMAP",
           color_cells_by = "drug_of_interest",
           label_cell_groups = F) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave("tahoe_pseudobulk_MrVI_UMAP_z-space_HDAC.png",
       dpi = 600,
       height = 4, width = 6)

for (gene in 1:30) {
  plot_cells(cds_temp,
             reduction_method = "UMAP",
             genes = paste0("X",gene), 
             label_cell_groups = F, 
             cell_size = 0.5
            ) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 8),
          legend.key.height = unit(0.25, "cm"),
          legend.key.width = unit(0.25, "cm"))
  ggsave(paste0("/Users/rossgiglio/Desktop/tahoe_100M/z_space/tahoe_pseudobulk_MrVI_UMAP_z-space_X",gene,".png"),
         dpi = 600,
         height = 2, width = 2.75)
}


z_space_organ <- z_space_for_print %>%
  separate(cell_drug, into = c("cell_name", "drug"), sep = "_") %>%
  left_join(cell_line_meta, by = c("cell_name"))

z_space_organ_mean <- z_space_organ %>%
  group_by(Organ) %>%
  summarise(mean_X4 = mean(X4), sd = sd(X4))

ggplot(z_space_organ_mean %>%
         mutate(Organ = fct_reorder(Organ, -mean_X4)), 
       aes(x = Organ, 
           y = mean_X4)) +
  geom_bar(stat = "identity",
           aes(fill = Organ), 
           show.legend = F) +
  scale_y_continuous(breaks = seq(-1,2,0.5),
                     limits = c(-0.5,2.5)) +
  ylab("Z-space X4") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title.x = element_blank()) +
  monocle3:::monocle_theme_opts()
ggsave(file = "tahoe_pseudobulk_MrVI_UMAP_z-space_X4_bar.png",
       dpi = 600,
       height = 2, width = 3)  

z_space_organ_mean <- z_space_organ %>%
  group_by(Organ) %>%
  summarise(mean_X26 = mean(X26), sd = sd(X26))
  
  
ggplot(z_space_organ_mean %>%
         mutate(Organ = fct_reorder(Organ, -mean_X26)), 
       aes(x = Organ, 
           y = mean_X26)) +
  geom_bar(stat = "identity",
           aes(fill = Organ), 
           show.legend = F) +
  # scale_y_continuous(breaks = seq(-1,2,0.5),
  #                    limits = c(-0.5,2.5)) +
  ylab("Z-space X26") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title.x = element_blank()) +
  monocle3:::monocle_theme_opts()
ggsave(file = "tahoe_pseudobulk_MrVI_UMAP_z-space_X26_bar.png",
       dpi = 600,
       height = 2, width = 3)  
  