library(ggcorrplot)
library(patchwork)
library(tidyverse)
library(viridis)

surface <- readxl::read_xlsx('GSM_Cortable.xlsx')
patch <- readxl::read_xlsx('PMM_Cor.Table.xlsx')

surface_matrix <- surface %>%
  select(-X__1) %>%
  as.matrix()

patch_matrix <- patch %>%
  select(-X__1) %>%
  as.matrix()

rownames(surface_matrix) <- colnames(surface_matrix)
rownames(patch_matrix) <- colnames(patch_matrix)

col <- viridis::viridis(n = 3)
col <- c('gray0', 'white', 'gray1')

gg_surface <- ggcorrplot(surface_matrix, title = 'Surface metrics', 
                         show.legend = FALSE, 
                         lab = TRUE, type = 'lower',
                         colors = col, lab_size = 3, lab_col = 'white')

gg_patch <- ggcorrplot(patch_matrix, title = 'Patch metrics',
                       lab = TRUE, type = 'lower', 
                       colors = col, lab_size = 3, lab_col = 'white')

gg_overall <- gg_surface + gg_patch

ggsave('Correlation_3.jpeg', plot = gg_overall, width = 25, units='cm')
