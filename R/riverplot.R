plot_river <- function(anno.df1, anno.df2, common.cells,...)
  {
    source("~/zizhen/My_R/map_river_plot.R")
    source("~/zizhen/My_R/sankey_functions.R")
    cols = c("sample_id","cluster_id","cluster_label","cluster_color")
    df1 = anno.df1 %>% filter(sample_id %in% common.cells) %>% select(cols)
    df2 = anno.df2 %>% filter(sample_id %in% common.cells) %>% select(cols)
    colnames(df2)[-1] = paste0("map_",     colnames(df2)[-1])
    df = left_join(df1, df2)
    g= river_plot(df, ...)
  }
