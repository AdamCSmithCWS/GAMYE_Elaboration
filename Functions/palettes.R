# palettes



breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
trend_map_labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
trend_map_labls = paste0(trend_map_labls, " %")

trend_plot_cats <- function(x, labls = trend_map_labls,
                            t_breaks = c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)){
  Tplot <- cut(x,breaks = c(-Inf, t_breaks, Inf),labels = labls)
  return(Tplot)
}


map_palette_v <- c("#fde725", "#dce319", "#b8de29", "#95d840", "#73d055", "#55c667",
                   "#238a8d", "#2d708e", "#39568c", "#453781", "#481567")

map_palette_s <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")


names(map_palette_v) <- trend_map_labls
names(map_palette_s) <- trend_map_labls


