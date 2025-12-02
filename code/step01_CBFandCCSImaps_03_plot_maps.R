library(ggseg)
library(ggsegGlasser)
library(ggplot2)
library(dplyr)

# define plot funtion
plot_bil_map <- function(df_plot, limb, limt, savename){
  # use ggplot + geom_brain 
  pdf(file = savename, width = 10, height = 5)
  
  p <- ggplot(df_plot, aes(fill = score)) +
    geom_brain(
      atlas = glasser,
      colour = "grey",
      hemi = NULL,        # hemis：'left', 'right'，default both
      side = NULL,        # view：'lateral', 'medial'，default both
      size = 0.05,
      show.legend = TRUE
    ) +
   scale_fill_gradient2(
      low = "#2166ac",      # dark blue
      mid = "white",        
      high = "#b2182b",     # dark red
      midpoint = (limb+limt)/2,
      limits = c(limb, limt),
      breaks = c(limb, limt),
      labels = c(limb, limt),
      oob = scales::squish,
      na.value = "#9d9d9d"
    ) +
    guides(fill = guide_colorbar(
      title = NULL,                 # colorbar title
      direction = "horizontal",     # horizontal colorbar
      title.position = "top",       # title above colorbar
      label.position = "bottom",    # label below colorbar
      barwidth = unit(2, "cm"),     # colorbar width
      barheight = unit(0.4, "cm"),  # colorbar height
      ticks = FALSE                 # tick
    )) +
    theme_void() +
    theme(
      legend.position = "bottom",         
      legend.justification = "center",   
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  print(p)
  dev.off()
}


plot_uni_map <- function(df_plot, limb, limt, savename){
  pdf(file = savename, width = 8, height = 5)
  
  p <- ggplot(df_plot, aes(fill = score)) +
    geom_brain(
      atlas = glasser,
      colour = "grey",
      hemi = 'left',       
      side = NULL,       
      size = 0.05,
      show.legend = TRUE
    ) +
    scale_fill_gradient2(
      low = "#2166ac",     
      mid = "white",        
      high = "#b2182b",    
      midpoint = (limb+limt)/2,
      limits = c(limb, limt),
      breaks = c(limb, limt),
      labels = c(limb, limt),
      oob = scales::squish,
      na.value = "grey"
    ) +
    guides(fill = guide_colorbar(
      title = NULL,                
      direction = "horizontal",     
      title.position = "top",      
      label.position = "bottom",   
      barwidth = unit(2, "cm"),    
      barheight = unit(0.4, "cm"), 
      ticks = FALSE               
    )) +
    theme_void() +
    theme(
      legend.position = "bottom",        
      legend.justification = "center",   
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  print(p)
  dev.off()
}



# set input and output path
InPath = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data/group/'
OutPath = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/Draft/figure/maps/'


# plot bil PScs results
filename = 'CBF_Mean_PScs_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -1, 1, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_Mean_PScs_qc_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -1, 1, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_Mean_PScs2_qc_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -1, 1, paste(OutPath, filename, '.pdf', sep = ""))




# plot lobe
filename = 'lobe'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, 1, 5, paste(OutPath, filename, '.pdf', sep = ""))




# plot bil CBF value PC1 results
filename = 'CBF_PC1_Val_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_PC2_Val_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_PC3_Val_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))




# plot bil CBF value PC1 results ---layer
filename = 'CBF_PC1_Val_L6_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_PC1_Val_L5_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_PC1_Val_L4_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_PC1_Val_L3_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_PC1_Val_L2_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CBF_PC1_Val_L1_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))



# plot bil mito results
filename = 'Mito_CI_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'Mito_CII_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'Mito_CIV_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'Mito_MitoD_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'Mito_MRC_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'Mito_TRC_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))



# plot Bigbrain MPC
filename = 'BB_G1_bil'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -0.1, 0.1, paste(OutPath, filename, '.pdf', sep = ""))
plot_bil_map(df_plot, -0.000001, 0.000001, paste(OutPath, 'BB_G1_bil_2part.pdf', sep = ""))




# plot function gradient
filename = 'FunctionGradient1'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -7, 7, paste(OutPath, filename, '.pdf', sep = ""))





# plot SNR and CNR-G1
filename = 'SNRmean'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, 1, 7, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'CNR_G1'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -5, 10, paste(OutPath, filename, '.pdf', sep = ""))






# plot cell maps
filename = 'Endo_uni'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'Per_uni'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'In6a_uni'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'In6b_uni'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))

filename = 'Oligo_uni'
df_plot = read.csv(paste(InPath, filename, '.csv', sep = ""))
df_plot$score[df_plot$score == 0] <- NA 
plot_bil_map(df_plot, -2, 2, paste(OutPath, filename, '.pdf', sep = ""))



