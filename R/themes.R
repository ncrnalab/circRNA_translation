MS_colors <- c(
   "signal" = "blue",
   "noise" = "grey",
   "Uniprot" = "#2166AC",
   "CDS" = "#2166AC",
   "circRNA" = "#B2182B",
   "COVID19" = "#9e46ab",
   "sno/miRNA" = "#9e46ab"
)


library(RColorBrewer)
colfunc <- colorRampPalette(c("dark blue", "orange"))
frame_colors <- colfunc(3)


riboseq_colors <- c(
   "frame 1" = frame_colors[1],
   "frame 2" = frame_colors[2],
   "frame 3" = frame_colors[3]
)



labels <- c("uc001qop.2" = "\u03B2-ACTIN",
           "uc003sot.4" = "GAPDH",
           "CDS" = "\u03B2-ACTIN",
           "circRNA" = "circRNA",
           "highq" = "High quality",
           "lowq" = "Low quality")




ggtheme <-
  theme_bw () + 
  theme (
    panel.background = element_rect(fill = "white", color="grey", size=1),
    text = element_text (size = 22, color="black"),
    #axis.title = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    strip.text = element_text(size = 20, color = "black"),
    #legend.title = element_text(size = 17, face = NULL, color = "black"),
    legend.text = element_text(size = 20, face = NULL, color = "black"),
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    strip.background = element_blank()
  ) 
