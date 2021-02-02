library (RCurl)
library (data.table)  
library (tidyverse)
library (usethis)


#source ("~/projects/SFPQ/DESeq2_common.R")
source ("R/riboseq_utils.R")
source ("R/themes.R")

riboseq_data <- data.frame (fread ("data/RiboSeq_featureCounts.txt"))



# add circRNA data
circ_data <- data.frame (fread ('zcat data/RiboSeq_circ.txt.gz'))

circ_data <-
   circ_data %>% filter (real_rel_pos > -8 & real_rel_pos < 6) %>% # only bsj-spanning reads
      group_by (srr, "read_length" = rc) %>%
      summarize (val=sum(reads), feature = "circRNA", n=n()) 

# set zeros for all circ/
srr <- unique (riboseq_data$srr)  
rl <- unique (riboseq_data$read_length)
   
circ_grid <- expand.grid (srr=srr, read_length=rl)
 
circ_grid <- merge (circ_grid, circ_data, by=c("srr", "read_length"), all.x=T)
circ_grid$val[is.na (circ_grid$val)] <- 0
circ_grid$feature <- "circRNA"


bind_data <- rbind (
   as.data.frame (riboseq_data[,colnames (riboseq_data)]),
   as.data.frame (circ_grid[,colnames (riboseq_data)]))

# expression

bind_data <- bind_data %>% group_by (srr) %>% mutate (myrpm = 1e6 * val / sum(val)) 


# add quality
#qual_data <- data.frame (fread(getURL("sftp://login.genome.au.dk/faststorage/home/tbhmbi/HOME/annotations/circRNA/hsa_ref_riboseq_quality.txt", userpwd=rawToChar(as.raw (userpwd_genomedk)), connecttimeout=120), fill=T))
qual_data <- data.frame (fread ("data/RiboSeq_qual.txt"))
qual_data$read_length <- as.integer (gsub ("size", "", qual_data$rc))


bind_data <- 
   bind_data %>% 
   left_join (qual_data, by=c("srr", "read_length")) 


# convert to rpm
#depth_data <- data.frame (fread(getURL("sftp://login.genome.au.dk/faststorage/home/tbhmbi/RiboSeq/hsa/read_counts2.txt", userpwd=rawToChar(as.raw (userpwd_genomedk)), connecttimeout=120), fill=T))

depth_data <- data.frame (fread("data/RiboSeq_mapstats.txt"))

colnames (depth_data) <- c("srr", "read_length", "count")
depth_data$srr <- get_srr (depth_data$srr)


r2c <- structure (depth_data$count, names=paste (depth_data$srr, depth_data$read_length, sep="_"))

bind_data$rpm <- bind_data$val * 1e6 / r2c[paste (bind_data$srr, bind_data$read_length, sep="_")]

gg <- list()

use_features <- c("CDS", "circRNA")


gg[["counts"]] <-
ggplot (bind_data %>% filter (feature %in% use_features) %>% group_by (read_length, feature) %>% summarize (n=sum(val)), aes (x=read_length, y=n, fill=feature)) +
   geom_col () + 
   scale_y_log10(labels=dropZero) +
   scale_x_continuous(breaks = c(25, 30, 35)) +
   scale_fill_manual(values = MS_colors) + 
   facet_wrap (~feature) + #, scales = "free_y") + 
   labs (x="Read length (nts)", y="#Reads", color="") + 
   ggtheme + 
   theme (legend.position = "none")

use_features <- c("CDS", "circRNA", "sno/miRNA")
   
gg[["counts_w_snomi"]] <-
ggplot (bind_data %>% filter (feature %in% use_features) %>% group_by (read_length, feature) %>% summarize (n=sum(val)), aes (x=read_length, y=n, fill=feature)) +
   geom_col () + 
   scale_y_log10(labels=dropZero) +
   scale_x_continuous(breaks = c(25, 30, 35)) +
   scale_fill_manual(values = MS_colors) + 
   facet_wrap (~feature) + #, scales = "free_y") + 
   labs (x="Read length (nts)", y="#Reads", color="") + 
   ggtheme + 
   theme (legend.position = "none")


use_features <- c("CDS", "circRNA", "sno/miRNA", "lncRNA")
   
gg[["counts_w_snomi_and_lncrna"]] <-
ggplot (bind_data %>% filter (feature %in% use_features) %>% group_by (read_length, feature) %>% summarize (n=sum(val)), aes (x=read_length, y=n, fill=feature)) +
   geom_col () + 
   scale_y_log10(labels=dropZero) +
   scale_x_continuous(breaks = c(25, 30, 35)) +
   scale_fill_manual(values = MS_colors) + 
   facet_wrap (~feature) + #, scales = "free_y") + 
   labs (x="Read length (nts)", y="#Reads", color="") + 
   ggtheme + 
   theme (legend.position = "none")


bind_data <- bind_data %>%
   group_by (read_length) %>%
   mutate (qq = findInterval (-log(p), quantile(-log(p), probs=0:4/5)))
   

use_features <- c("CDS", "circRNA", "sno/miRNA")

gg[["RPMMvsQuality"]] <-
   ggplot (bind_data %>% filter (myrpm > 0, feature %in% use_features), aes (x=qq, y=myrpm, color=feature, fill=feature)) +
      geom_smooth(method="lm", show.legend = FALSE) +
      scale_y_log10(labels = dropZero) +
      facet_wrap (~ read_length, nrow=1) + 
      scale_color_manual(values=MS_colors) +
      scale_fill_manual(values=MS_colors) +
      scale_x_continuous (labels = c("L", "", "", "", "H")) +
      labs (x="Quality", y="RPMM", color="") + 
      #scale_color_discrete(discrete = TRUE)+
      #scale_fill_discrete(discrete = TRUE)  +
      ggtheme

use_features <- c("CDS", "circRNA", "sno/miRNA", "lncRNA")

gg[["RPMMvsQuality_w_lncRNA"]] <-
   ggplot (bind_data %>% filter (myrpm > 0, feature %in% use_features), aes (x=qq, y=myrpm, color=feature, fill=feature)) +
      geom_smooth(method="lm", show.legend = FALSE) +
      scale_y_log10(labels = dropZero) +
      facet_wrap (~ read_length, nrow=1) + 
      scale_color_manual(values=MS_colors) +
      scale_fill_manual(values=MS_colors) +
      scale_x_continuous (labels = c("L", "", "", "", "H")) +
      labs (x="Quality", y="RPMM", color="") + 
      #scale_color_discrete(discrete = TRUE)+
      #scale_fill_discrete(discrete = TRUE)  +
      ggtheme


bind_data$feature <- factor (bind_data$feature, levels = c("CDS", "lncRNA", "sno/miRNA", "circRNA"))


gg[["RPMMvsQuality_legend"]] <-
   ggplot (bind_data %>% filter (myrpm > 0, feature %in% use_features), aes (x=qq, y=myrpm, color=feature)) +
      geom_line (size=1.2) +
      scale_y_log10(labels = dropZero) +
      scale_color_manual(values=MS_colors) +
      facet_wrap (~ read_length, nrow=1) +
      scale_x_continuous (labels = c("L", "", "", "", "H")) +
      labs (x="Quality", y="RPMM", color="") +
      ggtheme



cairo_pdf("figures/riboseq_correlations_W7.pdf", onefile=T, width=7)
   print (gg)
dev.off()


cairo_pdf("figures/riboseq_correlations_W10.pdf", onefile=T, width=10)
   print (gg)
dev.off()


