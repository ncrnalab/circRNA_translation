# Load libs
library (tidyverse)
library (Biostrings)
library (RCurl)
library (data.table)  

source ("R/colors.R") # overwrite ggtheme
source ("R/riboseq_utils.R") 


dfAll <- data.frame (fread ("zcat data/RiboSeq_ref.txt.gz"))
dfQual <- data.frame (fread ("data/RiboSeq_qual.txt"))


fa = readBStringSet("data/RiboSeq_ref.fa")
fa_names <- names (fa)
fa_names <- sapply (fa_names, function (x) unlist (strsplit (x, " "))[1])
names (fa) <- fa_names

# only 29nt
read_lengths <- list(
   "29" = list(min=29, max=29)
)


dfAll.fdr <- dfAll %>%
     group_by (srr, rc) %>%
     summarize (mfdr = mean (fdr)) %>%
     group_by (rc) %>%
     summarize (lowq = quantile (mfdr, 0.95),
                highq = quantile (mfdr, 0.05)) %>%
     as.data.frame ()


highq <- structure (dfAll.fdr$highq, names=dfAll.fdr$rc)
lowq <- structure (dfAll.fdr$lowq, names=dfAll.fdr$rc)

dfAll$highq <- dfAll$fdr <= highq[as.character(dfAll$rc)]
dfAll$lowq <- dfAll$fdr >= lowq[as.character(dfAll$rc)]

refseq <- lapply (names (fa), function (x) list(ref = x)) %>% as.list()

x <- refseq[[1]]

refseq <- map (refseq, function (x) {
   
   x$seq <- toString (fa[[x$ref]])
   
   x$orf <- FindOrf (x$seq) %>% mutate (length = end-start) %>% top_n (length, n=1) %>% as.list()
   
   x$alpha_frame <- x$orf$frame

   n <- names (read_lengths)[1]
   rl <- read_lengths[[n]]
   
   x$rl <- map2 (read_lengths, names (read_lengths), function (rl, n) {
     
      fdr_quality <- "highq"
      rl$quality <- map (c("highq","lowq"), function (fdr_quality) {

         lquality <- list (quality = fdr_quality)
         
         m.sub.append <- dfAll %>% 
            filter (header == x$ref) %>%
            filter (rc >= rl$min & rc <= rl$max) %>%
            as.data.frame () #l.gene.append[[header]]
       
         m.sub.append <- m.sub.append[m.sub.append[,fdr_quality],]
          
         if (nrow(m.sub.append) == 0) {
             return (NULL);
         }
           
         xstart <- x$orf$start
         xend <- x$orf$end
         frame_number <- x$orf$frame
         frame <- paste ("frame", frame_number)
           
         if (xend > nchar(x$seq)) {
             return (NULL);
         }
           
         dfFrame <- data.frame (header=x$ref, 
                                xstart=xstart, xend=xend, 
                                frame_number=frame_number)
                                
          
          m.sub.cds <- m.sub.append %>% filter (real_pos >= xstart, real_pos <= xend)
          lquality$reads <- data.frame (type = "CDS", name = x$ref, qual = fdr_quality, reads = sum (m.sub.cds$reads))
            
          # print
          
          farray <- c(1,2,3,1,2,3)
           
          m.sub.append$xheader <- x$ref
          m.sub.append$qual <- fdr_quality
          m.sub.append$frame_number <- ((m.sub.append$real_pos-1) %% 3) + 1
          m.sub.append$frame <- paste ("frame", farray[m.sub.append$frame_number-(x$alpha_frame-1)+3])
      
          
          
          gg <- ggplot (m.sub.append, aes (x=real_pos, y=reads, fill=frame, color=frame)) + 
             ggtheme
      
          gg <- gg + geom_bar (stat="identity") 
           
          yrange <- ggplot_build(gg)$layout$panel_scales_y[[1]]$range$range
          #yrange <- ggplot_build(gg)$layout$panel_ranges[[1]]$y.range
           
          dfFrame$y <- -(yrange[2]) * 0.015 * (dfFrame$frame_number+1) + (yrange[2]) * 0.009
          dfFrame$frame <- paste ("frame", farray[dfFrame$frame_number-(x$alpha_frame-1)+3])
           
          gg <- gg + geom_segment (data=dfFrame, aes(x = xstart, y = y, xend = xend, yend = y, color=frame), size=2.1) + 
             scale_color_manual(values=riboseq_colors, guide=F) +
             scale_fill_manual(values=riboseq_colors) + 
             facet_wrap (~qual, labeller = as_labeller(labels)) +
             labs (x="P-site position", y="#Reads", fill="")
           
          lquality$gg[[paste (x$ref, rl$min, rl$max, fdr_quality)]] <- gg
          lquality
      
      })
         
      rl
   })
  
   x
})

# get ggplots
print_gg <- map (refseq, function (x) {
  map (x$rl, function (rl) {
     map (rl$quality, function (lquality) {
        lquality$gg
     })
  })
})
    

#get mapped reads
dfReads.cds <- map (refseq, function (x) {
  map (x$rl, function (rl) {
     map (rl$quality, function (lquality) {
        lquality$reads
     }) %>% bind_rows ()
  }) %>% bind_rows ()
}) %>% bind_rows ()



dfAll <- fread ('zcat data/RiboSeq_circ.txt.gz')

dfAll$highq <- dfAll$fdr <= highq[as.character(dfAll$rc)]
dfAll$lowq <- dfAll$fdr >= lowq[as.character(dfAll$rc)]
 

rl <- read_lengths[[1]]

circ <- map (read_lengths, function (rl) {
  
   fdr_quality <- "highq"
   rl$quality <- map (c("highq","lowq"), function (fdr_quality) {
      
       lquality <- list (quality = fdr_quality)

       m.sub.append <- dfAll %>% 
                           filter (rc >= rl$min & rc <= rl$max) %>%
                           filter (real_rel_pos >= -8, real_rel_pos <=6) %>%
                           as.data.frame ()
      
       m.sub.append <- m.sub.append[m.sub.append[,fdr_quality],]
          
       lquality$reads <- data.frame (type = "circRNA", name = "top1000", qual = fdr_quality, reads = sum (m.sub.append$reads))

       lquality
   })
   
   rl
})


dfReads.circ <- 
   map (circ, function (rl) {
     map (rl$quality, function (lquality) {
        lquality$reads
     }) %>% bind_rows ()
   }) %>% bind_rows ()



df.bind <- rbind (dfReads.cds, dfReads.circ) %>% filter (name != "uc003sot.4")
df.bind$qual <- factor (df.bind$qual, levels= c("lowq", "highq"))

print_gg[["reads"]] <-
   ggplot (df.bind, aes (x=qual, y=reads, fill=type)) + 
   geom_col () +
   scale_x_discrete(labels = c("highq" = "High", "lowq" = "Low")) + 
   scale_fill_manual (values = MS_colors) +
   labs (x="Quality", y="#Reads") +
   facet_wrap (~type, scales = "free_y", labeller = as_labeller (labels)) +
   ggtheme + 
   theme (legend.position = "none")
   

cairo_pdf("figures/riboseq_reads.pdf", onefile=T)
   print (print_gg)
dev.off()


   
