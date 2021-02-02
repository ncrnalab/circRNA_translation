library (tidyverse)
library (data.table)
library (RCurl)

source ("R/themes.R")


MS_data <- list(
   "circRNA" = list (file="data/MS_circRNA.txt.gz", circ=T, grep="HUBNER"),
   "SARS-CoV-2" = list (file="data/MS_COVID19.txt.gz", circ=F, grep="SARS2"),
   "Uniprot" = list (file="data/MS_Uniprot_subsample.txt.gz", circ=F, grep="HUMAN")
   
)
   

MS_data <- map2 (MS_data, names (MS_data), function (msn, n) {

   dfMS <- data.frame (fread (paste ("zcat", msn$file)))
   
   print (nrow (dfMS))

   if (!"Matches" %in% colnames (dfMS) && "Unique..Proteins." %in% colnames (dfMS)) {
      dfMS$Matches = ifelse (dfMS$Unique..Proteins. == "yes", 1, 2)
      
   }

   if (!"Leading.proteins" %in% colnames (dfMS) && "Leading.razor.protein"  %in% colnames (dfMS)) {
      dfMS$Leading.proteins <- dfMS$Leading.razor.protein
   }
   
   
   
   dfMS_u <- 
      dfMS %>% filter (Unique..Proteins. == "yes", Matches <= 1)
   
   if (msn$grep != "") dfMS_u <- dfMS_u %>% filter (grepl (msn$grep, Leading.proteins))

   if (msn$circ) {

      # check if peptide is across BSJ
      lp <- dfMS_u$Leading.proteins
      lp <- gsub ("REV__", "", lp)

      postfix <- sapply (strsplit (lp, "@"), function (x) x[2])
      dfMS_u$bsj <- as.numeric (sapply (strsplit (postfix, "_"), function (x) x[1]))
      dfMS_u$circ <- ifelse (dfMS_u$Reverse == "+",
                             dfMS_u$Start.position <= dfMS_u$Protein.Length - dfMS_u$bsj,
                             dfMS_u$End.position >= dfMS_u$bsj)

   } else {
      dfMS_u$circ <- FALSE
      dfMS_u$bsj <- -1

   }

   dfMS_u$decoy <- dfMS_u$Reverse == "+"
   dfMS_u$decoy_text <- ifelse (dfMS_u$Reverse == "+", "noise", "signal")
   dfMS_u$decoy_text2 <- ifelse (dfMS_u$Reverse == "+", "noise", n)
   
   dfMS_signal <- dfMS_u %>% filter (!decoy)

   set.seed (1234)
   dfMS_u$random <- sample(1:nrow(dfMS_u), nrow(dfMS_u))

   dfMS_u$name = n
   
   dfMS_u$decoy_count <- ifelse (dfMS_u$decoy, dfMS_u$MS.MS.count, 0)
   dfMS_u$signal_count <- ifelse (!dfMS_u$decoy, dfMS_u$MS.MS.count, 0)
   
   # msn$dfFDR <-
   #    dfMS_u %>% 
   #     arrange (PEP, random) %>%
   #     mutate (TPR=cumsum(signal_count)/sum(signal_count), FPR=cumsum(decoy_count)/sum(decoy_count), FDR=2*cumsum(decoy_count)/(cumsum(2*decoy_count+signal_count)), hits = cumsum(signal_count), decoys = cumsum(decoy_count)) %>%
   #     arrange (-PEP, random) %>%
   #     mutate (Q=cummin (FDR), rank=rank(PEP))
   # 
    msn$dfFDR <-
      dfMS_u %>% 
       arrange (PEP, random) %>%
       mutate (#TPR=cumsum(signal_count)/sum(signal_count), 
               #FPR=cumsum(decoy_count)/sum(decoy_count), 
               FDR=cumsum(decoy_count)/(cumsum(decoy_count+signal_count)), 
               targets = cumsum(signal_count), 
               decoys  = cumsum(decoy_count)) %>%
       arrange (-PEP, random) %>%
       mutate (Q=cummin (FDR), rank=rank(PEP))

  
    df <- msn$dfFDR 
    
    print (nrow (dfMS_u))
    
    msn   
   
   
})



# combine data
dfFDR.ALL <- map2 (MS_data, names (MS_data), function (msn, n) {
  
   cols <- c("name", "PEP", "Q", "MS.MS.count", "decoy", "decoy_text", "decoy_text2", "circ")
   msn$dfFDR[,cols]
   
}) %>% bind_rows ()
   

dfFDR.ALL %>%
   group_by (name) %>%
   summarize (max (Q))



gg <- list()

ggtheme <- ggtheme + theme (axis.text.x = element_text (size=14))
ggtheme <- ggtheme + theme (axis.text.y = element_text (size=14))

dfFDR.ALL$name <- factor (dfFDR.ALL$name, levels = c("Uniprot", "circRNA", "SARS-CoV-2"))
dfFDR.ALL$decoy_text2 <- factor (dfFDR.ALL$decoy_text2, levels = c("Uniprot", "circRNA", "SARS-CoV-2", "noise"))


gg[['density']] <-
ggplot (dfFDR.ALL, aes (x=PEP, color=decoy_text2)) +
   geom_density (aes(y=..count..), size=1.2) +
   scale_color_manual (values = MS_colors) +
   scale_y_log10() + 
   scale_x_continuous(breaks = c(0.01, 0.05, 0.09)) +
   facet_wrap (~name, nrow=1) + 
   labs (color="", x="PEP", y="#Peptides") +
   ggtheme
   

dfFDR.sum <- dfFDR.ALL %>% group_by (decoy_text2, decoy, name) %>%
   summarize (n=sum (MS.MS.count), max_fdr = max (Q)) 



gg[['npeptides+decoys']] <-
   ggplot (dfFDR.sum) + 
   scale_y_log10() +
   geom_col (aes (x=decoy, y=n, fill=decoy_text2)) +
   geom_text (aes (x=decoy, y=n,label = n), vjust=0, size=5) +
   scale_fill_manual (values = MS_colors) +
   #scale_y_continuous(dropZero) + 
   #facet_wrap (~name, nrow=1, scales="free_y") +
   facet_wrap (~name, nrow=1) +
   labs (x="", y="#Peptides", fill="") +
   ggtheme + theme (axis.text.x = element_blank())

 


gg[['FDRbyPEP']] <-
   ggplot (dfFDR.ALL %>% filter (!decoy), aes (x=PEP, y=Q, color=name)) + 
     geom_line(size = 1.3) +  
     scale_color_manual(values = MS_colors) +
     labs (x="Ranked by score", y="FDR", color="") +
     facet_wrap (~name, nrow=1) +
     ggtheme # + # ylim (c(0,1)) + 




gg[['FDRbyPEP_facet']] <-
   ggplot (dfFDR.ALL %>% filter (!decoy), aes (x=PEP, y=Q, color=name)) + 
     geom_point(size = 1.3) +
     geom_line(size = 1.1) +  
     scale_color_manual(values = MS_colors) +
     labs (x="PEP", y="FDR", color="") +
     facet_wrap (~name) + #, scales = "free_x") +
     scale_x_continuous(breaks=c(0.01, 0.05, 0.09)) + 
     ggtheme # + # ylim (c(0,1)) + 



gg[['FDRbyPEP_wdecoy']] <-
ggplot (dfFDR.ALL, aes (x=PEP, y=Q, color=decoy_text2)) + 
  geom_point(size = 1.3) +  
  labs (x="PEP", y="FDR", color="") +
  scale_color_manual (values = MS_colors) +
  ggtheme + ylim (c(0,1)) + 
  scale_x_continuous(breaks=c(0.01, 0.05, 0.09)) + 
  facet_wrap (~name, nrow=1) + 
  theme (legend.position = c(.05, .95),
         legend.justification = c(0,1))
 


cairo_pdf("figures/MS.pdf", onefile=T, width=9)
   print (gg)
dev.off()


