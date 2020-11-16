


#write the function
dropZero <- function(l){
  lnew <- c()
  for(i in l) {
     
     if (!is.na (i) && as.numeric (i) == 0) {
       lnew <- c(lnew, "0")
     } else {
       sci <- formatC(as.numeric(i), format = "e", digits = 0)
       sci0 <- gsub("e\\+0", "e", sci, perl = TRUE)
       sci0 <- gsub("e\\-0", "e-", sci0, perl = TRUE)

       if (sci0 == "1e0") {
          sci0 <- "1"
       }
       
       lnew <- c(lnew, sci0)
    }
  }
  as.character(lnew)
}



get_srr <- function (v) {
  
  as.character (sapply (v, function (x) {
    
    ss <- strsplit (x, "\\.|_|/")[[1]]
    
    i <- grep ("SRR", ss)
    
    r <- ifelse (length (i) == 1,  ss[i], x)
    
    r
    
  }))
  
}



extract_from_append <- function (data, start, end) {
  
  if (!("real_pos" %in% colnames (data)) || !("reads" %in% colnames (data))) {
    
    print ("ERROR IN DATA (extract_from_append)")
    
  }
  
  dfPos <- subset (data, data$real_pos >= start & data$real_pos < end)
  
  dfA <- data.frame(real_pos=start:(end-1), reads=0)
  
  if (nrow (dfPos) > 0) {
    
    dfA <- aggregate (dfPos$reads, by=list(dfPos$real_pos), sum)
    colnames (dfA) <- c("real_pos", "reads")
    
  }
  
  dfReturn <- merge (data.frame (real_pos=start:(end-1)), dfA, by="real_pos", all.x=T)
  dfReturn[is.na(dfReturn)] <- 0
  
  return (dfReturn)
  
}

FindOrf <- function (regseq) {
  
  #orf <- gregexpr ("ATG(?:[ATGC]{3}){5,}(?:TAG|TAA|TGA)", regseq, ignore.case = T, perl=T)
  
  gr <- gregexpr ("ATG", regseq, ignore.case = T)
  dfStart <- data.frame (pos=gr[[1]], frame=(gr[[1]]-1) %% 3, type="START")
  gr <- gregexpr ("TAG|TAA|TGA", regseq, ignore.case = T)
  dfStop <- data.frame (pos=gr[[1]], frame=(gr[[1]]-1) %% 3, type="STOP")    
  
  dfStartStop <- rbind (dfStart, dfStop)
  dfStartStop <- dfStartStop[order(dfStartStop$pos),]
  
  dfOrf <- data.frame ()
  f <- 2
  for (f in 0:2) {
    
    atg <- -1
    
    dfF <- subset (dfStartStop, dfStartStop$frame==f)
    
    if (nrow(dfF) > 0) {
      
      for (i in 1:nrow(dfF)) {
        
        #print (dfF[i,])
        if (dfF[i,"type"] == "START" && atg == -1) {
          
          atg <- dfF[i,"pos"]
          #print (atg)
        }
        
        else if (dfF[i,"type"] == "STOP" && atg != -1) {
          
          orf_start <- atg
          orf_end <- dfF[i,"pos"]
          
          #print (orf_end)
          
          if (orf_end - orf_start >= 30){ # ONLY ORF > 30 nt (10 aa)
            
            dfOrf <- rbind (dfOrf,
                            data.frame (start=orf_start, end=orf_end, frame=f+1))
          }
          
          atg <- -1
          
        }
      }  
    }
  }
  
  return (dfOrf)
  
}

