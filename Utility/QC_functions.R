require(tidyverse)

runQC_Traces_singlechannel <- function(trace_df,QC_save_file,QC_reset = FALSE) {
  if (file.exists(QC_save_file) && QC_reset == FALSE) {
    pws <- read.csv(QC_save_file,stringsAsFactors = FALSE)
  }
  else {
    pws <- tibble(plate_well_sub = unique(trace_df$plate_well_sub),
                  keep_trace = -1) # -1 is unanalyzed flag
  }
  
  i <- 1
  unread_locs <- which(pws$keep_trace==-1)
  
  if (length(unread_locs) != 0) {
    i <- unread_locs[1] #go to first unread (allows for saving)
    
    
    while (i <= nrow(pws)) {
      traces <- trace_df %>% filter(plate_well_sub == pws$plate_well_sub[i])
      
      if (sum(is.na(traces$signal)) == 0) { #make sure you can plot the traces, if not throw out
        plot(traces$time,traces$signal,
             col = traces$peak_ID,pch = 16,
             main = paste("Trace",i))
        lines(traces$time,traces$signal,col = "black")
        
        
        print(paste0("i = ",i))
        ui <- readline()
        
        if (ui == 1) {
          pws$keep_trace[i] <- 1
          
          i <- i+1
        }
        if (ui == 2) {
          pws$keep_trace[i] <- 0
          
          i <- i+1
        }
        if (ui == 3) {
          pws$keep_trace[i] <- 2 # "maybe"
          
          i <- i+1
        }
        if (ui == 'q' || ui == 'Q') {
          break
        }
        if (ui == 'b' || ui == 'B') {
          i <- i-1
          while (pws$keep_trace[i] == -2) { #skip errors
            i <- i-1
          }
        }
      }
      else {
        print(paste("oopsie on ",i,pws$plate_well_sub[i]))
        pws$keep_trace[i] <- -2 #error flag
        
        i <- i+1
      }
    }
    write.csv(pws,QC_save_file,row.names = FALSE)
  }
  
  return(pws)
}



#run for two traces at once
runQC_Traces <- function(trace_df,QC_save_file,QC_reset = FALSE) {
  if (file.exists(QC_save_file) && QC_reset == FALSE) {
    pws <- read.csv(QC_save_file,stringsAsFactors = FALSE)
  }
  else {
    pws <- tibble(plate_well_sub = unique(trace_df$plate_well_sub),
                  keep_trace = -1) # -1 is unanalyzed flag
  }
  
  i <- 1
  unread_locs <- which(pws$keep_trace==-1)
  
  if (length(unread_locs) != 0) {
    i <- unread_locs[1] #go to first unread (allows for saving)
    
    
    while (i <= nrow(pws)) {
      traces <- trace_df %>% filter(plate_well_sub == pws$plate_well_sub[i])
      
      if (sum(is.na(traces$signal.calcium)) + sum(is.na(traces$signal.contractility)) == 0) { #make sure you can plot the traces, if not throw out
        par(mfrow= c(2,1))
        plot(traces$time,traces$signal.calcium,
             col = traces$peak_ID.calcium,pch = 16,
             main = paste("Calcium",i))
        plot(traces$time,traces$signal.contractility,
             col = traces$peak_ID.contractility,pch = 16,
             main = paste("Contractility",i))
        
        
        
        print(paste0("i = ",i))
        ui <- readline()
        
        if (ui == 1) {
          pws$keep_trace[i] <- 1
          
          i <- i+1
        }
        if (ui == 2) {
          pws$keep_trace[i] <- 0
          
          i <- i+1
        }
        if (ui == 3) {
          pws$keep_trace[i] <- 2 # "maybe"
          
          i <- i+1
        }
        if (ui == 'q' || ui == 'Q') {
          break
        }
        if (ui == 'b' || ui == 'B') {
          i <- i-1
          while (pws$keep_trace[i] == -2) { #skip errors
            i <- i-1
          }
        }
      }
      else {
        print(paste("oopsie on ",i,pws$plate_well_sub[i]))
        pws$keep_trace[i] <- -2 #error flag
        
        i <- i+1
      }
    }
    write.csv(pws,QC_save_file,row.names = FALSE)
  }
  
  return(pws)
}

