require(tidyverse)
require(progress)

# From a column, parse the well, well.letter and well.number.
# MODIFIED TO NOT REQUIRE "Well" IN THE MATCHED STRING, STILL REQUIRES "__"
parse.well <- function(csvfile,traceID=alias){
  # Inputs:
  # - csvfile: the result of performing read_csv on a contractility or KinetiCal file
  # - traceID: a columnname existent in csvfile from which to extract the well information.
  #             By default in contractility analysis and KinetiCal, the column with this info is called "alias"
  # Outputs:
  # - three (3) new columns in csvfile: well, well.letter, and well.number
  
  well.string <- csvfile %>% 
    select({{traceID}}) %>% 
    unlist() %>% 
    str_extract(., "__[A-Z]_\\d{3}") %>% 
    str_remove(.,"__")
  
  csvfile$well.letter <- substr(well.string,start = 1,stop=1)
  csvfile$well.number <- substr(well.string,start = 4,stop=5)
  csvfile$well <- paste0(csvfile$well.letter,csvfile$well.number)
  csvfile  
}

#DWS added
parse.subposition <- function(csvfile,traceID=alias) {
  # Inputs:
  # - csvfile: the result of performing read_csv on a contractility or KinetiCal file
  # - traceID: a columnname existent in csvfile from which to extract the well information.
  #             By default in contractility analysis and KinetiCal, the column with this info is called "alias"
  # Outputs:
  # - one (3) new columns in csvfile: subposition,subposition.row,subposition.col
  
  
  #r_0005_c_0005
  subposition.string <- csvfile %>% 
    select({{traceID}}) %>% 
    unlist() %>% 
    str_extract(., "r_\\d{4}_c_\\d{4}")
  
  csvfile$subposition.row <- substr(subposition.string,start = 5,stop=6)
  csvfile$subposition.col <- substr(subposition.string,start = 12,stop=13)
  csvfile$subposition <- paste0("r",csvfile$subposition.row,"c",csvfile$subposition.col)
  csvfile  
}

# The contractility analysis and KinetiCal may contain data from different plates,
# or different timepoints of the same plate. 
# We will use all information in traceID EXCEPT the well information, to create unique names for a plate/timepoint
make.plate.names <- function(csvfile,traceID=alias){
  # Inputs:
  # - csvfile: the result of performing read_csv on a contractility or KinetiCal file
  # - traceID: a columnname existent in csvfile from which to extract the plate information
  # Outputs:
  # - one (1) new column in csvfile: plate.name
  
  all.but.well.string <- csvfile %>% 
    select({{traceID}}) %>% 
    unlist() %>%
    str_remove(., "Well__[A-Z]_\\d{3}.*")
  csvfile$plate.name <- all.but.well.string
  csvfile
}

make.plate.names.nochannel <- function(csvfile,traceID=alias){
  # Inputs:
  # - csvfile: the result of performing read_csv on a contractility or KinetiCal file
  # - traceID: a columnname existent in csvfile from which to extract the plate information
  # Outputs:
  # - one (1) new column in csvfile: plate.name
  
  all.but.well.string <- csvfile %>% 
    select({{traceID}}) %>% 
    unlist() %>%
    str_remove(., "_Well__[A-Z]_\\d{3}.*")
  csvfile$plate.name <- all.but.well.string
  csvfile
}

# Standard loading of contractility analysis. 
# This function will return a list with the $metrics and $traces dataframes.
# They will share the 'alias' column in common.

loadContractility <- function(contractility.filepath, loadMetrics = TRUE, loadTraces = TRUE){
  contractility.csv <- read_csv(contractility.filepath,na = "NaN",show_col_types = FALSE)
  # Initialize NULL output
  output.list <- list("metrics" = NULL,
                      "traces"  = NULL)
  if (loadMetrics){
    message("Loading metrics", appendLF = F)
    metrics <- contractility.csv %>% 
      dplyr::select(-c(time,signal,peak_ID))
    output.list$metrics <- metrics
    message("  done")
    
  }else{
    message("Loading metrics  skipped")
  }
  
  if (loadTraces){
    b <- contractility.csv %>% 
      group_by(alias) %>% 
      dplyr::select(alias,time,signal,peak_ID) %>% 
      group_split()
    
    pb <- progress_bar$new(total=length(b),
                           format = "Loading traces [:bar] :percent eta :eta",
                           width = 60)
    
    bbb <- lapply(b,function(trace.df){
      pb$tick()
      alias <- trace.df %>% 
        dplyr::select(alias) %>% 
        unlist()
      
      time <- trace.df %>% 
        dplyr::select(time) %>% 
        parse_trace()
      
      signal <- trace.df %>%
        dplyr::select(signal) %>%
        parse_trace()
      
      peak_ID <- trace.df %>%
        dplyr::select(peak_ID) %>%
        parse_trace() %>%
        as.factor()
      
      df <- data.frame(alias, time, signal, peak_ID,row.names=NULL)
      df
    })
    traces <- bind_rows(bbb)
    output.list$traces <- traces
    message("Loading traces   done")
  }else{
    message("Loading traces   skipped")
  }
  return(output.list)
}

# As of today, KinetiCal's output follows the same structure as Contractility, so this function is a wrapper
# for loadContractility, just changing the names of the variables for the user.
loadKinetiCal <- function(kinetical.filepath, ...){
  loadContractility(contractility.filepath = kinetical.filepath,...)
} 

# The contractility analysis and KinetiCal encode trace data as a string wrapped by {} and separated by ;
# This function unwraps trace data from column names given by xCol and yCol and returns them as a tidy dataframe.
# This function is more flexible than loadContractilityTraces, as it only loads 1 column of trace 
# (for instance only signal, or only time).
tidyTraceFromCSV <- function(csvfile, traceID, traceCol){
  b <- csvfile %>% 
    group_by({{traceID}}) %>% 
    dplyr::select({{traceID}},{{traceCol}}) %>% 
    group_split()
  trace.df <- b[[1]]
  bbb <- lapply(b,function(trace.df){
    ID <- trace.df %>% 
      dplyr::select({{traceID}}) %>% 
      unlist()
    
    x <- trace.df %>% 
      dplyr::select({{traceCol}}) %>% 
      parse_trace()
    
    
    df <- data.frame(ID, x, row.names=NULL) %>%
      rename({{traceID}} := ID,
             {{xCol}}    := traceCol)
    df
  })
  bind_rows(bbb)
}

# Takes trace data in string format and converts it to numeric. 
# This is a helper function for tidyTraceFromCSV
parse_trace <- function(traceString){
  # Remove the {} 
  c<-gsub(x=traceString,pattern = "[{]",replacement = "")
  c<-gsub(x=c,pattern = "[}]",replacement = "")
  # Split by ;
  c <- strsplit(c,split=";",fixed=TRUE) %>% 
    unlist()
  c <- c %>% as.numeric()
  c
}

#assumes it is receiving a recording from one field of view
downstrokeCalciumValue <- function(combined_trace,pct = 0.5,contractilityTrace = signal.contractility,calciumTrace = signal.calcium) {
  dPeak <- max(combined_trace$contractilityTrace)
  
  target <- pct * dPeak
  #print(paste("target: ",target))
  thresholds <- numeric()
  
  combined_trace %>% select({{contractilityTrace}},{{calciumTrace}}) %>%
    
  
  for (i in 1:(nrow(combined_trace)-1)) {
    if (combined_trace[i,"contractility"]>target &&
        combined_trace[i+1,"contractility"]<=target) {#check if going down and crossing threshold
      thresholds <- c(thresholds,(approx(combined_trace[i:(i+1),"contractility"],combined_trace[i:(i+1),"calcium"],xout=target)$y)) #linear interpolation to find the matching value in sparse data
    }
  }
  
  return(mean(thresholds))
}

#assumes you have one continous trace (e.g. from one movie from one well/subposition)
calculateCalciumAtContractility <- function(testdata,pct = 0.5,conTrace = signal.contractility,calTrace = signal.calcium) {
  
  con <- testdata %>% select({{conTrace}}) %>% unlist(use.names = F)
  
  cal <- testdata %>% select({{calTrace}}) %>% unlist(use.names = F)
  
  
  dPeak <- max(con)
  dMin <- min(con)
  
  #target is pct * amplitude, so need to account for baseline
  target <- (dPeak-dMin) * pct + dMin
  
  #find points that are greater than the target
  a1 <- con > target
  
  #shift these by 1
  a2 <- a1 %>% shift(-1)
  
  #choose points that are discordant (one less and one greater)
  #and where the starting point was greater than the target
  a3 <- xor(a1,a2) & a1
  
  #get the indices
  ind <- which(a3)
  
  # do the math - y = m * (xTarget-x0)+y0 with xTarget being target
  # and (x0,y0) being the point at con[ind],cal[ind]
  calvalues <- ((cal[ind+1]-cal[ind])/(con[ind+1]-con[ind]))*(target-con[ind])+cal[ind]
  
  return(mean(calvalues))
}

calciumAtContractilityDataframe <- function(testdata) {
  return(data.frame(alias = unique(testdata$alias),cal50 = calculateCalciumAtContractility(testdata)))
}
