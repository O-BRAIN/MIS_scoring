###  Calculating the MIS score for Ospan task data
###  
###

###
###  Script by Mathis Lammert (lammert@cbs.mpg.de)
###
###  last edit: July 26, 2021 by Hendrik Hartmann (hehartmann@cbs.mpg.de)
###

# Side note
# Raw files of the Ospan task can come in all different shapes and file formats, depending on
# the script and software used to apply the task. Users of this script have to adjust it to their
# file structure accordingly. To provide a general idea of what the script does and what shape the
# data needs to have we provide to dummy data sets. We hope that by running the script with
# this dummy data users are able to follow the structure of this script and can adjust it to their needs.

# Preperation -------------------------------------------------------------

# install the pacman package first and then install not yet installed packages or load already
# installed packages necessary for this script

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, EnvStats, psych, car, pastecs, psych)

# Clean your environment
rm(list = ls())

# Set working directory to location of your Ospan files
setwd("path/to/files")

# load the dummy data sets and save as list elements in "data"
files <- list.files(pattern=".txt")
data <- lapply(files, 
                read.table, 
                sep=",",
                header = T,
                dec = ".")

# Add the subject IDs as name for the list elements of "data"
ids <- read.table(text = files, sep = "_")$V3 %>% parse_number(.) %>% as.numeric()
data <- Map(cbind, ID = ids, data)
names(data) = ids



# Extracting relevant data
ospan.sets <- function(x) 
{
  
  # line-nr of start of test-phase
  ltest <- grep("Start_Test", data[[as.character(x)]]$Code)
  
  # line-nr and contents of set-starts
  lsets <- grep("Test_Setgroesse", data[[as.character(x)]]$Code)
  contset <- as.character(data[[as.character(x)]]$Code[lsets])
  
  # Start result-array: ID, Group and Set-sizes (should also be of length 15)
  res <- data.frame(
    ID = x,
    Setsize = read.table(text = contset, sep = "_" )$V3
  )
  
  # line-nr, line-content and value (1 = correct, 0 = uncorrect) of calculations
  lcalc <- grep("Test_FB", data[[as.character(x)]]$Code)
  contcalc <- as.character(data[[as.character(x)]]$Code[lcalc])
  errorclasscalc <- substr(contcalc, 9, 12)
  corrcalc <- ifelse(substr(contcalc, 9, 12) == "corr", 1,0)
  
  # summing up correct calculation per set
  for(y in 1:length(res$Setsize)) {
    start <- sum(res$Setsize[1:y]) - res$Setsize[y] + 1
    end <-  sum(res$Setsize[1:y])
    res$Correct.Calc[y] <- sum(corrcalc[start:end])
  }
  
  # summing up correct calculation per set: ignoring the first calculation (needed for later applicatin)
  for(y in 1:length(res$Setsize)) {
    start <- sum(res$Setsize[1:y]) - res$Setsize[y] + 1 + 1
    end <-  sum(res$Setsize[1:y])
    res$Correct.Calc.MinusOne[y] <- sum(corrcalc[start:end])
  }
  
  # line-nr, line-content (1 = correct, 0 = uncorrect) of log-intern post-calculation-summary
  lseq1 <- grep("Sequenz", data[[as.character(x)]]$Code)[10:24]
  contseq1 <- as.character(data[[as.character(x)]]$Code[lseq1])
  temp2 <- read.table(text = contseq1, sep = "_" )
  
  # Extracting Set-Nr., correct (presentation-own) Calculation of correct Calculations and given Sequenz from summary
  res$Set.Nr <- temp2$V2
  res$Set.Seq.Nr <- temp2$V3
  res$Correct.Calc.intern <- temp2$V6
  res$Seq.given <- as.character(temp2$V10)
  
  # line-nr & line-content to get the sequence as stated by the subject
  for(i in 1:length(res$Setsize))
  {
    nseq2 <- grep(paste("Test_",res$Set.Nr[i],"_",res$Set.Seq.Nr[i],"_BS_seq",sep = ""), data[[as.character(x)]]$Code)
    ifelse(
      length(nseq2) == 0,         # if no attempt was made to state a sequence, return 0
      res$Seq.stated[i] <- "",
      {                           # otherwise only return the last attempt of each trial
        lseq2 <- tail(nseq2, n = 1)
        contseq2 <- as.character(data[[as.character(x)]]$Code[lseq2])
        temp <- read.table(text = contseq2, sep = "_", colClasses = c("character"))
        res$Seq.stated[i] <- as.character(temp$V6)
      }
    )
  }
  # Return result table
  return(res)
  
} # End function 

# Produce a list with all Ospan_sets to speed up further calls, making use of the above function.
data_sets <- lapply(ids,ospan.sets)  
names(data_sets) = ids


# Computation of the MIS score and math training score
comp.mis <- function (x, df = FALSE)
{
  # Stop function if x is not part of data (e.g. excluded subjects)
  if(!length(data_sets[[as.character(x)]]) > 0) stop("Error: requested ID is not part of the dataset")

  MIS_ov <- data.frame(data_sets[[as.character(x)]], t(mapply(
    function(CalcCor,CalcTot, Seq.given, Seq.stated) {
      
      # Ratio of correct over total results (the first calculation is ignored, as it does not influence the Working Memory)
      RatioCalc <- CalcCor / CalcTot
      
      # Ratio of correct over total letters (every letter from the given sequence is checked within the stated sequence)
      RatioItem <- sum(unlist(strsplit(Seq.given,"")) %in% unlist(strsplit(Seq.stated,""))) / nchar(Seq.given)
      
      # Ratio of the (sum of the) longest stated sequences in the right order over the longest possible sequence.
      # Only count sequences of at least 2 letters
      
      # a) first define start and end points for splitting the given sequence in all possible and needed combinations: 
      # e.g. for a 4-letter-sequence this: 1:2, 1:3, 1:4, 2:3; 3:4
      txx <- expand.grid(Start = 1:nchar(Seq.given), Stop = 1:nchar(Seq.given))
      txx <- txx[which(txx$Start < txx$Stop),]
      txx <- txx[with(txx, order(txx$Start)), ]
      
      # b) creating a list with all possible sequences
      MIS_seq <- data.frame(AllSubstr = mapply(function(SubstrStart, SubstrStop) {
        substr(Seq.given, start = SubstrStart, stop = SubstrStop)
      },
      SubstrStart = txx$Start, 
      SubstrStop = txx$Stop), 
      stringsAsFactors = FALSE)
      
      # c) looking for existance of each possible sequence within the stated Sequence
      MIS_seq$exist <- !sapply(sapply(MIS_seq$AllSubstr, grep, x = Seq.stated), identical, y = integer(0) )
      MIS_seq$length <- nchar(MIS_seq$AllSubstr)
      
      # d) delete all non matched sequences
      MIS_seq <- subset(MIS_seq, exist == TRUE)
      
      # e) Finding the longest sequence, than delete all sequences from the list, that contain letters from this sequence
      MIS_seq_res <- ifelse(nrow(MIS_seq) > 0,max(MIS_seq$length),0)  
      
      # f) calculation the Order ratio:
      RatioOrder <- MIS_seq_res / nchar(Seq.given)
      
      # MIS Formula, as described above
      res_MIS <- RatioCalc * ( RatioItem + RatioOrder ) / 2
      
      
      
      return(data.frame(RatioCalc, RatioItem, RatioOrder, res_MIS))
      
      
    },
    
    CalcCor <- data_sets[[as.character(x)]]$Correct.Calc.MinusOne,
    CalcTot <- data_sets[[as.character(x)]]$Setsize - 1,
    Seq.given <- data_sets[[as.character(x)]]$Seq.given,
    Seq.stated <- data_sets[[as.character(x)]]$Seq.stated
  )))
  
  math.train <- function(x)
  {
    # line-nr for start of test-phase
    lmath <- grep("Start_2.Training", data[[as.character(x)]]$Code)
    # line-nr, line-content and value (1 = correct, 0 = uncorrect) of calculations
    lcalc <- grep("T2_FB", data[[as.character(x)]]$Code)
    contcalc <- as.character(data[[as.character(x)]]$Code[lcalc])
    math_train <- ifelse(substr(contcalc, 7, 10) == "corr", 1,0) %>% sum(.)/length(lcalc)*100
    
    return(math_train)
  }
  math_train <-math.train(x)
  
  #return values
  temp <- data.frame(
    ID = x,
    MIS = sum(unlist(MIS_ov["res_MIS"])),
    MIS_math = sum(unlist(MIS_ov["RatioCalc"])),
    MIS_item = sum(unlist(MIS_ov["RatioItem"])),
    MIS_sequ = sum(unlist(MIS_ov["RatioOrder"])),
    math_accuracy = tail(data_sets[[as.character(x)]]$Correct.Calc.intern, n = 1),
    math_train <- math_train
  )
  
  ifelse(df == FALSE, return(temp), return(MIS_ov)) 
}

# apply function to all your subjects
data_mis <- lapply(ids, comp.mis)
names(data_mis) = ids

