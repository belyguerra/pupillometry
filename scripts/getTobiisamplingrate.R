######################
### FILE.DIRECTORY ###
######################
#file.directory <- "C:/Users/Maria/Dropbox/New Data/sequentialset/raw data/Filtered_and_interpolated_files/Outputfiles/Only EXP rows and rel.time and rel.pupil"

file.directory <- "~/EyeTracking/SET/data/raw/T1"

#create an empty list
allsubj.rowdurations <- vector()

## 1. read in files
subj.list <- vector()
filenames <- list.files(file.directory, pattern="*188-1.csv")
for (file.name in filenames) {
  file <- vector()
  file <- read.delim(paste(file.directory, "/", file.name, sep=""), sep='\t', header=TRUE)
  subj.list <- c(subj.list, file[1,'Subject'])

  row.times <- vector()
  for (row in 1:dim(file[2])) {
    row.time <- file$TimestampSec[row] + file$TimestampMicrosec[row]*0.000001
    row.times <- c(row.times, row.time)
  }
  
  row.durations <- vector()
  for (row in 1:(length(row.times) - 1)) {
    row.duration <- row.times[row + 1] - row.times[row]
    row.durations <- c(row.durations, row.duration)
  }
  
  subj.rowdurations <- cbind(row.durations, rep(file[1,'Subject'], length(row.durations)))
  colnames(subj.rowdurations) <- c("row.durations", "subj.id")
  
  allsubj.rowdurations <- rbind(allsubj.rowdurations, subj.rowdurations)
  allsubj.rowdurations <- as.data.frame(allsubj.rowdurations)
 }


summaries <- vector()
#for (subj in levels(as.factor(allsubj.rowdurations(subj.id)))) {
  #s <- summary(allsubj.rowdurations$row.durations[allsubj.rowdurations$subj.id == subj])
  s <- summary(allsubj.rowdurations$row.durations[allsubj.rowdurations$subj.id])
  summaries <- c(summaries, s[4])
#}

#summary.df <- cbind(summaries, subj.list)
summary.df <- cbind(summaries, filenames)
write.table(summary.df, file = "~/EyeTracking/SET/data/SET_true_samplingrate_LSAT_T1.csv", sep=" ", append=TRUE)


#############################################################################################
### MARIA RESULTS: THE TRUE SAMPLING RATE IS AT LEAST 18ms (FOR SUBJ 25 & 26 EVEN 26ms AND 29ms) ###
#############################################################################################

# > summary.df
# summaries subj.list
# Mean   0.01821        10
# Mean   0.02084        11
# Mean   0.01800        12
# Mean   0.01802        13
# Mean   0.01873        14
# Mean   0.01868        15
# Mean   0.01820        16
# Mean   0.01864        17
# Mean   0.01910        18
# Mean   0.01846        19
# Mean   0.01840        20
# Mean   0.01811        21
# Mean   0.01904        22
# Mean   0.01834        23
# Mean   0.01820        24
# Mean   #0.02630        25#
# Mean   #0.02927        26#
# Mean   0.01846        28
# Mean   0.01873        29
# Mean   0.01848         3
# Mean   0.01896        30
# Mean   0.01858         4
# Mean   0.01851         5
# Mean   0.01910         6
# Mean   0.01858         7
# Mean   0.01988         8
# Mean   0.01881         9

# > mean(summary.df[1])
# [1] 0.01821

#############################################################################################
### BELEN RESULTS: THE TRUE SAMPLING RATE IS AT LEAST 10ms 
#############################################################################################
# > summary.df
# summaries subj.list
# Mean   0.01049       101
# Mean   0.01046       102
# Mean   0.01101       103
# Mean   0.01014       104
# Mean   0.01000       105
# Mean   0.01091       106
# Mean   0.01063       107
# Mean   0.01035       108
# Mean   0.01122       109
# Mean   0.01053       110
# Mean   0.01028       111
# Mean   0.01050       112
# Mean   0.01071       113
# Mean   0.01041       114
# Mean   0.01141       115
# Mean   0.01043       116

#> mean(summary.df[1])
#[1] 0.01049
