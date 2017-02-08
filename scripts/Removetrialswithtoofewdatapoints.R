##################
### PARAMETERS ###
##################

file.directory <- "~/EyeTracking/SET/data/raw/T2/Filtered_and_interpolated_files"

trial.threshold <- 50    # Please indicate in percent how much valid each trial must at least have. For example,
                         # if trial.threshold = 50, all data points from trials that have more than 50% NAs will
                         # be replaced with NAs. That means, data from those trials will not contribute to any subsequent analyses.

############
### CODE ###
############

## 0 create new folders to put output files and plots
#dir.create(file.path(file.directory, "trials_removed"))   #for output files

## 1. read in files
filenames <- list.files(file.directory, pattern="207.filtered_and_interpolated.90%.csv")
for (file.name in filenames) {
  seqSET.file <- vector()
  seqSET.file <- read.csv(paste(file.directory, "/", file.name, sep=""))
  filename.nocsv <- unlist(strsplit(file.name, split='.csv', fixed=TRUE))[1]   #filename.nocsv is only the name of the file, without ending '.csv'
  
## 2. calculate how much valid data there is in each trial
  seqSET.file$remove_trial <- seqSET.file$interp_av_pupil  # create new column in data file that will have invalid trials removed; start with all the data that interp_av_pupil has
  remove.counter <- 0   #counts how many trials' data points have been removed
  trials <- levels(as.factor(seqSET.file$TrialId))
  for (trial in trials) {
    trial.dat <- seqSET.file$interp_av_pupil[seqSET.file$TrialId == trial & seqSET.file$TrainorExp == 'Exp' & seqSET.file$CurrentObject != 'Fixation']
    valid.points <- sum(!is.na(trial.dat)) / length(trial.dat) * 100

## 3. replace trials that have less than trial.threshold % of valid points with NA
    if (valid.points < trial.threshold) {
      seqSET.file$remove_trial[seqSET.file$TrialId == trial & seqSET.file$TrainorExp == 'Exp' & seqSET.file$CurrentObject != 'Fixation'] <- NA  # subsequently remove data from the column if trials don't have enough data
      cat('Subject ', seqSET.file[1, 'Subject'], ' does not have enough data in trial ', trial, '. Data points of this trial will be removed.\n', sep='')
      remove.counter <- remove.counter + 1
    }
  }

write.csv(seqSET.file, file=paste(file.directory, "/trials_removed/", filename.nocsv, '.', remove.counter, '_invalid_trials_removed.csv', sep=''))
}
