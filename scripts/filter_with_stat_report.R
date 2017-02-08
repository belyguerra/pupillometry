################################################################################
####### PLEASE SPECIFY THE PARAMETERS AND THE FILE DIRECTORY HERE: #############
################################################################################

file.directory <- "~/EyeTracking/SET/data/raw/T1"
                    # indicate the folder that contains all the files you want to filter and interpolate (must be .csv files)
                    # put the directory into " " and use '/' (instead of '\')

se <- 5             # number of standard errors around loess regression function where pupil values are still accepted
numPoints <- 80     # number of data points on both sides of the regressed point that are used to calculate the loess regression
max.interp <- 50                 # only gaps that are shorter than max.interp rows will be interpolated; gaps that are longer (where
                                  # more consecutive rows of data are missing) will not be interpolated. With a sampling rate of 8ms, 25 rows correspond to 200ms.

min.percent.of.valid.data <- 25   # indicates how many percent of data points in the raw file need to
                                  # be valid in order for the script to execute; errors are likely to occur
                                  # if less than 30% of the data are valid; good interpolation needs more than 50%
                                  # (but decision may be made on trail basis)


################################################################################
###### THE 'REAL' CODE STARTS HERE: ############################################
################################################################################

# create new folders to put output files and plots
#dir.create(file.path(file.directory, "Filtered_and_interpolated_files"))   #for output files
#dir.create(file.path(file.directory, "Plots"))                             #for output plots

# prepare data.frame that contains stats about how much data has been interpolated and filterd for each subject
filter_stats <- data.frame(Subject = NA, row.numb = NA, valid.data = NA, block = NA, filtered_numb = NA, filtered_perc = NA, interpolated_numb = NA, interpolated_perc = NA)
resting_stats <- data.frame(Subject = NA, resting = NA)

# read file into R
filenames <- list.files(file.directory, pattern = "*178-1.csv")
for (file.name in filenames) {
  raw.dat <- vector()
  raw.dat <- read.delim(paste(file.directory, "/", file.name, sep = ""), sep="\t", header=TRUE)
  #raw.dat <-read.csv(paste(file.directory,"/",file.name,sep=""),header=TRUE)
  
  filename.nocsv_v1 <- unlist(strsplit(file.name, split = '.', fixed = TRUE))[1]   #only the name of the file, without ending '.csv'
  filename.nocsv <- unlist(strsplit(filename.nocsv_v1, split = '-', fixed = TRUE))[2]

  # replace invalid pupil data with "NA"
  raw.dat$DiameterPupilLeftEye[raw.dat$DiameterPupilLeftEye == -1] <- NA
  raw.dat$DiameterPupilLeftEye[raw.dat$ValidityLeftEye == 4] <- NA
  raw.dat$DiameterPupilRightEye[raw.dat$DiameterPupilRightEye == -1] <- NA
  raw.dat$DiameterPupilRightEye[raw.dat$ValidityRightEye == 4] <- NA
  
  numb.of.val <- sum(!is.na(raw.dat$DiameterPupilLeftEye)) + sum(!is.na(raw.dat$DiameterPupilRightEye))
  perc.of.val <- numb.of.val*100 / (dim(raw.dat)[1]*2)
  
  cat('In file ', filename.nocsv, ', ', numb.of.val, ' out of ', dim(raw.dat)[1]*2, ' data points ',
      '(', round(perc.of.val, 2), '%)', ' have valid values.\n', sep = '')
  row.numb <- dim(raw.dat)[1]
  valid.data <- numb.of.val / 2
  
  # if less than 'min.percent.of.valid.data' data points have valid values, a copy of the file is made that indicates how many data points have valid values, but the two additional columns are not added; that means, the rest of the script does not execute
  if (perc.of.val < min.percent.of.valid.data) {
    write.csv(raw.dat, file = paste(file.directory, "/", file.name, '.filtered_and_interpolated.', round(perc.of.val), '%.csv', sep = ''))
    
  } else {
  
    # This function takes the average of diameters or whichever eye has data if only one eye has data
    get.values <- function(dat) {
      rVal <- (dat$DiameterPupilLeftEye + dat$DiameterPupilRightEye) / 2
      rVal[is.na(rVal) & !is.na(dat$DiameterPupilLeftEye)] <- dat$DiameterPupilLeftEye[is.na(rVal) & !is.na(dat$DiameterPupilLeftEye)]
      rVal[is.na(rVal) & !is.na(dat$DiameterPupilRightEye)] <- dat$DiameterPupilRightEye[is.na(rVal) & !is.na(dat$DiameterPupilRightEye)]
      return(rVal)
    }
    
    # collector for cleaned pupil data
    filtered_av_pupil <- vector()
    interp_av_pupil <- vector()
    
    #split file into pieces of 8010 rows (blocksize); with overlap of 80 rows (numPoints)
    blocksize <- 8000
    for (block in seq(1, (dim(raw.dat)[1]), blocksize)) {
      block.dat <- matrix()

      if ((block + blocksize) < dim(raw.dat)[1]) {
        block.dat <- raw.dat[block:(block + blocksize + 2 * numPoints),]
      } else {
        block.dat <- raw.dat[block:dim(raw.dat)[1],]
      }

      span <- numPoints / dim(block.dat)[1]
      
      # get data values (x = ID, y = pupildiameter)
      dat <- data.frame()
      dat <- data.frame(X = block.dat$ID, Y = get.values(block.dat))
      
      val.percent <- sum(!is.na(dat$Y))/dim(dat)[1]
      cat('In file ', filename.nocsv, ', block ', block, ': ', sum(!is.na(dat$Y)), ' out of ', dim(dat)[1],
          ' points (', round(val.percent * 100,2), '%) have valid values.\n', sep='')
      
      # fit loess to data
      dat.lo <- loess(Y ~ X, dat, span = span)
      dat.fit <- predict(dat.lo, se = 1)
      
      # store fit values as well as upper and lower bounds
      fit <- rep(NA, dim(dat)[1])
      fit[!is.na(dat$Y)] <- dat.fit$fit
      lb <- rep(NA, dim(dat)[1])
      lb[!is.na(dat$Y)] <- dat.fit$fit - (dat.fit$se.fit * se)
      ub <- rep(NA, dim(dat)[1])
      ub[!is.na(dat$Y)] <- dat.fit$fit + (dat.fit$se.fit * se)
      
      # diagnostic plots and message
      png(filename = paste(file.directory, "/Plots/", filename.nocsv, "loess.rows", block, ".png", sep=''))
      plot(dat$X, dat$Y, pch = '.', type = 'l')
      lines(dat$X, fit, pch = '.', type = 'l', col = 'blue')
      lines(dat$X, lb, pch = '.', type = 'l', col = 'red')
      lines(dat$X, ub, pch = '.', type = 'l', col = 'red')
      dev.off()
      
      passFilter <- vector()
      passFilter <- ifelse(is.na(dat$Y) | dat$Y < lb | dat$Y > ub, 0, 1)
      
      png(filename = paste(file.directory, "/Plots/", filename.nocsv, "filtered.rows", block, ".png", sep = ''))      
      plot(dat$X, dat$Y, pch = 20, col = ifelse(passFilter == 0, 'red', 'black'),
           main = paste('Subject', raw.dat[1, 'Subject'], ', data points', dat$X[1], ' - ', dat$X[length(dat$X)], sep = ' '))
      dev.off()
      plot(dat$X, dat$Y, pch = 20, col = ifelse(passFilter == 0, 'red', 'black'),
           main = paste('Subject', raw.dat[1, 'Subject'], ', data points', dat$X[1], ' - ', dat$X[length(dat$X)], sep = ' '))
      
      perc <- (sum(!is.na(dat$Y)) - sum(passFilter))/sum(!is.na(dat$Y))
      filtered_numb <- (sum(!is.na(dat$Y)) - sum(passFilter))
      cat('In file ', filename.nocsv, ', block ', block, ': ', filtered_numb, ' out of ',
          sum(!is.na(dat$Y)), ' points (', round(perc * 100, 2), '%) have been filtered.\n', sep = '')
      
      # replace filtered values with NA (passFilter is 0 when points should be filtered); safe filtered pupil data to filtered_av_pupil
      dat$Y[!passFilter] <- NA
      
      if (block == 1) {
        if (blocksize > dim(raw.dat)[1]) {
          filtered_av_pupil <- c(filtered_av_pupil, dat$Y[1:dim(raw.dat)[1]])
        } else {
          filtered_av_pupil <- c(filtered_av_pupil, dat$Y[1:(blocksize + numPoints)])
      }}
      else if ((block + blocksize) < dim(raw.dat)[1] & block > 1) {
        filtered_av_pupil <- c(filtered_av_pupil, dat$Y[(numPoints + 1):(blocksize + numPoints)])
      }
      else if ((block + blocksize) > dim(raw.dat)[1]) {
        filtered_av_pupil <- c(filtered_av_pupil, dat$Y[(numPoints + 1):(dim(block.dat)[1])])
      }
      
      # replace missing values with 'loess' fit points
      Missinglist <- vector()
      Missinglist <- dat$X[is.na(dat$Y) & !is.na(dat$X)]
      dat.lo <- loess(Y ~ X, dat, span = span)
      missing.predict <- vector()
      missing.predict <- predict(dat.lo, data.frame(X = Missinglist))
      dat$Y[is.na(dat$Y) & !is.na(dat$X)] <- missing.predict
      
      # replace interpolated data points for missing sequences > 200ms (max.interp = 25) with NA - if Missinglist contains more than 25 values at all
      if (length(Missinglist) > max.interp) {
        dat$drop <- rep(0, length(dat$X))           #vector collecting info about length of missing data
        for (i in 1:(length(Missinglist) - max.interp)) {
          if ((Missinglist[i] + max.interp) == (Missinglist[i + max.interp])) {
            for (j in 0:max.interp) {
              dat$drop[dat$X == Missinglist[i+j]] <- 1
            }
          }
        }
        dat$Y[dat$drop == 1] <- NA         #replace pupildata with NAs if to many missing values in a row
      }

      # plot interpolated points and message
      prozent <- (length(Missinglist) - sum(is.na(dat$Y)))/length(Missinglist)
      interpolated_numb <- (length(Missinglist) - sum(is.na(dat$Y)))
      cat('In file ', filename.nocsv, ', block ', block, ': ', interpolated_numb,
          ' out of ', length(Missinglist), ' missing data points (', round(prozent*100, 2),
          '%) have been interpolated.\n', sep = '')
      
      filter_stats_row <- data.frame(Subject = raw.dat$Subject[1], row.numb = row.numb, valid.data = valid.data, block = block, filtered_numb = filtered_numb, filtered_perc = perc * 100, interpolated_numb = interpolated_numb, interpolated_perc = prozent*100)
      filter_stats <- rbind(filter_stats, filter_stats_row)
      
      png(filename = paste(file.directory, "/Plots/", filename.nocsv, "interpolated.rows", block, ".png", sep=''))      
      plot(dat$X, dat$Y, pch=20, col = ifelse(dat$X[!is.na(dat$X)] %in% Missinglist, 'forestgreen', 'black'),
           main = paste('Subject', raw.dat[1, 'Subject'], ' (data points', dat$X[1], ' - ', dat$X[length(dat$X)],
                        '): interpolated data points', sep=' '))
      dev.off()
      plot(dat$X, dat$Y, pch=20, col = ifelse(dat$X[!is.na(dat$X)] %in% Missinglist, 'forestgreen', 'black'),
           main = paste('Subject', raw.dat[1, 'Subject'], ' (data points', dat$X[1], ' - ', dat$X[length(dat$X)],
                        '): interpolated data points', sep=' '))
      
      # safe interpolated data into variable interp_av_pupil; notice overlap of 50 rows (numPoints)
      if (block == 1) {
        if (blocksize > dim(raw.dat)[1]) {
          interp_av_pupil <- c(interp_av_pupil, dat$Y[1:dim(raw.dat)[1]])
        } else {
          interp_av_pupil <- c(interp_av_pupil, dat$Y[1:(blocksize + numPoints)])
        }}
      else if ((block + blocksize) < dim(raw.dat)[1] & block > 1) {
        interp_av_pupil <- c(interp_av_pupil, dat$Y[(numPoints + 1):(blocksize + numPoints)])
      }
      else if ((block + blocksize) > dim(raw.dat)[1]) {
        interp_av_pupil <- c(interp_av_pupil, dat$Y[(numPoints + 1):(dim(block.dat)[1])])
      }
    }
    #attach cleaned data as new column to data
    raw.dat$interp_av_pupil <- interp_av_pupil
    raw.dat$filtered_av_pupil <- filtered_av_pupil
    write.csv(raw.dat, file = paste(file.directory, "/Filtered_and_interpolated_files/", filename.nocsv,
                                  '.filtered_and_interpolated.', round(perc.of.val), '%.csv', sep = ''))
    
    cat('Subject', raw.dat$Subject[1], 'has an average resting pupil diameter of', mean(interp_av_pupil, na.rm = T), 'cm.\n')
    resting_stats.row <- data.frame(Subject = raw.dat$Subject[1], resting = mean(interp_av_pupil, na.rm = T))
    resting_stats <- rbind(resting_stats, resting_stats.row)
}}
filter_stats <- filter_stats[2:dim(filter_stats)[1],]

filter_stats$subj_filtered_row_numb <- NA
for (subj in levels(as.factor(filter_stats$Subject))) {
  filter_stats$subj_filtered_row_numb[filter_stats$Subject == subj] <- sum(filter_stats$filtered_numb[filter_stats$Subject == subj])
}

filter_stats$subj_filtered_row_perc <- filter_stats$subj_filtered_row_numb / filter_stats$valid.data
filter_stats$subj_filtered_row_perc_total <- filter_stats$subj_filtered_row_numb / filter_stats$row.numb
mean(filter_stats$subj_filtered_row_perc)
mean(filter_stats$subj_filtered_row_perc_total)

filter_stats$subj_interpolated_row_numb <- NA
for (subj in levels(as.factor(filter_stats$Subject))) {
  filter_stats$subj_interpolated_row_numb[filter_stats$Subject == subj] <- sum(filter_stats$interpolated_numb[filter_stats$Subject == subj])
}

filter_stats$subj_interpolated_row_perc <- filter_stats$subj_interpolated_row_numb / filter_stats$valid.data
filter_stats$subj_interpolated_row_perc_total <- filter_stats$subj_interpolated_row_numb / filter_stats$row.numb
mean(filter_stats$subj_interpolated_row_perc)
mean(filter_stats$subj_interpolated_row_perc_total)

resting_stats <- resting_stats[2:dim(resting_stats)[1],]
write.csv(resting_stats, '~/EyeTracking/SET/data/RestingStats_T1_v3.csv')
