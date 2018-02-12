__author__ = 'belyguerra'
import numpy as np
import pandas as pd
import glob
import os
import scipy
import matplotlib.pyplot as plt
import settings

########################################################################################################
########################################  Script Overview   ############################################
########################################################################################################
# Goal: calculate Task and Item Evoked Pupil Response for SET
# We calculate the change in pupil dilation that occurs after the presentation of an item relative to baseline
# Duration: Item presentation= 1000ms, Interstimulus Interval(ISI)= 500ms, Response duration set to 1000ms
# Baseline TEPR = 0-1500ms; Baseline IEPR = Previous item fixation cross (e.g. 1000-1500ms)
# Pupil change = pupil_diameter_Tn - mean_baseline_pupil_diameter
########################################################################################################


########################################################################################################
######################################  User defined parameters   ######################################
########################################################################################################
#AVG_SAME_Time creates time bins before any calculations
AVG_SAME_TIME = True
#TIME_ROUND_TO defines the size of the time bis in ms (don't set to 0)
TIME_ROUND_TO = 50

# if ROLLING then use rolling avg for pupil data, else uses plain avg of time buckets
ROLLING = False
WINDOW_SIZE = 100

#Set file path where the subject data is located (I use the first path for testing, so do not erase)
#filepath = 'C:\\Users\\belyguerra\\Documents\\ReasoningTraining\\SET\\test'
filepath = '/home/bunge/bguerra/EyeTracking/SET/data/raw/T1/Filtered_and_interpolated_files/trials_removed/'
########################################################################################################


########################################################################################################
######################################  Experiment & Data parameters   #################################
########################################################################################################

#Create a data frame, lists, and name of new columns that will be populated with input and manipulated data
dfSUMMARY = pd.DataFrame()
headers = []
all_lines = []
new_headers = ['Window', 'PupilAvg', 'PupilAvgRoll', 'TrialBaseline', 'TEPR', 'TEPR_fix', 'IEPR', 'WindowTimeNormalized']

#The columns we will be reading from the input data files
colSubject = 2
colTETTime = 5
colDiameterPupilLeftEye = 15
colValidityLeftEye = 17
colDiameterPupilRightEye = 22
colValidityRightEye = 24
colTrialId = 25
colACC = 28
colRT = 29
colCategory = 30
colSETornoSET = 31
colTrainorExp = 32
colPosition_false_shape = 33
colCurrentObject = 40
colInterpPupil = 41 #Pupil data: filtered and interpolated
colFilteredPupil = 42
colRemoveTrials = 43 #Pupil data: filtered, interpolated, and NAs in trials with less than 50% valid data

#Dictionary that has placeholders we can use to later separate fixation vs item presentation stages
#The keys are the value of the CurrentObject window, and we are just setting them equal to a number
mapWindow = {
    'fixation' : 'fixation',
    'firstitem' : '1',
    'seconditem' : '2',
    'thirditems' : '3',
    'fourthitem' : '4',
    'feedback' : 'feedback',
}

#Dictionary that defines the length of each object presentation
#The key is the object (e.g. item1) and the values are another dictionary
#The second dictionary contains the min value and length for that window
####This can probably just be a tupple - don't know if that makes things lighter computationally :)
windowNorm = {
    'item1' : {'min' : 0, 'len' : 1000},
    'fix1' : {'min' : 1000, 'len' : 500},
    'item2' : {'min' : 1500, 'len' : 1000},
    'fix2' : {'min' : 2500, 'len' : 500},
    'item3' : {'min' : 3000, 'len' : 1000},
    'fix3' : {'min' : 4000, 'len' : 500},
    'item4' : {'min' : 4500, 'len' : 1000},
    'fix4' : {'min' : 5500, 'len' : 500},
    'response' : {'min' : 6000, 'len' : 1000}
}
########################################################################################################

########################################################################################################
######################################  Definition of Functions   ######################################
########################################################################################################
#Converts string into float numbers; the data is read as all strings as default
#Since the columns containing 'NA' values do not have anything that we would average, we can set it to whatever
def str2float(strval):
    if not strval:
        return 0.0
    if len(strval) == 0:
        return 0.0
    if strval.upper() == 'NA':
        return 0.0 #this does not matter for averages
    if strval.upper() == 'NAN':
        return 0.0

    return float(strval)

#Converts string into integer numbers; the data is read as all strings as default
#Since the columns containing 'NA' values do not have anything that we would average, we can set it to whatever
def str2int(strval):
    if not strval:
        return 0
    if len(strval) == 0:
        return 0
    if strval.upper() == 'NA':
        return 0 #this does not matter for averages
    if strval.upper() == 'NAN':
        return 0

    return int(strval)

#Calculates the average pupil diameter for both eyes. This average is used for calculating TERP/IERP
#The function takes the values for only one eye if the other eye does not have valid data
#If both eyes have invalid data, the program later on ignores those rows from calculations
#This step is performed before creating time bins
def CalcAvgPupilDialation(validityRight, diamPupilRight, validityLeft, diamPupilLeft):
    diamAvg = -1
    if validityRight == 0 and validityLeft == 0:
        diamAvg = (diamPupilLeft + diamPupilRight) / 2
    elif validityLeft == 4 and validityRight == 0:
        diamAvg = diamPupilRight
    elif validityLeft == 0 and validityRight == 4:
        diamAvg = diamPupilLeft

    return diamAvg

#Checks if we have valid pupil data for a time point
def isValidPupilDialation(validityRight, validityLeft):
    return validityRight == 0 or validityLeft == 0

#Rounds a number (num) to nearest number we want (nearest). Define "nearest" top of file
def roundToNearest(num, nearest):
    return int(nearest * round(float(num)/nearest))

#Creates an empty dictionary for each so it can later take the values from the analysis
def setDefaultValuesForNewColumns(filelines, new_columns):
    row_number = -1
    for row in filelines:
        row_number += 1
        new_columns[row_number] = {}
        for new_header in new_headers:
            new_columns[row_number][new_header] = ''


#
def combineTimeBuckets(filelines, new_columns):
    windowStartTime = -1

    #adding new variable to restart the window time count at 0 in the new trial
    prevTrialId = None

    filelines_new = []
    sameTimePupils = []
    prevTime = None
    prevRow = None
    row_number = -1
    prevWindow = None
    prevCurrentObject = None

    new_window_list = []
    new_normalized_time_list = []
    new_pupil_rolling_avg_list = []

    pupil_rolling_avg_buffer = []

    row_number = -1
    num_buckets_passed = 0
    for row in filelines:
        row_number += 1
        rawTime = str2float(row[colTETTime])
        window = new_columns[row_number]['Window'] #accesses the window value for the row
        #adding new variable to restart the window time count at 0 in the new trial
        currentTrialId = str2int(row[colTrialId])
        currentObject = row[colCurrentObject]
        pupilValidOnly = str2float(row[colRemoveTrials])

        if len(window) == 0:
            continue

        #added
        if currentObject == "Fixation":
            continue

        #added
        if prevTrialId != currentTrialId:
            windowStartTime = rawTime
            pupil_rolling_avg_buffer = []

        time_from_start = rawTime - windowStartTime
        norm_time = roundToNearest(time_from_start, TIME_ROUND_TO)
        new_cols[row_number]['WindowTimeNormalized'] = norm_time
        avgDialation = pupilValidOnly

        if prevRow is None:
            prevRow = row
            prevTime = norm_time
            prevWindow = window
            prevCurrentObject = currentObject

        if norm_time != prevTime or prevWindow != window or prevCurrentObject != currentObject:
            num_buckets_passed += 1
            avgForBin = np.mean(sameTimePupils)
            prevRow[colRemoveTrials] = str(avgForBin)

            sameTimePupils = []
            filelines_new.append(prevRow)
            new_window_list.append(prevWindow)
            new_normalized_time_list.append(prevTime)
            prevRow = row

            if len(pupil_rolling_avg_buffer) > 0:
                first_time = pupil_rolling_avg_buffer[0][0]
                last_time = pupil_rolling_avg_buffer[-1][0]
                if last_time - first_time >= WINDOW_SIZE:
                    pupils_to_avg = filter(lambda d: last_time - d[0] <= WINDOW_SIZE, pupil_rolling_avg_buffer)
                    avg = np.mean([p[1] for p in pupils_to_avg])
                    for x in range(num_buckets_passed):
                        new_pupil_rolling_avg_list.append(str(avg))
                    num_buckets_passed = 0

                    # remove pupil data greater than 100 from latest
                    pupil_rolling_avg_buffer = pupils_to_avg

        sameTimePupils.append(avgDialation)
        pupil_rolling_avg_buffer.append((norm_time, avgDialation))

        prevWindow = window
        prevTrialId = currentTrialId
        prevCurrentObject = currentObject
        prevTime = norm_time

    if prevRow is not None:
        num_buckets_passed += 1
        avgForBin = np.mean(sameTimePupils)
        prevRow[colRemoveTrials] = str(avgForBin)
        filelines_new.append(prevRow)
        new_window_list.append(prevWindow)
        new_normalized_time_list.append(prevTime)

        if len(pupil_rolling_avg_buffer) > 0:
            last_time = pupil_rolling_avg_buffer[-1][0]
            pupils_to_avg = filter(lambda d: last_time - d[0] <= WINDOW_SIZE, pupil_rolling_avg_buffer)
            avg = np.mean([p[1] for p in pupils_to_avg])
            for x in range(num_buckets_passed):
                new_pupil_rolling_avg_list.append(str(avg))
            num_buckets_passed = 0

            # remove pupil data greater than 100 from latest
            pupil_rolling_avg_buffer = filter(lambda d: last_time - d[0] <= WINDOW_SIZE, pupil_rolling_avg_buffer)
        else:
            for x in range(num_buckets_passed):
                new_pupil_rolling_avg_list.append('NA')

    return filelines_new, new_window_list, new_normalized_time_list, new_pupil_rolling_avg_list

#Defines the new names of the stimuli seen by subjects (e.g. first item, first fixation, etc)
def calcWindowNames(filelines, new_columns):
    #timeMap [key = TrialID] [value = another dictionary{key = CurrentObject][value= {first time, end time Object}]}
    timeMap = {}

    #initializing values
    prevCurrentObject = ''
    prevTrialId = -1
    prevTimestamp = 0

    row_number = -1
    for row in filelines:
        row_number += 1
        currentObject = row[colCurrentObject]
        currentTrialId = str2int(row[colTrialId])
        currentTime = str2float(row[colTETTime])

        if currentObject.lower() in mapWindow:
            #we set the Window column according to the mapping we have, at first just using place holder (1,2,etc)
            new_columns[row_number]['Window']  = mapWindow[currentObject.lower()]

            #we check if currentObject changed from the last row
            if prevCurrentObject != currentObject or currentTrialId != prevTrialId:

                #Sets the 'end' time of the previous currentObject as soon as we encounter a new currentObject
                if prevTimestamp != 0:
                    timeMap[prevTrialId][prevCurrentObject]['end'] = prevTimestamp

                if currentTrialId not in timeMap:
                    timeMap[currentTrialId] = {}

                timeMap[currentTrialId][currentObject] = {
                    'start' : currentTime,
                    'end' : 0.0,
                }

            prevCurrentObject = currentObject
            prevTrialId = currentTrialId
            prevTimestamp = currentTime

    timeMap[prevTrialId][prevCurrentObject]['end'] = prevTimestamp

    row_number = -1
    for row in filelines:
        row_number += 1

        currentObject = row[colCurrentObject]
        currentTrialId = str2int(row[colTrialId])
        currentTime = str2float(row[colTETTime])

        if currentObject.lower() in mapWindow:

            if currentObject.lower() == 'fixation' or currentObject.lower() == 'feedback':
                continue

            start = timeMap[currentTrialId][currentObject]['start']
            end = timeMap[currentTrialId][currentObject]['end']

            itemNumber = int(mapWindow[currentObject.lower()])

            if currentTime - start < 1000:
                new_cols[row_number]['Window'] = 'item' + new_cols[row_number]['Window']
            elif (currentTime - start < 1500):
                new_cols[row_number]['Window'] = 'fix' + new_cols[row_number]['Window']
            elif itemNumber == 4:
                new_cols[row_number]['Window'] = 'response'
            else:
                new_cols[row_number]['Window'] = 'fix' + new_cols[row_number]['Window']


# main starts here
### glob is slow!!! ###
for f in glob.glob(filepath + '/*.csv'):
    filename = os.path.basename(f)
    print 'processing file:', filename

    filelines = []
    with open(f, 'r') as csvfile:
        first = True
        for line in csvfile:
            if first:
                first = False
                if len(headers) == 0:
                    headers_data = line.strip().split(',')

                    for h in headers_data:

                        if h[0] == '"' and h[-1] == '"':
                            headers.append(h[1:-1])
                        else:
                            headers.append(h)
                continue

            data_temp = line.strip().split(',')

            if len(data_temp) < 44:
                continue

            data = []
            for d in data_temp:
                if d[0] == '"' and d[-1] == '"':
                    data.append(d[1:-1])
                else:
                    data.append(d)

            filelines.append(data)

    if len(filelines) == 0:
        continue

    new_cols = {}

    #initializing values
    prevCurrentObject = ''
    prevTrialId = -1
    prevTimestamp = 0

    #baselineMap[key=TrialId] [value = dic object{count, total, average for dilations, where currentObject = first item}
    baselineMap = {}

    #fimap [key=TrialId] [value = dic{[key = 'fix#'][value= count, total, average for dilations]
    fixmap = {}

    #windowmap [key=TrialId] [value = dic{[key = window name][value = count, total, average, start time, end time]
    windowmap = {}

    #key is TrialId, value is dic for every field in summary dataframe
    trialDataMap = {}
    biasMap = {}
    encodingStrategyMap = {}

    biasMap['SET'] = []
    biasMap['noSET'] = []

    encodingStrategyMap['span1 or span2'] = []
    encodingStrategyMap['span3'] = []

    ### get rid of practice ###
    filelines = [row for row in filelines if row[colTrainorExp] == 'Exp']
    filelines = [row for row in filelines if len(row[colRemoveTrials]) > 0 and row[colRemoveTrials].upper() != 'NA']

    ### initialize values for new columns ###
    setDefaultValuesForNewColumns(filelines, new_cols)

    ### set window values ###
    calcWindowNames(filelines, new_cols)

    ### if set to avg time bins, then aggregate rows with same time bin ###
    if AVG_SAME_TIME:
        filelines, new_window_list, new_normalized_times, new_pupil_rolling_avg_list = combineTimeBuckets(filelines, new_cols)

        ### since number of rows changed, have to reset new columns ###
        new_cols = {}
        setDefaultValuesForNewColumns(filelines, new_cols)
        row_number = -1
        for row in filelines:
            row_number += 1
            new_cols[row_number]['Window'] = new_window_list[row_number]
            new_cols[row_number]['WindowTimeNormalized'] = new_normalized_times[row_number]
            new_cols[row_number]['PupilAvgRoll'] = new_pupil_rolling_avg_list[row_number]

    ##set baseline map
    row_number = -1
    for row in filelines:
        row_number += 1

        currentObject = row[colCurrentObject]
        currentTrialId = str2int(row[colTrialId])
        currentTime = str2float(row[colTETTime])
        validityLeft = str2int(row[colValidityLeftEye])
        validityRight = str2int(row[colValidityRightEye])
        diamPupilLeft = str2float(row[colDiameterPupilLeftEye])
        diamPupilRight = str2float(row[colDiameterPupilRightEye])
        subject = row[colSubject]
        acc = str2int(row[colACC])
        rt = str2int(row[colRT])
        category = row[colCategory]
        setOrNoSet = row[colSETornoSET]
        positionFalseShape = row[colPosition_false_shape]
        pupilInterp = str2float(row[colInterpPupil])
        pupilFilter = str2float(row[colFilteredPupil])
        pupilValidOnly = str2float(row[colRemoveTrials])

        if len(new_cols[row_number]['PupilAvgRoll']) == 0 or new_cols[row_number]['PupilAvgRoll'].upper() == 'NA':
            continue
        pupilAvgRolling = str2float(new_cols[row_number]['PupilAvgRoll'])

        #initialize average pupil diameter to -1 because it's an invalid value that is filtered out later
        diamAvg = -1

        new_cols[row_number]['PupilAvg'] = pupilValidOnly

        if currentObject.lower() in mapWindow:

            #set constant data for trial
            if currentTrialId not in trialDataMap:
                trialDataMap[currentTrialId] = {}
                trialDataMap[currentTrialId]['Subject'] = subject
                trialDataMap[currentTrialId]['Trial'] = currentTrialId
                trialDataMap[currentTrialId]['ACC'] = acc
                trialDataMap[currentTrialId]['RT'] = rt
                trialDataMap[currentTrialId]['Category'] = category
                trialDataMap[currentTrialId]['SETornoSET'] = setOrNoSet
                trialDataMap[currentTrialId]['Position_false_shape'] = positionFalseShape

            #we check if currentObject changed from the last row
            if prevCurrentObject != currentObject or currentTrialId != prevTrialId:

                #we check if currentObject is equal to first item. If it is we initialize the baselineMap
                if currentObject.lower() == 'firstitem':
                    baselineMap[currentTrialId] = {
                        'count' : 0,
                        'total' : 0.0,
                        'avg' : 0.0,
                    }

            if currentObject.lower() == 'firstitem':
                diamAvg = pupilAvgRolling if ROLLING else pupilValidOnly

                if diamAvg > 0:
                    baselineMap[currentTrialId]['count'] += 1
                    baselineMap[currentTrialId]['total'] += diamAvg

            prevCurrentObject = currentObject
            prevTrialId = currentTrialId
            prevTimestamp = currentTime

    ## calculate avg for each trial id in baseline map
    for trialId in baselineMap.keys():
        if baselineMap[trialId]['count'] > 0:
            baselineMap[trialId]['avg'] = baselineMap[trialId]['total'] / baselineMap[trialId]['count']

    print 'DONE WITH FIRST PASS!'

    ## filling in fix map
    row_number = -1
    for row in filelines:
        row_number += 1

        currentObject = row[colCurrentObject]
        currentTrialId = str2int(row[colTrialId])
        currentTime = str2float(row[colTETTime])
        validityLeft = str2int(row[colValidityLeftEye])
        validityRight = str2int(row[colValidityRightEye])
        diamPupilLeft = str2float(row[colDiameterPupilLeftEye])
        diamPupilRight = str2float(row[colDiameterPupilRightEye])
        subject = row[colSubject]
        acc = str2int(row[colACC])
        rt = str2int(row[colRT])
        category = row[colCategory]
        setOrNoSet = row[colSETornoSET]
        positionFalseShape = row[colPosition_false_shape]
        window = new_cols[row_number]['Window']
        pupilInterp = str2float(row[colInterpPupil])
        pupilFilter = str2float(row[colFilteredPupil])
        pupilValidOnly = str2float(row[colRemoveTrials])
        pupilAvgRolling = str2float(new_cols[row_number]['PupilAvgRoll'])

        diamAvg = -1

        if not window.startswith('fix'):
            continue

        ### d'oh fixation starts with fix ###
        if currentObject.lower() == 'fixation':
            continue

        if currentTrialId not in fixmap:
            fixmap[currentTrialId] = {}

        itemNumber = int(mapWindow[currentObject.lower()])

        if itemNumber not in fixmap[currentTrialId]:
            fixmap[currentTrialId][itemNumber] = {
                'count' : 0,
                'total' : 0.0,
                'avg' : 0.0
            }

        diamAvg = pupilAvgRolling if ROLLING else pupilValidOnly

        if diamAvg > 0:
            fixmap[currentTrialId][itemNumber]['count'] += 1
            fixmap[currentTrialId][itemNumber]['total'] += diamAvg


    for trialId in fixmap.keys():
        for fixNum in fixmap[trialId].keys():
            if fixmap[trialId][fixNum]['count'] > 0:
                fixmap[trialId][fixNum]['avg'] = fixmap[trialId][fixNum]['total'] / fixmap[trialId][fixNum]['count']

    print 'DONE PASS 2!'

    ###pass 3 fill in windowmap based on window values
    prevWindow = -1
    prevCurrentTrialId = -1
    prevTime = -1
    row_number = -1
    for row in filelines:
        row_number += 1

        currentObject = row[colCurrentObject]
        currentTrialId = str2int(row[colTrialId])
        currentTime = str2float(row[colTETTime])
        validityLeft = str2int(row[colValidityLeftEye])
        validityRight = str2int(row[colValidityRightEye])
        diamPupilLeft = str2float(row[colDiameterPupilLeftEye])
        diamPupilRight = str2float(row[colDiameterPupilRightEye])
        subject = row[colSubject]
        acc = str2int(row[colACC])
        rt = str2int(row[colRT])
        category = row[colCategory]
        setOrNoSet = row[colSETornoSET]
        positionFalseShape = row[colPosition_false_shape]
        window = new_cols[row_number]['Window']
        pupilInterp = str2float(row[colInterpPupil])
        pupilFilter = str2float(row[colFilteredPupil])
        pupilValidOnly = str2float(row[colRemoveTrials])
        pupilAvgRolling = str2float(new_cols[row_number]['PupilAvgRoll'])

        diamAvg = -1

        if currentObject.lower() in mapWindow:

            if prevWindow != window or currentTrialId != prevCurrentTrialId:

                if currentTrialId not in windowmap:
                    windowmap[currentTrialId] = {}

                windowmap[currentTrialId][window] = {
                    'start' : currentTime,
                    'end' : 0.0,
                    'count' : 0,
                    'total' : 0.0,
                    'avg' : 0.0
                }

                if prevWindow != -1:
                    windowmap[prevCurrentTrialId][prevWindow]['end'] = prevTime

            diamAvg = pupilAvgRolling if ROLLING else pupilValidOnly

            if diamAvg > 0:
                windowmap[currentTrialId][window]['count'] += 1
                windowmap[currentTrialId][window]['total'] += diamAvg

                ### set biasMap and encodingStrategyMap values
                if acc == 1:
                    if window == 'fix3':
                        if setOrNoSet == 'SET':
                            biasMap['SET'].append(diamAvg)
                        elif positionFalseShape == '3':
                            biasMap['noSET'].append(diamAvg)

                        if setOrNoSet == 'SET' and category == 'span3':
                            encodingStrategyMap['span3'].append(diamAvg)
                        elif setOrNoSet == 'SET':
                            encodingStrategyMap['span1 or span2'].append(diamAvg)

                    elif window == 'fix2':
                        if setOrNoSet == 'SET' and category == 'span3':
                            encodingStrategyMap['span3'].append(diamAvg)
                        elif setOrNoSet == 'SET':
                            encodingStrategyMap['span1 or span2'].append(diamAvg)


            prevWindow = window
            prevCurrentTrialId = currentTrialId
            prevTime = currentTime

    windowmap[prevCurrentTrialId][prevWindow]['end'] = prevTime

    # set default values
    for trialId in windowmap.keys():
        trialDataMap[trialId]['fix1_pupil_av'] = 'NA'
        trialDataMap[trialId]['fix2_pupil_av'] = 'NA'
        trialDataMap[trialId]['fix3_pupil_av'] = 'NA'
        trialDataMap[trialId]['fix4_pupil_av'] = 'NA'

    for trialId in windowmap.keys():
        for window in windowmap[trialId].keys():
            if windowmap[trialId][window]['count'] > 0:
                windowmap[trialId][window]['avg'] = windowmap[trialId][window]['total'] / windowmap[trialId][window]['count']
            if trialId in trialDataMap and window == 'fix1':
                trialDataMap[trialId]['fix1_pupil_av'] = windowmap[trialId][window]['avg']
            elif trialId in trialDataMap and window == 'fix2':
                trialDataMap[trialId]['fix2_pupil_av'] = windowmap[trialId][window]['avg']
            elif trialId in trialDataMap and window == 'fix3':
                trialDataMap[trialId]['fix3_pupil_av'] = windowmap[trialId][window]['avg']
            elif trialId in trialDataMap and window == 'fix4':
                trialDataMap[trialId]['fix4_pupil_av'] = windowmap[trialId][window]['avg']

    ##calc TEPR_fix2fix3avg
    for trialId in windowmap.keys():
        if trialId in trialDataMap:

            fix2Total = 0
            fix2Count = 0
            fix3Total = 0
            fix3Count = 0

            if 'fix2' in windowmap[trialId]:
                fix2Count = windowmap[trialId]['fix2']['count']
                fix2Total = windowmap[trialId]['fix2']['total']
            if 'fix3' in windowmap[trialId]:
                fix3Count = windowmap[trialId]['fix3']['count']
                fix3Total = windowmap[trialId]['fix3']['total']

            total = fix2Total + fix3Total
            count = fix2Count + fix3Count
            if count > 0:
                trialDataMap[trialId]['TEPR_Fix2Fix3avg'] = total / count
            else:
                trialDataMap[trialId]['TEPR_Fix2Fix3avg'] = 'NA';

    print 'DONE PASS 3!'

    ##pass 4 used to set new column values based on baseline and fix maps
    row_number = -1
    for row in filelines:
        row_number += 1

        currentObject = row[colCurrentObject]
        currentTrialId = str2int(row[colTrialId])
        currentTime = str2float(row[colTETTime])
        validityLeft = str2int(row[colValidityLeftEye])
        validityRight = str2int(row[colValidityRightEye])
        diamPupilLeft = str2float(row[colDiameterPupilLeftEye])
        diamPupilRight = str2float(row[colDiameterPupilRightEye])
        subject = row[colSubject]
        acc = str2int(row[colACC])
        rt = str2int(row[colRT])
        category = row[colCategory]
        setOrNoSet = row[colSETornoSET]
        positionFalseShape = row[colPosition_false_shape]
        window = new_cols[row_number]['Window']
        pupilInterp = str2float(row[colInterpPupil])
        pupilFilter = str2float(row[colFilteredPupil])
        pupilValidOnly = str2float(row[colRemoveTrials])
        pupilAvgRolling = str2float(new_cols[row_number]['PupilAvgRoll'])

        diamAvg = -1

        if currentObject.lower() in mapWindow:

            if currentObject.lower() == 'fixation' or currentObject.lower() == 'feedback':
                continue

            ## calculate normalized window time
            #if window in windowNorm:

                ### old logic, normalizing to specific window
                #minNorm = windowNorm[window]['min']
                #lenNorm = windowNorm[window]['len']
                #minTime = windowmap[currentTrialId][window]['start']
                #maxTime = windowmap[currentTrialId][window]['end']

                # if maxTime != minTime:
                #     normTime = minNorm + lenNorm * (currentTime-minTime)/(maxTime - minTime)
                # else:
                #     normTime = 0

                # new_cols[row_number]['WindowTimeNormalized'] = round(normTime, -1) #used to be -1

                ### new logic, round to nearest 20 ms ###
                # if 'item1' in windowmap[currentTrialId]:
                #     minTime = windowmap[currentTrialId]['item1']['start']
                #
                #     roundTime = roundToNearest(currentTime - minTime, TIME_ROUND_TO)
                #     if roundTime < 0:
                #         roundTime = 0
                #
                #     new_cols[row_number]['WindowTimeNormalized'] = roundTime

            itemNumber = int(mapWindow[currentObject.lower()])
            if currentTrialId in baselineMap:
                new_cols[row_number]['TrialBaseline'] = baselineMap[currentTrialId]['avg']

            diamAvg = pupilAvgRolling if ROLLING else pupilValidOnly

            if diamAvg > 0:
                if currentTrialId in baselineMap:
                    new_cols[row_number]['TEPR'] = diamAvg - baselineMap[currentTrialId]['avg']

                if currentTrialId in fixmap:
                    if 1 in fixmap[currentTrialId]:
                        new_cols[row_number]['TEPR_fix'] = diamAvg - fixmap[currentTrialId][1]['avg']

                    if itemNumber - 1 in fixmap[currentTrialId]:
                        avg = fixmap[currentTrialId][itemNumber - 1]['avg']
                        new_cols[row_number]['IEPR'] = diamAvg - avg
                    elif itemNumber in fixmap[currentTrialId]:
                        avg = fixmap[currentTrialId][itemNumber]['avg']
                        new_cols[row_number]['IEPR'] = diamAvg - avg

    ### append to full set of rows
    row_number = -1
    for row in filelines:
        row_number += 1
        data_in_row = []

        ### filter out columns we don't want ###
        cnt = 0
        for data in row:
            if cnt in settings.columns:
                data_in_row.append(data)
            cnt += 1

        for new_header in new_headers:
            data_in_row += [str(new_cols[row_number][new_header])]

        all_lines.append(','.join(data_in_row))

    #append summary data
    tempSummaryDF = pd.DataFrame(trialDataMap.values())
    dfSUMMARY = dfSUMMARY.append(tempSummaryDF, ignore_index=True)

    print 'DONE!!!'

### print out csv###

###filter out columns we don't want in headers ###
cnt = 0
temp_headers = []
for header in headers:
    if cnt in settings.columns:
        temp_headers.append(header)
    cnt += 1
headers = temp_headers
### add new headers
headers += new_headers

with open('SET_PupilTS_LSAT_T1.csv', 'w') as csvfile:
    ### write headers
    csvfile.write('%s\n' % ','.join(headers))
    for line in all_lines:
        csvfile.write('%s\n' % line)

dfSUMMARY.to_csv('SET_PupilAvg_LSAT_T1.csv')

print 'REALLY DONE'
########################################################################################################
