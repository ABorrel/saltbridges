import time
import formatCharacter
import repertory


def initAction(actionName):
    """Initialisation of log file
    in: action
    out: file open and date Start"""
    
    dateStart = time.strftime('%X %x')
    dateStartFileName = formatCharacter.date(dateStart)
    print dateStart, ": ", actionName
    begin = time.clock()
    fileLog = open(repertory.logFile() + dateStartFileName + actionName.replace (" ", ""), "w")
    fileLog.write(actionName + "\n")
    fileLog.write(dateStart + "\n")

    return begin, fileLog


def endAction(actionName, startTime, logFile):
    """End log file
    in: name action, time start, log file
    out: close file"""
    
    dateEnd = time.strftime('%X %x')
    end = time.clock()
    timeExecution = end - startTime
    hour = int(timeExecution / 3600)
    timeRest = timeExecution - (hour * 3600)
    minute = int(timeRest / 60)
    second = float(timeRest - (minute * 60))
    print
    print "Time execution :", hour, "h", minute, "min", second, "s"
    logFile.write("Time execution :" + str(hour) + "h" + str(minute) + "min" + str(second) + "s" + "\n")
    print dateEnd, ": End ", actionName
    logFile.write(dateEnd + ": End " + actionName)
    logFile.close()
