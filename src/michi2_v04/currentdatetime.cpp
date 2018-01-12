#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include "currentdatetime.h"

const std::string currentDateTime() {
    // Get current date/time, format is YYYY-MM-DD.HH:mm:ss
    // http://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X %Z", &tstruct);
    return buf;
};




