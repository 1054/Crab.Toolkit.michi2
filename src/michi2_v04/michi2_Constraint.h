#ifndef H_michi2Constraint
#define H_michi2Constraint

#include <stdio.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>    // std::transform
#include <string>
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw
#include "michi2_DataClass.h"

using namespace std;


struct ConstraintStructure {
    // OBSOLETE SINCE 20180110
    int fromLIB=0; int toLIB=0; // if value>0 then it's LIB id; if =-1 then use full SED.
    int fromPAR=0; int toPAR=0; // if value>0 then it's PAR id; if =-1 then use integration; if=0 then use LIB index.
    double fromFactor=1; double toFactor=1; // the factor to multiply to from/to data.
    double fromLowerX=0; double toLowerX=0; // the affected X range of from/to data.
    double fromUpperX=0; double toUpperX=0; // the affected X range of from/to data.
    int OperatorType=0; std::string OperatorTypeStr=""; // OperatorType: 0 "="; 1 ">="; -1 "<="; 2 ">"; -2 "<";
    // for example
    // 1 1 -1 2 1 means lib1 par1 <= lib2 par1, which is the case for LVG two component fitting, first component should have colder temperature.
    // 4 -1 0 -1 -1 means lib4(2e5,2.002e5) integration*2.5 = full(8,1000) integration
    // note that the left hand variable is the TO variable
};



class michi2Constraint {
public:
    int fromLIB, toLIB; // if value>0 then it's LIB id; if =-1 then use full SED.
    int fromPAR, toPAR; // if value>0 then it's PAR id; if =-1 then use integration; if=0 then use LIB index.
    double fromFactor; double toFactor; // the factor to multiply to from/to data.
    double fromMathPower; double toMathPower; // the factor to multiply to from/to data.
    double fromMultiplication; double toMultiplication; // the factor to multiply to from/to data.
    double fromAddition; double toAddition; // the value to add to from/to data (after multiplication).
    double fromLowerX; double toLowerX; // the affected X range of from/to data.
    double fromUpperX; double toUpperX; // the affected X range of from/to data.
    int OperatorType; std::string OperatorTypeStr; // OperatorType: 0 "="; 1 ">="; -1 "<="; 2 ">"; -2 "<";
    // for example
    // 1 1 -1 2 1 means lib1 par1 <= lib2 par1, which is the case for LVG two component fitting, first component should have colder temperature.
    // 4 -1 0 -1 -1 means lib4(2e5,2.002e5) integration*2.5 = full(8,1000) integration
    // note that the left hand variable is the TO variable
    
    michi2Constraint(std::vector<std::string> input_arg);
    
    michi2Constraint(const char *input_arg1, const char *input_arg2, const char *input_arg3, const char *input_arg4, const char *input_arg5);
    
    ~michi2Constraint();
    
    void clear();
    
    
    // Extract double numbers from a given string
    // loop to extract all possible double numbers
    std::vector<double> extractStringDouble(std::string InputStr);
    
    int parse(std::vector<std::string> input_args, int verbose = 1);
    
    bool check(struct mnchi2parallelParams * pParams, int debug = 0);
    
    bool check(std::vector<michi2DataClass *> SDLIBS, int debug = 0);
    
    bool check(struct mnchi2parallelParams * pParams, std::vector<michi2DataClass *> SDLIBS, int debug = 0);
};





#endif
