#ifndef H_michi2DataClass
#define H_michi2DataClass

#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <regex> gcc >=4.9 (CentOS >=8)
//#include <iterator>
//#include "CrabStringClean.cpp"
//#include "CrabStringReadColumn.cpp"
//#include "CrabTableReadColumn.cpp"
//#include "CrabTableReadInfo.cpp"
//#include "CrabTableGetLineCount.cpp"

using namespace std;

//extern std::string CrabStringTrim(const string &t, int trimflag);



/* Struct Data Class */
/*
 -> XStr    the X array in string format, full size
 -> YStr    the Y array in string format, full size
 -> X       the X array that matches obs, the same size as obs data
 -> Y       the Y array that matches obs, the same size as obs data
 -> CVAR    the column number of variable X and Y
 -> CPAR    the column number of parameter 1,2,3,...
 -> NVAR    the number count of variable, should be 2 because there are only two variable X and Y
 -> NPAR    the number count of parameters, can be any >=1, e.g. 2 if we have T_kin and n_H_2 when doing LVG modelling
 */
class michi2DataClass {
public:
    std::vector<double> X;    std::vector<double> Y;    std::vector<double>  YErr;     std::vector<int>  Matched;
    std::vector<string> XStr; std::vector<string> YStr; std::vector<string>  YErrStr;  long XCol; long YCol; long YErrCol; long XNum; long YNum;
    std::vector<std::vector<double> > FVAR;             std::vector<long>    NVAR;     std::vector<long> CVAR;
    std::vector<std::vector<double> > FPAR;             std::vector<long>    NPAR;     std::vector<long> CPAR;   std::vector<string> TPAR;
    std::ifstream                     FileStream;       std::string          FilePath;
    // FVAR[0] = X,    FVAR[1] = Y
    // FPAR[0] = PAR1, FPAR[1] = PAR2, ...
    std::vector<string> FilterCurveFilePath; // <added><20171001>
    michi2DataClass(const char * InputFile, int verbose = 1);
    ~michi2DataClass();
    const char *michi2sprint(const char* wsPrefix, long wi, const char* wsSuffix);
    int michi2stoi(string str);
    int michi2wstoi(wstring wstr);
    std::vector<double> michi2stod(std::vector<string> strVec);
    std::vector<double> michi2wstod(std::vector<wstring> wstrVec);
    std::vector<std::string> getDataBlock(long lineNumber, long lineCount, int debug = 0); // this function load data block into this->X, this->Y
    std::vector<std::string> readFilterCurveListFile(const char * InputFile);
};



/* Struct For Parallel */
struct mnchi2parallelParams {
    std::vector<michi2DataClass *> SDOBSList;
    std::vector<michi2DataClass *> SDLIBList;
    //std::vector<std::ofstream *> SDOUTList;
    std::vector<std::string> InputLibList;
    std::vector<std::string> OutputTableList;
    long  ibOBS;     //    begin index of OBS list
    long  idOBS;     //          index of OBS list
    long  ieOBS;     //      end index of OBS list
    bool  irOBS;     // rounding index of OBS list 就是计算是否进位
    std::vector<long> ibLIBList; //    begin index of LIB list
    std::vector<long> idLIBList; //          index of LIB list
    std::vector<long> ieLIBList; //      end index of LIB list
    std::vector<bool> irLIBList; // rounding index of LIB list 就是计算是否进位
    long iBegin; long iEnd; long i;
    long nObs; long nLib; long nRow;
};





#endif
