#ifndef H_michi2Constraint
#define H_michi2Constraint

#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>    // std::transform
#include <string>
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw

using namespace std;


struct ConstraintStructure {
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
    double fromMultiplication; double toMultiplication; // the factor to multiply to from/to data.
    double fromAddition; double toAddition; // the value to add to from/to data (after multiplication).
    double fromLowerX; double toLowerX; // the affected X range of from/to data.
    double fromUpperX; double toUpperX; // the affected X range of from/to data.
    int OperatorType; std::string OperatorTypeStr; // OperatorType: 0 "="; 1 ">="; -1 "<="; 2 ">"; -2 "<";
    // for example
    // 1 1 -1 2 1 means lib1 par1 <= lib2 par1, which is the case for LVG two component fitting, first component should have colder temperature.
    // 4 -1 0 -1 -1 means lib4(2e5,2.002e5) integration*2.5 = full(8,1000) integration
    // note that the left hand variable is the TO variable
    
    michi2Constraint(std::vector<std::string> input_arg) {
        this->clear();
        this->parse(input_arg);
    };
    
    michi2Constraint(const char *input_arg1, const char *input_arg2, const char *input_arg3, const char *input_arg4, const char *input_arg5) {
        this->clear();
        std::vector<std::string> input_args;
        input_args.push_back(input_arg1);
        input_args.push_back(input_arg2);
        input_args.push_back(input_arg3);
        input_args.push_back(input_arg4);
        input_args.push_back(input_arg5);
        this->parse(input_args);
    };
    
    ~michi2Constraint() {};
    
    void clear() {
        fromLIB=0; toLIB=0; // if value>0 then it's LIB id; if =-1 then use full SED.
        fromPAR=0; toPAR=0; // if value>0 then it's PAR id; if =-1 then use integration; if=0 then use LIB index.
        fromFactor=1; toFactor=1; // the factor to multiply to from/to data (obsolete).
        fromMultiplication=1.0; toMultiplication=1.0; // the factor to multiply to from/to data.
        fromAddition=0.0; toAddition=0.0; // the value to add to from/to data (after multiplication).
        fromLowerX=0; toLowerX=0; // the affected X range of from/to data.
        fromUpperX=0; toUpperX=0; // the affected X range of from/to data.
        OperatorType=0; OperatorTypeStr=""; // OperatorType: 0 "="; 1 ">="; -1 "<="; 2 ">"; -2 "<";
        // for example
        // 1 1 -1 2 1 means lib1 par1 <= lib2 par1, which is the case for LVG two component fitting, first component should have colder temperature.
        // 4 -1 0 -1 -1 means lib4(2e5,2.002e5) integration*2.5 = full(8,1000) integration
        // note that the left hand variable is the TO variable
    };
    
    
    // Extract double numbers from a given string
    // loop to extract all possible double numbers
    std::vector<double> extractStringDouble(std::string InputStr) {
        double TempDbl = 0.0;
        std::string TempStr = InputStr+" "; // in case the string got empty
        std::string ExctStr = ""; // extracted str
        std::vector<double> OutputDbl;
        while (std::string::npos != TempStr.find_first_of(".+-Ee0123456789")) {
            TempStr = TempStr.substr(TempStr.find_first_of(".+-Ee0123456789"));
            size_t TempPos = TempStr.find_first_not_of(".+-Ee0123456789");
            if(std::string::npos != TempPos) {
                ExctStr = TempStr.substr(0,TempPos);
                TempStr = TempStr.substr(TempPos);
            } else {
                ExctStr = TempStr;
            }
            if(ExctStr!="") {
                if(0 == ExctStr.find_first_of(".+-0123456789")) {
                    //std::cout << "DEBUG: ExctStr = " << ExctStr << std::endl;
                    TempDbl = std::stod(ExctStr.c_str());
                    OutputDbl.push_back(TempDbl);
                }
                ExctStr="";
            }
        }
        return OutputDbl;
    }
    
    int parse(std::vector<std::string> input_args, int verbose = 1) {
        clear();
        std::string TempConstraintStr;
        std::string PrevConstraintStr;
        if(verbose>=1) {
            std::cout << "Input Constraint:";
        }
        int i=-1;
        for(int j=1; j<=5; j++) {
            TempConstraintStr = input_args[i+j];
            //
            // convert string to upper case
            std::transform(TempConstraintStr.begin(),TempConstraintStr.end(),TempConstraintStr.begin(),::toupper);
            //
            // store previous string
            if(j==1) {PrevConstraintStr=TempConstraintStr;} // PrevConstraintStr is used as previous container.
            //
            // if input "lib3" LIB id, then set constraint from/to LIB 3 (number starting from 1).
            // if input "full" then set constraint from/to full SED.
            //
            // if input "index" then set constraint on index.
            // if input "integration" then set constraint on integration.
            // if input "par3" LIB PAR id, then set constraint on LIB PAR.
            //
            if (TempConstraintStr=="SED" || (std::string::npos!=TempConstraintStr.find("FULL")) ) {
                //
                // This means the constraint is applied to the whole SED instead of each LIB
                // e.g. LIB3 NORM(2.14e5) EQ SED INT(8,1000)
                //
                if(j==4) {this->fromLIB = -1;}
                else if(j==1) {this->toLIB = -1;}
                
            } else if ( TempConstraintStr=="L" || (std::string::npos!=TempConstraintStr.find("LIR")) || (std::string::npos!=TempConstraintStr.find("VLV")) ) {
                //
                // This means the constraint calculates IR_LUMINOSITY instead of using PARAMETER of the target LIB
                // e.g. LIB3 NORM(2.14e5) EQ SED LIR(8,1000)
                //
                if(j==5) {this->fromPAR = -3; this->fromFactor = 1.0; this->fromLowerX = 0.0; this->fromUpperX = 0.0;}
                else if(j==2) {this->toPAR = -3; this->toFactor = 1.0; this->toLowerX = 0.0; this->toUpperX = 0.0;}
                //
                // read INT X range
                std::vector<double> TempConstraintDbl = extractStringDouble(TempConstraintStr);
                if(TempConstraintDbl.size()>1) {
                    if(j==5) {this->fromLowerX = TempConstraintDbl[0]; this->fromUpperX = TempConstraintDbl[1];}
                    else if(j==2) {this->toLowerX = TempConstraintDbl[0]; this->toUpperX = TempConstraintDbl[1];}
                }
                //
                // read INT multiplication factor if available
                if(TempConstraintDbl.size()>2) {
                    if(j==5) {this->fromMultiplication = TempConstraintDbl[2]; }
                    else if(j==2) {this->toMultiplication = TempConstraintDbl[2]; }
                }
                //
                // read INT addition factor if available
                if(TempConstraintDbl.size()>3) {
                    if(j==5) {this->fromAddition = TempConstraintDbl[3]; }
                    else if(j==2) {this->toAddition = TempConstraintDbl[3]; }
                }
                //
                // clear
                TempConstraintDbl.clear();
                
            } else if ( (std::string::npos!=TempConstraintStr.find("INT")) || (std::string::npos!=TempConstraintStr.find("INTEGR"))) {
                //
                // This means the constraint calculates INTEGRATION instead of using PARAMETER of the target LIB
                // e.g. LIB3 NORM(2.14e5) EQ SED INT(8,1000)
                //
                if(j==5) {this->fromPAR = -2; this->fromFactor = 1.0; this->fromLowerX = 0.0; this->fromUpperX = 0.0;}
                else if(j==2) {this->toPAR = -2; this->toFactor = 1.0; this->toLowerX = 0.0; this->toUpperX = 0.0;}
                //
                // read INT X range
                std::vector<double> TempConstraintDbl = extractStringDouble(TempConstraintStr);
                if(TempConstraintDbl.size()>1) {
                    if(j==5) {this->fromLowerX = TempConstraintDbl[0]; this->fromUpperX = TempConstraintDbl[1];}
                    else if(j==2) {this->toLowerX = TempConstraintDbl[0]; this->toUpperX = TempConstraintDbl[1];}
                }
                //
                // read INT multiplication factor if available
                if(TempConstraintDbl.size()>2) {
                    if(j==5) {this->fromMultiplication = TempConstraintDbl[2]; }
                    else if(j==2) {this->toMultiplication = TempConstraintDbl[2]; }
                }
                //
                // read INT addition factor if available
                if(TempConstraintDbl.size()>3) {
                    if(j==5) {this->fromAddition = TempConstraintDbl[3]; }
                    else if(j==2) {this->toAddition = TempConstraintDbl[3]; }
                }
                //
                // clear
                TempConstraintDbl.clear();
                
            } else if (TempConstraintStr=="F" || (std::string::npos!=TempConstraintStr.find("FLUX")) || (std::string::npos!=TempConstraintStr.find("NORM"))) {
                //
                // This means the constraint calculates NORMALIZATION instead of using PARAMETER of the target LIB
                // e.g. LIB3 NORM(2.14e5) EQ SED INT(8,1000)
                //
                if(j==5) {this->fromPAR = -1; this->fromFactor = 1.0; this->fromLowerX = 0.0; this->fromUpperX = 0.0;}
                else if(j==2) {this->toPAR = -1; this->toFactor = 1.0; this->toLowerX = 0.0; this->toUpperX = 0.0;}
                //
                // read NORM X pos
                std::vector<double> TempConstraintDbl = extractStringDouble(TempConstraintStr);
                if(TempConstraintDbl.size()>0) {
                    if(j==5) {this->fromLowerX = TempConstraintDbl[0]; this->fromUpperX = TempConstraintDbl[0];}
                    else if(j==2) {this->toLowerX = TempConstraintDbl[0]; this->toUpperX = TempConstraintDbl[0];}
                }
                //
                // read NORM multiplication factor if available
                if(TempConstraintDbl.size()>1) {
                    if(j==5) {this->fromMultiplication = TempConstraintDbl[1]; }
                    else if(j==2) {this->toMultiplication = TempConstraintDbl[1]; }
                }
                //
                // read NORM addition factor if available
                if(TempConstraintDbl.size()>2) {
                    if(j==5) {this->fromAddition = TempConstraintDbl[2]; }
                    else if(j==2) {this->toAddition = TempConstraintDbl[2]; }
                }
                //
                // clear
                TempConstraintDbl.clear();
                
            } else if (TempConstraintStr=="I" || (std::string::npos!=TempConstraintStr.find("INDEX"))) {
                //
                // This means the constraint uses INDEX instead of PARAMETER of the target LIB
                // e.g. LIB2 INDEX EQ LIB3 INDEX
                //
                if(j==5) {this->fromPAR = 0;}
                else if(j==2) {this->toPAR = 0;}
                
            } else if(std::string::npos != TempConstraintStr.find_first_of("0123456789")) {
                //
                // it's normal mode, set constraint on some PAR of some LIB
                // e.g. LIB2 PAR2 GE LIB3 PAR2
                //
                std::string TempConstraintStr2 = TempConstraintStr.substr(TempConstraintStr.find_first_of("0123456789"));
                if(j==1){this->toLIB = std::stoi(TempConstraintStr2);}
                else if(j==2){this->toPAR = std::stoi(TempConstraintStr2);}
                else if(j==4){this->fromLIB = std::stoi(TempConstraintStr2);}
                else if(j==5){this->fromPAR = std::stoi(TempConstraintStr2);}
                
            } else if(j==3) {
                //
                // j==3 is the operator
                //
                if     (std::string::npos != TempConstraintStr.find("=")  ) {this->OperatorType=0;  this->OperatorTypeStr="EQ";}
                else if(std::string::npos != TempConstraintStr.find(">=") ) {this->OperatorType=1;  this->OperatorTypeStr="GE";}
                else if(std::string::npos != TempConstraintStr.find("<=") ) {this->OperatorType=-1; this->OperatorTypeStr="LE";}
                else if(std::string::npos != TempConstraintStr.find(">")  ) {this->OperatorType=2;  this->OperatorTypeStr="GT";}
                else if(std::string::npos != TempConstraintStr.find("<")  ) {this->OperatorType=-2; this->OperatorTypeStr="LT";}
                else if(std::string::npos != TempConstraintStr.find("EQ") ) {this->OperatorType=0;  this->OperatorTypeStr="EQ";}
                else if(std::string::npos != TempConstraintStr.find("GE") ) {this->OperatorType=1;  this->OperatorTypeStr="GE";}
                else if(std::string::npos != TempConstraintStr.find("LE") ) {this->OperatorType=-1; this->OperatorTypeStr="LE";}
                else if(std::string::npos != TempConstraintStr.find("GT") ) {this->OperatorType=2;  this->OperatorTypeStr="GT";}
                else if(std::string::npos != TempConstraintStr.find("LT") ) {this->OperatorType=-2; this->OperatorTypeStr="LT";}
                else {
                    std::cout << std::endl;
                    std::cout << "Error! The input constraint operator " << input_args[i+j] << " could not be understood!" << std::endl;
                    std::cout << std::endl;
                    return -1;
                }
            } else {
                std::cout << std::endl;
                std::cout << "Error! The input constraint argument " << input_args[i+j] << " could not be understood!" << std::endl;
                std::cout << std::endl;
                return -1;
            }
            PrevConstraintStr = TempConstraintStr;
        }
        if(verbose>=1) {
            std::cout << std::endl;
            std::cout << "\t";
            //
            if(this->toLIB==-1) {
                std::cout << "SED";
            } else {
                std::cout << "LIB" << this->toLIB;
            }
            std::cout << " ";
            if(this->toPAR==-3) {
                std::cout << "vLv(" << this->toLowerX << "," << this->fromUpperX << ")";
            } else if(this->toPAR==-2) {
                std::cout << "INT(" << this->toLowerX << "," << this->toUpperX << ")";
            } else if(this->toPAR==-1) {
                std::cout << "NORM(" << this->toLowerX << ")";
            } else if(this->toPAR==0) {
                std::cout << "INDEX" << this->toPAR;
            } else {
                std::cout << "PAR" << this->toPAR;
            }
            if(this->toMultiplication!=1.0) {
                // TOOD float equal-comparision is inaccurate
                std::cout << "*" << this->toMultiplication;
            }
            if(this->toAddition!=0.0) {
                // TOOD float equal-comparision is inaccurate
                std::cout << "*" << this->toAddition;
            }
            //
            //std::cout << " ";
            std::cout << " " << this->OperatorTypeStr << " ";
            //std::cout << " ";
            //
            if(this->fromLIB==-1) {
                std::cout << "SED";
            } else {
                std::cout << "LIB" << this->fromLIB;
            }
            std::cout << " ";
            if(this->fromPAR==-3) {
                std::cout << "vLv(" << this->fromLowerX << "," << this->fromUpperX << ")";
            } else if(this->fromPAR==-2) {
                    std::cout << "INT(" << this->fromLowerX << "," << this->fromUpperX << ")";
            } else if(this->fromPAR==-1) {
                std::cout << "NORM(" << this->fromLowerX << ")";
            } else if(this->fromPAR==0) {
                std::cout << "INDEX" << this->fromPAR;
            } else {
                std::cout << "PAR" << this->fromPAR;
            }
            if(this->fromMultiplication!=1.0) {
                // TOOD float equal-comparision is inaccurate
                std::cout << "*" << this->fromMultiplication;
            }
            if(this->fromAddition!=0.0) {
                // TOOD float equal-comparision is inaccurate
                std::cout << "*" << this->fromAddition;
            }
            //
            std::cout << std::endl;
            std::cout << std::endl;
        }
        return 0;
    };
};

#endif
