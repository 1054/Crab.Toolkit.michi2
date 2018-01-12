#include "michi2_Constraint.h"

using namespace std;



michi2Constraint::michi2Constraint(std::vector<std::string> input_arg) {
    this->clear();
    int parse_ok = this->parse(input_arg);
    if(parse_ok<0) { exit (EXIT_FAILURE); } // 2018-01-10 exit if failed to understand the constraint
};



michi2Constraint::michi2Constraint(const char *input_arg1, const char *input_arg2, const char *input_arg3, const char *input_arg4, const char *input_arg5) {
    this->clear();
    std::vector<std::string> input_args;
    input_args.push_back(input_arg1);
    input_args.push_back(input_arg2);
    input_args.push_back(input_arg3);
    input_args.push_back(input_arg4);
    input_args.push_back(input_arg5);
    int parse_ok = this->parse(input_args);
    if(parse_ok<0) { exit (EXIT_FAILURE); } // 2018-01-10 exit if failed to understand the constraint
};



michi2Constraint::~michi2Constraint() {
    this->clear();
};



void michi2Constraint::clear() {
    fromLIB=0; toLIB=0; // if value>0 then it's LIB id; if =-1 then use full SED.
    fromPAR=0; toPAR=0; // if value>0 then it's PAR id; if =-1 then use integration; if=0 then use LIB index.
    fromFactor=1; toFactor=1; // the factor to multiply to from/to data (obsolete).
    fromMathPower=std::nan(""); toMathPower=std::nan(""); // the factor to the power of the from/to data. Added since 2018-01-10.
    fromMultiplication=std::nan(""); toMultiplication=std::nan(""); // the factor to multiply to from/to data.
    fromAddition=std::nan(""); toAddition=std::nan(""); // the value to add to from/to data (after multiplication).
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
std::vector<double> michi2Constraint::extractStringDouble(std::string InputStr) {
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



int michi2Constraint::parse(std::vector<std::string> input_args, int verbose) {
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
        if(j==1) {
            PrevConstraintStr = TempConstraintStr; // PrevConstraintStr is used as the previous string container.
        }
        //
        // if input "lib3" LIB id, then set constraint from/to LIB 3 (number starting from 1).
        // if input "full" then set constraint from/to full SED.
        //
        // if input "index" then set constraint on index.
        // if input "integration" then set constraint on integration.
        // if input "par3" LIB PAR id, then set constraint on LIB PAR.
        //
        if(j==3) {
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
        } else if ( TempConstraintStr=="SED" || (std::string::npos!=TempConstraintStr.find("FULL")) ) {
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
            
        } else if ( (TempConstraintStr=="F") || (std::string::npos!=TempConstraintStr.find("FLUX")) || (std::string::npos!=TempConstraintStr.find("NORM"))) {
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
            
        } else if ( (TempConstraintStr=="I") || (std::string::npos!=TempConstraintStr.find("INDEX")) ) {
            //
            // This means the constraint uses INDEX instead of PARAMETER of the target LIB
            // e.g. LIB2 INDEX EQ LIB3 INDEX
            //
            if(j==5) {this->fromPAR = 0;}
            else if(j==2) {this->toPAR = 0;}
            
        } else if ( (TempConstraintStr=="V") || (std::string::npos!=TempConstraintStr.find("VALUE")) ) {
            //
            // 2018-01-10
            // This means the constraint uses VALUE instead of PARAMETER of the target LIB
            // e.g. LIB2 PAR2 EQ VALUE 3.0
            //
            if(j==4) {this->fromLIB = -2; }
            else if(j==1) {this->toLIB = -2; }
            
        }  else if ( (PrevConstraintStr=="VALUE") && (TempConstraintStr.find_first_of("0123456789.+-")==0) ) {
            //
            // 2018-01-10
            // This means the constraint uses VALUE instead of PARAMETER of the target LIB
            // e.g. LIB2 PAR2 EQ VALUE 3.0
            //
            if(j==5) {this->fromAddition = std::stod(TempConstraintStr); this->fromPAR = -99; }
            else if(j==2) {this->toAddition = std::stod(TempConstraintStr); this->toPAR = -99; }
            
        } else if( (std::string::npos!=TempConstraintStr.find_first_of("0123456789")) ) {
            //
            // it's normal mode, set constraint on some PAR of some LIB
            // e.g. LIB2 PAR2 GE LIB3 PAR2
            //
            size_t TempPosNumber = TempConstraintStr.find_first_of("0123456789");
            std::string TempConstraintStr2 = TempConstraintStr.substr(TempPosNumber);
            if(j==1){this->toLIB = std::stoi(TempConstraintStr2);}
            else if(j==2){this->toPAR = std::stoi(TempConstraintStr2);}
            else if(j==4){this->fromLIB = std::stoi(TempConstraintStr2);}
            else if(j==5){this->fromPAR = std::stoi(TempConstraintStr2);}
            //
            // in addition, if power is given
            // e.g. LIB2 PAR2 GE LIB3 PAR2^-0.25*2.0-50
            // -- NOTE that we must put "PAR2" at the beginning (TODO).
            //
            size_t TempPosMathPower, TempPosMathPowerEnd;
            TempPosMathPower = TempConstraintStr.find_first_of("^",TempPosNumber+1); // try to find MathPower sign
            if(std::string::npos!=TempPosMathPower) {
                TempPosMathPowerEnd = TempConstraintStr.find_first_not_of("0123456789Ee.+-",TempPosMathPower+1);
                if(std::string::npos!=TempPosMathPowerEnd) {
                    std::string TempConstraintStr3 = TempConstraintStr.substr(TempPosMathPower+1,TempPosMathPowerEnd-TempPosMathPower-1+1);
                    if(j==5) {this->fromMathPower = std::stod(TempConstraintStr3); }
                    else if(j==2) {this->toMathPower = std::stod(TempConstraintStr3); }
                }
            }
            //
            // in addition, if multiplication is given
            // e.g. LIB2 PAR2 GE LIB3 PAR2*2.0-50
            // -- NOTE that we can not give something like "2.0*PAR2-50", we must put "PAR2" at the beginning (TODO).
            // -- NOTE that we can not give something like "2.0/PAR2-50", we must put "PAR2" at the beginning (TODO).
            //
            size_t TempPosMultiplication, TempPosMultiplicationEnd;
            TempPosMultiplication = TempConstraintStr.find_first_of("*",TempPosNumber+1); // try to find multiplication sign
            if(std::string::npos!=TempPosMultiplication) {
                TempPosMultiplicationEnd = TempConstraintStr.find_first_not_of("0123456789Ee.+-",TempPosMultiplication+1);
                if(std::string::npos!=TempPosMultiplicationEnd) {
                    std::string TempConstraintStr3 = TempConstraintStr.substr(TempPosMultiplication+1,TempPosMultiplicationEnd-TempPosMultiplication-1+1);
                    if(j==5) {this->fromMultiplication = std::stod(TempConstraintStr3); }
                    else if(j==2) {this->toMultiplication = std::stod(TempConstraintStr3); }
                }
            } else {
                TempPosMultiplication = TempConstraintStr.find_first_of("/",TempPosNumber+1); // otherwise, try to find division sign
                if(std::string::npos!=TempPosMultiplication) {
                    TempPosMultiplicationEnd = TempConstraintStr.find_first_not_of("0123456789Ee.+-",TempPosMultiplication+1);
                    if(std::string::npos!=TempPosMultiplicationEnd) {
                        std::string TempConstraintStr3 = TempConstraintStr.substr(TempPosMultiplication+1,TempPosMultiplicationEnd-TempPosMultiplication-1+1);
                        if(j==5) {this->fromMultiplication = 1.0/std::stod(TempConstraintStr3); }
                        else if(j==2) {this->toMultiplication = 1.0/std::stod(TempConstraintStr3); }
                    }
                }
            }
            //
            // in addition, if addition is given
            // e.g. LIB2 PAR2 GE LIB3 PAR2/2.0-50
            // -- NOTE that we can not give something like "-50+PAR2/2.0", we must put "PAR2" at the beginning.
            // -- NOTE that we can not give something like "PAR2/2.0-100+50", we can only set one value "-50".
            //
            size_t TempPosAddition, TempPosAdditionEnd;
            TempPosAddition = TempConstraintStr.find_first_of("+",TempPosNumber+1); // try to find plus sign
            if(std::string::npos!=TempPosAddition) {
                TempPosAdditionEnd = TempConstraintStr.find_first_not_of("0123456789Ee.+-",TempPosAddition+1);
                if(std::string::npos!=TempPosAdditionEnd) {
                    std::string TempConstraintStr3 = TempConstraintStr.substr(TempPosAddition+1,TempPosAdditionEnd-TempPosAddition-1+1);
                    if(j==5) {this->fromAddition = std::stod(TempConstraintStr3); }
                    else if(j==2) {this->toAddition = std::stod(TempConstraintStr3); }
                }
            } else {
                TempPosAddition = TempConstraintStr.find_first_of("-",TempPosNumber+1); // otherwise, try to find minus sign
                if(std::string::npos!=TempPosAddition) {
                    TempPosAdditionEnd = TempConstraintStr.find_first_not_of("0123456789Ee.+-",TempPosAddition+1);
                    if(std::string::npos!=TempPosAdditionEnd) {
                        std::string TempConstraintStr3 = TempConstraintStr.substr(TempPosAddition+1,TempPosAdditionEnd-TempPosAddition-1+1);
                        if(j==5) {this->fromAddition = -std::stod(TempConstraintStr3); }
                        else if(j==2) {this->toAddition = -std::stod(TempConstraintStr3); }
                    }
                }
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
        if(this->toLIB==-2) {
            std::cout << "VALUE"; // 2018-01-10 allow to set constant value
        } else if(this->toLIB==-1) {
            std::cout << "SED";
        } else {
            std::cout << "LIB" << this->toLIB;
        }
        std::cout << " ";
        if(this->toPAR==-99) {
            std::cout << ""; // 2018-01-10 allow to set constant value
        } else if(this->toPAR==-3) {
            std::cout << "vLv(" << this->toLowerX << "," << this->fromUpperX << ")";
        } else if(this->toPAR==-2) {
            std::cout << "INT(" << this->toLowerX << "," << this->toUpperX << ")";
        } else if(this->toPAR==-1) {
            //std::cout << "NORM(" << this->toLowerX << ")"; // TODO currently we can only set constraint to the normalization of the templates, i.e., aCOE in michi2_MinPack.cpp, but could not set to the normalization at a specific X position. For example, if we have a SED template normalized to total integrated flux density of 1.0 and it is LIB5, and we set a constraint of "-constraint LIB5 NORM EQ SED(8,1000)", then it means the total integrated flux density of LIB5 will be always fixed to the integration of the full SED from LowerX=8 to UpperX=1000. We can not yet set something like "-constraint LIB5 NORM(21500) EQ SED(8,1000)" to make the LIB5 Y value at X=21500 fixed to SED(8,1000).
            std::cout << "NORM";
        } else if(this->toPAR==0) {
            std::cout << "INDEX";
        } else {
            std::cout << "PAR" << this->toPAR;
        }
        if(!std::isnan(this->toMathPower)) {
            std::cout << "^" << this->toMathPower;
        }
        if(!std::isnan(this->toMultiplication)) {
            std::cout << "*" << this->toMultiplication;
        }
        if(!std::isnan(this->toAddition)) {
            std::cout << std::showpos << this->toAddition << std::noshowpos;
        }
        //
        //std::cout << " ";
        std::cout << " " << this->OperatorTypeStr << " ";
        //std::cout << " ";
        //
        if(this->fromLIB==-2) {
            std::cout << "VALUE"; // 2018-01-10 allow to set constant value
        } else if(this->fromLIB==-1) {
            std::cout << "SED";
        } else {
            std::cout << "LIB" << this->fromLIB;
        }
        std::cout << " ";
        if(this->fromPAR==-99) {
            std::cout << ""; // 2018-01-10 allow to set constant value
        } else if(this->fromPAR==-3) {
            std::cout << "vLv(" << this->fromLowerX << "," << this->fromUpperX << ")";
        } else if(this->fromPAR==-2) {
            std::cout << "INT(" << this->fromLowerX << "," << this->fromUpperX << ")";
        } else if(this->fromPAR==-1) {
            //std::cout << "NORM(" << this->fromLowerX << ")"; // see the note above at "else if(this->toPAR==-1)"
            std::cout << "NORM";
        } else if(this->fromPAR==0) {
            std::cout << "INDEX";
        } else {
            std::cout << "PAR" << this->fromPAR;
        }
        if(!std::isnan(this->fromMathPower)) {
            std::cout << "^" << this->fromMathPower;
        }
        if(!std::isnan(this->fromMultiplication)) {
            std::cout << "*" << this->fromMultiplication;
        }
        if(!std::isnan(this->fromAddition)) {
            std::cout << std::showpos << this->fromAddition << std::noshowpos;
        }
        //
        std::cout << std::endl;
        std::cout << std::endl;
    }
    return 0;
};




bool michi2Constraint::check(struct mnchi2parallelParams * pParams, int debug) {
    std::vector<michi2DataClass *> SDLIBS;
    return this->check(pParams,SDLIBS,debug);
}

bool michi2Constraint::check(std::vector<michi2DataClass *> SDLIBS, int debug) {
    return this->check(NULL,SDLIBS,debug);
}

bool michi2Constraint::check(struct mnchi2parallelParams * pParams, std::vector<michi2DataClass *> SDLIBS, int debug) {
    bool ConstraintOK = true;
    michi2Constraint *TempConstraint = this;
    if( SDLIBS.size()>0 && TempConstraint->fromLIB>0 && TempConstraint->toLIB>0 && TempConstraint->fromPAR>0 && TempConstraint->toPAR>0 ) {
        // set constraint: toLIB.toPAR \Operator\ fromLIB.fromPAR
        //           e.g., LIB4  PAR3   GT        LIB3    PAR3
        //
        // 2018-01-10: now allow simple ^ */ +- operations on PAR
        //           e.g., LIB4  PAR3^-0.5*2.0+50  GT  LIB3  PAR3^-0.5-50  -- no bracket please
        //
        double TempPAR1 = SDLIBS[TempConstraint->toLIB-1]->FPAR[TempConstraint->toPAR-1][0];
        double TempPAR2 = SDLIBS[TempConstraint->fromLIB-1]->FPAR[TempConstraint->fromPAR-1][0];
        if(debug>=1) {
            std::cout << "michi2Constraint::check: Setting constraint: LIB" << TempConstraint->toLIB << " PAR" << TempConstraint->toPAR << " (" << TempPAR1 << ") " << TempConstraint->OperatorTypeStr << " LIB" << TempConstraint->fromLIB << " PAR" << TempConstraint->toPAR << " (" << TempPAR2 <<  ")";
        }
        if(!std::isnan(TempConstraint->toMathPower)) { TempPAR1 = pow(TempPAR1,TempConstraint->toMathPower); }
        if(!std::isnan(TempConstraint->fromMathPower)) { TempPAR2 = pow(TempPAR2,TempConstraint->fromMathPower); }
        if(!std::isnan(TempConstraint->toMultiplication)) { TempPAR1 = TempPAR1 * (TempConstraint->toMultiplication); }
        if(!std::isnan(TempConstraint->fromMultiplication)) { TempPAR2 = TempPAR2 * (TempConstraint->fromMultiplication); }
        if(!std::isnan(TempConstraint->toAddition)) { TempPAR1 = TempPAR1 + (TempConstraint->toAddition); }
        if(!std::isnan(TempConstraint->fromAddition)) { TempPAR2 = TempPAR2 + (TempConstraint->fromAddition); }
        if(TempConstraint->OperatorType == 0) {
            ConstraintOK = (TempPAR1 == TempPAR2); // set constraint: toLIB.toPAR EQ fromLIB.fromPAR
        } else if(TempConstraint->OperatorType == 1) {
            ConstraintOK = (TempPAR1 >= TempPAR2); // set constraint: toLIB.toPAR GE fromLIB.fromPAR
        } else if(TempConstraint->OperatorType == -1) {
            ConstraintOK = (TempPAR1 <= TempPAR2); // set constraint: toLIB.toPAR LE fromLIB.fromPAR
        } else if(TempConstraint->OperatorType == 2) {
            ConstraintOK = (TempPAR1 > TempPAR2); // set constraint: toLIB.toPAR GT fromLIB.fromPAR
        } else if(TempConstraint->OperatorType == -2) {
            ConstraintOK = (TempPAR1 < TempPAR2); // set constraint: toLIB.toPAR LT fromLIB.fromPAR
        }
        if(debug>=1) {
            std::cout << " OK? " << ConstraintOK << std::endl;
        }
    } else if( pParams && TempConstraint->fromLIB>0 && TempConstraint->toLIB>0 && TempConstraint->fromPAR==0 && TempConstraint->toPAR==0 ) {
        // set constraint: toLIB.toINDEX \Operator\ toINDEX.fromINDEX
        //           e.g., LIB4  INDEX    EQ        LIB3    INDEX
        //
        // 2018-01-10: do not allow any ^ */ +- operation on INDEX
        //
        double TempINDEX1 = pParams->idLIBList.at(TempConstraint->toLIB-1);
        double TempINDEX2 = pParams->idLIBList.at(TempConstraint->fromLIB-1);
        if(debug>=1) {
            std::cout << "michi2Constraint::check: Setting constraint: LIB" << TempConstraint->toLIB << " INDEX (" << TempINDEX1 << ") EQ LIB" << TempConstraint->fromLIB << " INDEX (" << TempINDEX2 <<  ")";
        }
        if(TempConstraint->OperatorType == 0) {
            ConstraintOK = (TempINDEX1 == TempINDEX2); // set constraint: toLIB.toINDEX EQ fromLIB.fromINDEX
        } else if(TempConstraint->OperatorType == 1) {
            ConstraintOK = (TempINDEX1 >= TempINDEX2); // set constraint: toLIB.toINDEX GE fromLIB.fromINDEX
        } else if(TempConstraint->OperatorType == -1) {
            ConstraintOK = (TempINDEX1 <= TempINDEX2); // set constraint: toLIB.toINDEX LE fromLIB.fromINDEX
        } else if(TempConstraint->OperatorType == 2) {
            ConstraintOK = (TempINDEX1 > TempINDEX2); // set constraint: toLIB.toINDEX GT fromLIB.fromINDEX
        } else if(TempConstraint->OperatorType == -2) {
            ConstraintOK = (TempINDEX1 < TempINDEX2); // set constraint: toLIB.toINDEX LT fromLIB.fromINDEX
        }
        if(debug>=1) {
            std::cout << " OK? " << ConstraintOK << std::endl;
        }
    } else if( SDLIBS.size()>0 && TempConstraint->fromLIB==-2 && TempConstraint->toLIB>0 && TempConstraint->fromPAR==-99 && TempConstraint->toPAR>0 ) {
        // 2018-01-10
        // set constraint: toLIB.toPAR \Operator\ VALUE.fromAdd
        //           e.g., LIB4  PAR3   GT        VALUE 0.5
        //           e.g., LIB4  PAR3   LT        VALUE 1.0
        //
        // 2018-01-10: if the fromLIB==-2 and fromPAR==-99, then constrain toLIB.toPAR by the input constant value TempConstraint->fromAddition.
        //
        double TempPAR1 = SDLIBS[TempConstraint->toLIB-1]->FPAR[TempConstraint->toPAR-1][0];
        double TempVALUE2 = TempConstraint->fromAddition;
        if(debug>=1) {
            std::cout << "michi2Constraint::check: Setting constraint: LIB" << TempConstraint->toLIB << " PAR" << TempConstraint->toPAR << " (" << TempPAR1 << ") " << TempConstraint->OperatorTypeStr << " VALUE (" << TempVALUE2 <<  ")";
        }
        if(TempConstraint->OperatorType == 0) {
            ConstraintOK = (TempPAR1 == TempVALUE2); // set constraint: toLIB.toPAR EQ VALUE.fromAddition
        } else if(TempConstraint->OperatorType == 1) {
            ConstraintOK = (TempPAR1 >= TempVALUE2); // set constraint: toLIB.toPAR GE VALUE.fromAddition
        } else if(TempConstraint->OperatorType == -1) {
            ConstraintOK = (TempPAR1 <= TempVALUE2); // set constraint: toLIB.toPAR LE VALUE.fromAddition
        } else if(TempConstraint->OperatorType == 2) {
            ConstraintOK = (TempPAR1 > TempVALUE2); // set constraint: toLIB.toPAR GT VALUE.fromAddition
        } else if(TempConstraint->OperatorType == -2) {
            ConstraintOK = (TempPAR1 < TempVALUE2); // set constraint: toLIB.toPAR LT VALUE.fromAddition
        }
        if(debug>=1) {
            std::cout << " OK? " << ConstraintOK << std::endl;
        }
    } else {
        // set constraint: toLIB[LowerX:UpperX]*Factor = toLIB[LowerX:UpperX]*Factor
        // <TODO><20160719><dzliu>
        // ITS TOO HARD TO IMPLEMENT SO MANY THINGS
        // Constraining the normalization of a LIB is not what should be constrained here.
        // This code should only be a tool to compute minimum chi2 from OBS data with multiple LIB data.
        // We can set constraints on LIB parameters, LIB index, but not LIB normalization!
    }
    return ConstraintOK;
}













