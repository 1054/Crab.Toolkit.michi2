/*
 compile:
     g++ -pthread main.cpp -o michi2_v02_test
     g++ -pthread main.cpp -o michi2_v02_201507
     g++ -pthread main.cpp -o michi2_v02_201509
     g++ -pthread main.cpp -o michi2_v02_201509_DL07
     clang++ -std=c++11 -pthread main.cpp -o michi2_v04_201607
     clang++ -std=c++11 -pthread main.cpp -o michi2_v04_20160727
 deploy:
     \cp a.out ../../06_LVG_Synthesis/03_Minimize_Chi2/michi2_v04
 test: 
     ./michi2_v04 -obs flux-co-ngc1068.dat -lib lib.lvg lib.lvg lib.lvg -out out.out
     ./michi2_v04_201607  -obs test/flux-co-ngc1068.dat  -lib test/lib_big_co.lvg test/lib_big_co.lvg   -out test/out.out   -constraint "LIB1" "PAR3" "<" LIB2 par2
     ./michi2_v04_201607  -obs test/flux-co-ngc1068.dat  -lib test/lib_big_co.lvg test/lib_big_co.lvg   -out test/out.out   -constraint LIB1 INDEX EQ LIB2 INDEX
*/

#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>    // std::transform
#include <string>
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw
#include <time.h>
#include "michi2_Constraint.h"
#include "michi2_v04.cpp"

using namespace std;




int main(int argc, char **argv)
{
    /**
    // You’ll notice that argv[0] is the path and name of the program itself.
    std::cout << "CrabTableReadInfo: Input " << argc-1 << " arguments." << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << argv[i] << std::endl;
    }
    **/
    if(argc>1)
    {
        // michi2(argv[1],argv[2]);
        // michi2("flux-co.dat","lib.lvg");
        // m2chi2("flux-co.dat","lib.lvg","lib.lvg");
        // m2chi2("flux-co.dat","flower.lvg","flower.lvg");
        std::vector<std::string> FileObsList; // Obs Flux & FlueErr List
        std::vector<std::string> FileLibList; // Lib
        std::vector<std::string> FileOutList; // Out
        std::vector<std::string> InputFilterCurveList; // Obs Filter Curve List
        Constraints.clear();
        std::cout << "Welcome!" << std::endl << std::endl;
        int i=0;
        while(i<argc) {
            if(0==strncmp(argv[i],"-obs",4))
            {
                i++;
                while(i<argc) {
                    if(0!=strncmp(argv[i],"-",1)) {
                        FileObsList.push_back(argv[i]);
                    } else {break;}
                    i++;
                }
            }
            else if(0==strncmp(argv[i],"-filter",7))
            {
                i++;
                while(i<argc) {
                    if(0!=strncmp(argv[i],"-",1)) {
                        InputFilterCurveList.push_back(argv[i]);
                    } else {break;}
                    i++;
                }
            }
            else if(0==strncmp(argv[i],"-lib",4))
            {
                i++;
                while(i<argc) {
                    if(0!=strncmp(argv[i],"-",1)) {
                        FileLibList.push_back(argv[i]);
                    } else {break;}
                    i++;
                }
            }
            else if(0==strncmp(argv[i],"-out",4))
            {
                i++;
                while(i<argc) {
                    if(0!=strncmp(argv[i],"-",1)) {
                        FileOutList.push_back(argv[i]);
                    } else {break;}
                    i++;
                }
            }
            else if(0==strncmp(argv[i],"-redshift",9) || 0==strcmp(argv[i],"-z"))
            {
                i++;
                if(i<argc) {
                    if(strspn(argv[i],"+-.eE0123456789")==strlen(argv[i])) {
                        InfoRedshift = argv[i]; i++;
                    } else {
                        std::cout << "Error! Input redshift is invalid!" << std::endl;
                    }
                }
            }
            else if(0==strncmp(argv[i],"-parallel",9) || 0==strcmp(argv[i],"-p"))
            {
                i++;
                if(i<argc) {
                    if(strspn(argv[i],"0123456789")==strlen(argv[i])) {
                        NumbParallel = atoi(argv[i]); i++;
                    } else {
                        std::cout << "Error! Input parallel is invalid!" << std::endl;
                    }
                }
            }
            else if(0==strncmp(argv[i],"-constrain",10) || 0==strncmp(argv[i],"-con",4) || 0==strcmp(argv[i],"-s"))
            {
                // i++;
                // a valid constraint requires 5 arguments:
                // LIB1 PAR1 OPERATOR LIB2 PAR2
                // for example -constraint LIB3-DL07 PAR1-UMIN EQ LIB4-DL07 PAR1-UMIN
                if((i+5)<argc) {
                    michi2Constraint *TempConstraint = new michi2Constraint(argv[i+1], argv[i+2], argv[i+3], argv[i+4], argv[i+5]);
                    Constraints.push_back(TempConstraint);
                    i+=5;
                } else {
                    std::cout << std::endl;
                    std::cout << "Error! The input constraint should have 5 arguments, e.g. lib3 par3 < lib2 par3." << std::endl;
                    std::cout << "Examples:" << std::endl;
                    std::cout << "    -constraint lib2 par3 eq lib3 par3 # e.g. DL07 cold+warm SED, cold Umin should be equal to warm Umin." << std::endl;
                    std::cout << "    -constraint lib1 par1 lt lib2 par1 # e.g. LVG two-component CO SLED, cold T_kin should be less than warm T_kin." << std::endl;
                    std::cout << "Notes:" << std::endl;
                    std::cout << "    We can set several constraints for example" << std::endl;
                    std::cout << "    -constraint lib1 par1 le lib2 par1 -constraint lib1 par2 le lib2 par2" << std::endl;
                    std::cout << "    This is for LVG two-component CO SLED, both T_kin and n_H2 are set to be smaller in colder component." << std::endl;
                    // <TODO><20160719><dzliu> std::cout << "    -constraint \"lib4[200000]\" \"flux\" eq \"full[8:1000]\" \"integration*0.0005\" # e.g. set radio 200000um flux be the 8-1000 integration times 0.0005." << std::endl;
                    // <TODO><20160719><dzliu>
                    // ITS TOO HARD TO IMPLEMENT SO MANY THINGS
                    // See notes in michi2_v04.cpp
                    //
                    std::cout << std::endl;
                    return -1;
                }
                i++;
            }
            else
            {
                i++;
            }
        }
        std::cout << "Please make sure that you already have these files under current directory:" << std::endl;
        for(int iObs = 0; iObs < FileObsList.size(); iObs++) {
            std::cout << "\t" << FileObsList[iObs];
        }
        std::cout << std::endl;
        for(int iLib = 0; iLib < FileLibList.size(); iLib++) {
            std::cout << "\t" << FileLibList[iLib];
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "We will output to:" << std::endl;
        for(int iOut = 0; iOut < FileOutList.size(); iOut++) {
            std::cout << "\t" << FileOutList[iOut];
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "OK, let's go! Begin at " << currentDateTime() << std::endl;
        std::cout << std::endl;
        
        // if(0==strlen(FileLi4)) {
        //     if(0==strlen(FileLi3)) {
        //         if(0==strlen(FileLi2)) {
        //             michi2(FileObs,FileLib,FileOut);
        //         } else {
        //             m2chi2(FileObs,FileLib,FileLi2,FileOut);
        //         }
        //     } else {
        //         m3chi2(FileObs,FileLib,FileLi2,FileLi3,FileOut);
        //     }
        // } else {
        //     m4chi2(FileObs,FileLib,FileLi2,FileLi3,FileLi4,FileOut);
        // }

        //<TODO><UNCOMMNET>//
        mnchi2(FileObsList,FileLibList,FileOutList,InputFilterCurveList);
        //<TODO><UNCOMMNET>//
        
//        if(argc>=2) {
//            michi2("flux-co.dat","lib.lvg",argv[1]);
//            // <TODO> // m2chi2("flux-co.dat","lib.lvg","lib.lvg",argv[1]);
//        } else {
//            michi2("flux-co.dat","lib.lvg");
//            // <TODO> // m2chi2("flux-co.dat","lib.lvg","lib.lvg");
//        }
        
        std::cout << std::endl;
        std::cout << "Finally! End at " << currentDateTime() << std::endl;
        std::cout << std::endl;
    } else
    {
        
        std::cout << std::endl;
        std::cout << "Usage: \n";
        std::cout << "       michi2_v04 -obs flux-co.dat -lib lib.lvg -out output.csv\n";
        std::cout << "       \n";
        std::cout << "       michi2_v04 -obs flux-co.dat -lib lib.lvg lib.lvg -out output.csv\n";
        std::cout << "       \n";
        std::cout << "       michi2_v04 -obs flux-co.dat -lib Star.SED DL07Lo.SED DL07Hi.SED -out output.csv\n";
        std::cout << "       \n";
        // std::cout << "Version: \n\t michi2_v04 " << "2014-08-22 Orme des Merisiers" << " copyleft " << std::endl;
        // std::cout << "Version: \n\t michi2_v04 " << "2015-04-09 Orme des Merisiers" << std::endl;
        // std::cout << "Version: \n\t michi2_v04 " << "2016-07-14 Nanjing" << std::endl;
        std::cout << "       michi2_v04 -obs flux-obsframe.dat \\\n";
        std::cout << "                  -redshift 6.3 \\\n";
        std::cout << "                  -lib Star.SED AGN.SED DL07Hi.SED DL07Lo.SED Radio.SED \\\n";
        std::cout << "                  -out output.dat \\\n";
        std::cout << "                  -constrain LIB3 INDEX EQ LIB4 INDEX \\\n";
        std::cout << "                  -constrain LIB5 NORM EQ SED \"vLv(8,1000)*1.061619121e-06\" \\\n";
        std::cout << "                  -filter filter.list\n";
        std::cout << "                  # (*) The first constraint means that we lock the index of DL07Hi.SED and\n";
        std::cout << "                  #     DL07Lo.SED to be the same, i.e. same Umin, qPAH, etc.\n";
        std::cout << "                  # (*) The second constraint means that we lock the normalization of Radio.SED\n";
        std::cout << "                  #     to total IR luminosity vLv(8,1000) following the IR-radio correlation,\n";
        std::cout << "                  #     i.e. with q_IR = 2.4, S_(1.4GHz)/S_(IR,8-1000) = 1.0/3750/10^2.4 = 1.061619121e-06\n";
        std::cout << "                  # (*) We can also read filter curves according to the input \"filter.list\" file,\n";
        std::cout << "                  #     which should have two columns, wavelength and filter curve file path without white space,\n";
        std::cout << "                  #     and rows should exactly correspond to \"flux-obsframe.dat\".\n";
        std::cout << "                  #     For wavelength without applicable filter curve, just put \"none\" in the second column.\n";
        std::cout << "                  #     The filter curve file path should refer to a two-column ASCII file which contains\n";
        std::cout << "                  #     obs-frame wavelength and filter transmission value normalized to 1.\n";
        std::cout << "       \n";
        std::cout << "Version: \n";
        std::cout << "         michi2_v04 " << "2017-10-01 Heidelberg" << std::endl;
        std::cout << std::endl;
        
        /*
         
        std::cout << "michi2" << std::endl;
        std::cout << "  aim: this small code computes the chi-square between observed dataset and library dataset. Datasets must contain two independent variables, the first Var1 is used to match the obs and lib, the second Var2 is used to calculated chi-square. " << std::endl;
        std::cout << "  use: michi2 -obs \"flux-co.dat\" -lib \"flower.lvg\"" << std::endl;
        std::cout << "  use: m2chi2 -obs \"flux-co.dat\" -lib \"flower.lvg\" \"flower.lvg\"" << std::endl;
        std::cout << std::endl;
         
         */
    }
    return 0;
}

