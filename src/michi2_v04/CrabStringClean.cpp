/*
 
 We must be careful about codecs. 
 All files input should be utf8. 
 
 */

#ifndef H_CrabStringClean
#define H_CrabStringClean
#include <iostream>
#include <string>
#include <wchar.h>

using namespace std;

std::string CrabStringTrim(const string &t, int trimflag = 0);

std::wstring CrabStringTrim(const wstring &t, int trimflag = 0);



std::string CrabStringTrim(const string &t, int trimflag)
{
    std::string str = t;
    std::string strspace = " "; // trim space
    size_t found;
    // trim leading space
    if(0==trimflag || 2==trimflag) {
        found = str.find_first_not_of(strspace);
        if (found != string::npos) {
            if(found>0) {
                // std::cout << found << std::endl;
                str.erase(0,found);
            }
        }
        else
            str.clear(); // str is all whitespace
    }
    // trim trailing space
    if(0==trimflag || 2==trimflag) {
        found = str.find_last_not_of(strspace);
        if (found != string::npos) {
            if(found<str.length()-1) {
                // std::cout << found << std::endl;
                str.erase(found+1);
            }
        }
        else
            str.clear(); // str is all whitespace
    }
    return str;
}

std::wstring CrabStringTrim(const wstring &t, int trimflag)
{
    std::wstring str = t;
    std::wstring ws = L" "; // trim space
    size_t found;
    // trim leading space
    if(0==trimflag || 2==trimflag) {
        found = str.find_first_not_of(ws);
        if (found != string::npos) {
            if(found>0) {
                // std::cout << found << std::endl;
                str.erase(0,found);
            }
        }
        else
            str.clear(); // str is all whitespace
    }
    // trim trailing space
    if(0==trimflag || 2==trimflag) {
        found = str.find_last_not_of(ws);
        if (found != string::npos) {
            if(found<str.length()-1) {
                // std::cout << found << std::endl;
                str.erase(found+1);
            }
        }
        else
            str.clear(); // str is all whitespace
    }
    return str;
}

#endif