#ifndef H_CrabTableReadInfo
#define H_CrabTableReadInfo
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "CrabStringClean.cpp"
#include "CrabStringUnicode.cpp"
#include "CrabStringFindWholeWord.cpp"

using namespace std;



std::string CrabTableReadInfo(const char *InputFile, const char *KeyName, const char *CommentMark = "#");



/*

std::wstring ctristrtrim(const wstring &t, int trimflag = 0);

std::wstring ctristrtrim(const wstring &t, int trimflag)
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

*/



std::string CrabTableReadInfo(const char *InputFile, const char *InputKeyName, const char *CommentMark)
{
    // http://www.cplusplus.com/forum/beginner/70805/
    std::string OutputText;
    int debug = 0;
    if(strlen(InputKeyName)>0 && strlen(InputFile)>0) {
        std::wstring line;
        std::wstring keyname = CrabStringWideString(InputKeyName); // From utf8 char * to utf16 wstring
        // std::wstring commark(CommentMark,CommentMark+strlen(CommentMark)); // From char * to wstring: http://stackoverflow.com/questions/2573834/c-convert-string-or-char-to-wstring-or-wchar-t
        std::wstring commark = CrabStringWideString(CommentMark); // From utf8 char * to utf16 wstring
        std::wstring keyvalue;
        std::wstring keycomment;
        
        if(debug) { std::cout << "DEBUG:InputFile==" << InputFile << std::endl; }
        if(debug) { std::cout << "DEBUG:InputKeyName==" << InputKeyName << std::endl; }
        if(debug) { std::cout << "DEBUG:CommentMark==" << CommentMark << std::endl; }
        
        std::ifstream backstory (InputFile); std::string backline; // Do not use wstring to read file! <TODO>
        if (backstory.is_open()) {
            while (backstory.good()) {
                getline(backstory,backline);
                // if(debug) { std::wcout << "DEBUG:" << backline.length() << std::endl; }
                line = CrabStringWideString(backline); // From utf8 char * to utf16 wstring
                // if(debug) { std::wcout << "DEBUG:" << line.length() << std::endl; }
                line = CrabStringTrim(line);
                // if(debug) { std::wcout << "DEBUG:" << line.length() << std::endl; }
                // <Corrected><20150201><DzLIU> Read Galfit COMP_11 COMP_1 Problem
                size_t foundpos1 = CrabStringFindWholeWord(line,keyname,L" ");
                size_t foundpos2 = CrabStringFindWholeWord(line,keyname,L"=");
                size_t foundpos = string::npos;
                if(foundpos1!=string::npos) {foundpos = foundpos1;}
                if(foundpos2!=string::npos) {foundpos = foundpos2;}
                if(foundpos!=string::npos) {
                    if(debug) { std::cout << "DEBUG:"; csu_print_wstr(line); std::wcout << std::endl; }
                    // size_t foundpos = line.find_first_not_of(keyname); // find keyname as a whole word // <Corrected><20150201><DzLIU> Read Galfit COMP_11 COMP_1 Problem
                    // wstring trimline = line;                           // find keyname as a whole word // <Corrected><20150201><DzLIU> Read Galfit COMP_11 COMP_1 Problem
                    // trimline.erase(0,foundpos);                        // find keyname as a whole word // <Corrected><20150201><DzLIU> Read Galfit COMP_11 COMP_1 Problem
                                                                          // <TODO><20150201><DzLIU> Still can not deal with this: "Happy Birthday = 2010" but keyname is "Happy".
                    wstring trimline = line;
                    trimline.erase(0,foundpos+keyname.length());
                    if(debug) { std::cout << "DEBUG:"; csu_print_wstr(trimline); std::wcout << std::endl; } // <Debug><20150201><DzLIU> Read Galfit COMP_11 COMP_1 Problem
                    trimline = CrabStringTrim(trimline);
                    if(debug) { std::cout << "DEBUG:"; csu_print_wstr(trimline); std::wcout << std::endl; } // <Debug><20150201><DzLIU> Read Galfit COMP_11 COMP_1 Problem
                    if(0==trimline.find(L"=")) { // find keyname as a whole word
                        trimline.erase(0,1);
                        foundpos = trimline.find(commark); // find comment
                        if(string::npos!=foundpos) { // no need to say 0<=foundpos because foundpos is unsigned.
                            keycomment = trimline;
                            keycomment.erase(0,foundpos);
                            trimline.erase(foundpos);
                        }
                        keyvalue = CrabStringTrim(trimline);
                        OutputText = csu_wstr_to_utf8(keyvalue);
                        // [show full information] std::cout << keyname << " = " << keyvalue << std::endl;
                        // std::wcout << keyvalue << std::endl; // <TODO> Verbose
                    }
                    if(debug) { std::cout << "DEBUG:ReadKeyValue==" << OutputText << std::endl; }
                }
            }
            backstory.close();
        }
        else {
            std::cout << "CrabTableReadInfo: Unable to open file! Please check " << InputFile << std::endl;
        }
    }
    return OutputText;
}


#endif