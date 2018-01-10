#ifndef H_CrabTableGetLineCount
#define H_CrabTableGetLineCount
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;



string ctglcstrtrim(const string &t, int trimflag = 0);

long CrabTableGetLineCount(const char *InputFile, const char * CommentMark = "#");





string ctglcstrtrim(const string &t, int trimflag)
{
    string str = t;
    string strspace = " "; // trim space
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



long CrabTableGetLineCount(const char *InputFile, const char * CommentMark)
{
    //
    long linco = 0;
    //
    if(strlen(InputFile)>0) {
        std::string line;
        std::string trline;
        std::string commark(CommentMark); // From char * to string: http://stackoverflow.com/questions/2573834/c-convert-string-or-char-to-string-or-wchar-t
        std::ifstream backstory (InputFile);
        if (backstory.is_open()) {
            while (backstory.good()) {
                // get this line
                getline(backstory,line);
                // get it trimmed
                trline = ctrcstrtrim(line);
                // whether this is empty line?
                if(trline.empty()) {
                    if(linco>0) {break;} else {continue;}
                }
                // whether this is comment line?
                if(0==trline.find(commark)) {
                    continue;
                }
                // whether this is content line?
                linco++;
            }
            backstory.close();
        }
        else {
            std::cout << "CrabTableGetLineCount: InputFile is invalid! Please check " << InputFile << std::endl;
        }
    }
    return linco;
}



#endif