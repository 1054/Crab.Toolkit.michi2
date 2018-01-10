/*
 
 last update:
 
 2015-06-15 add argument BreakupMark in functions
 
 */

#ifndef H_CrabStringReadColumn
#define H_CrabStringReadColumn
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;



string csrcstrtrim(const string &t, int trimflag = 0);

long csrcstrfindwholeword(const string &t, const string &strsearch, const string &BreakupMark = " ");

string csrcstrfindwholeword(const string &t, size_t pos, const string &BreakupMark = " ");

std::vector<std::string> CrabStringReadColumn(std::vector<std::string> InputText, const char *InputColHead, const char *CommentMark = "#", const char *BreakupMark = " ");

std::vector<std::string> CrabStringReadColumn(std::vector<std::string> InputText, int InputColNumber, const char *CommentMark = "#", const char *BreakupMark = " ");






string csrcstrtrim(const string &t, int trimflag)
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


long csrcstrfindwholeword(const string &t, const string &strsearch, const string &BreakupMark)
{
    size_t foundpos = 0;
    // std::cout << strsearch << " " << strsearch.length() << " " << t.find(strsearch) << std::endl;
    while(string::npos!=foundpos) {
        foundpos = t.find(strsearch,foundpos);
        // std::cout << t << "\n" << strsearch << "\n" << foundpos << t.substr(foundpos-1,BreakupMark.length()) << t.substr(foundpos+strsearch.length(),BreakupMark.length()) << std::endl;
        // std::cout << (BreakupMark==t.substr(foundpos-1,BreakupMark.length())) << std::endl;
        // std::cout << (BreakupMark==t.substr(foundpos+strsearch.length(),BreakupMark.length())) << std::endl;
        // sleep(1.0);
        // find keyname as a whole word
        if(string::npos!=foundpos) {
            if(BreakupMark==t.substr(foundpos-1,BreakupMark.length()) &&
               BreakupMark==t.substr(foundpos+strsearch.length(),BreakupMark.length())) { break; }
        } else {
            break;
        }
        foundpos++;
    }
    return foundpos;
}


string csrcstrfindwholeword(const string &t, size_t pos, const string &BreakupMark)
{
    string foundstr;
    if(0==pos) {
        size_t foundposL = 0;
        size_t foundposR = t.find(BreakupMark)-1;
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    } else if(t.length()-1==pos) {
        size_t foundposL = t.rfind(BreakupMark)+1;
        size_t foundposR = t.length()-1;
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    } else if(0<pos && t.length()-1>pos) {
        string strL = t.substr(0,pos);
        string strR = t.substr(pos);
        size_t foundposL = strL.rfind(BreakupMark)+1;
        size_t foundposR = strR.find(BreakupMark)-1+pos;
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    }
    return foundstr;
}






std::vector<std::string> CrabStringReadColumn(std::vector<std::string> InputText, const char *InputColHead, const char * CommentMark, const char *BreakupMark)
{
    //
    std::vector<std::string> OutputArray;
    // OutputArray.push_back(superman);
    //
    if(strlen(InputColHead)>0 && InputText.size()>0) {
        std::string line;
        std::string trline;
        std::string colname(InputColHead); // From char * to string: http://stackoverflow.com/questions/2573834/c-convert-string-or-char-to-string-or-wchar-t
        std::string colvalue;
        std::string commark(CommentMark); // From char * to string: http://stackoverflow.com/questions/2573834/c-convert-string-or-char-to-string-or-wchar-t
        std::string breakup(BreakupMark); // <added><20150615><dzliu>
        long colpos = -1;
        if (InputText.size()>0) {
            for (int i=0; i<InputText.size(); i++) {
                // get this line
                line = InputText[i];
                // get it trimmed
                trline = csrcstrtrim(line);
                // whether this is empty line?
                if(trline.empty()) {
                    if(OutputArray.size()>0) {break;} else {continue;}
                }
                // whether this is comment line?
                if(0==trline.find(commark)) {
                    continue;
                }
                // whether this is header line? (we take the first non-comment line as the header line.)
                if(-1==colpos) {
                    colpos = csrcstrfindwholeword(line,colname,breakup); // <updated><20150615><dzliu>,breakup
                    if(string::npos==colpos) {
                        std::cout << "CrabStringReadColumn: Column Header not found!" << std::endl;
                        break;
                    }
                }
                // whether this is content line? (we read String content by matching the vertical position!)
                if(colpos>=0) {
                    for(int j=0; j<colname.length(); j++) {
                        colvalue = csrcstrfindwholeword(line,colpos+j,breakup); // <updated><20150615><dzliu>,breakup
                        if(!colvalue.empty()){
                            colvalue = csrcstrtrim(colvalue);
                            OutputArray.push_back(colvalue);
                            // std::cout << colvalue << std::endl; // <TODO> Verbose
                            // std::cout << colvalue << "(" << colpos << ":" << colpos+colname.length()-1 << ")" << std::endl;
                            break;
                        }
                    }
                }
            }
        }
        else {
            std::cout << "CrabStringReadColumn: Unable to read text!" << std::endl;
        }
    }
    return OutputArray;
}



std::vector<std::string> CrabStringReadColumn(std::vector<std::string> InputText, int InputColNumber, const char * CommentMark, const char *BreakupMark)
{
    //
    std::vector<std::string> OutputArray;
    // OutputArray.push_back(superman);
    //
    if(InputColNumber>0 && InputText.size()>0) {
        std::string line;
        std::string trline;
        std::string colname;
        std::string colvalue;
        std::string commark(CommentMark); // From char * to string: http://stackoverflow.com/questions/2573834/c-convert-string-or-char-to-string-or-wchar-t
        std::string breakup(BreakupMark); // <added><20150615><dzliu>
        long colpos = -1;
        if (InputText.size()>0) {
            for (int i=0; i<InputText.size(); i++) {
                // get this line
                line = InputText[i];
                // get it trimmed
                trline = csrcstrtrim(line);
                // whether this is empty line?
                if(trline.empty()) {
                    if(OutputArray.size()>0) {break;} else {continue;}
                }
                // whether this is comment line?
                if(0==trline.find(commark)) {
                    continue;
                }
                // whether this is content line? (we read String content by counting the column number!)
                if(1==1) {
                    colpos = 0; 
                    for(long j=0; j<line.length(); j++) {
                        colvalue = csrcstrfindwholeword(line,j,breakup); // <updated><20150615><dzliu>,breakup
                        if(!colvalue.empty()){
                            colpos--;
                            if((-colpos)==InputColNumber) {
                                colpos = j; break;
                            } else {
                                j+=colvalue.length();
                            }
                        }
                    }
                    if(colpos<0) {
                        std::cout << "CrabStringReadColumn: Column Number is too large!" << std::endl;
                 /* } else if(colpos==0) {
                        std::cout << "CrabStringReadColumn: Empty line! <TODO> Stop or not?" << std::endl; // <Modified><20150201><DzLIU>
                    } else if(colpos>0) {
                        OutputArray.push_back(colvalue);
                        // std::cout << colvalue << std::endl; // <TODO> Verbose
                        // std::cout << colvalue << "(" << colpos << ":" << colpos+colname.length()-1 << ")" << std::endl;
                    }
                 */
                    } else if(colpos>=0) {
                        OutputArray.push_back(colvalue);
                        // std::cout << colvalue << std::endl; // <TODO> Verbose
                        // std::cout << colvalue << "(" << colpos << ":" << colpos+colname.length()-1 << ")" << std::endl;
                    }
                }
            }
        }
        else {
            std::cout << "CrabStringReadColumn: InputText is invalid or InputNumber is not positive!" << std::endl;
        }
    }
    return OutputArray;
}


#endif
