/*
 
 We must be careful about codecs.
 All files input should be utf8. 
 
 last update:
 
     2015-03-16 fixed problem of reading whole word
                "if(string::npos==foundposR)"
                "if(string::npos==foundposL)"
                "if(wstring::npos==foundposR)"
                "if(wstring::npos==foundposL)"
 
     2015-06-15 add BreakupMarkOnceMore
 
 */

#ifndef H_CrabStringFindWholeWord
#define H_CrabStringFindWholeWord
#include <iostream>
#include <iomanip> // std::setfill, std::setw
#include <cstdlib>
#include <string>
#include <locale>
#include <wchar.h>

using namespace std;

long CrabStringFindWholeWord(const string &t, const string &strsearch, const string &BreakupMark = " ", const string &BreakupMarkOnceMore = "", int debug = 0);

long CrabStringFindWholeWord(const wstring &t, const wstring &strsearch, const wstring &BreakupMark = L" ", const wstring &BreakupMarkOnceMore = L"", int debug = 0);

string CrabStringFindWholeWord(const string &t, size_t pos, const string &BreakupMark = " ", const string &BreakupMarkOnceMore = "", int debug = 0);

wstring CrabStringFindWholeWord(const wstring &t, size_t pos, const wstring &BreakupMark = L" ", const wstring &BreakupMarkOnceMore = L"", int debug = 0);



long CrabStringFindWholeWord(const string &t, const string &strsearch, const string &BreakupMark, const string &BreakupMarkOnceMore, int debug)
{
    size_t foundpos = 0;
    while(string::npos!=foundpos) {
        foundpos = t.find(strsearch,foundpos);
        if(string::npos!=foundpos) {
            if(foundpos>=BreakupMark.length() && foundpos+strsearch.length()<t.length()-BreakupMark.length()) {
                if(BreakupMark==t.substr(foundpos-1,BreakupMark.length()) &&
                   BreakupMark==t.substr(foundpos+strsearch.length(),BreakupMark.length())) {break;}
                if(BreakupMarkOnceMore.length()>0){
                    if(BreakupMarkOnceMore==t.substr(foundpos-1,BreakupMarkOnceMore.length()) &&
                       BreakupMarkOnceMore==t.substr(foundpos+strsearch.length(),BreakupMarkOnceMore.length())) {break;}
                }
            }
            if(foundpos==0 && foundpos+strsearch.length()<t.length()-BreakupMark.length()) {
                if(BreakupMark==t.substr(foundpos+strsearch.length(),BreakupMark.length())) {break;}
                if(BreakupMarkOnceMore.length()>0){
                    if(BreakupMarkOnceMore==t.substr(foundpos+strsearch.length(),BreakupMarkOnceMore.length())) {break;}
                }
            }
            if(foundpos>=BreakupMark.length() && foundpos+strsearch.length()==t.length()-BreakupMark.length()-1) {
                if(BreakupMark==t.substr(foundpos-1,BreakupMark.length())) {break;}
                if(BreakupMarkOnceMore.length()>0){
                    if(BreakupMarkOnceMore==t.substr(foundpos-1,BreakupMarkOnceMore.length())) {break;}
                }
            }
            // <Corrected><20150316><dzliu> when the header is at the end of the line
            if(foundpos+strsearch.length()-1==t.length()-1) {
                break;
            }
        } else {
            break;
        }
        foundpos++;
    }
    return foundpos;
}



long CrabStringFindWholeWord(const wstring &t, const wstring &strsearch, const wstring &BreakupMark, const wstring &BreakupMarkOnceMore, int debug)
{
    size_t foundpos = 0;
    while(wstring::npos!=foundpos) {
        foundpos = t.find(strsearch,foundpos);
        if(debug) {
            std::wcout << "CrabStringFindWholeWord: Debuging: Searching \"" << strsearch << "\" we got foundpos = " << foundpos << std::endl;
            std::wcout << "CrabStringFindWholeWord: Debuging: BreakupMarkOnceMore = \"" << BreakupMarkOnceMore << "\"" << std::endl; // <debug><20170324><dzliu>
        }
        if(wstring::npos!=foundpos) {
            if(foundpos>=BreakupMark.length() && foundpos+strsearch.length()<t.length()-BreakupMark.length()) {
                if(BreakupMark==t.substr(foundpos-1,BreakupMark.length()) &&
                   BreakupMark==t.substr(foundpos+strsearch.length(),BreakupMark.length())) {break;}
                if(BreakupMarkOnceMore.length()>0){
                    if(BreakupMarkOnceMore==t.substr(foundpos-1,BreakupMarkOnceMore.length()) &&
                       BreakupMarkOnceMore==t.substr(foundpos+strsearch.length(),BreakupMarkOnceMore.length())) {break;}
                }
            }
            if(foundpos==0 && foundpos+strsearch.length()<t.length()-BreakupMark.length()) {
                if(BreakupMark==t.substr(foundpos+strsearch.length(),BreakupMark.length())) {break;}
                if(BreakupMarkOnceMore.length()>0){
                    if(BreakupMarkOnceMore==t.substr(foundpos+strsearch.length(),BreakupMarkOnceMore.length())) {break;}
                }
            }
            if(foundpos>=BreakupMark.length() && foundpos+strsearch.length()==t.length()-BreakupMark.length()-1) {
                if(BreakupMark==t.substr(foundpos-1,BreakupMark.length())) {break;}
                if(BreakupMarkOnceMore.length()>0){
                    if(BreakupMarkOnceMore==t.substr(foundpos-1,BreakupMarkOnceMore.length())) {break;}
                }
            }
            // <Corrected><20150316><dzliu> when the header is at the end of the line
            if(foundpos+strsearch.length()-1==t.length()-1) {
                break;
            }
        } else {
            break;
        }
        foundpos++;
    }
    return foundpos;
}



string CrabStringFindWholeWord(const string &t, size_t pos, const string &BreakupMark, const string &BreakupMarkOnceMore, int debug)
{
    string foundstr;
    if(0==pos) {
        string strR = t;
        size_t foundposL = 0;
        size_t foundposR = t.find(BreakupMark);
        if(debug) {
            std::cout << "CrabStringFindWholeWord: Debuging: BreakupMarkOnceMore = \"" << BreakupMarkOnceMore << "\"" << std::endl; // <debug><20170324><dzliu>
        }
        if(BreakupMarkOnceMore.length()>0){ // <added><20150615><dzliu> use another BreakupMarkOnceMore
            if(string::npos==foundposR){foundposR = strR.find(BreakupMarkOnceMore);}else{
                if(strR.find(BreakupMarkOnceMore)<foundposR){foundposR=strR.find(BreakupMarkOnceMore);}
            }
        }
        if(string::npos==foundposR){foundposR=t.length()-1;}else{
            foundposR=foundposR-1;
        }
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    } else if(t.length()-1==pos) {
        string strL = t;
        size_t foundposL = t.rfind(BreakupMark);
        if(BreakupMarkOnceMore.length()>0){ // <added><20150615><dzliu> use another BreakupMarkOnceMore
            if(string::npos==foundposL){foundposL = strL.rfind(BreakupMarkOnceMore);}else{
                if(strL.rfind(BreakupMarkOnceMore)>foundposL){foundposL=strL.rfind(BreakupMarkOnceMore);}
            }
        }
        if(string::npos==foundposL){foundposL=0;}else{
            foundposL=foundposL+1;
        }
        size_t foundposR = t.length()-1;
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    } else if(0<pos && t.length()-BreakupMark.length()-1>pos) {
        string strL = t.substr(0,pos);
        string strR = t.substr(pos);
        size_t foundposL = strL.rfind(BreakupMark);
        if(BreakupMarkOnceMore.length()>0){ // <added><20150615><dzliu> use another BreakupMarkOnceMore
            if(string::npos==foundposL){foundposL = strL.rfind(BreakupMarkOnceMore);}else{
                if(strL.rfind(BreakupMarkOnceMore)>foundposL){foundposL=strL.rfind(BreakupMarkOnceMore);}
            }
        }
        if(string::npos==foundposL){foundposL=0;}else{
            if(BreakupMarkOnceMore.length()>0){if(strL.rfind(BreakupMarkOnceMore)>foundposL){foundposL=strL.rfind(BreakupMarkOnceMore);}}
            foundposL=foundposL+1;
        }
        size_t foundposR = strR.find(BreakupMark);
        if(BreakupMarkOnceMore.length()>0){ // <added><20150615><dzliu> use another BreakupMarkOnceMore
            if(string::npos==foundposR){foundposR = strR.find(BreakupMarkOnceMore);}else{
                if(strR.find(BreakupMarkOnceMore)<foundposR){foundposR=strR.find(BreakupMarkOnceMore);}
            }
        }
        if(string::npos==foundposR){foundposR=t.length()-1;}else{
            foundposR=foundposR-1+pos;
        }
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    }
    return foundstr;
}



wstring CrabStringFindWholeWord(const wstring &t, size_t pos, const wstring &BreakupMark, const wstring &BreakupMarkOnceMore, int debug)
{
    wstring foundstr;
    if(0==pos) {
        wstring strR = t;
        size_t foundposL = 0;
        size_t foundposR = t.find(BreakupMark);
        if(BreakupMarkOnceMore.length()>0){ // <added><20150615><dzliu> use another BreakupMarkOnceMore
            if(wstring::npos==foundposR){foundposR = strR.find(BreakupMarkOnceMore);}else{
                if(strR.find(BreakupMarkOnceMore)<foundposR){foundposR=strR.find(BreakupMarkOnceMore);}
            }
        }
        if(wstring::npos==foundposR){foundposR=t.length()-1;}else{
            foundposR=foundposR-1;
        }
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    } else if(t.length()-1==pos) {
        wstring strL = t;
        size_t foundposL = t.rfind(BreakupMark);
        if(BreakupMarkOnceMore.length()>0){ // <added><20150615><dzliu> use another BreakupMarkOnceMore
            if(wstring::npos==foundposL){foundposL = strL.rfind(BreakupMarkOnceMore);}else{
                if(strL.rfind(BreakupMarkOnceMore)>foundposL){foundposL=strL.rfind(BreakupMarkOnceMore);}
            }
        }
        if(wstring::npos==foundposL){foundposL=0;}else{
            foundposL=foundposL+1;
        }
        size_t foundposR = t.length()-1;
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    } else if(0<pos && t.length()-1>pos) {
        wstring strL = t.substr(0,pos);
        wstring strR = t.substr(pos);
        size_t foundposL = strL.rfind(BreakupMark);
        if(BreakupMarkOnceMore.length()>0){ // <added><20150615><dzliu> use another BreakupMarkOnceMore
            if(wstring::npos==foundposL){foundposL = strL.rfind(BreakupMarkOnceMore);}else{
                if(strL.rfind(BreakupMarkOnceMore)>foundposL){foundposL=strL.rfind(BreakupMarkOnceMore);}
            }
        }
        if(wstring::npos==foundposL){foundposL=0;}else{
            foundposL=foundposL+1;
        }
        size_t foundposR = strR.find(BreakupMark);
        if(BreakupMarkOnceMore.length()>0){ // <added><20150615><dzliu> use another BreakupMarkOnceMore
            if(wstring::npos==foundposR){foundposR = strR.find(BreakupMarkOnceMore);}else{
                if(strR.find(BreakupMarkOnceMore)<foundposR){foundposR=strR.find(BreakupMarkOnceMore);}
            }
        }
        if(wstring::npos==foundposR){foundposR=t.length()-1;}else{
            foundposR=foundposR-1+pos;
        }
        foundstr = t.substr(foundposL,foundposR-foundposL+1);
    }
    return foundstr;
}



#endif
