/*
 
 We must be careful about codecs. 
 All files input should be utf8. 
 
 */

#ifndef H_CrabStringUnicode
#define H_CrabStringUnicode 
#include <iostream>
#include <iomanip> // std::setfill, std::setw
#include <cstdlib>
#include <string>
#include <vector>
#include <clocale> // <20160725><dzliu> modified from <locale> to <clocale> for mac.
#include <wchar.h>

using namespace std;

std::wstring CrabStringWideString(std::string str1, const char *Locale=NULL);
std::wstring CrabStringWideString(const char *str1, const char *Locale=NULL);

std::string csu_wstring_to_string(std::wstring str1);

std::wstring csu_utf8_to_wstr(std::string& src);
std::string csu_wstr_to_utf8(std::wstring& src);
void csu_utf8_to_wstr(std::string &src, std::wstring &dest);
void csu_wstr_to_utf8(std::wstring &src, std::string &dest);
void csu_print_str_chars(std::string str);
void csu_print_wstr_chars(std::wstring wstr);
void csu_print_str(std::string str, const char *Locale="zh_CN.UTF-8");
void csu_print_wstr(std::wstring wstr, const char *Locale="zh_CN.UTF-8");




std::wstring CrabStringWideString(std::string stdsS1, const char *Locale)
{
    return CrabStringWideString(stdsS1.c_str(),Locale);
}


std::wstring CrabStringWideString(const char *charS1, const char *Locale)
{
    //
    int debug =0;
    //
    if(debug) {
        charS1 = " 你好啊 !"; // u8" 你好啊 !" (since C++11) 
    }
    //
    std::size_t byteS1 = strlen(charS1);
    //
    if(debug) {
        std::cout << "DEBUG:csu_wstring() charS1=" << charS1 << " byteS1=" << strlen(charS1) << std::endl;
    }
    // 
    std::string stdmS1(charS1,byteS1); // std::string multi-byte String
    std::wstring stdwS1 = csu_utf8_to_wstr(stdmS1); // std::string wide-char String
    //
    if(debug) {
        std::cout << "DEBUG:csu_wstring() stdmS1="; csu_print_str(stdmS1); std::cout << " lenS1=" << strlen(stdmS1.c_str());
        std::cout << " stdmS1.find(好)=" << stdmS1.find("好"); csu_print_str_chars(stdmS1); std::cout << std::endl;
        std::wcout << L"DEBUG:csu_wstring() stdwS1="; csu_print_wstr(stdwS1); std::cout << " lenWS1=" << wcslen(stdwS1.c_str());
        std::wcout << L" stdwS1.find(好)=" << stdwS1.find(L"好"); csu_print_wstr_chars(stdwS1); std::cout << std::endl;
    }
    //
    // std::wstring stdwS2 = L" 你好";
    // if(debug) {
    //     std::cout << "DEBUG:csu_wstring() stdwS2?" << stdwS2.find(L"好") << " bytWS2=" << wcslen(stdwS2.c_str()) << std::endl;
    //     std::cout << "DEBUG:csu_wstring() stdwS2?" << stdwS2.find(L"好"); csu_print_wstr_chars(stdwS2);
    // }
    //
    return stdwS1;
}


std::string csu_wstring_to_string(std::wstring str1)
{
    //
    unsigned int bytlen1 = sizeof str1; // L"你好" is 8 bytes on x64.
    unsigned int bytlen2 = (unsigned int)(bytlen1/sizeof(wchar_t)*sizeof(char));
    char *chars2 = new char[bytlen2+1];
    memset(chars2,'\0',bytlen2+1);
    std::cout << "DEBUG:csu_wstring_to_string() bytlen1=" << bytlen1 << std::endl;
    std::size_t bytnew2 = std::wcstombs(chars2,str1.c_str(),bytlen1);
    std::string str2(chars2);
    return str2;
}






std::wstring csu_utf8_to_wstr(std::string& src)
{
    std::wstring dest;
    csu_utf8_to_wstr(src,dest);
    return dest;
}

std::string csu_wstr_to_utf8(std::wstring& src)
{
    std::string dest;
    csu_wstr_to_utf8(src,dest);
    return dest;
}

void csu_utf8_to_wstr(std::string &src, std::wstring &dest)
{
    // http://www.linuxquestions.org/questions/programming-9/wstring-utf8-conversion-in-pure-c-701084/
	dest.clear();
	wchar_t w = 0;
	int bytes = 0;
	wchar_t err = L'�';
	for (size_t i = 0; i < src.size(); i++){
		unsigned char c = (unsigned char)src[i];
		if (c <= 0x7f){//first byte
			if (bytes){
				dest.push_back(err);
				bytes = 0;
			}
			dest.push_back((wchar_t)c);
		}
		else if (c <= 0xbf){//second/third/etc byte
			if (bytes){
				w = ((w << 6)|(c & 0x3f));
				bytes--;
				if (bytes == 0)
					dest.push_back(w);
			}
			else
				dest.push_back(err);
		}
		else if (c <= 0xdf){//2byte sequence start
			bytes = 1;
			w = c & 0x1f;
		}
		else if (c <= 0xef){//3byte sequence start
			bytes = 2;
			w = c & 0x0f;
		}
		else if (c <= 0xf7){//3byte sequence start
			bytes = 3;
			w = c & 0x07;
		}
		else{
			dest.push_back(err);
			bytes = 0;
		}
	}
	if (bytes)
		dest.push_back(err);
}

void csu_wstr_to_utf8(std::wstring& src, std::string& dest)
{
    // http://www.linuxquestions.org/questions/programming-9/wstring-utf8-conversion-in-pure-c-701084/
	dest.clear();
	for (size_t i = 0; i < src.size(); i++){
		wchar_t w = src[i];
		if (w <= 0x7f)
			dest.push_back((char)w);
		else if (w <= 0x7ff){
			dest.push_back(0xc0 | ((w >> 6)& 0x1f));
			dest.push_back(0x80| (w & 0x3f));
		}
		else if (w <= 0xffff){
			dest.push_back(0xe0 | ((w >> 12)& 0x0f));
			dest.push_back(0x80| ((w >> 6) & 0x3f));
			dest.push_back(0x80| (w & 0x3f));
		}
		else if (w <= 0x10ffff){
			dest.push_back(0xf0 | ((w >> 18)& 0x07));
			dest.push_back(0x80| ((w >> 12) & 0x3f));
			dest.push_back(0x80| ((w >> 6) & 0x3f));
			dest.push_back(0x80| (w & 0x3f));
		}
		else
			dest.push_back('?');
	}
}

void csu_print_str_chars(std::string str)
{
    //
    unsigned char ch;
    const char *ptr = str.c_str();
    for(int k=0; k<=strlen(ptr); k++) {
        ch = ptr[k]; std::cout << " ["<<k<<"]=" << setw(2*sizeof(char)) << setfill('0') << hex << uppercase << (int)(ch) << std::flush;
    }
    // std::cout << std::endl;
}

void csu_print_wstr_chars(std::wstring wstr)
{
    //
    unsigned wchar_t ch;
    const wchar_t *wptr = wstr.c_str();
    for(int k=0; k<=wcslen(wptr); k++) {
        ch = wptr[k]; std::cout << " ["<<k<<"]=" << setw(sizeof(wchar_t)) << setfill('0') << hex << uppercase << (int)(ch) << std::flush;
    }
    // std::cout << std::endl;
}

void csu_print_str(std::string str, const char *Locale)
{
    //
    std::cout << str << std::flush;
}

void csu_print_wstr(std::wstring wstr, const char *Locale)
{
    //
    std::setlocale(LC_ALL,Locale);
    std::wcout << wstr << std::flush;
    std::setlocale(LC_ALL,"");
    //
    // how to get current locale? use setlocale(LC_ALL,"").
    // http://stackoverflow.com/questions/12170488/how-to-get-current-locale-of-my-environment
    // cout << "DEBUG:csu_wstring() LC_ALL=" << setlocale(LC_ALL,NULL) << endl;
    // cout << "DEBUG:csu_wstring() LC_CTYPE=" << setlocale(LC_CTYPE,NULL) << endl;
}

#endif
