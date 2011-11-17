#ifndef _PRN_H_
#define _PRN_H_

#include "_tt.hh"
#include "_inp.hh"
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <map>
#include <strstream>
#include <cstring>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <functional>

namespace Print_table {
   using namespace std;
   const int WC = 28;  // width of configuration column
   const int WT =  2;  // width of term
   const int WJ =  5;  // width of J value
   const int WE = 15;  // total energy
   const int WL = 11;  // Levels
   const int WS =  9;  // Splitting
   const int WM = 10;  // Lifetimes
   extern const char* atom[];
   extern const char* roman[];
   extern const int sa, sr;
   extern int _NEL_;
   extern double  _RZ_;
   extern int _Z_;
   extern string* _tr_type_;
   void Print_files(Dt::Base_table*);
   void Print_files(Dt::Base_table*,string&, string&);
   void file_ZEL_N(string& F,int& Z);
}

#endif
