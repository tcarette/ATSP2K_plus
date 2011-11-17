#ifndef _REP_LIB_H
#define _REP_LIB_H
#include <string>
#include <map>
#include "repr.hh"
#include <iostream>
#include <fstream>
#include <strstream>
#include <string>
#include <cmath>

int const sz = 256;

string& ret_str(string&);

int ret_ind_LS(string&, map<string,int>&);

string& u_format(string&, double, double);

#endif

