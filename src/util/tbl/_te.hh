#ifndef _TE_H_
#define _TE_H_
#include "_tt.hh"
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


using namespace Dt;
using namespace std;

class Error_lsj {
public:
   string _lsj1;
   string _lsj0;
   string _lsj2;
   string _lsj3;
   string _lsj4;
   string _lsj5;
   string _lsj6;
   string _lsj7;
   ofstream FLOG;
   Error_lsj();

   ~Error_lsj() {FLOG.close();}

   void print_ok() {
      cerr << _lsj1 << _lsj0 << _lsj2 << _lsj3 <<
              _lsj4 << _lsj5 << _lsj6 << _lsj0 << _lsj7;
   }
   virtual void print_error(ST* t)  { 
      cerr << "Incorrect .lsj format\n";
      t->print_ST_error();        
      cerr << "The expected format is the next 5 consequitive lines:\n"; 
      print_ok();
   }
};

class Error_lsj_empty : public Error_lsj {
public:
   virtual void print_error(ST* t) {
      cerr << " Error: incorrect .lsj format\n";
      Error_lsj::print_error(t);
   }
};

class Error_lsj_S : public Error_lsj {
public:
   virtual void print_error(ST* t) {
      cerr << " Error: The characters \"S =\" in "
              "the wrong place or missing\n";
      Error_lsj::print_error(t);
   }
};

class Error_lsj_negatives : public Error_lsj{
public:
   virtual void print_error(ST* t)  {
      cerr << "Error: There were negative values:\n";
      Error_lsj::print_error(t);
   }
};

class Error_lsj_line1: public Error_lsj {
public:
   virtual void print_error(ST* t) {
      cerr << "Error: The line below has incorrect format:\n";
      Error_lsj::print_error(t);
      cerr << "Expected: ";
      cerr << _lsj2 << endl;
   }
};

class Error_lsj_line2: public Error_lsj {
public:
   virtual void print_error(ST* t) {
      cerr << "Error: The line below has incorrect format:\n";
      Error_lsj::print_error(t);
      cerr << "Expected: ";
      cerr << _lsj3 << endl;
   }
};
class Error_lsj_line3: public Error_lsj {
public:
   virtual void print_error(ST* t) {
      cerr << "Error: The line below has incorrect format:\n";
      Error_lsj::print_error(t);
      cerr << "Expected: ";
      cerr << _lsj4 << endl;
   }
};
class Error_lsj_line4: public Error_lsj {
public:
   virtual void print_error(ST* t) {
      cerr << "Error: The line below has incorrect format:\n";
      Error_lsj::print_error(t);
      cerr << "Expected: ";
      cerr << _lsj5 << endl;
   }
};
class Error_lsj_line5: public Error_lsj {
public:
   virtual void print_error(ST* t) {
      cerr << "Error: The line below has incorrect format:\n";
      Error_lsj::print_error(t);
      cerr << "Expected: ";
      cerr << _lsj6 << endl;
   }
};

#endif
 
