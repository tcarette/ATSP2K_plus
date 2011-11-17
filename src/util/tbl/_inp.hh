#ifndef _INP_H_
#define _INP_H_
#include "_tt.hh"
#include <iostream>
#include <vector>
#include <cctype>
#include <strstream>
#include <algorithm>
#include <string>


namespace Input_ {
   class Init;
   using namespace std;
   typedef vector<int>::const_iterator VCI;
   class Io_obj;
   void prog_header();
   void check_input_files(std::string&);
   void check_arg(int, char*);
   void e_levels(Dt::Base_table&);
   class Filter_inp;
   int inp_loop (Filter_inp*); 
   class Filter_lev;
   void set_add_lev( Dt::Base_table& t, Input_::Filter_lev& f);
   void select_inp();
   extern std::string file_in, file_j, file_lsj, file_lev, file_tr, file_new;
   extern double Z;
  // template<class T> class Inp_ini; 
   class Inp_ini;
   class Inp_lsjini;
   class inp_lini;
   extern std::string* _filebase_;
   extern std::string* _file_in_;
   extern bool _Y_lsj;
   extern bool _Y_j;
   void help_msg();

}

class Input_::Init {
public:
   int static count;
public:
   Init();
   ~Init();
};

namespace {
  Input_::Init FI;
}

class Input_::Io_obj {
public:
   virtual Io_obj* clone() const = 0;
   virtual ~Io_obj() {}
};


// make the two below a template if another 
//namespace Input_ {
//   template<class T> class Inp_ini : public T{
//}
class Input_::Inp_ini : 
public Dt::ETable_jini, public Input_::Io_obj {
public:
   Io_obj* clone() const {return new Inp_ini(*this);}
   Inp_ini(const std::string& fn) : ETable_jini(fn) {}
   static Io_obj* if_j(const std::string& fn) { return new Inp_ini(fn);}
};

class Input_::Inp_lsjini :   
public Dt::ETable_lsjini, public Input_::Io_obj {
public:
   Io_obj* clone() const {return new Inp_lsjini(*this);}
   Inp_lsjini(const std::string& fn) : ETable_lsjini(fn) {}
   static Io_obj* if_lsj(const std::string& fn) { return new Inp_lsjini(fn);}
};


//general, input "frame"
class Input_::Filter_inp {
public:
   class Retry {
   public:
      virtual const char* message() {return 0;}
  };
  virtual void start() {}
  virtual int read() = 0;
  virtual void write() {}
  virtual void compute() {}
  virtual int result() = 0;
  virtual int retry(Retry& m) {cerr << m.message() << '\n'; return 2;}
  virtual ~Filter_inp() {}
};
 
//specialized to read unphysical levels
class Input_::Filter_lev : public Input_::Filter_inp {
   static const int sz = 256;
   istream& is;
   ostream& os;
   vector<int> ve;
   char buf[sz];
   int i; 
public:
   bool cutoff;
   int read();
   void compute();
   int result();
  friend void Input_::set_add_lev(Dt::Base_table&, Input_::Filter_lev&);
   Filter_lev(istream& ii, ostream& oo) : is(ii), os(oo), i(0),
              cutoff(false)  {}
   ~Filter_lev() {}
};

#endif



