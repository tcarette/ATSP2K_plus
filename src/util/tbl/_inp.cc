#include "_inp.hh"
#include "_prn.hh"
#include <iostream>
#include <vector>
#include <cctype>
#include <strstream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <string>


int Input_::Init::count = 0;
std::string* Input_::_file_in_;
std::string* Input_::_filebase_;

std::string& ini_obj() {
   static std::string s(""); 
//   std::cout << s << "::" << &s << std::endl;
   return s;
}

bool Input_::_Y_lsj;
bool Input_::_Y_j;

Input_::Init::Init() {
  if (count++ == 0) {
     _filebase_ = new string("");
     _file_in_ = new string("");
     Print_table::_tr_type_ = new string("");
  }
}

Input_::Init::~Init() {
   if (count == 0) {
     delete _filebase_;
     delete _file_in_;
   }
}

void Input_::prog_header() {
   cout << "\n\n             tables: Utility for Computing Levels and"
           " Lifetimes\n\n";
}

void Input_::help_msg() {
   string sf(80,'#');
   cout << sf << endl;
   cout << "Help message about the features of tables:\n\n";
   cout << " . This program takes either a *.l, *.j or *.lsj as argument.\n"
           "   If you do not provide any argument the program prints this\n"
           "   help message and prompts for the input file\n\n";

   cout << ". In the case of *.lsj the user will need to provide Z and NEL\n\n";
         
   cout << " . tables lists all levels and shows the duplicates (unphysical)" 
              " that will\n"
           "   be erased by default. The user has the following options:\n"
           "   - select additional levels to be erased by typing their" 
                 " numbers;\n"
           "   - select a cut-off point by typing the level followed "
                 "by a \"+\" sign,\n"
           "     all levels beyond that point will be erased;\n"
           "   - from the levels listed as duplicates select levels to keep;\n";

   cout<<"\n . tables creates the following files (input ATOM.j or ATOM.lsj):\n"
           "   - ATOM-lev.dat.lt : levels and lifetimes\n"
           "   - ATOM-lev.dat.gj : levels and g_J factors\n"
           "   - ATOM-lin.dat.rt : transition energies, S, gf, aki\n"
           "   - ATOM-lin.dat.wl : transition energies, wavelengths\n"
           "   - ATOM.lsj.new   : ATOM.lsj is cleaned and sorted\n\n";
   cout << " . ATOM.lsj remains unmodified\n\n";
   cout << " . Sorting of transitions and *.lsj is lower->upper\n";
   cout << sf << endl;
}

void Input_::check_input_files(std::string& p) { 
   int sz = p.size();
   std::string s;

   if (p.size() > 3) {
      s = p.substr(sz - 4, sz);
      if (s == ".lsj") {
         _filebase_ = new string(p.substr(0,sz - 4));
         cout << 
            "\n Enter Z (atomic number) and NEL (number of electrons):\n";
         char line[80];
         cin.getline(line,80,'\n'); 
         for (int p = 0; p < 80; p++) {
            if (!(isdigit(line[p]) || 
                isspace(line[p] || line[p] == '\n'))) line[p] = ' ';
            if (line[p] == '\n') break;
         }

         istrstream inb(line);
         int iZ, iNEL;
         inb >> iZ >> iNEL;
         if (iZ < 1 || iZ > 100) 
            {cout << " Z must be between 1:100..\n"; exit(1);}
         else if (iNEL < 1 || iNEL > 40) 
            {cout << " NEL must be between 1:40..\n"; exit(1);}
         else if (iNEL > iZ) 
            {cout << " Z < NEL: invalid..\n"; exit(1);}
         Print_table::_Z_ = iZ;
         Print_table::_NEL_ = iNEL; 
         return;
      }
   }
//   if (p.size() > 8) {
//       s = p.substr(sz - 8, sz); 
//      if (s == "-lev.dat") {
//        *_filebase_ = p.substr(0,sz - 7);
//          return;
//       }
//   }
   s = p.substr(sz - 2, sz);
   if (s == ".j" || s == ".l") {
      _filebase_ = new string(p.substr(0,sz - 2)); 
      return;
   }

   cout << "\n\n\"" << p << "\" not recognized, "
           " only *.j, *.l, *.lsj accepted! Exitting...\n\n\n"; 
   exit(1);
}

void Input_::check_arg(int argc, char* argv) {
   prog_header();
   std::string file_in;
   if (argc < 2 ) {
      Input_::help_msg();
      std::cout << "\nEither a *.j, *.l, *.lsj "
      " files is expected as an argument,\nto continue, enter "
      "a filename:\n";
       getline(cin,file_in,'\n');
//      std::cin >> file_in; 
//      cin.ignore();
   } else {
     string S(argv);
     file_in = S;
   }
cerr << " Input file is \"" << file_in << "\"" << endl;

   {
      ifstream in(file_in.c_str()); 
      if (!in)  {
         cout << "\n\"" << file_in << "\" not found! Quitting...\n\n";
         exit(1);
      }
   } 
   check_input_files(file_in);
   {
       std::string if_lsj = *_filebase_ + ".lsj";
       ifstream in_test(if_lsj.c_str());
       if (in_test) _Y_lsj = true;
       std::string if_j  = *_filebase_ + ".j";
       ifstream in_j(if_j.c_str());
       if (in_j) _Y_j = true;
    }
   _file_in_ =  new string(file_in);
}

// prompt for input
void Input_::e_levels(Dt::Base_table& t) {
   std::cout << "\nThe following duplicate states will be erased:\n\n";
   int i = 0, j = 0;
   for (Dt::LCI cit = t.lv.begin(); cit != t.lv.end(); cit++) {
   ++i;
      if (!(*cit)->unph) {
         std::cout << std::setw(3) << i << " ";
         if (!(++j%15)) std::cout << std::endl;
      }
   }
   if(!j) {
       std::cout << "No duplicate states are present!..\n";
       std::cout << "\nSelect levels to be erased.";
   }
   else {
      std::cout << "\n\nSelect more unphyisical,";
   }
   std::cout << " use a \"+\" at the end for a cut-off point. Example:\n"
        "3 8 +, will additionally remove 3, 8 and higher.\n";
}

//function reading input 
int Input_::inp_loop (Filter_inp* p) {
   for (;;) {
      try {
          p->start();
          while (p->read()) { }
          p->compute();
          p->write();
          return p->result();
      }
      catch (Filter_inp::Retry& m) { if (int i = p->retry(m)) return i; }
      catch(...) { cerr << "Fatal filter error\n"; return 1; }
   }
}
   
int Input_::Filter_lev::read() {
   char c;
   is.get(c);
   if (c == '+') cutoff = true;
   if (c == '+' || c == '\n') { buf[i++] = '\0'; return 0;}
   if (isdigit(c) || isspace(c)) { buf[i++] = c;} 
   return is.good();
}

void Input_::Filter_lev::compute() { 
   istrstream istrm(buf,sz);
   i = 0; int n = 0;
   while (istrm >> n) { 
      ve.push_back(n); 
      i++; }
}

int Input_::Filter_lev::result() { 
    sort(ve.begin(),ve.end());
    vector<int>::iterator pc = unique(ve.begin(),ve.end());
    ve.erase(pc,ve.end()); 
    return 0;
}

void Input_::set_add_lev(Dt::Base_table& t, Input_::Filter_lev& f) {
//  adding or "unerasing" levels
   using namespace Dt; 
   int static ic = 0;
   int sz = f.ve.size();
   if (sz == 0) return;
   int max = f.ve[sz-1];
   for (Dt::LI itl = t.lv.begin(); itl != t.lv.end(); itl++) {
      ++ic;
      for (VCI itv = f.ve.begin(); itv != f.ve.end(); itv++) {
         if (ic == *itv)  
            (*itl)->unph = ((*itl)->unph) ? false : true; 
         if (f.cutoff && (ic >= max)) {
            (*itl)->unph = false;
         }
         //cout << (*itl)->unph << endl;
      }  
   }
}



