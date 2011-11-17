//#include "_inp.hh"
#include "_tt.hh"
#include "_prn.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <sstream>
using namespace std;

Input_::Io_obj* get_obj(const string& fn) {
   typedef Input_::Io_obj* (*PF) (const string&);
   map<string,PF> io_map;
   string s1 = *Input_::_filebase_ + ".j";
   string s2 = *Input_::_filebase_ + ".lsj";
   string s3 = *Input_::_filebase_ + ".l";
   io_map[s1] = &Input_::Inp_ini::if_j;
   io_map[s2] = &Input_::Inp_lsjini::if_lsj;
   io_map[s3] = &Input_::Inp_ini::if_j;
   PF f = io_map[fn]; 
   return f(fn); 
}

int main(int argc, char* argv[]) {
//  check the arguments and get the filename
    Input_::check_arg(argc, argv[1]);
    string fn = *Input_::_file_in_;
    
//  make a object as a function ot the file type(lsj,j,dat)
    Input_::Io_obj* p = get_obj(fn);

    cout << " \nCompute lifetimes from this .lsj file y/n: ";
    string LT;
    getline(cin,LT,'\n');
     
    if (LT[0] != 'y') {
        cout << " You have entered " << LT 
             << " lifetimes will not be computed\n";
    }
    bool lifetimes_ = false;
    if (LT[0] == 'y') lifetimes_ = true; 
    
//  convert pt to base 
    Dt::Base_table* pt = dynamic_cast<Dt::Base_table*>(p); 

//  get input from user (unphysical)
    Input_::e_levels(*pt);
    Input_::Filter_lev fobj(cin,cout);
    Input_::inp_loop(&fobj);
    set_add_lev(*pt, fobj);
//    pt->print_lev(cout,false);
    cout << pt->lv.size() << " levels initially:"<< endl;
    pt->LV_rm_unph();
    cout << pt->lv.size() << " levels after erasing unphysical" << endl;
//  set term and configuration representation for the table levles
//    pt->rep_TC();

// if a *.lsj file is present compute transition data

    if (Input_::_Y_lsj) pt->JTR_rmv(); 
    if (Input_::_Y_lsj && lifetimes_) pt->JTR_lifetimes();

    pt->LV_print(cout);

    bool db_yes = false;
    if (Input_::_Y_lsj) {
       cin.ignore();
       cout << "\n Save database files (y/n): ";
       string save_db;
       getline(cin,save_db,'\n');
       if (save_db[0] == 'y') {
          db_yes = true;
          string REF = "";
          cout << "   Enter a reference number: ";
          getline(cin,REF,'\n'); 
          if (!REF.size()) {
             cout << " Error: Must enter a reference number, "
                       "exit(1) \n";
             exit(1);
          }
          for (int i = 0; i < REF.size(); i++) {
            if (!(isdigit(REF[i]))) { 
               cout << " Error: Wrong input, must be a number, "
                       "exit(1): " << REF << "\n"; 
               exit(1);
            }
          }
          cout << "   Enter up to 10 characters file extension: ";
          string str_ext;
          getline(cin,str_ext,'\n');
          cout << "...saving in database format:";
          string EXT("");
          if (str_ext.size() > 10) EXT = str_ext.substr(0,10);
          else EXT = str_ext;
          for (int i = 0; i < EXT.size(); i++) {
            if (isspace(EXT[i])) EXT[i] = '_'; 
          }
          Print_table::Print_files(pt, REF, EXT);
       } else {
         cerr << "\n...not saving in database format:";
         Print_table::Print_files(pt);
       }
    }

    if (!Input_::_Y_lsj) Print_table::Print_files(pt);
}

