#include "_prn.hh"
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
#include <sstream>

using namespace std;

int Print_table::_Z_;
int Print_table::_NEL_;
std::string* Print_table::_tr_type_;
double Print_table::_RZ_;

void Print_table::file_ZEL_N(string& F,int& Z) {
    map<int,string> MP;
    MP[1] ="H"; MP[2]= "He";MP[3]= "Li";MP[4]= "Be";MP[5]= "B"; MP[6]= "C";
    MP[7] ="N"; MP[8]= "O" ;MP[9]= "F"; MP[10]="Ne";MP[11]="Na";MP[12]="Mg";
    MP[13]="Al";MP[14]="Si";MP[15]="P"; MP[16]="S" ;MP[17]="Cl";MP[18]="Ar";
    MP[19]="K"; MP[20]="Ca";MP[21]="Sc";MP[22]="Ti";MP[23]="V"; MP[24]="Cr";
    MP[25]="Mn";MP[26]="Fe";MP[27]="Co";MP[28]="Ni";MP[29]="Cu";MP[30]="Zn";
    MP[31]="Ga";MP[32]="Ge";MP[33]="As";MP[34]="Se";MP[35]="Br";MP[36]="Kr";
    MP[37]="Rb";MP[38]="Sr";MP[39]="Y"; MP[40]="Zr";MP[41]="Nb";MP[42]="Mo";
    MP[43]="Tc";MP[44]="Ru";MP[45]="Rh";MP[46]="Pd";MP[47]="Ag";MP[48]="Cd";
    MP[49]="In";MP[50]="Sn";MP[51]="Sb";MP[52]="Te";MP[53]="I"; MP[54]="Xe";
    MP[55]="Cs";MP[56]="Ba";MP[57]="La";MP[58]="Ce";MP[59]="Pr";MP[60]="Nd";
    MP[61]="Pm";MP[62]="Sm";MP[63]="Eu";MP[64]="Gd";MP[65]="Tb";MP[66]="Dy";
    MP[67]="Ho";MP[68]="Er";MP[69]="Tm";MP[70]="Yb";MP[71]="Lu";MP[72]="Hf";
    MP[73]="Ta";MP[74]="W"; MP[75]="Re";MP[76]="Os";MP[77]="Ir";MP[78]="Pt";
    MP[79]="Au";MP[80]="Hg";MP[81]="Tl";MP[82]="Pb";MP[83]="Bi";MP[84]="Po";
    MP[85]="At";MP[86]="Rn";MP[87]="Fr";MP[88]="Ra";MP[89]="Ac";MP[90]="Th";
    MP[91]="Pm";MP[92]="U"; MP[93]="Np";MP[94]="Pu";MP[95]="Am";MP[96]="Cm";
    MP[96]="Bk";MP[98]="Cf";MP[99]="Es";MP[100]="Fm";
    MP[101]="Md";MP[102]="No";MP[103]="Lr";

    F=MP[Z];
}

void Print_table::Print_files(Dt::Base_table* pt,
                  string& REF, string& EXT) {

   istringstream ist(REF);
   int ref_N;
   ist >> ref_N;

   int NEL = Print_table::_NEL_;
   int Z = Print_table::_Z_;
   string str_NEL;
   file_ZEL_N(str_NEL,NEL);
   
   string str_Z;
   ostringstream ost(str_Z);
   ost << Z;
   str_Z = ost.str();

   string str_refN;
   ostringstream ostr(str_refN);
   ostr << ref_N;
   str_refN = ostr.str();
   
   string FB = str_NEL + "_" + str_Z + "." + str_refN + 
             + "." + EXT;
   
   string file_lsj,file_lsj_bin,file_lev,file_lin,file_tr_aux,file_lv_aux;
   if(Input_::_Y_lsj) {
      file_lsj_bin = FB + ".lsj.bin";
      ofstream os(file_lsj_bin.c_str(),ios_base::out|ios_base::binary);
      pt->LSJ_print_bin(os,ref_N);
      os.close();
   }

   if(Input_::_Y_lsj) {
       file_lsj = FB + ".lsj.db";
       cerr << " To modify the reference energy:\n   -open file: " 
            << file_lsj << ", edit the first line\n   -rename"
            << " the file to .lsj"
            << " and rerun tables on the modified file";
       ofstream os(file_lsj.c_str());
       pt->LSJ_print(os);
       os.close();
   }

   {
       file_lev = FB + "-lev.lt.db";
       ofstream os(file_lev.c_str());
       pt->LV_print(os);
       os.close();
   }

   if (Input_::_Y_j) { 
       file_lv_aux = FB + "-lev.gj.db";
       ofstream os(file_lv_aux.c_str());
       pt->LV_print_aux(os);
       os.close();
   }

    if(Input_::_Y_lsj) {
       file_lin = FB +"-lin.dat.rt";
       ofstream os(file_lin.c_str());
       pt->JTR_print(os);
       os.close();
   }

    if(Input_::_Y_lsj)  {
       file_tr_aux = FB +"-lin.dat.mp";
       ofstream os(file_tr_aux.c_str());
       pt->JTR_print_MP(os);
       os.close();
   }

   cerr << "\n All saved in database format:\n";
   cerr << "   reference number : " << ref_N << endl;
   cerr << "   extension        : " << EXT << endl;
   cerr << "   lsj binary data  : " << file_lsj_bin << endl;
   cerr << "   lsj ascii data   : " << file_lsj << endl;
   cerr << "   levels           : " << file_lev << endl;
   if (Input_::_Y_j) cerr << "   Zeeman g_J       : " << file_lv_aux << endl;
   cerr << "   rates(multiplet) : " << file_tr_aux << endl; 
   cerr << "Place the files in the directory " << 
        "hf5:/home/wwwdata/data/" << str_NEL << endl 
        << "END" <<endl;
}

void Print_table::Print_files(Dt::Base_table* pt) {
   string FB = *Input_::_filebase_;
   string file_lsj,file_lev,file_lv_aux,file_lin,file_tr_aux; 
   if(Input_::_Y_lsj) { 
       file_lsj = FB + ".lsj.db";
       ofstream os(file_lsj.c_str());
       pt->LSJ_print(os);
       os.close();
   }
 
   { 
       file_lev = FB + "-lev.lt.db";
       ofstream os(file_lev.c_str());
       pt->LV_print(os);
       os.close();
   }

   if (Input_::_Y_j) {
       file_lv_aux = FB + "-lev.gj.db";
       ofstream os(file_lv_aux.c_str());
       pt->LV_print_aux(os);
       os.close();
   }

    if(Input_::_Y_lsj) {
       file_lin = FB +"-lin.dat.rt";
       ofstream os(file_lin.c_str());
       pt->JTR_print(os);
       os.close();
   }
   
    if(Input_::_Y_lsj)  {
       file_tr_aux = FB +"-lin.dat.mp";
       ofstream os(file_tr_aux.c_str());
       pt->JTR_print_MP(os);
       os.close();
   }
   cerr << " The data is saved in:\n";
   if (Input_::_Y_lsj) cerr << "   lsj data      : " << file_lsj << endl;
   cerr << "   levels        : " << file_lev << endl;
   if (Input_::_Y_j)   cerr << "   g_J factors   : " << file_lv_aux << endl;
   if (Input_::_Y_lsj) cerr << "   rates(sorted) : " << file_lin << endl;
   if (Input_::_Y_lsj) cerr << "   rates(mult)   : " << file_tr_aux << endl
                            << "END" << endl;
}

const char* Print_table::atom[] = {
  "H",  "He", "Li", "Be",  "B",  "C",  "N", "O",  "F",  "Ne",
  "Na", "Mg", "Al", "Si",  "P",  "S", "Cl", "Ar", "K",  "Ca",
  "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",  "Y", "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
  "Sb", "Te",  "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
  "Lu", "HF", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ar", "Ac", "Th",
  "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
  "Md", "No", "Lr", "Unq", "Unp", "Unh", "Uns"};

const int Print_table::sa = sizeof(atom)/sizeof(*atom);

const char* Print_table::roman[] = { 
   "I", "II", "III", "IV", "V", "VI", "VII", "VIII",
   "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII",
   "XIX", "XX", "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI", "XXVII",
   "XXVIII", "XXIX", "XXX", "XXXI", "XXXII", "XXXIII", "XXXIV", "XXXV",
   "XXXVI", "XXXVII", "XXXVIII", "XXXIX","39+","40+","41+","42+","43+",
    "44+","45+","46+","47+","48+","49+","50+","51+","52+","53+","54+",
    "55+","56+","57+","58+","59+","60+","61+","62+","63+","65+",
    "65+","66+","67+","68+","69+","70+","71+","72+","73+","74+",
    "75+","76+","77+","78+","79+","80+","81+","82+","83+","84+",
    "85+","86+","87+","88+","89+","90+","91+","92+","93+","94+",
    "95+","96+","97+","98+","99+","100+","101+","102+","103+","104+",
    "105+","106+","107+"};

const int Print_table::sr = sizeof(roman)/sizeof(*roman);

