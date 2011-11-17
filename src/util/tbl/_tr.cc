#include "_te.hh"
#include "_tt.hh"
#include "_prn.hh"
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
#include <cmath>
#include <sstream>
#include <ctime>


void appendOdd(string&);

//----------------------------------------------------
//-----  non-members included in namespace transition   
//-----  also some of them called from members         
//----------------------------------------------------
void Dt::ST::convto_LSJ(list<JTR_tr*>& T, double ge) {

typedef vector<string>::reverse_iterator RVS;
   int i = 0;
   string SS[5] = {""};
   for (RVS p = VS.rbegin(); p != VS.rend(); p++) {
       SS[i++] = (*p);
       if (i > 4) break;
   }
   
   //if (!check_empty()) throw;

   double d2(0),d4(0),d6(0);
//cout << "SS[0]" << SS[0] << endl;
   if (SS[0].size() > 0)  {
      // read velocity values
      string::size_type i = 0;
      i = SS[0].find("D");
      if (i != string::npos) {
         SS[0].replace(i,1,"e");
         i = SS[0].find("D",i+1);
         if (i != string::npos) {
             SS[0].replace(i,1,"e");
             i = SS[0].find("D",i+1);
             if (i != string::npos) SS[0].replace(i,1,"e");
          }
      }
      if (i == string::npos) {
         i = SS[0].find("E",5);
         if (i != string::npos) {
           SS[0].replace(i,1,"e");
           i = SS[0].find("E",i+1);
           if (i != string::npos) {
               SS[0].replace(i,1,"e");
               i = SS[0].find("E",i+1);
               if (i != string::npos) SS[0].replace(i,1,"e");
            }
         }
      }
      istringstream ist(SS[0]);
      ist >> d2 >> d4 >> d6;
   }

   string t;
   double d1(0),d3(0),d5(0);
//cout << "SS[1]" << SS[1] << endl;
   if (SS[1].size() > 0) {
     // read length values;
      string::size_type i(0),i1(0);
      i = SS[1].find("D");
      if (i != string::npos) {
         SS[1].replace(i,1,"e");
         i = SS[1].find("D",i+1);
         if (i != string::npos) {
            SS[1].replace(i,1,"e");
            i = SS[1].find("D",i+1);
            if (i != string::npos) SS[1].replace(i,1,"e");
         }
      }
      if (i == string::npos) {
         i = SS[1].find("E",5);
         if (i != string::npos) {
            SS[1].replace(i,1,"e");
            i = SS[1].find("E",i+1);
            if (i != string::npos) {
                SS[1].replace(i,1,"e");
                i = SS[1].find("E",i+1);
                if (i != string::npos) SS[1].replace(i,1,"e");
             }
         }
      }
      istringstream ist(SS[1]);
      string _s;
      ist >> t >> _s >> _s >> d1 >> _s >> _s >> d3 >> _s >> _s >> d5;
   } else return;

   double d7(0),d8(0),d9(0);
//cout << "SS[2]" << SS[2] << endl;
   if (SS[2].size() > 0) {
     // read WL;
      istringstream ist(SS[2]);
      string _s;
      ist >> d9 >> _s >> d8 >> _s >> d7;
   } else return;

   int i1;
   double d10(0);
   string c1;
//cout << "SS[3]" << SS[3] << endl;
   if (SS[3].size() > 0) {
      istringstream ist(SS[3]);
      string _s;
      ist >> i1 >> d10 >> c1;
   } else return;

   int i2;
   double d11(0);
   string c2;
//cout << "SS[4]" << SS[4] << endl;
   if (SS[4].size() > 0) {
      istringstream ist(SS[4]);
      string _s;
      ist >> i2 >> d11 >> c2;
   } else return;

   if (d11 > d10) {
cerr << " Warning! Lines:\n" << SS[4] << SS[3] << endl;
cerr << "An interchanged tranisiton (energies)"
        ", do you want to continue enter \"y\": ";
string cont;
getline(cin,cont,'\n');
if (cont[0] == 'y') {
    cerr << " You have entered " << cont
         << " the program continues\n";
} else {
   cerr << " tables exited on lines:\n";
   cerr << SS[4] << SS[3] << endl;
   exit(1);
}

      double dtmp = d10;
      string stmp = c1;
      int itmp = i1;
      i1 = i2;
      c1=c2;
      d10=d11;
      d11=dtmp;
      c2=stmp;
      i2 = itmp;
   }

d8 = fabs(d8);
d7 = fabs(d7);
d1 = fabs(d1);
d2 = fabs(d2);
d3 = fabs(d3);
d4 = fabs(d4);
d5 = fabs(d5);
d6 = fabs(d6);

   JTR_tr* J = new JTR_tr(d10,d11,d9,c1,c2,i1,i2,
                    d8,d7,t,d1,d2,d3,d4,d5,d6,ge);
   if (d5) T.push_back(J);
//cout << "T.size() : " << T.size() << endl;
}

void expected_format(string& L) {
cerr << " ... Found string ionization, but the format is incorrect, the line reads:n";
cerr << L << endl;
cerr << " expected format is :\n";
cerr << " expected format.........\n";
}
  
//  reading of *.lsj files, this function calls a nonmember rdstr_lsj
void Dt::Base_table::read_lsj(const string& fn) {
   string s = *Input_::_filebase_ + ".lsj";  // added

   appendOdd(s);
   ofstream FLOG("tables.log");
   ifstream in_lsj(s.c_str());  
   if(in_lsj) cout << "...reading " << s << "...\n"; 
   else {
      cout << s << " not found ..." << endl; 
      return;
   }
//
   const int sz = buf;
   static char t_read[sz];

   double GE_read(0);
   if(in_lsj.getline(t_read,sz)) {
      string G(t_read);
      
      if (G.substr(0,16) == " Ground Energy =") {
         istringstream ist(G.substr(16,21));
         ist>>GE_read;
         Base_table::E_min = GE_read; // set the Min energy to what is in the .lsj;
cout <<  " ... found Ground Energy = " << GE_read  << endl;
      }
   }

   IE_limit = 0;
   list<string*> SL;
   int no_transitions = 0;
   bool _M_ = false;
   while(in_lsj.getline(t_read,sz)) {
      string READ_LINE(t_read);

      string::size_type a = READ_LINE.find("Ionization");
      if (a != string::npos) {
          string::size_type b1, b2, b3;
          b1 = 6 + READ_LINE.find("Term ="); // +7 positions
          b2 = 3 + READ_LINE.find("J =");    // +4 positions
          b3 = 8 + READ_LINE.find("Energy =");  // +9 postions

          if (b1 != string::npos && b2 != string::npos && 
                  b3 != string::npos ) {
             if (b1 < b2 && b2 < b3 && b3 < READ_LINE.size()) { 
                string sT = READ_LINE.substr(b1+1,b2-b1-1);
                string sJ = READ_LINE.substr(b2+1,b3-b2-1);
                string sE = READ_LINE.substr(b3+1);
                istringstream ist1(sT);
                string::size_type f1 = sT.find("NA");
                ist1 >> termstr_LIMIT;
                istringstream ist2(sJ);
                string::size_type f2 = sJ.find("NA");
                ist2 >> jstr_LIMIT;
                istringstream ist3(sE);
                string::size_type f3 = sE.find("NA");
                if (f3 == string::npos) {
                   double d(0);
                   ist3 >> IE_limit;
                } else {
                   cerr << "... The energy for Ionization LIMIT"
                           " was \"NA\", set it to proper value, replacing \"NA\":\n";
                }
             } else expected_format(READ_LINE);
          }  else expected_format(READ_LINE);
      }

      if (_M_) {
        string* s = new string(READ_LINE);
        if (!s->size()) {
//cout << ">>" << *s << "<<" << endl;
           SL.push_back(s);
           _M_ = false;
        }
      }
      string* A = new string(READ_LINE);
//cout << "::" << READ_LINE << "::" << endl;
      SL.push_back(A);
      if (READ_LINE.size() > 6) {
         if (READ_LINE[5] == 'S') {
            no_transitions++;
            _M_ = true;
         }
      }
      READ_LINE = "";
   }
 
   if (!IE_limit)  {
      cerr << " No Ionization Limit (IL) was found:\n"
              " To include an IL, add the following line after the line with Ground Energy =...\n"
              "   Ionization Limit = ......\n\n";
   }
               
   if (no_transitions) {

   } else if (SL.size()) {

   } else return;

   typedef list<string*>::const_iterator LSI;

   ST* t = new ST;

   bool NT=false;
   bool inserted = false;
   for (LSI p = SL.begin(); p != SL.end(); p++) {
     string AA(*(*p));
     inserted = false;

     if (NT) {
       t->VS.push_back(AA);
       t->convto_LSJ(this->lsj,GE_read);
       t->VS.clear();
       NT = false;
       inserted = true;
     }

     if (AA.size() && !inserted) {
       t->VS.push_back(AA);
       if (AA[5] == 'S') NT = true; 
     }
    
   }

//   for (LSI p = SL.begin(); p != SL.end(); p++) {
//     string AA(*(*p));
//
//     t->VS.push_back(AA);
//     if (AA[5] == 'S') {
//        ++p;
//        if (p != SL.end()) {
//           AA = *(*p);
//        } else {
//           AA = "";
//        } 
//        t->VS.push_back(AA);
//        try {
//           t->convto_LSJ(this->lsj,GE_read);
//        } catch (Error_lsj& E) {
//          E.print_error(t);
//        }
//        t->VS.clear();
//        if (p == SL.end()) break;
//     } 
//   }
//
//   if (no_transitions == 0) {
//
//   }

   delete t;
   SL.clear();

}


