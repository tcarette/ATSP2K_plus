#include<sys/types.h>
#include<dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include "repr.hh"
#include "rep_lib.hh"
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <strstream>
#include <algorithm>
#include <list>
#include <cmath>
#include <functional>
#include <iomanip>


double Rep::RZ(int& Z) {
//compute Rydberg constant
   double ZZ(Z);
   double zmu = 0;
   if (ZZ == 1) zmu = 1;
   else if ( ZZ > 10) {zmu = 2 * ZZ + 1 + (ZZ - 11)/2;}
   else if (!(int(ZZ)%2) || (ZZ == 7)) {zmu = 2 * ZZ;}
   else {zmu = 2 * ZZ + 1;}
   return RZN/(1+RZD/zmu);
}

float conv_J (string& s) {
   float J;  
   char* t = &s[0];
   std::istrstream str(t);
   str >> J;
   string::size_type i = s.find_first_of("/");
   if (i != string::npos) J /= 2;
   return J;
}

Rep::LSJ_Rep::LSJ_Rep(string& sm) {
   int a1 = sm.find_first_of("&");
   string tmp = sm.substr(0,a1-1);
   conf_mchf = ret_str(tmp);
   int a2 = sm.find_first_of("&",a1+1);
   tmp = sm.substr(a1+1,a2-a1-1);
   term_mchf = ret_str(tmp);
   int a3 = sm.find_first_of("&",a2+1);
   tmp = sm.substr(a2+1,a3-a2-1);
   j_string_mchf = ret_str(tmp);
   JV = conv_J(j_string_mchf);
   int a4 = sm.find_first_of("&",a3+1);
   tmp = sm.substr(a3+1,a4-a3-1);
   parity_mchf = ret_str(tmp);
   int a5 = sm.find_first_of("&",a4+1);
   tmp = sm.substr(a4+1,a5-a4-1);
   conf_nist = ret_str(tmp);
   int a6 = sm.find_first_of("&",a5+1);
   tmp = sm.substr(a5+1,a6-a5-1);
   term_nist = ret_str(tmp);
   int a7 = sm.find_first_of("&",a6+1);
   tmp = sm.substr(a6+1,a7-a6-1); 
   j_string_nist = ret_str(tmp);
   int a8 = sm.find_first_of("&",a7+1);
   tmp = sm.substr(a7+1,a8-a7-1);
   conf_tex = ret_str(tmp);
   tmp = sm.substr(a8+1);
   term_tex = ret_str(tmp);
//   cout << "::" << conf_mchf <<"::"<<term_mchf<<"::"<<j_string_mchf<<
//        "::"<<parity_mchf<<"::"<<conf_nist<<"::"<<term_nist<<"::"<<
//            j_string_nist<<"::"<<conf_tex <<"::"<< term_tex<<"::"<<endl;
}

Rep::LS_Rep::LS_Rep(string& st){
   int T = 3;
   Ugroup = st.substr(0,2);
   int a1 = st.find_first_of("&");
   string tmp = st.substr(T,a1-1-T);
   conf = ret_str(tmp);
   int a2 = st.find_first_of("&",a1+1);
   tmp = st.substr(a1+1,a2-a1-1);
   term = ret_str(tmp); 
   int a3 = st.find_first_of("&",a2+1);
   tmp = st.substr(a2+1,a3-a2-1);
   parity = ret_str(tmp);
   int a4 = st.find_first_of("&",a3+1);
   tmp = st.substr(a3+1,a4-a3-1);
   tex_C = ret_str(tmp);
   tmp = st.substr(a4+1);
   tex_T = ret_str(tmp);
//   cout << "::"<< conf << "::" << term << "::" << tex_C << "::" << tex_T << "::" <<endl;
}

Rep::LS_lev::LS_lev(string& s1, int ZZ, int nn) : n(nn), Z(ZZ) {
//               10        20        30        40        50        60
//   ::0123456789_123456789_123456789_123456789_123456789_123456789_1234567890
//s0 :: Z =  10 n =  4
//s1 ::   7 -121.84352565  2s(2).2p(2)3P2_3P.3p_4D
//s2 ::   7 -121.58145918  2s(2).2p(2)3P2_3P.3d_4D
//s3 ::   57515.36 CM-1      1738.67 ANGS(VAC)      1738.67 ANGS(AIR)
//s4 :: E1  length:   S =  1.02663D+01   GF =  1.79359D+00   AKI =  1.97881D+08
//s5 ::    velocity:  S =  8.50566D+00   GF =  1.48599D+00   AKI =  1.63944D+08

// read line s1

   string tmp = s1.substr(20);;
   Conf = ret_str(tmp);

   char* t = &s1[4];
   istrstream is1(t);
   is1 >> E;
//   t = &s1[20];
//   istrstream is2(t);
//   is2 >> tmp;
}

Rep::Rep(const char* argv) : RZN(109737.31534),RZD(548.697e-06) {

   map<int,map<LS_lev*,double> > VM_lev;

//#
    DIR* dd = opendir(".");
    cout << "...reading all .ls files in the current directory\n";
    if (!dd) {
        perror("Error opening current directory \"./\"!\n");
        exit(1); //directory is not readable
    }

// read all entries
    while (dirent* dp = readdir(dd)) {
       string read_file(dp->d_name);
       int fnl = read_file.size();

       if (fnl > 3) {  //read only files with at least 4 char
          struct stat statv;
          stat(read_file.c_str(), &statv);
          if (&statv && statv.st_size > 100) { // do not read files < 100 B
             if (read_file.substr(fnl-3) == ".ls") { // read only .ls files
                char t_read[sz];
                ifstream inf(read_file.c_str());
                while(inf.getline(t_read,sz)) {
                   if(t_read[1] == 'Z') {
                      string s0(t_read);
                      char* t = &s0[5];
                      istrstream is1(t);
                      int ZZ;
                      is1 >> ZZ;
                      t = &s0[13];
                      istrstream is2(t);
                      int nn;
                      is2 >> nn;
             
                      inf.getline(t_read,sz);
                      string s1(t_read);
                      LS_lev* A = new LS_lev(s1,ZZ,nn);
                      Map_LEV[ZZ][nn].push_back(new Lev(A));
                    
                      VI p1 = find_if 
                           (Map_LZ[ZZ].begin(),Map_LZ[ZZ].end(),match_Lev(A));
                      if (p1==Map_LZ[ZZ].end())
                         { 
                            Map_LZ[ZZ].push_back(A); 
                         }
                      else if ((*p1)->E > A->E) {
                         (*p1)->E = A->E;
                         (*p1)->n = A->n;
                      }
             
                      inf.getline(t_read,sz);
                      string s2(t_read);
                      LS_lev* B = new LS_lev(s2,ZZ,nn);
                      Map_LEV[ZZ][nn].push_back(new Lev(B));
             
                      VI p2 = find_if
                           (Map_LZ[ZZ].begin(),Map_LZ[ZZ].end(),match_Lev(B));
                      if (p2==Map_LZ[ZZ].end())
                         { 
                            Map_LZ[ZZ].push_back(B); 
                         }
                      else if ((*p2)->E > B->E) {
                          (*p2)->E = B->E;
                          (*p2)->n = B->n;
                      }
             //cout << ZZ << " " << nn << s1 << endl;
             //cout << ZZ << " " << nn << s2 << endl;
                   }
                }
                cout << "...reading file " << read_file << endl;
             }
          }
       }
    }

    closedir(dd);
//#


   for (MI p = Map_LZ.begin(); p != Map_LZ.end(); p++) {
      int i = 0;
//      cout << "size of = " << (*p).second.size() << "::" << i<< endl;
      sort((*p).second.begin(),(*p).second.end(),Comp_LS_lev());
      for (VI pi = (*p).second.begin(); pi != (*p).second.end(); pi++) {
         (*pi)->ind = i++;
//         cout << (*pi)->ind << "::" << (*pi)->Z <<"::"<< (*pi)->n << (*pi)-> E << "::" << (*pi)->Conf << "::" << endl;
      }
   }
   cout << "...printing spectra in _A.log\n";
   comp_spec();
}


bool Lev_low (const Rep::Lev* A, const Rep::Lev* B) {
   return A->e < B->e;
}

void printH1(ostream& os, int Z, int n) {
   os << endl << endl << "Table Z = ";
   os << Z;
   os << ";  n = ";
   os << n << endl;
   string f(80,'-');
   os << f << endl;
//  30. 2s(2).2p(6).3p_2P.4p_1D   -1068.18597279         -1947121.54
//  31. 2s(2).2p(6).3p_2P.4p_1S   -1068.05105764         -1976731.71
   os << "    " << endl;
   os << "      Level" << setw(40);
   os << "Total Energy" << setw(20);
   os << endl << f << endl;
   os << endl ;
}

void printF(ostream& os) {
   string f(80,'-');
   os << endl << f << endl << endl;
}

void printH0(ostream& os) {
   os << "Table of LS energies present in the .ls file"
         " for each Z and n\n";
   os << endl << endl;
}

ostream& operator<<(ostream& os, Rep::Lev* A) {
   int cfs = 40 - A->c.size();
   os.setf(ios::left,ios::adjustfield);
   os << A->c << setw(cfs);
   os.setf(ios::fixed, ios::floatfield);
   os.setf(ios::right,ios::adjustfield);
   os.precision(8);
   os << A->e << setw(20);
   os.precision(2);
   os << A->ee << setw(15);
   return os;
}

void Rep::comp_spec(){
   typedef map<int,map<int,vector<Rep::Lev*> > >::iterator MZ;
   typedef map<int,vector<Rep::Lev*> >::iterator MN;
   typedef vector<Rep::Lev*>::iterator ML;

   ofstream FLOG("_A.log");
   cout << "...Look in file _A.log for printouts of the spectra\n";
   printH0(FLOG);
   for (MZ pz = Map_LEV.begin(); pz != Map_LEV.end(); pz++) {
      for (MN  pn = pz->second.begin();  
               pn != pz->second.end(); pn++) {
         // sort within Z and n ranges 
         sort(pn->second.begin(),pn->second.end(),CompLev());         
         vector<Lev*>::iterator 
               p = unique(pn->second.begin(),pn->second.end(),LevEq());
         pn->second.erase(p,pn->second.end());
//get the ground energy
         double GE = pn->second[0]->e;
//get Z
         int Z(pz->first);
         int i = 0;
         printH1(FLOG,pz->first,pn->first);
         for (ML pl = pn->second.begin(); 
                  pl != pn->second.end(); pl++) {
// compute the level
            (*pl)->ee = -(GE - (*pl)->e)*2*RZ(Z); 
            FLOG << setw(4) << i++ << ". " << *pl << endl;
//            cout << i++ <<"::" << (*pl)->z 
//                        <<"::" << (*pl)->n 
//                        <<"::" << (*pl)->c 
//                        <<"::" << (*pl)->e 
//                        <<"::" << (*pl)->ee << endl;
         }
         printF(FLOG);
      }
   }
}



Rep::Rep() {
   ifstream in("LEVELS.ref");
   static char t_read[sz];
   while (t_read[0] != 'M') {
      in.getline(t_read,sz);
      string s(t_read);
   }

   int i = 0;  // read MCHF  rep
   while (in.getline(t_read,sz)) {
      string s(t_read);
      if (t_read[0] == 'L') break;
      if(s.size() >= 50) {
         LSJ_Rep* m = new LSJ_Rep(s);
         m->ind = i++;
         vec_LSJ_Rep.push_back(m);
      }
   }

   i = 0; // read LS rep
   while (in.getline(t_read,sz)) {
      string s(t_read);
      if(s.size() > 50) {
         LS_Rep* o = new LS_Rep(s);
         o->ind = i++;
         vec_LS_Rep.push_back(o);
      }
   }
}

void Rep::print() {
   CI_LSJ_Rep pj;
   CI_LS_Rep ps;
   cout << "Size of LSJ_Rep = " << vec_LSJ_Rep.size() << endl;
   cout << "Size of LS_Rep = " << vec_LS_Rep.size() << endl;
   for (pj = vec_LSJ_Rep.begin(); pj !=vec_LSJ_Rep.end(); pj++) {
     cout << (*pj)->ind << "::" << (*pj)->JV << "::" <<
             (*pj)->conf_mchf << "::" << (*pj)->term_mchf << "::" <<
             (*pj)->j_string_mchf << "::" << (*pj)->parity_mchf << "::" <<
             (*pj)->conf_nist << "::" << (*pj)->term_nist << "::" <<
             (*pj)->j_string_nist << "::" << (*pj)->conf_tex << "::" <<
             (*pj)->term_tex << "::" << endl;
 }
   for (ps = vec_LS_Rep.begin(); ps !=vec_LS_Rep.end(); ps++) {
//     cout << (*ps)->ind << "::" << (*ps)->conf << "::" << 
//             (*ps)->term << "::" << (*ps)->parity << "::" <<
//             (*ps)->tex_C << "::" << (*ps)->tex_T << "::" << endl;
   }
}

