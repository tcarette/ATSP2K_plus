#include "lev.hh"
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



float conv_J (string& s) {
   float J;  
   char* t = &s[0];
   std::istrstream str(t);
   str >> J;
   string::size_type i = s.find_first_of("/");
   if (i != string::npos) J /= 2;
   return J;
}

Lev::LS_lev::LS_lev(string& s1, int nn, int ZZ) : n(nn), Z(ZZ) {
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
   Conf_L = ret_str(tmp);
   char* t = &s1[5];
   istrstream is1(t);
   is1 >> E;
//   t = &s1[20];
//   istrstream is2(t);
//   is2 >> tmp;
}


Lev::LSJ_Lev::LSJ_Lev(string& sm) {
   int T = 3;
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

Lev::LS_Lev::LS_Lev(string& st){
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

Lev::Lev() {
   ifstream inf("LEVELS.ref");
   map<int,map<LS_lev*,double>> VM_lev;
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

         inf.getline(t_read,sz);
         string s2(t_read);
         LS_lev* A = new LS_lev(s2,ZZ,nn);

      }
   }
}

void Lev::print() {
   CI_LSJ_Lev pj;
   CI_LS_Lev ps;
   cout << "Size of LSJ_Lev = " << vec_LSJ_Lev.size() << endl;
   cout << "Size of LS_Lev = " << vec_LS_Lev.size() << endl;
   for (pj = vec_LSJ_Lev.begin(); pj !=vec_LSJ_Lev.end(); pj++) {
     cout << (*pj)->ind << "::" << (*pj)->JV << "::" <<
             (*pj)->conf_mchf << "::" << (*pj)->term_mchf << "::" <<
             (*pj)->j_string_mchf << "::" << (*pj)->parity_mchf << "::" <<
             (*pj)->conf_nist << "::" << (*pj)->term_nist << "::" <<
             (*pj)->j_string_nist << "::" << (*pj)->conf_tex << "::" <<
             (*pj)->term_tex << "::" << endl;
 }
   for (ps = vec_LS_Lev.begin(); ps !=vec_LS_Lev.end(); ps++) {
//     cout << (*ps)->ind << "::" << (*ps)->conf << "::" << 
//             (*ps)->term << "::" << (*ps)->parity << "::" <<
//             (*ps)->tex_C << "::" << (*ps)->tex_T << "::" << endl;
   }
}

