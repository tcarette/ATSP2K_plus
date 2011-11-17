#include "ls.hh"
#include "repr.hh"
#include "rep_lib.hh"
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <strstream>
#include <algorithm>
#include <functional>
#include <list>
#include <cmath>
#include<iomanip>
#include<sys/types.h>
#include<dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>

#include "dat.hh"

const string sA(80,'-');
const string sB1(" Z  n       EL            EU        ");
const string sB2("S(L)      S(V)      gf(L)    gf(V)   Error");
const string LS_dat_hdr = sA + '\n' + sB1 + sB2 + '\n' + sA + '\n';
const string LS_dat_ft = sA + '\n';

LS::LS_tr::LS_tr(string& s1,string& s2, string& s3, 
                 string& s4,string& s5, int ZZ, int nn, 
                            map<string,int>& LS_ind) : Z(ZZ), n(nn) {
//               10        20        30        40        50        60
//   ::0123456789_123456789_123456789_123456789_123456789_123456789_1234567890
//s0 :: Z =  10 n =  4
//s1 ::   7 -121.84352565  2s(2).2p(2)3P2_3P.3p_4D
//s2 ::   7 -121.58145918  2s(2).2p(2)3P2_3P.3d_4D
//s3 ::   57515.36 CM-1      1738.67 ANGS(VAC)      1738.67 ANGS(AIR)
//s4 :: E1  length:   S =  1.02663D+01   GF =  1.79359D+00   AKI =  1.97881D+08
//s5 ::    velocity:  S =  8.50566D+00   GF =  1.48599D+00   AKI =  1.63944D+08

   s4[27] = 'e';
   s4[47] = 'e';
   s4[68] = 'e';
   s5[27] = 'e';
   s5[47] = 'e';
   s5[68] = 'e';

// read line s1

   string tmp = s1.substr(20);;
   Conf_L = ret_str(tmp);
   char* t = &s1[4];
   istrstream is1(t);
   is1 >> EL;
//   t = &s1[20];
//   istrstream is2(t);
//   is2 >> tmp;
map<string,int>::iterator posI = LS_ind.find(Conf_L);
if (posI != LS_ind.end()) ind_i = LS_ind[Conf_L];  //use map
//   map<string,int>::iterator posI = LS_ii[Z].find(Conf_L);
//   if (posI != LS_ii[Z].end()) ind_i = LS_ii[Z][Conf_L];  //use map
   else ind_i = -1;

//read line s2
   tmp = s2.substr(20);
   Conf_U = ret_str(tmp);
   t = &s2[4];
   istrstream is3(t);
   is3 >> EU;
//   t = &s2[20];
//   istrstream is4(t);
//   is4 >> tmp;
map<string,int>::iterator posF = LS_ind.find(Conf_U);
if(posF != LS_ind.end()) ind_f = LS_ind[Conf_U];
//   map<string,int>::iterator posF = LS_ii[Z].find(Conf_U);
//   if(posF != LS_ii[Z].end()) ind_f = LS_ii[Z][Conf_U];
   else ind_f = -1;
 
//read line s3
   t = &s3[0];
   istrstream is5(t);
   is5 >> LS_trE;
   t = &s3[16];
   istrstream is6(t);
   is6 >> w_VAC;
   t = &s3[39];
   istrstream is7(t);
   is7 >> w_AIR;

//read line s4
   TR_type = s4.substr(1,2); 
   t = &s4[20];
   istrstream is8(t);
   is8 >> SL_L;
   t = &s4[40];
   istrstream is9(t);
   is9 >> GF_L;
   t = &s4[61];
   istrstream is10(t);
   is10 >> AKI_L;

//readline s5
   t = &s5[20];
   istrstream is11(t);
   is11 >> SL_V;
   t = &s5[40];
   istrstream is12(t);
   is12 >> GF_V;
   t = &s5[61];
   istrstream is13(t);
   is13 >> AKI_V;
   double Max_S = (SL_L <= SL_V) ? SL_V : SL_L;
   if (SL_L-SL_V) perc_LS_error = fabs((SL_L-SL_V)/Max_S);
   else perc_LS_error = 0;

//cout << "::"<<n<<"::"<<Conf_L<<"::"<<Conf_U<<"::"<<EL<<"::"<<EU<<"::"
//           <<LS_trE<<"::"<<w_VAC<<"::"<<w_AIR<<"::"<<TR_type<<"::"
//           <<SL_L<<"::"<<GF_L<<"::"<<AKI_L<<"::"<<SL_V<<"::"
//           <<GF_V<<"::"<<AKI_V<<"::"<<perc_LS_error<<"::"<<Z<<"::" 
//           <<ind_i<<"::"<<ind_f<<"::"<<endl;

}

LS::LS(const char* argv) : R(new Rep(argv)) {
     set_SEQ();
   //ofstream FLOG("_A.log");
   //ifstream inf("A.LS");
   //ifstream inf(argv);
   //if (!inf) {
  //    cout << argv << " not found, exiting...\n";
  //    exit (2);
  // }

   map<int,vector<LS_tr*> > ALL_Z_LS;
   map<string,int> LS_ind;

   map<int,map<string,int> >LS_ii;

   for (MI mi = R->Map_LZ.begin(); mi != R->Map_LZ.end(); mi++) {
      for (unsigned int i = 0; i < (*mi).second.size(); i++) {
         string SC = (*mi).second[i]->Conf;
         int iz = (*mi).first;
         int ii = (*mi).second[i]->ind;
         LS_ii[iz][SC] = ii;
      }
   }
   cout << "...rereading .ls" << endl;
// open the current directory
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
                      inf.getline(t_read,sz);
                      string s1(t_read);
                      inf.getline(t_read,sz);
                      string s2(t_read);
                      inf.getline(t_read,sz);
                      string s3(t_read);
                      inf.getline(t_read,sz);
                      string s4(t_read);
                      inf.getline(t_read,sz);
                      string s5(t_read);
                      char* t = &s0[5];
                      istrstream is1(t);
                      int ZZ;
                      is1 >> ZZ;
                      t = &s0[13];
                      istrstream is2(t);
                      int nn;
                      is2 >> nn;
             
                      LS_tr* A = new LS_tr(s1,s2,s3,s4,s5,ZZ,nn,LS_ii[ZZ]);
             
                      if (A->ind_i != -1 && A->ind_f != -1) {
                         // using map of Z and vector transitions
                         ALL_Z_LS[ZZ].push_back(A);  
                      }
                      else if (A->ind_i == -1) {
                      }
                      else if (A->ind_f == -1)  {
                      }
                   }
                }
             }
          }
       }
    }

   typedef map<int,vector<LS_tr*> >::iterator ZLS;

   cout << "...Sorting LS" << endl;
   for (ZLS pi = ALL_Z_LS.begin(); pi != ALL_Z_LS.end(); pi++) {
       sort(pi->second.begin(), pi->second.end(), Comp_LS_ind());
       for (unsigned int i = 0; i < pi->second.size(); i++) {
//cout << "::"<<pi->second[i]->Z<< "::"<<pi->second[i]->n<<"::";
//cout << pi->second[i]->Conf_L << "=" << pi->second[i]->EL << "::"; 
//cout << pi->second[i]->Conf_U << "=" << pi->second[i]->EU <<"::"<< endl;
       }
//cout << pi->first << endl;
   }
   nZ = ALL_Z_LS.size();
   ZLS pi = ALL_Z_LS.begin();
   Zdiff = pi->first;
   cout << "...LS spectrum for each n and Z in file _A.log\n";
   cout << "...LS errors in files: " << argv << "_Z_nn" << endl;

   if (nZ) {
      LS_print_dat(ALL_Z_LS,argv);
   } 
   else {
     cout << "...no files with LS data found!... The program uses this "
          "format:\n\n";
cout<<"  Transition between files:\n";
cout<<"  E   # initial state     \n";
cout<<"  O   # final state       \n";
cout<< endl;
cout<<" Z =  12 n =  4\n";
cout<<"   4 -199.79148383  2s(2).2p(6).3s_2S.3p_3P                         \n";
cout<<"   2 -199.46542750  2s(2).2p(6).3s_2S.4s_3S                         \n"; 
cout<<"   71559.46 CM-1      1397.44 ANGS(VAC)      1397.44 ANGS(AIR)\n";
cout<<" E1  length:   S =  2.48666D+01   GF =  5.40515D+00 "
      "  AKI =  6.15408D+09\n";
cout<<"    velocity:  S =  1.40851D+00   GF =  3.06162D-01 "
      "  AKI =  3.48583D+08\n\n";
    }

//   cout << "...computing LS errors " << endl;
//   LS_mk_Err(ALL_Z_LS);
   //cout << "...LS errors in .tex in \"_A.ls_acc.tex\""<< endl;
//   cout << "...LS errors in ascii in \"_A.log\"" << endl; 
   //LS_print_tex();
   delete R;
}

void LS::LS_print_dat(map<int,vector<LS_tr*> >& ALL_Z_LS, const char* argv) {
   typedef map<int,vector<LS_tr*> >::iterator ZLS;
   for (ZLS pi = ALL_Z_LS.begin(); pi != ALL_Z_LS.end(); pi++) {
       sort(pi->second.begin(), pi->second.end(), Comp_LS_ind());
// Leave only unique transitions: Z, n, Eupper, Elower, ConfU, ConfL diff.

       vector<LS_tr*>::iterator
               p = unique(pi->second.begin(),pi->second.end(),LS_trEq());
       pi->second.erase(p,pi->second.end());

       char buf[15];
       ostrstream ost(buf,sizeof(buf));
       int zs = pi->first;
       string sA;
       ost << argv << "_Z_" << zs << '\0';
       //ofstream out(buf);
   
       ostringstream OA(sA); 
       OA << "LS_" << SEQ << "_Z_" << zs;
       string fn = OA.str();
       ofstream out(fn.c_str());

       //streambuf* SB = cout.rdbuf();
       //cout.rdbuf(out.rdbuf());
       string si;
       string sf;
       out << LS_dat_hdr;
       cout << LS_dat_hdr;
       for (int i = 0; i < pi->second.size(); i++) {
          if (si != pi->second[i]->Conf_L || sf != pi->second[i]->Conf_U) {
              si = pi->second[i]->Conf_L;
              sf = pi->second[i]->Conf_U;
              int W1 = 30 - si.size();
              int W2 = 30 - sf.size();
              out << endl << setw(W1) << si << " " << setw(W2) 
                          <<  sf << " " << endl;
              cout << endl << setw(W1) << si << " " << setw(W2)
                          <<  sf << " " << endl;

          }
          cout.precision(7);
          cout.setf(ios::showpoint);
          cout.setf(ios::fixed, ios::floatfield);
          cout << " "<<pi->second[i]->Z<< " "<<pi->second[i]->n<<" ";
          cout << pi->second[i]->EL << " " << pi->second[i]->EU <<" ";
          cout.precision(3);
          cout.setf(ios::showpoint);
          cout.setf(ios::scientific, ios::floatfield);
          cout << pi->second[i]->SL_L << " " << pi->second[i]->SL_V << " ";
          cout << pi->second[i]->GF_L << " " << pi->second[i]->GF_V << " ";
          cout.precision(3);
          cout.setf(ios::showpoint);
          cout.setf(ios::fixed,ios::floatfield);
          cout << pi->second[i]->perc_LS_error << endl;

          out.precision(7);
          out.setf(ios::showpoint);
          out.setf(ios::fixed, ios::floatfield);
          out << " "<<pi->second[i]->Z<< " "<<pi->second[i]->n<<" ";
          out << pi->second[i]->EL << " " << pi->second[i]->EU <<" ";
          out.precision(3);
          out.setf(ios::showpoint);
          out.setf(ios::scientific, ios::floatfield);
          out << pi->second[i]->SL_L << " " << pi->second[i]->SL_V << " ";
          out << pi->second[i]->GF_L << " " << pi->second[i]->GF_V << " ";
          out.precision(3);
          out.setf(ios::showpoint);
          out.setf(ios::fixed,ios::floatfield); 
          out << pi->second[i]->perc_LS_error << endl;
       } 
     out << LS_dat_ft << pi->second.size() 
          << " LS transitions.\n" << LS_dat_ft;
     cout << LS_dat_ft << pi->second.size()
          << " LS transitions.\n" << LS_dat_ft;

//      cout << pi->first << endl;
   }
}

int LS::ret_ind(string& si) {
   for (int i = 0; i < R->vec_LS_Rep.size(); i++) {
      string str_ref = R->vec_LS_Rep[i]->conf + '_' + R->vec_LS_Rep[i]->term;
      if ( si == str_ref) return R->vec_LS_Rep[i]->ind;
   }
   return -1;
}

void LS::LS_mk_Err(map<int,vector<LS_tr*> >& ALL_Z_LS) {
   nZ = ALL_Z_LS.size();
   int ZZ = nZ;
   typedef map<int,vector<LS_tr*> >::iterator ZLS;
   ZLS pn = ALL_Z_LS.begin();
   for (ZLS pi = ALL_Z_LS.begin(); pi != ALL_Z_LS.end(); pi++) {
      int indv = pi->first - Zdiff;  // vector index = Z - Z neutral 
      for (int i = 0; i < pi->second.size(); i++) {
         int ind_I = ret_ind(pi->second[i]->Conf_L);   // return index of conf Lower
         int ind_F = ret_ind(pi->second[i]->Conf_U);   // return ind of conf Upper
         if (ind_I != -1 || ind_F != -1) { 
            int NN = pi->second[i]->n;
            CI_LSerr pe = find_if(LS_ERR.begin(),LS_ERR.end(),Find_ind_LS_tr(ind_I,ind_F)); 
            if (pe == LS_ERR.end()) {  // insert a new transition
               LS_err* A = new LS_err(pn->first,ZZ,pi->second[i],R);
               A->vec_Z[indv] = pi->first;
               A->vec_N[indv] = pi->second[i]->n;
               A->vec_Perc_err[indv] = pi->second[i]->perc_LS_error;
               A->vec_GF[indv] = pi->second[i]->GF_L;
//       cout << A << endl;  // 
               LS_ERR.push_back(A);
            } else {                   // old transition;
               if(NN > (*pe)->vec_N[indv]) {
                  (*pe)->vec_N[indv] = NN;
                  (*pe)->vec_Perc_err[indv] = pi->second[i]->perc_LS_error; 
                  (*pe)->vec_GF[indv] = pi->second[i]->GF_L;
                  int n1 = (*pe)->vec_N[indv];
                  int n2 = pi->second[i]->n;
                  float f1 = (*pe)->vec_Perc_err[indv];
                  float f2 = pi->second[i]->perc_LS_error;
//       cout << *pe << endl;  //
               } else if(NN == (*pe)->vec_N[indv]) {
                  if ((*pe)->vec_Perc_err[indv] < pi->second[i]->perc_LS_error) {
                     (*pe)->vec_Perc_err[indv] = pi->second[i]->perc_LS_error;           
                     (*pe)->vec_GF[indv] = pi->second[i]->GF_L;
                     int n1 = (*pe)->vec_N[indv];
                     float f1 = (*pe)->vec_Perc_err[indv];
                     float f2 = pi->second[i]->perc_LS_error;
//       cout << *pe << endl;  //
                  }
               }
            }
         } 
      }
   }
   cout << "... " << LS_ERR.size() << " LS transitions." << endl;
} 

LS::LS_err::LS_err(int Z,int zz,LS::LS_tr* tr,Rep* R) : ZZ(zz) {
   UgroupU = tr->UgroupU;
   UgroupL = tr->UgroupL;
   ind_i = tr->ind_i;
   ind_f = tr->ind_f;
   int i = 0;
   while (i < R->vec_LS_Rep.size() && R->vec_LS_Rep[i]->ind != ind_i) { ++i; }
   if (i < R->vec_LS_Rep.size()) {
      Lconf = R->vec_LS_Rep[i]->conf;
      Lterm = R->vec_LS_Rep[i]->term;
      Lconf_tex = R->vec_LS_Rep[i]->tex_C;
      Lterm_tex = R->vec_LS_Rep[i]->tex_T;
      UgroupL = R->vec_LS_Rep[i]->Ugroup;
   } else {
      cout << "error at:\n";
   } 
   i = 0;
   while (i < R->vec_LS_Rep.size() && R->vec_LS_Rep[i]->ind != ind_f) { ++i; }
   if (i < R->vec_LS_Rep.size()) {
      Uconf = R->vec_LS_Rep[i]->conf;
      Uterm = R->vec_LS_Rep[i]->term;
      Uconf_tex = R->vec_LS_Rep[i]->tex_C;
      Uterm_tex = R->vec_LS_Rep[i]->tex_T;
      UgroupU = R->vec_LS_Rep[i]->Ugroup;
   } else {
      cout << "error at:\n";
   }
   vec_Z.resize(ZZ,0);
   vec_GF.resize(ZZ,0);
   vec_N.resize(ZZ,0);
   vec_Perc_err.resize(ZZ,0);
}

//void LS::LS_print_tex() {
//   ofstream FLOG("_A.log", std::ios::out | std::ios::app);
//   FLOG << "Table 1:Table LS accuracies. (gf-values on the second line)\n";
//   string P(61+nZ*10,'-');
//   FLOG << P << endl;
//   FLOG << setw(60) << " ";
//   for (int i = 0; i < nZ; i++) {
//      if (Zdiff + i < 10) FLOG << setw(9) << "Z=" <<  Zdiff + i;
//      if (Zdiff + i >= 10) FLOG << setw(8) << "Z=" <<  Zdiff + i;
//   }
//   FLOG << endl;
//   FLOG << P << endl;
//   sort(LS_ERR.begin(), LS_ERR.end(), Comp_LS_err());
//   for (CI_LSerr pp = LS_ERR.begin(); pp != LS_ERR.end(); pp++) {
//      FLOG << *pp << endl;
//   }
//   FLOG << P << endl;
//
//   ofstream out("_A.ls_acc.tex");
//   out << THeader_LS << endl;
//   int cnt = 0;
//
//   string CL,CU,TL,TU;
//   for (CI_LSerr pp = LS_ERR.begin(); pp != LS_ERR.end(); pp++) {
//      cnt++;
//      string spL, spU;
//
//
//      if (CL == (*pp)->Lconf_tex && CU == (*pp)->Uconf_tex && (cnt-1)%40) {
//         spL = "       & " + (*pp)->Lterm_tex;
//         spU = "       & " + (*pp)->Uterm_tex;
//      } else if (CL == (*pp)->Lconf_tex && CU != (*pp)->Uconf_tex && (cnt-1)%40) {
//         CU = (*pp)->Uconf_tex;
//         spL = "          & " + (*pp)->Lterm_tex;
//         spU = (*pp)->Uconf_tex + " & " + (*pp)->Uterm_tex;
//      } else {
//         CL = (*pp)->Lconf_tex;
//         CU = (*pp)->Uconf_tex;
//         spL = (*pp)->Lconf_tex + " & " + (*pp)->Lterm_tex;
//         spU = (*pp)->Uconf_tex + " & " + (*pp)->Uterm_tex;
//      }
///*
//      if (CL != (*pp)->Lconf_tex || !cnt%40) {
//         CL = (*pp)->Lconf_tex;
//         CU = (*pp)->Uconf_tex;
//         spL = (*pp)->Lconf_tex + " & " + (*pp)->Lterm_tex; 
//         spU = (*pp)->Uconf_tex + " & " + (*pp)->Uterm_tex;
//      } else if (CL == (*pp)->Lconf_tex && CU != (*pp)->Uconf_tex) {
//         CU = (*pp)->Uconf_tex;
//         spL = "          & " + (*pp)->Lterm_tex;
//         spU = (*pp)->Uconf_tex + " & " + (*pp)->Uterm_tex;
//      } else (CL == (*pp)->Lconf_tex && CU == (*pp)->Uconf_tex ) {
//         spL = "       & " + (*pp)->Lterm_tex;
////         spU = "       & " + (*pp)->Uterm_tex;
//      } 
//*/
//      out << spL << " & " << spU << " & " << endl;
//      int nZ = (*pp)->vec_Z.size();
 //     for (int ii = 0; ii < (*pp)->ZZ; ii++) {
// //        out.precision(3);
//         out.setf(ios::showpoint);
//         out.setf(ios::fixed, ios::floatfield);
//         float B = (*pp)->vec_Perc_err[ii]*100;
//         if (B < 0.01) out.precision(3);
//         if (B > 0.01) out.precision(2);
////         if (B > 1) out.precision(3);
//         if (B > 10) out.precision(4);
//
//         if (B > 10e-3) out << setw(6) << B;
//         else if(B < 10e-3 && B > 0) out << setw(6) << "0.000 ";
//       else out << setw(6)<< "" ;
//         if (B > 0.001)  out << setw(6) << B ;
//         else out << setw(6) << "";
//         if (ii == (*pp)->ZZ - 1 ) {}
//         else out << "&";
//      }
//      out << " \\\\" <<endl;
//      if (!(cnt%40)) out << TFooter_LS << endl << endl;
//      if (!(cnt%40)) out << THeader_LS;
//   }
//   out << TFooter_LS << endl;
//}

ostream& operator<<(std::ostream& os, LS::LS_err* A) {
      string ssi = A->Lconf + "_" + A->Lterm.substr(0,2);
      string ssf = A->Uconf + "_" + A->Uterm.substr(0,2);
      int W1 = 30 - ssi.size();
      int W2 = 30 - ssf.size();
      os << setw(W1) << ssi << " - ";
      os << setw(W2) << ssf << "  " ;
      int nZ = A->vec_Z.size();
      for (int ii = 0; ii < A->ZZ; ii++) {
         os.precision(3);
         os.setf(ios::showpoint);
         os.setf(ios::fixed, ios::floatfield);
         float B = A->vec_Perc_err[ii]*100;
         if (B > 0.001)  os << setw(10) << B << "";
         else if (B<10e-3 && B>0) os << setw(10) << "0.000";
         else os << setw(10) << "-" <<"";
      }
      os << endl;
      os << setw(30) << " ";
      os << setw(30) << " " ;
      for (int ii = 0; ii < A->ZZ; ii++) {
         os.precision(3);
         os.setf(ios::showpoint);
         os.setf(ios::scientific, ios::floatfield);
         float B = A->vec_GF[ii];
         if (B > 0)  os << setw(10) << B << "";
         else os << setw(10) << "-" <<"";
      }
      return os;
}


