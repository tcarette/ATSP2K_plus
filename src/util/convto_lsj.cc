#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <iomanip>
#include <sstream>
#include <map>
#include <cmath>

using namespace std;

/* 
EE = 2 * CONST * (E-GE)
E = (EE/2*CONST) + GE
*/

void fZ(double& ZZ, double& dRZ) {
   const double RZN = 109737.31534;
   const double RZD = 548.697e-06;
   double zmu = 0;
   if (ZZ == 1) zmu = 1;
   else if ( ZZ > 10) zmu = 2 * ZZ + 1 + (ZZ - 11)/2;
   else if (!(int(ZZ)%2) || (ZZ == 7)) zmu = 2 * ZZ;
   else zmu = 2 * ZZ + 1;
   dRZ = RZN/(1+RZD/zmu);
}

struct JTR_tr {
   double GE_read;
   bool unph;
   double lev_i, lev_f;
   double ei, ef, te;
   string ci, cf, ti, tf;
   int ji, jf;
   double angs_v, angs_a;
   string tr_type;
   double _sL, _sV, _gfL, _gfV, _akiL, _akiV;
   double _f_ik;
   JTR_tr(double& d1, double& d2, double& d3, string& s1, string& s2, int& i1,
          int& i2, double& d4, double& d5, string& s3, double& d6, double& d7,
          double& d8, double& d9, double& d10, double& d11,double& ge) :
      ei(d1), ef(d2), te(fabs(d3)), ci(s1), cf(s2), ji(i1), jf(i2), angs_v(fabs(d4)),
      angs_a(fabs(d5)), tr_type(s3), _sV(fabs(d6)), _sL(fabs(d7)), _gfV(fabs(d8)), _gfL(fabs(d9)),
      _akiV(fabs(d10)), _akiL(fabs(d11)), GE_read(ge) {
         int jj = (d1 <= d2) ? i1 : i2;
         _f_ik = _gfL/(1*jj + 1);
         lev_i = 0;
         lev_f = 0;
       }
};

typedef list<JTR_tr*>::const_iterator JCI;

ostream& operator<<(ostream& os, JCI jp) {
      os.setf(ios::right,ios::adjustfield);
      os << setw(4) << (*jp)->ji;
      os.setf(ios::fixed,ios::floatfield);
      os.precision(8);
      os << setw(14) << (*jp)->ei << "  ";
      os.setf(ios::left,ios::adjustfield);
      os << setw(59) <<  (*jp)->ci << endl;
      os.setf(ios::right,ios::adjustfield);
      os << setw(4) << (*jp)->jf;
      os.setf(ios::fixed,ios::floatfield);
      os.precision(8);
      os << setw(14) << (*jp)->ef << "  ";
      os.setf(ios::left,ios::adjustfield);
      os << setw(59) << (*jp)->cf << endl;
      os.precision(2);
      os.setf(ios::fixed,ios::floatfield);
      os.setf(ios::right,ios::adjustfield);
      os << setw(11) << (*jp)->te << " CM-1";
      os << setw(13) << (*jp)->angs_v << " ANGS(VAC)";
      os << setw(13) << (*jp)->angs_a << " ANGS(AIR)" << endl;
      os.precision(5);
      os.setf(ios::scientific,ios::floatfield);
      os << " " << (*jp)->tr_type;
         os << "  S =  " << (*jp)->_sL;
         os << "   GF =  " << (*jp)->_gfL;
         os << "   AKI =  " << (*jp)->_akiL << endl;
      if ((*jp)->_akiV) {
         os << "          " << (*jp)->_sV;
         os << "         " << (*jp)->_gfV;
         os << "          " << (*jp)->_akiV << endl;
      }
      os << "\n\n";
}

typedef map<int,map<int,map<char,map<string,double> > > >mapLevels;
typedef map<int,map<int,map<char,string> > > ML;
typedef map<string,map<int,double> > ME;


void print_lsj(list<JTR_tr*>& T) {
   typedef list<JTR_tr*>::iterator TI;
   ofstream OF("grasp.lsj");
   OF << endl << endl << endl;
   for (TI t = T.begin(); t != T.end(); t++) {
      OF << t << endl;
   }
}

template<char v>
struct WhichLev {enum {value = (int)v};};

void print_eachJ(ofstream& OF, double E,string& C) {
   OF << "   Ssms=   0.0000000000    g_J=   0.0000000000  g_JLS=   0.0000000000" << endl;
   OF.setf(ios_base::fixed,ios::floatfield);
   OF.precision(9);
   OF << setw(6) << 1 
      << setw(16) << E
      << "  " << C << endl;
   OF << " 1.00000000" << endl;
}

void print_J(ofstream& OF, int J, int n) {
   OF << setw(7) << "2*J =" << setw(5) << J << setw(10) << "NUMBER =" << setw(4) << n << endl;
}

void print_Z(ofstream& OF, int NEL, int Z) {

   OF.setf(ios_base::fixed,ios::floatfield);
   OF.precision(1);
   double ZD = double(Z);
   OF << "  Z=" << setw(2) << Z 
      << setw(7) << "Z =" << setw(6) << setprecision(1) << ZD 
      << setw(7) << "NEL =" << setw(4) << setprecision(1) << NEL 
      << setw(9) << "NCFG =" << setw(4) << 1 << endl << endl << endl;
}

struct Lev {
   string S;
   int J;
   double E;
   Lev(string& s, double e, int j) :S(s), E(e), J(j) {}
};

void print_j(list<JTR_tr*>& T,int NEL, int Z) {
   typedef list<JTR_tr*>::iterator TI;
   ofstream OF("grasp.j");

   map<double,map<int,string> > MJ;
   typedef map<double,map<int,string> >::iterator PJ;

   for (TI t = T.begin(); t != T.end(); t++) {
      double Ei = (*t)->ei;
      int Ji = (*t)->ji;
      string Si = (*t)->ci;
      MJ[Ei][Ji] = Si;
      double Ef = (*t)->ef;
      int Jf = (*t)->jf;
      string Sf = (*t)->cf;
      MJ[Ef][Jf] = Sf;
   }
   
   map<int,map<double,Lev*> >ML;
   for (PJ p = MJ.begin(); p != MJ.end(); p++) {
      double E = p->first;
      map<int,string>::iterator pi = p->second.begin();
      int J = pi->first;
      string S = pi->second;
      Lev *A = new Lev(S,E,J);
      ML[J][E] = A;
   }

   typedef map<int,map<double,Lev*> >::reverse_iterator MR;
   typedef map<double,Lev*>::iterator MI;
  
   print_Z(OF,NEL,Z);
   for (MR r = ML.rbegin(); r != ML.rend(); r++) {
      int J = r->first;
      int sz = r->second.size();
      print_J(OF,J,sz);
      for (MI pi = r->second.begin(); pi != r->second.end(); pi++) {
         double E = pi->first;
         string S = pi->second->S;
         print_eachJ(OF,E,S);
      }
      OF << endl << endl;
   }
   OF << " END" << endl;
}

void print_lev_help() {

cout << endl << endl;

cout << "                This program converts graspVU data to a .lsj and a .j files\n";
cout << endl;
cout << "  The program needs the following files:\n";
cout << "      case.e\n";
cout << "      case.ct\n";
cout << "  case.e is also a reference map between labels (.lsj notation) and "<< endl;
cout << "                  \"Pos J Parity\" (graspVU notation), the format is shown below\n";

cout << " Energy levels for ..." << endl;
cout << " Rydberg constant is   109737.31534" << endl;
cout << " No - Serial number of the state; Pos - Position of the state within the " << endl;
cout << " J/P block; Splitting is the energy difference with the lower neighbor" << endl;
cout << " -------------------------------------------------------------------------" << endl;
cout << " No Pos  J Parity Energy Total    Levels     Splitting  Configuration" << endl;
cout << "                      (a.u.)      (cm^-1)     (cm^-1)" << endl;
cout << " -------------------------------------------------------------------------" << endl;
cout << "  1  1  1/2 +    -241.4318363                           3s_2S" << endl;
cout << "  2  1  1/2 -    -241.1868119    53776.62    53776.62   3p_2P" << endl;
cout << "  3  1  3/2 -    -241.1857460    54010.58      233.96   3p_2P" << endl;
cout << "  4  1  5/2 +    -240.8999698   116731.20    62720.62   3d_2D" << endl;
cout << "  5  1  3/2 +    -240.8999575   116733.88        2.69   3d_2D" << endl;
cout << "  6  2  1/2 +    -240.8542360   126768.60    10034.71   4s_2S" << endl;
cout << " -------------------------------------------------------------------------" << endl;
cout << endl;
cout << " To create file .e:\n";
cout << "  1. Run \"rlevels *\" and create a file \"case.e\"" << endl;
cout << "  2. Append at the end of each line the labels\n";
cout << "  3. The combination Pos & J & Parity must be unique to allow correct\n"
        " mapping, or the program will exit(1)\n" << endl;
cout << " To create a file .ct, concatenate all .ct files\n";
cout << "\n Please make sure the labeling is consistent with what is on the Web\n";

cout << "\n The output will be saved in a file grasp.lsj\n";
}

bool check_files(string& E, string& CT) {

   ifstream i1(E.c_str());
   if (!i1) {
      cerr << " Error: " << E << " missing, or can not be open, exit(1)" << endl;
      exit(1);
   }

   ifstream i2(CT.c_str());
   if (!i2) {
      cerr << " Error: " << CT << " missing, or can not be open, exit(1)" << endl;
      exit(1);
   }
}

void get_NELZ(int& NEL, int& Z) {

   cerr << "\n Enter atomic Number Z : ";
   string _nel, _z;
   getline(cin,_z,'\n');
   istringstream istz(_z);
   istz >> Z;
   if (Z < 1 || Z > 110) 
      cerr << "Wrong input for Z: \"" << Z << "\", exit(1)." << endl;
      
   cerr << " Enter Number of Electrons, NEL : ";
   getline(cin,_nel,'\n');
   istringstream istn(_nel);
   istn >> NEL;
   if (NEL < 1 || NEL > 110) 
      cerr << "Wrong input for Z: \"" << Z << "\", exit(1)." << endl << endl;
}

void get_input(string& E, string& CT) {

   cerr << "Enter a case (.ct file and .e files must be present): ";
   string casename;
   getline(cin,casename,'\n');

   CT = casename + ".ct";
   E = casename + ".e";

   check_files(E,CT); 
}

bool read_e(mapLevels& M, string& E) {

   bool not_unique_levels(true);
   ifstream inf(E.c_str());
   if (!inf) {
      cerr << " Error, Can't open or read " << E << ", exit(1)" << endl;
      exit(1);
   }

   int sz = 256;
   char t_read[sz];
   bool line1(false);

   while (inf.getline(t_read,sz)) {
      string rd(t_read);
     
      if (rd.size() > 28) {
         if (!line1 && 
              rd.substr(0,29) == " No Pos  J Parity Energy Tota") {
            line1 = true;
         } 
         if (line1 && rd.substr(0,28) == "                      (a.u.)") break; 
      }
   } 
 
   int _i(0),i1(0),j2(0);
   char p1;
   double d1(0);
   string L1,jstr("     ");
 
   string f_(80,'-');
   cout << " Found the following levels: " << endl;
   cout<< f_ << endl;
   cout << setw(8) << "Pos" << setw(8) << "2*J" << setw(8) << "Parity" 
        << setw(20) << "Total E" << "   Level " << endl;
   cout << f_ << endl;
   while (inf.getline(t_read,sz)) {
      string rd(t_read);  
      if (rd.size() > 28 && rd.substr(0,10) != " ---------") {
         string s1 = rd.substr(0,29); 
         string s2;
         if (rd.size() > 54) {
            s2 = rd.substr(55);
            istringstream is2(s2);
            is2 >> L1;
         } else {
            cerr << " Error, line " << rd << " has missing label " << endl;
            exit(1);
         }

         istringstream ist(s1);
         ist >> _i >> i1 >> jstr >> p1 >> d1; 
        
         istringstream is1(jstr);
         is1 >> j2;

         if (jstr.find("/") == string::npos) j2 *= 2; // multiply by 2 if not "/"

         if (!M[i1][j2][p1][L1]) M[i1][j2][p1][L1] = d1;
         else {
           cerr << "\n\nError, the combination of Pos & J & Parity"
                   " is not unique, the program will exit!" << endl;
           cerr << "Error in Line::" << rd << endl;
           not_unique_levels = false; // set to false so that the program will exit
         }
   cout.precision(7);
   cout.setf(ios::fixed,ios::floatfield);
   cout << setw(8) << i1 << setw(8) << j2 << setw(8) << p1
        << setw(20) << d1 << "   " << L1 << endl;
//cout << "::" << i1 << "::" << j2 << "::" << p1 
//     << "::" << d1 << "::" << L1 << "::" << endl;
      } 
   } 
   cout << f_ << endl;
   return not_unique_levels;
}

void repl_charDE(string& S) {
   string::iterator p;
   for (p = S.begin(); p != S.end(); p++) {
     char A = *p;
     if (A == 'E' || A == 'D') *p = 'e';
   }
   return;
}

bool read_ct(mapLevels& M, string& CT, list<JTR_tr*>& TL) {
   ifstream inf(CT.c_str());
   if (!inf) {
     cerr << " Exiting! No such file: " << CT << endl;
     exit(1);
   }
  
   int sz = 256;
   char t_read[sz];
   int i = 0;

   vector<string>VS;

   while (inf.getline(t_read,sz)) { 
      string R(t_read);
      VS.push_back(R);
   }

   typedef vector<string>::iterator VI;
   
   string tr_type;
   string en_tr;

   int lcount = 0;
   int count_OK(0);
   int count_All(0);
   for (VI p = VS.begin(); p != VS.end(); p++,lcount++) {
// cout << *p << ":: before" << endl; 
//     repl_charDE(*p); 
//cout << *p << ":: after" << endl << endl;;
     if ((*p).size() == 34 && (*p).substr(16,18) == ")-pole transitions") {
       string o((*p).substr(14,2));
       istringstream ist(o);
       string A(" "); 
       A[0] = (*p)[1];
       string so;
       ist >> so;
       tr_type = A + so;
     }

     string B(" Lev  J P   Lev  J P"); //       E (Kays)         A (s-1)");

     if ((*p).substr(0,20) == B) {
        while (lcount++ < VS.size()-1 && (*p++).size() > 0) {
           repl_charDE(*p);
           int l1(0),j1(0),l2(0),j2(0);
           char c1,c2;
           double d1(0),d2(0);
           char T;
           double d3(0),d4(0);
           char C1;
           string jstr1, jstr2;
           double D1(0),D2(0),D3(0);

           if ((*p).size() >= 68) {
//               (*p)[33] = 'e';
//               (*p)[50] = 'e';
//               (*p)[65] = 'e';
//               if ((*p).size() >= 80) (*p)[80] = 'e';
               istringstream ist1((*p).substr(0,54));
               istringstream ist2((*p).substr(55));
               ist1 >> l1 >> jstr1 >> c1 >> l2 >> jstr2 >> c2 >> d1 >> T >> d2;
//cout << (*p) << endl;
//cout << l1 << jstr1 << c1 << l2 << jstr2 << c2 << d1 << T << d2 << endl;
               if (d2) {
                  if ((*p).size() >= 84) ist2 >> d3 >> d4;
                  else ist2 >> d3;
               }
               if (T=='M') {
                  D1 = d2; d2 = 0;
                  D2 = d3; d3 = 0;
                  D3 = d4; d4 = 0;
               } 
               if (T == 'C') {
                  lcount++;
                  if (lcount < VS.size()) {
                     p++;
                     repl_charDE(*p);
                     if ((*p).size() >= 68) {
//                        (*p)[50] = 'e';
//                        (*p)[65] = 'e';
//                        if ((*p).size() >= 80) (*p)[80] = 'e';
                        istringstream q1((*p).substr(0,54));
                        q1 >> C1 >> D1;
                        if (D1) {
                           istringstream q2((*p).substr(55));
                           if ((*p).size() >= 84) q2 >> D2 >> D3;
                           else q2 >> D2;
                        }
                     }
                  }
               }
if (d1) {
             count_All++;
             istringstream js1(jstr1);
             istringstream js2(jstr2);
             js1 >> j1;
             js2 >> j2;
             int j21(0),j22(0);

             if (jstr1.find("/") == string::npos) j21 = 2*j1;
             else j21 = j1;
             if (jstr2.find("/") == string::npos) j22 = 2*j2;
             else j22 = j2;

             string S1, S2;
             bool OK_transiton(true);
             typedef map<string,double>::iterator ISD;
             ISD p1 = M[l1][j21][c1].begin();
             if (p1 != M[l1][j21][c1].end()) S1 = p1->first;
             else OK_transiton = false;
             ISD p2 = M[l2][j22][c2].begin();
             if (p2 != M[l2][j22][c2].end()) S2 = p2->first;
             else OK_transiton = false; 
             //string S1 = M[l1][j1][c1];
             //string S2 = M[l2][j2][c2];

             double E1(0), E2(0);

             if (M[l1][j21][c1][S1]) E1 = M[l1][j21][c1][S1];
             else OK_transiton = false;

             if (M[l1][j21][c1][S1]) E2 = M[l2][j22][c2][S2];
             else OK_transiton = false;

             double a1(0),a2(0),a3(0);
             double da(fabs(d1)); 
             a1 = 100000000/(da);
if (d1 < 0) {
   std::swap(E1,E2);
   std::swap(j21,j22);
   std::swap(S1,S2);
} else {
//cout << j21 << "::" << j22 << "::" << d2 << "::" << D1 << endl;
    d2 *= (double(j21)+1)/(double(j22) + 1);
    D1 *= (double(j21)+1)/(double(j22) + 1);
//if (d2 == 0 || D1 == 0) cout << "::::::::" << j21 << "::" << j22 << "::" << d2 << "::" << D1 << endl;
}

JTR_tr* A = new JTR_tr(E1,E2,da,S1,S2,j21,j22,a1,a2,tr_type,d4,D3,d3,D2,d2,D1,a3);
   if (OK_transiton)  {
      count_OK++;
//      cout << "New transition " << endl;
      TL.push_back(A);
//if (d2 == 0 || D1 == 0) {
//      cout << S1 << "::"<<S2 << "::"<<E1 << "::"<<E2; 
//      cout << "::"<<l1 << "::"<<j1 << "::"<<c1 << "::"<<l2 << "::"<<j2 << "::"<<c2 << "::"<<d1 << "::"<<T << "::"<<d2;
//      cout << "::"<<d3 << "::"<<d4;
//      cout << "::"<<C1 << "::"<<D1 << "::"<<D2 << "::"<<D3 << endl;
//}
   } else {
      cerr << "  Error for line : " << *p << endl;
   }
}
            }
         }
      }
   }

   if (count_All != count_OK) 
      cout << " There were errors: All transitions = " 
           << count_All << ", only " << count_OK << " transitions in .lsj file\n";
   return true;
}

int main() {

   print_lev_help();

   list<JTR_tr*> T;

   ML m;
   ML m1; 
   ME me;
   string E,CT;
   int NEL, Z;
   get_input(E,CT);
   get_NELZ(NEL,Z);

   mapLevels M;
   if (!read_e(M,E)) {
      cerr << "\n Wrong input! Some combinations "
              "Pos J Parity were duplicated, exit(1)\n" << endl;;
      cerr << " One way to avoid this error is to prepare several"
              " groups of .e files and .ct files with unique Pos & J & Parity\n";
      exit(1);
   }

   read_ct(M,CT,T);

   map<int,string> MS;
//   if (!read_ct_old(T,m,m1,me,ct_file)) {
//      cerr << "\n Enter output from .cm file (energy data): ";
//      string cm_file;
//      getline(cin,cm_file,'\n');
//      read_ct(T,m,m1,ct_file);
//      read_cm(T,MS,cm_file); 
//   }
   cerr << " ...writing " << T.size() << " transitions in grasp.lsj\n";
   print_lsj(T);
   cerr << " ... writing levels in grasp.j for Z = " << Z << " and NEL = " << NEL << endl;
   print_j(T,NEL,Z);
   cerr << "\n\n The converted files are grasp.lsj and grasp.j" << endl;  
}
