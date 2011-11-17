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

using namespace std;

//----------------------------------------------------
//-----  non-members included in namespace transition   
//-----  also some of them called from members         
//----------------------------------------------------

void Dt::fZ(double& ZZ, double& dRZ) {
   double zmu = 0;
   if (ZZ == 1) zmu = 1;
   else if ( ZZ > 10) {zmu = 2 * ZZ + 1 + (ZZ - 11)/2;}
   else if (!(int(ZZ)%2) || (ZZ == 7)) {zmu = 2 * ZZ;}
   else {zmu = 2 * ZZ + 1;}
   dRZ = RZN/(1+RZD/zmu);
}


// compare if term tb duplicates ta
bool Dt::Un_t (LCI ta, LCI tb) {
   return (((*ta)->config == (*tb)->config) &&
           ((*ta)->term == (*tb)->term) &&
           ((*ta)->J == (*tb)->J));
}

void Dt::ETable_jini::set_ZNEL() {
   Print_table::_Z_ = int(Z);
   Print_table::_NEL_ = NEL;
   Print_table::_RZ_ = RZ;
}

void Dt::set_tr_type(string& s) {
   *Print_table::_tr_type_ = s;
}

void  Dt::conv_J(string& B, const int j) {
   string S;
   std::ostringstream ost(S);
   if (j%2) {
     if (j > 10) ost <<  j << "/2";
     else ost << " " << j << "/2";
   }
   else {
      if (j >= 20) ost << "  " << j/2;
      else ost << "   " << j/2;
   }
   B = ost.str();
}

//---------------------------------------------------------------
//------- comparing two transitions using initial-final energies  
//---------------------------------------------------------------
class Dt::Comp_lsjLU {
public:
   bool operator()(const JTR_tr* ta, const JTR_tr* tb) const {
      //return true if ta < tb
      if(ta->ef == tb->ef) return ta->ei <= tb->ei;
      return ta->ef <= tb->ef;
   }
};

class Dt::Comp_lsjUL {
public:
   bool operator()(const JTR_tr* ta, const JTR_tr* tb) const {
      //return true if ta < tb
      if(ta->ei == tb->ei) return ta->ef <= tb->ef;
      return ta->ei < tb->ei;
   }
};

class Dt::Comp_mp {
public:
   bool operator()(const JTR_mp* ta, const JTR_mp* tb) const {
      //if(ta->E_low != tb->E_low) return ta->E_low >= tb->E_low;
      //if(ta->E_high != tb->E_high) return ta->E_low >= tb->E_high;
      return ta->E_tr >= tb->E_tr;
   } 
};   

class Dt::Comp_lsjE {
public:
   bool operator()(const JTR_tr* ta, const JTR_tr* tb) const {
//      if(ta->te == tb->te) return ta->te <= tb->te;
      if(ta->tr_type != tb->tr_type) 
          return ta->tr_type <= tb->tr_type;
      return ta->te <= tb->te;
   }
};

class Dt::Comp_lsjgF {
public:
   bool operator()(const JTR_tr* ta, const JTR_tr* tb) const {
      if(ta->_f_ik == tb->_f_ik) return ta->te <= tb->te;
      return ta->_f_ik > tb->_f_ik;
   }
};

class Dt::Comp_lsjSL {
public:
   bool operator()(const JTR_tr* ta, const JTR_tr* tb) const {
      if(ta->_sL == tb->_sL) return ta->te <= tb->te;
      return ta->_sL > tb->_sL;
   }
};

class Dt::Comp_lsjWL {
public:
   bool operator()(const JTR_tr* ta, const JTR_tr* tb) const {
      if(ta->angs_v == tb->angs_v) return ta->angs_v <= tb->angs_v;
      return ta->angs_v <= tb->angs_v;
   }
};


//--------------------------------------------------------------
//-----     operator for sorting lsj tbl within the final term
//--------------------------------------------------------------
class Dt::Comp_lsj_int { 

public:
   bool operator()(const JTR_tr* ta, const JTR_tr* tb) const {
      //return true if ta < tb
      if(ta->cf == tb->cf) return ta->ef <= tb->ef;
      return ta->ei < tb->ei;
   }
};

//----------------------------------------------------
// ----member functions of Dt::Level--
//---------------------------------------------------

Dt::Level::Level(std::string sc = "",  int jj = 0,
      double d1 = 0, double d2 = 0, double d3 = 0,
      double d4 = 0, double d5 = 0, double d6 = 0) :
      tmp(sc), J(jj), totalE(d1), lev(d2), lifetime(d3),
      Ssms(d4), g_J(d5), g_JLS(d6), config_rep(" "),
      term_rep(" "), unph(true) , Splitting(0) { mk_term(); }
void mk_term();

Dt::JTR_tr::JTR_tr(double& d1, double& d2, double& d3, string& s1, 
      string& s2, int& i1, int& i2, double& d4, double& d5, string& s3, 
      double& d6, double& d7, double& d8, double& d9, double& d10, 
      double& d11, double& ge) :
      ei(d1), ef(d2), te(d3), ci(s1), cf(s2), ji(i1), jf(i2), angs_v(d4),
      angs_a(d5), tr_type(s3), _sL(d6), _sV(d7), _gfL(d8), _gfV(d9),
      _akiL(d10), _akiV(d11), GE_read(ge) {
         int jj = (d1 <= d2) ? i1 : i2;
         _f_ik = _gfL/(1*jj + 1);
         lev_i = 0;
         lev_f = 0;
      }

void Dt::Level::mk_term() {
   int ss = tmp.size();
   if (ss < 2) return;
   if (isdigit(tmp[ss-1]) || tmp[ss-1] == '*') {
      ss--;
      config = tmp.substr(0, ss - 3);
      term = tmp.substr(ss - 2, 3);
   } else {
      config = tmp.substr(0, ss - 3);
      term = tmp.substr(ss - 2, 2);
   }
}

// ---->>>member functions of Dt::JTR_tr<<<-----
// void get_terms();      // not implemented

//----------------------------------------------------------
// --- member functions of Abstract Dt::Base_table  
//--------------------------------------------------------
  

// sort the energies from  low  :->Dt::Base_table::sort_e() 
void Dt::Base_table::sort_e() {
   std::cout << "...sorting energies: ...\n\n";
//   sort(lv.begin(),lv.end(),Comp());
   lv.sort(Comp());
   LCI p = lv.begin();
   if (p == lv.end()) return;
   (*p)->config_rep = (*p)->config;
   (*p)->term_rep = (*p)->term;
   string C(" ");
   string T(" ");
   string curr_C(" ");
   string curr_T(" ");
   for( ; p != lv.end(); p++) {
      bool NC = false;
      if (curr_T != (*p)->term) {
         curr_C = (*p)->config;
         curr_T = (*p)->term;
         (*p)->config_rep = curr_C;
         (*p)->term_rep = curr_T;
         NC = true;
      }
      else (*p)->config_rep = " ";

      if (!NC) {
         if (curr_T != (*p)->term) {
            curr_T = (*p)->term;
            (*p)->term_rep = curr_T;
         }
         else (*p)->term_rep = " ";
      }

//      if (curr_SLC != (*pi)->Lconf_tex) {
//         curr_SLC = (*pi)->Lconf_tex;
//         sLC = curr_SLC;
//      } else { sLC = ""; }
//      if ( ++ptmp != lv.end()) {
//        if ((*ptmp)->term != (*p)->term) {
//           (*ptmp)->term_rep = (*ptmp)->term;
//        }
//        if ((*ptmp)->config != (*p)->config) {
//           (*ptmp)->config_rep = (*ptmp)->config;
//           (*ptmp)->term_rep = (*ptmp)->term;
//        }
//           (*ptmp)->term_rep = (*ptmp)->term;
//           (*ptmp)->config_rep = (*ptmp)->config;
//      }
   }
}

// sort levels according to terms and check for duplicates
// :->Dt::Base_table::sort_t() 
void Dt::Base_table::sort_t() {
   std::cout << "...arranging terms: ...\n\n";
   LI p1 = lv.begin();
   LI p2 = lv.begin();
   LI ptmp = lv.begin();
   while (p1 != lv.end()) {
      p1 = ptmp;
      p2 = ++ptmp;
      if (p2 == lv.end()) break;
      while (p2 != lv.end()) {
         if (Un_t(p1, p2)) {
//            (*p2)->config_rep = "dupl. . . . . . . . . . . .->";
//            (*p2)->term_rep = (*p2)->term;
            (*p2)->unph = false;
            if ((*ptmp)->unph) {
               lv.insert(ptmp,(*p2));
               p2 = lv.erase(p2);
            } else if(++p2 != lv.end()) {}
         } else if(++p2 != lv.end()) {}
      }
   }
//@@@@@@@@@@@@@2
   string C(" ");
   string T(" ");
   string curr_C(" ");
   string curr_T(" ");

  for( LI p = lv.begin(); p != lv.end(); p++) {
      bool NC = false;
      if (curr_C != (*p)->config) {
         curr_C = (*p)->config;
         curr_T = (*p)->term;
         (*p)->config_rep = curr_C;
         (*p)->term_rep = curr_T;
         NC = true;
      }
      else (*p)->config_rep = " ";

      if (!NC) {
         if (curr_T != (*p)->term) {
            curr_T = (*p)->term;
            (*p)->term_rep = curr_T;
         }
         else (*p)->term_rep = " ";
      }
   }
//@@@@@@@@@@@@@@@
}

void Dt::Base_table::rep_TC() {
   for (LI pl = lv.begin(); pl != lv.end(); pl++) {
         if (!((*pl)->unph)) {
            (*pl)->config_rep = "dupl. . . . . . . . . . . .->";
            (*pl)->term_rep = (*pl)->term;
         }   
   }
}

// :->remove unphysical levels using CpRmTbLv
void Dt::Base_table::LV_rm_unph() {
   lv.remove_if(CpRmTbLv());
//   splitt();
}

// After sorting compute splitting, and sort
void Dt::Base_table::splitt(){

   LI spb = lv.begin();
   LI spe = lv.begin();
   LI sptmp = lv.begin();
   while(spb != lv.end()) {
      spb = sptmp;
      do {
         if (++sptmp == lv.end()) break; 
         if (((*sptmp)->config == (*spb)->config) && 
              ((*sptmp)->term == (*spb)->term)) (*sptmp)->term_rep = " ";
      }
      while((((*sptmp)->config == (*spb)->config) &&
         ((*sptmp)->term == (*spb)->term)) && (sptmp != lv.end())); 
      spe = sptmp;
      if(++spe == lv.end() || sptmp == lv.end()) break;
      else {
         while (spe != lv.end()) {
            if (((*spe)->config == (*spb)->config) && 
                        ((*spe)->term == (*spb)->term)) {
               (*spe)->term_rep = " ";
               lv.insert(sptmp,(*spe));
               spe = lv.erase(spe);
            } else if(++spe != lv.end()) {}
         }
      }
   }
//@@@@@@@@@@@@@2
   string C(" ");
   string T(" ");
   string curr_C(" ");
   string curr_T(" ");

  for( LI p = lv.begin(); p != lv.end(); p++) {
      bool NC = false;
      if (curr_C != (*p)->config) {
         curr_C = (*p)->config;
         curr_T = (*p)->term;
         (*p)->config_rep = curr_C;
         (*p)->term_rep = curr_T;
         NC = true;
      }
      else (*p)->config_rep = " ";

      if (!NC) {
         if (curr_T != (*p)->term) {
            curr_T = (*p)->term;
            (*p)->term_rep = curr_T;
         }
         else (*p)->term_rep = " ";
      }
   }
//@@@@@@@@@@@@@@@

// end sort, compute splitting 
   LI pa = lv.begin();
   LI pb = lv.begin();
   while (pa != lv.end()) {
      while(++pb != lv.end()) { 
         if ((*pb)->config==(*pa)->config && (*pa)->term==(*pb)->term)
            {(*pb)->Splitting = (*pb)->lev - (*pa)->lev;}
         else break;
      }
      pa = pb;
   }
 
}

// This function does: 1. removes unphysical, 2. sorts *.lsj files, 
// 3.sorts again for terms, 4.computes lifetimes 5. prints the trasnitions.
// remove unphysical
void Dt::Base_table::JTR_rmv() {
   int A = lsj.size();
   cout << lsj.size() << " lsj transitions initially\n";
   for(JCI jp = lsj.begin(); jp != lsj.end(); jp++) {
      (*jp)->unph = JTR_unph(jp);
   }
// remove all lsj not present in table levels
   lsj.remove_if(CmpRmJTR());
   cout << "...sorting .lsj transitions\n" << endl;
   JTR_sortLSJ();
   JCI jr;
    cout << "...checking .lsj for duplicates" << endl;
//   jr = unique(lsj.begin(),lsj.end(),RmvDupl()); //remove any duplicate .lsj 
   lsj.erase(unique(lsj.begin(),lsj.end(),RmvDupl()),lsj.end());
   cout << lsj.size() << " transitions after removing of " 
        << A - lsj.size() << " duplicated transitions\n";
   if (lsj.size() == 0 ) 
        cerr << " There are 0 transitions, exit(1)\n", exit(1);
// sort the list of lsj 
   cout << "\nEnter a number to sort the transitions by:\n";
   cout << " 1. Transition energies.\n";
   cout << " 2. Lower-Upper.\n";
   cout << " 3. Upper-Lower.\n";
   cout << " 4. f value.\n";
   cout << " 5. Line strength.\n";
   cout << " 6. Wavelength(air).\n";
   int s;
   cin >> s;

   switch (s) {
   case 1:
      JTR_sortE();
      JTR_sortT();
      break;
   case 2:
      JTR_sortELU(); 
      JTR_sortT();
      break;
   case 3:
      JTR_sortEUL();
      JTR_sortT();
      break;
   case 4:
      JTR_sortgF();
      JTR_sortT();
      break;
   case 5:
      JTR_sortSL();
      JTR_sortT();
      break;
   case 6:
      JTR_sortWL(); 
      JTR_sortT();
      break;
   default:
      JTR_sortE(); 
      JTR_sortT();
   }
   JTR_sortELU();
   JTR_multiplet();
}

void Dt::Base_table::JTR_multiplet() {
   cout << " ...arranging by multiplets\n";

   JCI jt = lsj.begin();

   double E_ground (0);
   if (Base_table::E_min) 
      E_ground = Base_table::E_min;
   else 
      E_ground = (*jt)->ef;

   for(JCI jp = lsj.begin(); jp != lsj.end(); jp++) {
      bool ins = false;
      for(JMP jm = jmp.begin(); jm != jmp.end(); jm++) {
         if((*jp)->cf == (*jm)->C_low + '_' + (*jm)->T_low  &&
            (*jp)->ci == (*jm)->C_high + '_' + (*jm)->T_high ) ins = true;
         if (ins) (*jm)->MP_add(*jp);
         if (ins) break;
      }
      if (!ins) {  // insert a new multiplet
         JTR_mp* A = new JTR_mp(jp);
         //cout << (*jp)->cf << (*jp)->tf << (*jp)->ci << (*jp)->ti << endl;
         A->MP_add(*jp);
         jmp.push_back(A);
      }
   }
  
   Level_MP lmp(E_ground,this->lv);
   for(JMP jm = jmp.begin(); jm != jmp.end(); jm++) {
      (*jm)->lsjm.sort(Comp_lsjE());
      (*jm)->MP_computeE(E_ground,lmp);
      (*jm)->MP_compute_Aki();
//      cout << (*jm)->C_low<<"_"<<(*jm)->T_low
//           <<"-"<<(*jm)->C_high<<"_"<<(*jm)->T_high<<endl;
//      cout << (*jm)->lsjm.size() << "=(*jm)->lsjm.size()\n";
   }
}

class Dt::MP_Lev_comp {
public:
   bool operator()(const Dt::JTR_mp::lev_map a,
                   const Dt::JTR_mp::lev_map b) const {
     if (a.E != b.E) return a.E > b.E;
     else return a.j2 < b.j2;
   }
};

void Dt::JTR_mp::MP_add(JTR_tr* j) {
  JTR_mkterms(T_low,C_low,j->cf);
  JTR_mkterms(T_high,C_high,j->ci);
  lsjm.push_back(j);
}

Dt::Level_MP::Level_MP(double E_ground,std::list<Level*>& lv) {

   map<const char,int> Map_LS;
   Map_LS['S'] = 0;
   Map_LS['P'] = 1;
   Map_LS['D'] = 2;
   Map_LS['F'] = 3;
   Map_LS['G'] = 4;
   Map_LS['H'] = 5;
   Map_LS['I'] = 6;
   typedef map<string,int>::iterator SSI;
   map<string,int>m_g;
   for (LCI p = lv.begin(); p != lv.end(); p++) {
      int _IT_U = int((*p)->term[0] - 48);
      int _L_U = int(Map_LS[(*p)->term[1]]);
      int gg = _IT_U*(2*_L_U+1);
      string C = (*p)->config + "_" + (*p)->term;
      lev_MPG[C] = gg;
      int jj = (*p)->J;
      m_g[C] += jj+1; 
      double Z(Print_table::_Z_);
      double RZ(0);
      fZ(Z,RZ);
      double E = ((*p)->totalE-E_ground)*RZ*2*(jj+1)/gg;
      lev_MPE[C] += E;
   } 

   for (SSI p = lev_MPG.begin(); p != lev_MPG.end(); p++) {
      string C = (*p).first; 
      if (m_g[C] != lev_MPG[C]) lev_MPG[C] = 0;
   }
}

void Dt::JTR_mp::MP_computeE(double E_ground, Level_MP& lmp) {
   string SH = C_high+"_"+T_high;
   string SL = C_low+"_"+T_low;
   MP_gk = lmp.lev_MPG[SH];
   MP_gi = lmp.lev_MPG[SL]; 
   E_low = lmp.lev_MPE[SL];
   E_high = lmp.lev_MPE[SH];
   E_tr =  E_high - E_low;
   MP_WL = 1e8/E_tr;
}

void Dt::JTR_mp::MP_compute_Aki() {
   bool e1 = false;
   for (JCI p = lsjm.begin(); p != lsjm.end(); p++) {
      double a1, a2, a3;
      a1 = (*p)->_akiL;
      a2 = 1e8/(*p)->te;
      a3 = (*p)->_f_ik;
      double a4((*p)->ji+1);
      double a5((*p)->jf+1);
      double a6 = pow(a2,3);
      if ((*p)->tr_type == "E1") {
         MP_AKI += a1*a4*a6;
         MP_f_ik += a2*a3*a5;
         MP_S += (*p)->_sL;
         e1 = true;
      }
   }
   if (e1) MP_AKI /= (MP_gk*pow(MP_WL,3));
   if (e1) MP_f_ik /= (MP_WL*MP_gi);
//   cout << "::" << MP_AKI << "::" << MP_f_ik << "::" << MP_S <<"::"<< endl;
}

void Dt::Base_table::JTR_sortMP() {
   jmp.sort(Comp_MP());
}

//  a function comparing the levels, used in sorting *.lsj
bool Dt::Base_table::JTR_unph(JCI& jp) {
   for (LCI pi = lv.begin(); pi != lv.end(); pi++) {
      if(fabs((*pi)->totalE - (*jp)->ei) < 1e-8) {
         for(LCI pf = lv.begin(); pf != lv.end(); pf++)
            if(fabs((*pf)->totalE - (*jp)->ef) < 1e-8) return  true;
      }
   }
   for (LCI pi = lv.begin(); pi != lv.end(); pi++) {
      if(fabs((*pi)->totalE - (*jp)->ei) < 1e-7) {
         for(LCI pf = lv.begin(); pf != lv.end(); pf++)
            if(fabs((*pf)->totalE - (*jp)->ef) < 1e-7) {
                return  true;
            } 
      }
   }
   return false;
}

// all energies from *.lsj are sorted form low->high

void Dt::Base_table::JTR_sortE() {
   lsj.sort(Comp_lsjE());
}

void Dt::Base_table::JTR_sortELU() {
   lsj.sort(Comp_lsjLU());
}

void Dt::Base_table::JTR_sortEUL() {
   lsj.sort(Comp_lsjUL());
}

void Dt::Base_table::JTR_sortWL() {
   lsj.sort(Comp_lsjWL());
}

void Dt::Base_table::JTR_sortgF() {
   lsj.sort(Comp_lsjgF());
}

void Dt::Base_table::JTR_sortSL() {
   lsj.sort(Comp_lsjSL());
}

void Dt::Base_table::JTR_sortLSJ() {
   lsj.sort(Comp_LSJ()); 
}

// more specialized sorting on sections of the *.lsj (alike)
void Dt::Base_table::JTR_sortT() {
   JI jp = lsj.begin();
   JI ji = lsj.begin();
   while(jp != lsj.end()) {
      for(ji = jp; ji != lsj.end(); ++ji) { 
         if ((*ji)->ci != (*jp)->ci) { 
//             cout << "< >" << jp << endl;
//             cout << "<" <<++i << ">" << ji << endl;
             break;
         }
      }       
//  jp:  begining of the section
//->ji: end
//-------------------------------------
      JI pb = jp;
      JI pe = jp;
      JI ptmp = jp;
      while(pb != ji) {
         pb = ptmp;
         do {if (++ptmp == ji) break;}
         while (((*ptmp)->cf == (*pb)->cf) && ptmp != ji);
         pe = ptmp;
         if(++pe == ji || ptmp == ji) break;
         else {
            while (pe != ji) {
               if ((*pe)->cf == (*pb)->cf) {
                  if ((*pe)->te > (*pb)->te) lsj.insert(pb,(*pe));
                  else lsj.insert(ptmp,(*pe));
                  lsj.erase(pe); 
                  pe = ptmp;
               } else if(++pe != ji) {}
            } 
         } 
      }
//--------------------------------------
// = -> advance to the next section
      jp = ji;
   } 
}

// Computing the lifetimes

void Dt::Base_table::JTR_lifetimes() {

//   map<double,double> lt;
   typedef map<double,double>::const_iterator MCI;

   cout << " ...Computing Liefetimes...\n" << endl;

   typedef map<string,map<int,map<double,double> > >MMM;
   MMM mmm;
   for (JCI jb = lsj.begin(); jb != lsj.end(); jb++) {
      mmm[(*jb)->ci][(*jb)->ji][(*jb)->ei] += (*jb)->_akiL;
//cout << (*jb)->ci << (*jb)->ji << (*jb)->ei << " == " << (*jb)->_akiL << endl;
   }

   for (LI lp = lv.begin(); lp != lv.end(); lp++) {
     string _c = (*lp)->config + "_" + (*lp)->term;
     int _j = (*lp)->J;
     double _e = (*lp)->totalE;
     double _e_lev = mmm[_c][_j].begin()->first;
     double tr = mmm[_c][_j].begin()->second;
//cout << _c << _j << _e <<"::"<<_e_lev<< " == " << tr << endl;
     if (tr && fabs(_e - _e_lev) < 10e-8) {
        double _ti = 1/tr;
        (*lp)->lifetime = _ti;
     }
   }
}

 // separate terms, caled from JTR_print
void Dt::JTR_mkterms(string& term_i, string& config_i, string& sc_i) {
   int szi = sc_i.size();
   if (isdigit(sc_i[szi-1]) || sc_i[szi-1] == '*') {
      szi--;
      config_i = sc_i.substr(0,szi - 3);
      term_i = sc_i.substr(szi - 2, 4);
   } else {
      config_i = sc_i.substr(0,szi - 3);
      term_i = sc_i.substr(szi - 2, 3);
   }
} 

//##############################################
//##############################################
// printout .lsj.bin file 

// not implemented for now
void Dt::Base_table::CpLT() {}

//----------------------------------------------------
// ---  member functions of Dt::ETable_lsjini--
//----------------------------------------------------

Dt::ETable_lsjini::ETable_lsjini(const std::string& ss) :
file_name(ss), E_min(0) {

   Base_table::E_min = 0;
// read *.lsj file
   read_lsj(ss);
// compute levels   
   get_lev();
   std::cout << "\n... " << file_name
       << " reading OK...\n";
// sort for energies
   sort_e();
// sort terms
   sort_t();
// set conf. and term representations;
   rep_TC();
// compute RZ constant
   fZ(Z,RZ);
// set global _Z_ _NEL_
//   set_ZNEL();
   comp_lev();
   splitt();
   print_lev(std::cout);
   splitt();
//   print_lev(std::cout,false);
}

class Dt::Lev_map_comp {
public:
   bool operator()(const tmp_lv* a, const tmp_lv* b) const {
     if (a->e != b->e) return a->e > b->e;
     else {
        if (a->s != b->s) return a->s > b->s;
           else return a->i < b->i;
     }
   }
};

void Dt::ETable_lsjini::get_lev() { 
//map all levels found in *.lsj onto an array of energies

   map<tmp_lv*,Level*,Lev_map_comp> mpE;
   list<Level*> lv_del;
   for(JI jp = lsj.begin(); jp != lsj.end(); jp++) {
// for initial state:
      
       double E = (*jp)->ei;
       int J = (*jp)->ji;
       string cs = (*jp)->ci;
       Level* ni = new Level(cs,J,E);
       lv_del.push_back(ni);
// for the final state
       E = (*jp)->ef;
       J = (*jp)->jf;
       cs = (*jp)->cf;
       Level* nf = new Level(cs,J,E);
       lv_del.push_back(nf);
   }
   for (LI lp = lv_del.begin(); lp != lv_del.end(); lp++) {
      //double E = (*lp)->totalE;
     tmp_lv* E = new tmp_lv;
     E->e = (*lp)->totalE;
     E->i = (*lp)->J;
     E->s =  (*lp)->config + "_" + (*lp)->term;
     mpE[E] = *lp;
   }
   typedef map<tmp_lv*,Level*,Lev_map_comp>::const_iterator MCI;
   for (MCI mp = mpE.begin(); mp != mpE.end(); mp++) {
      Level* lp = new Level(mp->second);
      lv.push_back(lp); 
   }
}

void Dt::ETable_lsjini::comp_lev() { 
   LI p = lv.begin();
   if (!Base_table::E_min) E_min = (*p)->totalE;
   else E_min = Base_table::E_min;
   for (p = lv.begin(); p != lv.end(); p++)
      (*p)->lev = ((*p)->totalE - E_min) * 2 * RZ;
}

void Dt::ETable_lsjini::fZ(double& ZZ, double& RZ) { 
   ZZ = Print_table::_Z_;
   NEL = Print_table::_NEL_;
   Base_table::fZ(Z,RZ);
   set_ZNEL();
}
void Dt::ETable_lsjini::set_ZNEL() { 
   Print_table::_RZ_ = RZ;
}


//----------------------------------------------------
// ---  member functions of Dt::ETable_jini --
//----------------------------------------------------


// copy constructor -
Dt::ETable_jini::ETable_jini(const ETable_jini& t) 
           : RZ(t.RZ), Z(t.Z), dE(t.dE), NEL(t.NEL), iJ(t.iJ) ,
             file_name(t.file_name) {}
             
// construct table of transitions from file.j
// the constructor does it all for levels
Dt::ETable_jini::ETable_jini(const std::string& ss) : 
file_name(ss), E_min(0) {
   std::cout << "\n...reading " << file_name << "...\n\n";
   get_lev();
   std::cout << "\n... " << file_name
             << " reading OK...\n\n";
   sort_e();
   sort_t();
   rep_TC();
// set global _Z_ _NEL_
   fZ(Z,RZ);
//   set_ZNEL();
   comp_lev();
   splitt();
   print_lev(std::cout);
   read_lsj(ss);
}

// static integer counting how many levels are initially,
// used to specify what to remove
int Dt::ETable_jini::n_lev = 0;

// this function starts reading the *.j file, calls other members
void Dt::ETable_jini::get_lev() {
   std::ifstream j_in(file_name.c_str());
   std::string sw1("_123456789_123456789_123456789_123456789");
   std::string sw2("_123456789_123456789_123456789_123456789");
   std::cout << setw(80) << sw1 + sw2 << endl;
   const int sz = buf;
   static char t_read[sz];
   while(j_in.getline(t_read,sz)) {
      string test_end(t_read);
      bool end_j(false);
      string::size_type szt_e = test_end.find("END");
      string::size_type szt_s = test_end.find("*");
      if (szt_s > 2) szt_s = string::npos;
      if (szt_e != string::npos || szt_s != string::npos) end_j = true;
      if (end_j) { j_in.getline(t_read,sz); }
      if (read_a(t_read))
         for (int i = 0; i < 3; i++) j_in.getline(t_read,sz);
      else { j_in.getline(t_read,sz); j_in.getline(t_read,sz); }
      int nr = read_b(t_read);
      for (int i = 0; i < nr; i++) {
         j_in.getline(t_read,sz);
         read_c(t_read);
         j_in.getline(t_read,sz);
         read_d(t_read);
         int ii = ncfg/7;
         int ad = ii * 78;
         j_in.seekg(ad, std::ios::cur);
         if(ncfg%7) { j_in.getline(t_read,sz); }
      }
   }
}

// read the line with Z=...
bool Dt::ETable_jini::read_a(const char* p) {
   char p1[8], p2[4], p3[6], p4[8], a[2];
   static char rd[buf];
   for (int i = 0; i < buf; i++) rd[i] = p[i];
   std::istrstream in(rd);
   in >> p1 >> p2 >> a >> Z >> p3 >> a >> NEL >> p4 >> a >> ncfg;
   if (p[10] == 'Z') return true;
   return false;
}

// read the line with J, etc
int Dt::ETable_jini::read_b(char* p) {
   static char rd[buf];
   for (int i = 0; i < buf; i++){rd[i] = p[i];}
   int nJ(0);
   char p1[7], p2[9], a[2];
   std::istrstream in(rd);
   in >> p1 >> a >> iJ >> p2 >> a >> nJ;
   cout << "\n\n Labels for 2J = " << iJ << "\n";
   return nJ;
}

//read the line with Ssms, g_JTR, etc
void Dt::ETable_jini::read_c(char* p) { 
   static char rd[buf];
   for (int i = 0; i < buf; i++) rd[i] = p[i];
   char p1[6], p2[6], p3[6];
   std::istrstream in(rd);
   in >> p1 >> dSsms >> p2 >> dg_J >> p3 >> dg_JLS;
}

//read the line with the energy and configuration

std::string getOdd(std::string& conf);

void Dt::ETable_jini::read_d(char* p) {
   n_lev++;
   static char rd[buf];
   for (int i = 0; i < buf; i++) rd[i] = p[i];
   std::string configuration;
   int nn;
   std::istrstream in(rd); 
   in >> nn >> dE  >> configuration;
   configuration = getOdd(configuration);
   if (dE < E_min) E_min = dE;
   cout.setf(ios::fixed, ios::floatfield);
   cout.setf(ios_base::left);
   cout.precision(8);
   cout << " " << setw(30) << configuration << dE << endl;
   Level* nu = new Level(configuration,iJ,dE,0,0,dSsms,dg_J,dg_JLS);
   lv.push_back(nu);
//   std::cout << nu << endl;
}

void Dt::Base_table::fZ(double& ZZ, double& dRZ) {
//   using namespace Dt;
   double zmu = 0;  
   if (ZZ == 1) zmu = 1;
   else if ( ZZ > 10) {zmu = 2 * ZZ + 1 + (ZZ - 11)/2;}
   else if (!(int(ZZ)%2) || (ZZ == 7)) {zmu = 2 * ZZ;}
   else {zmu = 2 * ZZ + 1;}
   dRZ = RZN/(1+RZD/zmu);
}

// compute RZ as funcion of Z
void Dt::ETable_jini::fZ(double& ZZ, double& RZ) {
//   using namespace Dt;
   Base_table::fZ(Z,RZ);
   set_ZNEL();
}

// compute levels using the factor RZ and energy diff
void Dt::ETable_jini::comp_lev() {
   for (LI p = lv.begin(); p != lv.end(); p++)
      (*p)->lev = ((*p)->totalE - E_min) * 2 * RZ;
}


