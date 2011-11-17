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



void Dt::Base_table::paux_lv(ostream& os, LCI& lp) {
   int cfs = Print_table::WC - (*lp)->config_rep.size();
   os.setf(ios::left,ios::adjustfield);
   string j_str; conv_J(j_str,(*lp)->J);

   os << setw(28) << (*lp)->config_rep;
   os << setw(4)  << (*lp)->term_rep;
   os << setw(5)  << j_str << setw(0);

//...
   os.precision(2);
   os.setf(ios::fixed, ios::floatfield);
   os.setf(ios::right,ios::adjustfield);
   double dcmp = 1;
   if ((*lp)->lev > dcmp)
      os  << setw(Print_table::WL-1) << (*lp)->lev << " ";
   else {os << setw(Print_table::WL) << " ";}
  
   if ((*lp)->unph) {
     os.precision(5);
     os.setf(ios::scientific, ios::floatfield);
     if (fabs((*lp)->g_J) > 10e-8) {
        os << setw(13) << (*lp)->g_JLS << " ";
     }
     if (fabs((*lp)->g_J) > 10e-8) {
        os << setw(13) << (*lp)->g_J;
     }
     else  {os << " " ;}
   }
   os << endl;
}

// << operator for output of a single level
std::ostream&
Dt::operator<<(std::ostream& os, const Dt::Level* lv) {
   int cfs = Print_table::WC - lv->config_rep.size();
   os.setf(ios::left,ios::adjustfield);
   string j_str; conv_J(j_str,lv->J);
   int tw = lv->term_rep.size();

   os << setw(28) << lv->config_rep;
   os << setw(4)  << lv->term_rep;  
   os << setw(5)  << j_str << setw(0);
   os.precision(8);
   if ((lv->totalE+10) > 0) os.precision(10);
   else if ((lv->totalE+100) > 0 ) os.precision(9);
   else if ((lv->totalE+1000) > 0) os.precision(8);
   else os.precision(7);

   os.setf(ios::fixed, ios::floatfield);
   os.setf(ios::right,ios::adjustfield);
   os << lv->totalE;

//...
   os.precision(2);
   if (lv->lev >= 100000) os.precision(1);
   if (lv->lev >= 1000000) os.precision(0);
   
   double dcmp = 1;
   if (lv->lev > dcmp) 
      os  << setw(Print_table::WL-1) << lv->lev << " ";
   else {os << setw(Print_table::WL) << " ";}
   
   if (lv->unph) {
     os.precision(2);
     if (lv->Splitting >= 10000) os.precision(1);
     if (lv->Splitting >= 100000) os.precision(0);
     dcmp = 0.01;
     if (lv->Splitting > dcmp) os << setw(Print_table::WS-1) 
                                   << lv->Splitting << " ";
     else  {os << setw(Print_table::WS) << " ";}
     os.precision(4);
     os.setf(ios::scientific, ios::floatfield);
     if (lv->lifetime > 0.0) os << setw(Print_table::WM-1) << lv->lifetime;
     else  {os << " " ;}
   }
   return os;
}

// operator for output of asingle transition
void Dt::Base_table::paux_jtr(std::ostream& os, JCI& jp) {

/*
Z = 8  O I : O-like  (8 electrons).
------------------------------------------------------------------------------------------------------------------
Multiplet
Terms  gi gk Type         Ei          Ek     E(cm-1)       L(vac)     S        f_ik       A_ki
------------------------------------------------------------------------------------------------------------------


2s(2).2p(4)3P2                             2s(2).2p(3)4S3_4S.3d
3P  5D  9 15 E1           75.05         74.60      97194.83      1028.86  1.200e-05  3.938e-07  1.4887e+03
3F  3G  9 59 E2
        9  7 M2   1234526.61  3450821.39   2216292.97        45.12 1.636e-02 4.422e-11 1.8630e+02
        9  7 E1   1234526.61  3450821.39   2216292.97        45.12 6.045e-07 4.522e-07 1.9049e+06

*/

   string jstr_i; conv_J(jstr_i,(*jp)->jf);
   string jstr_f; conv_J(jstr_f,(*jp)->ji);
   os.setf(ios::right,ios::adjustfield);
   os << " " << jstr_i << " " <<  jstr_f;
   os.setf(ios::fixed, ios::floatfield);
   os.setf(ios::right,ios::adjustfield);
   os.precision(2);
   os << setw(12) <<  (*jp)->te;
   os << setw(15) <<  (*jp)->angs_v;
   os << setw(15) <<  (*jp)->angs_a;
   os << endl;
}

// output a single JTR transition
std::ostream& Dt::operator<<(std::ostream& os, JCI jp) {

   string si;  conv_J(si,(*jp)->jf);
   string sf;  conv_J(sf,(*jp)->ji);

   string jstr_i = si.substr(0,4);
   string jstr_f = sf.substr(0,4);
   os.setf(ios::right,ios::adjustfield);
   os << jstr_i << " " <<  jstr_f;
   os << "  " << (*jp)->tr_type;
   os.setf(ios::fixed, ios::floatfield);
   os.setf(ios::right,ios::adjustfield);
   os.precision(2);
   os << setw(14) <<  (*jp)->lev_f;
   os << setw(14) <<  (*jp)->lev_i;
   os.precision(2);
   os << setw(14) << (*jp)->te << setw(13) << (*jp)->angs_v;
   os.setf(ios::scientific, ios::floatfield);
   os.precision(3);
   os << setw(10) <<  (*jp)->_sL;  
//   os << setw(10) <<  (*jp)->_gfL << " " << (*jp)->_f_ik ;
   os << setw(10) <<  (*jp)->_f_ik ;
   os.setf(ios::scientific, ios::floatfield);
   os.precision(4);
   os << setw(11) <<  (*jp)->_akiL; 
//   os.precision(5);
//   os.setf(ios::fixed, ios::floatfield);
//   os.setf(ios::right,ios::adjustfield);
//   double dd = fabs(2*((*jp)->_sL-(*jp)->_sV)/((*jp)->_sV+(*jp)->_gfV));
//   os << setw(8) << dd; 
   return os;
}

//printing a multiplet
void Dt::JTR_mp::MP_print(ostream& os, bool printconf,Level_MP& lmp) {

// formatting
   if (printconf) {
      os.setf(ios::left,ios::adjustfield);
      os << endl << C_low << " - ";
      os.setf(ios::right,ios::adjustfield);
      os.setf(ios::left,ios::adjustfield);
      os << C_high << endl;
   } 
   const int WT = 3;
   int f1 = WT - T_low.size();
   int f2 = WT - T_high.size();
   os.setf(ios::right,ios::adjustfield);
   os << setw(f1) << T_low << setw(0) << " ";
   os << setw(f2) << T_high << setw(0) <<" ";

/*
Z = 8  O I : O-like  (8 electrons).
------------------------------------------------------------------------------------------------------------------
Multiplet
Terms  gi gk Type         Ei          Ek     E(cm-1)       L(vac)     S        f_ik       A_ki
------------------------------------------------------------------------------------------------------------------


2s(2).2p(4)3P2                             2s(2).2p(3)4S3_4S.3d
3P  5D  9 15 E1           75.05         74.60      97194.83      1028.86  1.200e-05  3.938e-07  1.4887e+03
3F  3G  9 59 E2
        9  7 M2   1234526.61  3450821.39   2216292.97        45.12 1.636e-02 4.422e-11 1.8630e+02
        9  7 E1   1234526.61  3450821.39   2216292.97        45.12 6.045e-07 4.522e-07 1.9049e+06
Z = 12  Mg I : Mg-like  (12 electrons).
2s(2).2p(6).3s(2) -- 2s(2).2p(6).3p_2P.3d
1S 3D  1  3            -0.00    1253422.27    1253422.27        79.78 8.083e-07 3.077e-06 1.0750e+06

       1  3 E1          0.00    1253422.27    1253438.47        79.78 8.083e-07 3.077e-06 1.0750e+06
       1  5 M2          0.00    1253422.27    1260286.95        79.35 8.810e-03 3.942e-11 8.3529e+00


2s(2).2p(6).3s(2) -- 2s(2).2p(6).3p_2P.3d
1S 3P  1  3            -0.00    1285376.03    1285376.03        77.80 1.156e-08 4.515e-08 1.6587e+04

              1  3 E1          0.00    1285376.03    1285392.64        77.80 1.156e-08 4.515e-08 1.6588e+04
 
*/

   string t1 = C_low+"_"+T_low;
   string t2 = C_high+"_"+T_high;

   if (MP_S+MP_f_ik+MP_AKI > 0 && 
          lmp.lev_MPG[t1] && lmp.lev_MPG[t2] && (T_low[0] == T_high[0])) {
      
      os.setf(ios::right,ios::adjustfield);
      os << setw(3) << MP_gi << setw(0) << " "
         << setw(2) << MP_gk << setw(0);
      os << "  E1";
      os.setf(ios::fixed, ios::floatfield);
      os.setf(ios::right,ios::adjustfield);
      os.precision(2);
      os << setw(14) <<  E_low;
      os << setw(14) <<  E_high;
      os.precision(2);
      os << setw(14) << E_tr << setw(13) << MP_WL;
      os.setf(ios::scientific, ios::floatfield);
      os.precision(3);
      os << setw(10) << MP_S ;
      os << setw(10) << MP_f_ik;
      os.setf(ios::scientific, ios::floatfield);
      os.precision(4);
      os << setw(11) << MP_AKI << setw(0);
   }
   
   os << endl; 
   for(JCI jp = lsjm.begin(); jp != lsjm.end(); jp++) {
      os << "        ";

      os.setf(ios::right,ios::adjustfield);
      os << setw(2) << (*jp)->jf+1 << setw(0) << " " <<
            setw(2) << (*jp)->ji+1 << setw(0);
      os << "  " << (*jp)->tr_type;
      os.setf(ios::fixed, ios::floatfield);
      os.setf(ios::right,ios::adjustfield);
      os.precision(2);
      os << setw(14) <<  (*jp)->lev_f;
      os << setw(14) <<  (*jp)->lev_i;
      os.precision(2);
      os << setw(14) << (*jp)->te << setw(13) << (*jp)->angs_v;
      os.setf(ios::scientific, ios::floatfield);
      os.precision(3);
      os << setw(10) <<  (*jp)->_sL;
//      os << setw(11) <<  (*jp)->_gfL << " " << (*jp)->_f_ik ;
      os << setw(10) <<  (*jp)->_f_ik ;
      os.setf(ios::scientific, ios::floatfield);
      os.precision(4);
      os << setw(11) <<  (*jp)->_akiL;
      os << endl;
//     os.precision(5);
//     os.setf(ios::fixed, ios::floatfield);
//     os.setf(ios::right,ios::adjustfield);
//     double dd = fabs(2*((*jp)->_sL-(*jp)->_sV)/((*jp)->_sV+(*jp)->_gfV));
//     os << setw(8) << dd;
   }
}

void Dt::Base_table::JTR_print_MP(std::ostream& os) {
   JCI jt = lsj.begin();

   double E_ground(0);

   if (Base_table::E_min) 
      E_ground = Base_table::E_min;
   else 
      E_ground = (*jt)->ef;

   Level_MP lmp(E_ground,this->lv);

   int n = int(Print_table::_Z_) - 1;
   int diff = n - Print_table::_NEL_;
   if (diff > Print_table::sr - 1) diff = Print_table::sr;
   int el = Print_table::_NEL_;

   if (n >106) // when *.lsj file used, or out of names
      os << n+1 << " = Z! check the atomic number!"  << endl;
   else { // there are atom names
      os << "Z = " << Print_table::_Z_ << "  ";
      if (el) os << Print_table::atom[n]<< " " ;
      if (el) os<< Print_table::roman[diff+1] ;
      if (el) os << " : " <<  Print_table::atom[el-1] << "-like ";
      if (el) os << " (" << el  << " electrons).";
      os << std::endl;
   }

   string f__(103,'-');
   os << f__ << endl; 
   os << "       Multiplet" << endl; 

/*
-----------------------------------------------------------------------------------------------------------------
       Multiplet
Terms  gi gk Type           Ei            Ek      E(cm-1)      L(vac)      S        f_ik       A_ki
-----------------------------------------------------------------------------------------------------------------

2s(2).2p(4)3P2                             2s(2).2p(3)4S3_4S.3d
3P  5D  9 15 E1           75.05         74.60      97194.83      1028.86  1.200e-05  3.938e-07  1.4887e+03
3F  3G  9 59 E2
        9  7 M2   1234526.61  3450821.39   2216292.97        45.12 1.636e-02 4.422e-11 1.8630e+02
        9  7 E1   1234526.61  3450821.39   2216292.97        45.12 6.045e-07 4.522e-07 1.9049e+06
*/


   os << "Terms  g_i g_k Type         E_i          E_k"
         "      E(cm-1)      L(vac)      S        f_ik       A_ki" << endl;
   os << f__ << endl << endl;

   JTR_sortMP();
   bool printconf = false;
   string conf_a, conf_b;
   for(JMP jm = jmp.begin(); jm != jmp.end(); jm++) {

   }
   for(JMP jm = jmp.begin(); jm != jmp.end(); jm++) {
      if (conf_a != (*jm)->C_low || conf_b != (*jm)->C_high) {
         conf_a = (*jm)->C_low;
         conf_b = (*jm)->C_high;
         printconf = true;
      }
      (*jm)->MP_print(os,printconf,lmp);
      printconf = false; 
      os << endl;
   }
   os << f__ << endl;
}

void Dt::Base_table::LSJ_print_bin(ostream& os, int& REF) {
   int NEL = Print_table::_NEL_;
   int Z = int(Print_table::_Z_);
   for (JCI jp = lsj.begin(); jp != lsj.end(); jp++) {
      int i = (*jp)->tr_type.size();
      os.write((char*)&i,sizeof(int));
      os.write((char*)(&(*jp)->tr_type[0]),i);
      os.write((char*)(&NEL),sizeof(int));
      os.write((char*)(&Z),sizeof(int));
      os.write((char*)(&(*jp)->_sL),sizeof(double));
      os.write((char*)(&(*jp)->_sV),sizeof(double));
      os.write((char*)(&(*jp)->_gfL),sizeof(double));
      os.write((char*)(&(*jp)->_gfV),sizeof(double));
      os.write((char*)(&(*jp)->_akiL),sizeof(double));
      os.write((char*)(&(*jp)->_akiV),sizeof(double));
      os.write((char*)(&(*jp)->angs_a),sizeof(double));
      os.write((char*)(&(*jp)->angs_v),sizeof(double));
      os.write((char*)(&(*jp)->te),sizeof(double));
      os.write((char*)(&(*jp)->ef),sizeof(double));
      os.write((char*)(&(*jp)->ei),sizeof(double));
      os.write((char*)(&(*jp)->jf),sizeof(int));
      os.write((char*)(&(*jp)->ji),sizeof(int));
      os.write((char*)(&(*jp)->GE_read),sizeof(double));
      os.write((char*)(&REF),sizeof(int));
      i = (*jp)->cf.size();
      os.write((char*)&i,sizeof(int));
      os.write((char*)(&(*jp)->cf[0]),i);
      i = (*jp)->ci.size();
      os.write((char*)&i,sizeof(int));
      os.write((char*)(&(*jp)->ci[0]),i);
   }
}

void Dt::Base_table::LSJ_print(ostream& os) {
   JCI jtmp = lsj.begin();
   if (jtmp == lsj.end()) return;
// get  ground energy
   double Eground(0);
   cerr.setf(ios::fixed, ios::floatfield);
   cerr.setf(ios::right,ios::adjustfield);
   cerr.precision(8);

   if (!(*jtmp)->GE_read) {
      Eground = (*jtmp)->ef;
        cerr << "\n::Reference energy set to the lowest found: ";
        cerr << setprecision(8) << Eground << endl << endl;
   }
   else {
      Eground = (*jtmp)->GE_read;
      cerr << "\n::Reference energy from the .lsj file: ";
      cerr << setprecision(8) << Eground << endl << endl;
   }

   os.setf(ios::fixed,ios::floatfield);
   os.precision(8);
   os << " Ground Energy = " << setw(20) << Eground;

   time_t t;
   time(&t);
   char* T = ctime(&t);
   os << " Last edited: " << T;

    os << " Ionization LIMIT: Term = ";
    if (termstr_LIMIT != "" || termstr_LIMIT != "NA") {
       os << termstr_LIMIT;
    } else {
       os << "NA";
    }
 
    os << " J = "; 
    if (jstr_LIMIT != "" || termstr_LIMIT != "NA") {
       os << jstr_LIMIT;
    } else {
       os << "NA";
    }

    os << " Energy = ";
    if (IE_limit) {
       cerr << " ... printing IE limit in .lsj...= " << IE_limit << endl;
       os.setf(ios::fixed,ios::floatfield);
       os.precision(8);
       os << IE_limit << endl;
    } else {
       cerr << " ... No Ionization Limit (IL) found, written \"NA\" in .lsj...\n";
       os << "NA" << endl;
    }

   int n = int(Print_table::_Z_) - 1;
   int diff = n - Print_table::_NEL_;
   if (diff > Print_table::sr - 1) diff = Print_table::sr;
   int el = Print_table::_NEL_;
   
   if (n >106) // when *.lsj file used, or out of names 
      os << n+1 << " = Z! check the atomic number!"  << endl;
   else { // there are atom names
      os << " Transitions of ";
      if (el) os << "Z = " << Print_table::_Z_ 
         << "  " << Print_table::atom[n]<< " " ;
      if (!el) os << "Z = " << Print_table::_Z_ ;
      if (el) os<< Print_table::roman[diff+1] ; 
      if (el) os << " : " <<  Print_table::atom[el-1] << "-like ";
      os << " (" << el  << " electrons). " << lsj.size() << std::endl;
   } 
   os << "\n\n";
   
   for (JCI jp = lsj.begin(); jp != lsj.end(); jp++) {
      os.setf(ios::right,ios::adjustfield);
      os << setw(4) << (*jp)->jf; 
      os.setf(ios::fixed,ios::floatfield);
      os.precision(8);
      os << setw(14) << (*jp)->ef << "  "; 
      os.setf(ios::left,ios::adjustfield);
      os << setw(59) <<  (*jp)->cf << endl;
      os.setf(ios::right,ios::adjustfield);
      os << setw(4) << (*jp)->ji;
      os.setf(ios::fixed,ios::floatfield);
      os.precision(8);
      os << setw(14) << (*jp)->ei << "  ";
      os.setf(ios::left,ios::adjustfield);
      os << setw(59) << (*jp)->ci << endl; 
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
}

void Dt::Base_table::JTR_print_aux(ostream& os) {
   int n = int(Print_table::_Z_) - 1;
   int diff = n - Print_table::_NEL_;
   if (diff > Print_table::sr - 1) diff = Print_table::sr;
   int el = Print_table::_NEL_;
   
   if (n >106) // when *.lsj file used, or out of names 
      os << n+1 << " = Z! check the atomic number!"  << endl;
   else { // 
      if (el) os << "Wavelengths of " << Print_table::atom[n]<< " " ;
      if (el) os<< Print_table::roman[diff+1] ; 
      if (el) os << " : " <<  Print_table::atom[el-1] << "-like ";
      if (el) os << " (" << el  << " electrons)." ;
      os << std::endl;
   } 

   string f__(80,'-');
   os << f__ << endl;
   os << "   Multiplet" << endl;
   os << "  J_i  -  J_k      E(cm-1)      "
         "ANGS(VAC)      ANGS(AIR)" << endl;
   os << f__ << endl << endl;
    
   JCI jtmp = lsj.begin();
   if (jtmp == lsj.end()) return;
  
   for(JCI jp = lsj.begin(); jp != lsj.end(); jp++) {
      if (((*jp)->cf != (*jtmp)->cf || 
              (*jp)->ci != (*jtmp)->ci) || jp == lsj.begin()) {   
         jtmp = jp;
         string sc, st;
         JTR_mkterms(st,sc,(*jp)->cf);
// formatting 
         const int WC = 29, WT = 3;
         volatile int f1 = WC - sc.size();
         volatile int f2 = WT - st.size();
         os.setf(ios::left,ios::adjustfield);
         os << endl << sc << setw(f1) << "";
         os.setf(ios::right,ios::adjustfield);
         os << st << setw(f2) << "  -   ";
         JTR_mkterms(st,sc,(*jp)->ci);
         f1 = WC - sc.size();
         f2 = WT - st.size();
         os.setf(ios::left,ios::adjustfield);
         os << sc << setw(f1) << "";
         os.setf(ios::right,ios::adjustfield);
         os << st << setw(f2) <<"" << endl;
//         os << sc << "   " << st << endl;
      }
//>>>>>>>>>> this prints each :jp: transition on :os:
      paux_jtr(os,jp); 
//>>>>>>>>>>
   }
  
   os << endl << f__ << endl;

}

void Dt::Base_table::LV_print_aux(ostream& os) {
// first line:
   int n = int(Print_table::_Z_) - 1;
   int el = Print_table::_NEL_;
   int diff = n - el;
   if (diff > Print_table::sr - 1) diff = Print_table::_NEL_; // out of range
   os << "Z = " << Print_table::_Z_ <<" " 
      << "  Energy levels and Zeeman factors, ";

   if (n >106) // when *.lsj file used, or out of names 
      os << n+1 << " = Z! check the atomic number!"  << endl;
   else { // there are atom names
      os << "Z = " << Print_table::_Z_ <<" " ;
      if(el) os << Print_table::atom[n]<< " " ;
      if(el) os<< Print_table::roman[diff+1] ;
      if(el) if (el) os << " : " <<  Print_table::atom[el-1] << "-like ";
      if(el) os << " (" << el  << " el).";
      os << std::endl;
 }
// next few lines 
   os.setf(ios::fixed, ios::floatfield);
   os << "  Rydberg constant is = " 
      << setprecision(4) << Print_table::_RZ_ << endl;
   string f__(80,'-');
   os << f__ <<  endl;
   os << "Configuration              " << " Term" << " J "
      << "    Levels " << "      g_JLS    " << "       g_J" << endl;
   os << "                              " << "     " << "   "
      << " cm^-1 \n";
   os << f__ << endl;

// print levels
   std::string comp_t = "";
   std::string comp_c = "";
   for(LCI p = lv.begin(); p != lv.end(); p++) {
      bool b1 = false;
      bool b2 = false;
      if (comp_t != ((*p)->term)) {comp_t = (*p)->term; b1 = true;}
      if (comp_c != ((*p)->config)) {comp_c = (*p)->config; b2 = true;}
      if ((*p)->unph) {
          if( b1 | b2) os << endl;
/// <<<<<<<<<< print :p: on :os:
          paux_lv(os,p);
/// <<<<<<<<
      }
   }

   os <<endl << f__ << endl;

}
//#####################################
//#####################################

// function printing *.lsj file 
void Dt::Base_table::JTR_print(ostream& os) {
   JCI jtmp = lsj.begin();
   if (jtmp == lsj.end()) return;
// get  ground energy
   double Eground(0);
  
   if (Base_table::E_min) 
      Eground = Base_table::E_min;
   else   
      Eground = (*jtmp)->ef;

   os.precision(8);
   
   int n = int(Print_table::_Z_) - 1;
   int diff = n - Print_table::_NEL_;
   if (diff > Print_table::sr - 1) diff = Print_table::sr;
   int el = Print_table::_NEL_;
   
   if (n >106) // when *.lsj file used, or out of names 
      os << n+1 << " = Z! check the atomic number!"  << endl;
   else { // there are atom names
      os << "Z = " << Print_table::_Z_ << "  ";
      if (el) os << Print_table::atom[n]<< " " ;
      if (el) os<< Print_table::roman[diff+1] ; 
      if (el) os << " : " <<  Print_table::atom[el-1] << "-like ";
      if (el) os << " (" << el  << " electrons).";
      os << std::endl;
   } 

   string f__(103,'-');
   os << f__ << endl;
   os << "Multiplet" << endl;
/*
Z = 12  Mg I : Mg-like  (12 electrons).
-------------------------------------------------------------------------------------------------------
Multiplet
Ji Jk Type          Ei            Ek     E(cm-1)       L(air)      S        f_ik       A_ki
-------------------------------------------------------------------------------------------------------


2s(2).2p(6).3s(2)            1S  -   2s(2).2p(6).3s_2S.3p         3P
 0  1 E1          0.00     303866.51     303870.44       329.09 8.995e-03 8.303e-03 1.7046e+08
 0  2 M2          0.00     303866.51     334529.05       298.93 7.147e+00 5.980e-10 8.9282e+00
*/


   os << "Ji Jk Type          E_i           E_k    E(cm-1)"
         "       L(air)      S        f_ik       A_ki" << endl;
   os << f__ << endl << endl;
    
// get Z;
   double Z(Print_table::_Z_);
   double RZ(0);
// compute the rydberg constatn
   fZ(Z,RZ);
   for(JCI jp = lsj.begin(); jp != lsj.end(); jp++) {
      if (((*jp)->cf != (*jtmp)->cf || 
              (*jp)->ci != (*jtmp)->ci) || jp == lsj.begin()) {   
         jtmp = jp;
         string sc, st;
         JTR_mkterms(st,sc,(*jp)->cf);
// formatting 
         const int WC = 29, WT = 3;
         volatile int f1 = WC - sc.size();
         volatile int f2 = WT - st.size();
         os.setf(ios::left,ios::adjustfield);
         os << endl << sc << setw(f1) << "";
         os.setf(ios::right,ios::adjustfield);
         os << st << setw(f2) << "   -   ";
         JTR_mkterms(st,sc,(*jp)->ci);
         f1 = WC - sc.size();
         f2 = WT - st.size();
         os.setf(ios::left,ios::adjustfield);
         os << sc << setw(f1) << "";
         os.setf(ios::right,ios::adjustfield);
         os << st << setw(f2) <<"" << endl;
//         os << sc << "   " << st << endl;
      }
      (*jp)->lev_i = ((*jp)->ei - Eground)*RZ*2;
      (*jp)->lev_f = ((*jp)->ef - Eground)*RZ*2;
      (*jp)->te = (*jp)->lev_i - (*jp)->lev_f;
      os <<  jp << endl; 
   }
  
   os << endl << f__ << endl;
}

void Dt::Base_table::LV_print(ostream& os) {
// first line:
   int n = int(Print_table::_Z_) - 1;
   int el = Print_table::_NEL_;
   int diff = n - el;
   if (diff > Print_table::sr - 1) diff = Print_table::_NEL_; // out of range
   os << endl << "Z = " << Print_table::_Z_
              << "  Energy levels and lifetimes for ";

   if (n >106) // when *.lsj file used, or out of names 
      os << n+1 << " = Z! check the atomic number!"  << endl;
   else { // there are atom names
      if (el) os << Print_table::atom[n]<< " " ;
      if (el) os<< Print_table::roman[diff+1] ;
      if (el) os << " : " <<  Print_table::atom[el-1] << "-like ";
      if (el) os << " (" << el  << " electrons).";
      os << std::endl;
 }
// next few lines 
   os.setf(ios::fixed, ios::floatfield);
   double IE_cm = 2*Print_table::_RZ_*(IE_limit - E_min);
   os << "  Rydberg constant is = " 
      << setprecision(4) << Print_table::_RZ_ << endl;
   string f__(80,'-');

   os << f__ <<  endl;

if (E_min || IE_limit) {
   if (E_min) {
      os << "  Ground Energy:   "
         <<  setprecision(8) << E_min << " a.u." << endl;
   }
   if (IE_limit) {
      os << "  Ionization LIMIT: " << termstr_LIMIT << " "  << jstr_LIMIT << " ";
         os << setprecision(8) << IE_limit << " a.u. ("
            << setprecision(2) << IE_cm << " cm-1)" << endl;
   }
   os << f__ <<  endl;
}

   os << "Configuration                " << " Term" << " J "
      << "  Energy Total" << "  Levels " << "Splitting "
      << " Lifetimes" << endl;
   os << "                              " << "     " << "   "
      << "    (a.u.)    " << "  cm-1 " << "  cm-1   "
      << "     s    " << endl;
   os << f__ << endl;

// print levels
   string comp_t = "";
   string comp_c = "";
   for(LCI p = lv.begin(); p != lv.end(); p++) {
      bool b1 = false;
      bool b2 = false;
      if (comp_t != ((*p)->term)) {comp_t = (*p)->term; b1 = true;}
      if (comp_c != ((*p)->config)) {comp_c = (*p)->config; b2 = true;}
      if ((*p)->unph) {
          if( b1 | b2) os << endl;
          os << *p << endl;
      }
   }

   os <<endl << f__ << endl;
}

// print all levels, also the header of the table
void Dt::Base_table::print_lev(ostream& os) {
   string f__(80,'-');
   os << f__ <<  endl << "     Configuration         "
      << "        Term" << "  J " << " Total Energy "
      << " Levels "  << endl << f__;

   static std::string comp_t = "";
   static std::string comp_c = "";
   int n_lev = 0;
   for(LCI p = lv.begin(); p != lv.end(); p++) {
      bool b1 = false;
      bool b2 = false;
      if (comp_t != ((*p)->term)) {comp_t = (*p)->term; b1 = true;}
      if (comp_c != ((*p)->config)) {comp_c = (*p)->config;b2 = true;}
      if( b1 | b2) os << endl;
      os << setw(3) << ++n_lev << ". " << *p << endl;
   }
   os << f__ << endl;
}


