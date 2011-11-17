#ifndef _TT_H_
#define _TT_H_
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

namespace Dt {
   using namespace std;
   class Lev_map_comp;
   class MP_Lev_comp;
   struct Level;
   struct Level_mp;
   struct JTR_tr;
   struct JTR_mp;
   struct Base_table;
   struct ETable_jini;
   struct ETable_lsjini;
   struct ETable_levini;
   typedef std::list<Level*>::const_iterator LCI;
   typedef std::list<Level*>::iterator LI;
   typedef std::list<JTR_tr*>::const_iterator JCI;
   typedef std::list<JTR_tr*>::iterator JI;
   typedef std::list<JTR_mp*>::iterator JMP;
   struct tmp_lv;
   class Comp;
   class Comp_mp;
   class Comp_MP;
   class Comp_lsjLU;
   class Comp_lsjUL;
   class Comp_lsjE;
   class Comp_lsjWL;
   class Comp_lsjgF;
   class Comp_lsjSL;
   class Comp_lsj_int;
   class CmpRmJTR;
   class CpRmTbLv;
   class RmvDupl; 
   class Comp_LSJ;
   struct Level_MP;
   bool Un_t(LCI, LCI);
   static const int buf(256);
   std::ostream& operator<<(ostream&, const Level*);
   std::ostream& operator<<(ostream&, JCI);
   const double RZN = 109737.31534;
   const double RZD = 548.697e-06;
   void conv_J(string&,const int);
   //void rdstr_lsj(string&, double&, double&, double&);
   void JTR_mkterms(string&, string&, string&);
   void set_tr_type(string&);
   void fZ(double&,double&);
   struct ST;
}

struct Dt::ST {
   vector<string> VS;
   ST(){}
   void convto_LSJ(list<JTR_tr*>&,double);
   bool check_empty() {
     typedef vector<string>::const_iterator SI;
     for (SI p = VS.begin(); p != VS.end(); p++) {
       if (!(*p).size()) return false;
     }
     return true;
   }
   void print_ST_error() {
     typedef vector<string>::const_iterator SI;
     for (SI p = VS.begin(); p != VS.end(); p++) {
       cerr << (*p) << endl;
     }
   }
   ~ST() { VS.clear();}
};

struct Dt::tmp_lv {
   string s;
   int i;
   double e;
};

struct Dt::Level {
   double GE_read;
   std::string config, config_rep;
   std::string term, term_rep;
   std::string tmp;
   int J, ss;
   bool unph;
   double totalE, lev, lifetime, Ssms, g_J, g_JLS, Splitting;
   Level(): lifetime(0),lev(0),Splitting(0){}
   Level(Level* cp) : config(cp->config), config_rep(cp->config_rep),
                  term(cp->term), term_rep(cp->term_rep),
                  J(cp->J), unph(cp->unph), totalE(cp->totalE),
                  Splitting(cp->Splitting), lifetime(0),lev(0){
   }
   Level(string, int, double, double, double, double, 
         double, double);
/*
   Level(std::string sc = "",  int jj = 0,
         double d1 = 0, double d2 = 0, double d3 = 0,
         double d4 = 0, double d5 = 0, double d6 = 0);
*/
   void mk_term();
};

struct Dt::Level_MP {
   map<string,double>lev_MPE;
   map<string,int>lev_MPG;
   Level_MP(double,std::list<Level*>&);
}; 

struct Dt::JTR_tr {
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
   JTR_tr(double&, double&, double&, string&, string&, int&,
          int&, double&, double&, string&, double&, double&,
          double&, double&, double&, double&,double&);
   void get_terms();
};

struct Dt::JTR_mp {
   struct lev_map {
     string c;
     string t;
     int j2;
     double E;
     lev_map (string sc, string st, int ii, double ee) :
         c(sc), t(st), j2(ii), E(ee) {}
   };
   double E_low, E_high, E_tr;
   string C_low, C_high, T_low, T_high;
   int MP_gi, MP_gk;
   double MP_WL;
   double MP_AKI;
   double MP_f_ik;
   double MP_S;
   string tr_type;
   list<JTR_tr*>lsjm;
   JTR_mp(JCI& p) : C_low((*p)->ci), C_high((*p)->cf),
                    T_low((*p)->ti), T_high((*p)->tf),
                    E_low(0), E_high(0),MP_gi(0),MP_gk(0),
                    MP_WL(0),tr_type("E1"),E_tr(0),
                    MP_f_ik(0), MP_S(0), MP_AKI(0){}
   void MP_add(JTR_tr*);
   void MP_computeE(double,Level_MP&);
   void MP_compute_Aki();
   void MP_print(ostream&,bool,Level_MP&);
};


// abstract Base_table : contains sorting functions
struct Dt::Base_table {
   double E_min;
   double IE_limit;
   string jstr_LIMIT;
   string termstr_LIMIT;
  
// each derived class provides reading and printing
   virtual void get_lev() = 0;
   virtual void comp_lev() = 0;
   virtual void set_ZNEL() = 0;
// read *.lsj
   virtual void read_lsj(const string&);
   std::list<JTR_tr*> lsj; 
   std::list<JTR_mp*> jmp;
// sorting functions in the base class
   std::list<Level*>lv;
   virtual void fZ(double&, double&);
   virtual void sort_e();
   virtual void sort_t();
   virtual void rep_TC();
   virtual void LV_rm_unph();
   virtual void splitt();
   virtual void print_lev(ostream& os);
   virtual void LV_print(ostream&);
   virtual void JTR_rmv();
   virtual bool JTR_unph(JCI&);
   virtual void JTR_sortE();
   virtual void JTR_sortELU();
   virtual void JTR_sortEUL();
   virtual void JTR_sortWL();
   virtual void JTR_sortgF();
   virtual void JTR_sortSL();
   virtual void JTR_sortT();
   virtual void JTR_sortLSJ();
   virtual void JTR_sortMP();
   virtual void JTR_lifetimes();
   virtual void JTR_print(ostream&);
   virtual void LSJ_print(ostream&);
   virtual void LSJ_print_bin(ostream&,int&);
   virtual void JTR_print_aux(ostream&);
   virtual void JTR_multiplet();
   virtual void JTR_print_MP(ostream&);
   virtual void LV_print_aux(ostream&);
   virtual void CpLT();
   virtual void paux_lv(ostream&, LCI&);
   virtual void paux_jtr(ostream&, JCI&);
};

struct Dt::ETable_jini : public Dt::Base_table {
   double RZ, E_min;
   static int n_lev;
   double Z, dE, dSsms,dg_J,dg_JLS;
   int NEL, ncfg, iJ;
//   std::list<Level*>lv;  // in the base class
   const std::string file_name;
   struct LS_term {string cfg, trm;};
// specialized copy contructor, removing unphisical!!
   ETable_jini() {}
   ETable_jini(const ETable_jini&);
//   ETable_jini* clone() const {return new ETable_jini(*this);}
   ETable_jini(const std::string&);
   void get_lev();
   bool read_a(const char*);
   int read_b(char*);
   void read_c(char*);
   void read_d(char*);
   void comp_lev();
//   void fZ();
   void fZ(double&, double&);
   double get_RZ() {return RZ;}
   void set_ZNEL();
};

struct Dt::ETable_lsjini : public Dt::Base_table {
   double RZ, E_min;
   static int n_lev;
   double Z, dE, dSsms,dg_J,dg_JLS;
   int NEL, ncfg, iJ;
   const std::string file_name;
   struct LS_term {string cfg, trm;};
// specialized copy contructor, removing unphisical!!
   ETable_lsjini(const std::string&);
   void get_lev();
   void comp_lev();
//   void fZ();
   void fZ(double&, double&);
   double get_RZ() {return RZ;}
   void set_ZNEL();
};

class Dt::Comp_MP{
public:
   bool operator()(const JTR_mp* ta, JTR_mp* tb) const {
   if (ta->C_low != tb->C_low)
      return (ta->E_low <=tb->E_low);
   if (ta->T_low != tb->T_low )
      return (ta->E_low <=tb->E_low);
   if (ta->C_high != tb->C_high)
      return (ta->E_high <=tb->E_high);
   if (ta->T_high != tb->T_high)
      return (ta->E_high <=tb->E_high);
   return (ta->E_tr >= tb->E_tr);
   }
};

class Dt::Comp_LSJ{
public:
   bool operator()(const JTR_tr* ta,const JTR_tr* tb) const {
   if (ta->cf != tb->cf) return (ta->cf <= tb->cf);
   if (ta->ci != tb->ci) return (ta->ci <= tb->ci);
   if (ta->tf != tb->tf) return (ta->tf <= tb->tf);
   if (ta->ti != tb->ti) return (ta->ti <= tb->ti);
   if (ta->ji != tb->ji) return (ta->ji <= tb->ji);
   if (ta->jf != tb->jf) return (ta->jf <= tb->jf);
   if (ta->ef != tb->ef) return (ta->ef <= tb->ef);
   if (ta->ei != tb->ei) return (ta->ei <= tb->ei);
   if (ta->tr_type != tb->tr_type) return (ta->tr_type <= tb->tr_type);
   return (ta->te <= tb->te);
   }
};

class Dt::RmvDupl {
public:
   bool operator()(const JTR_tr* ta,const JTR_tr* tb) const {
      return (ta->cf == tb->cf &&
              ta->ci == tb->ci &&
              ta->tf == tb->tf &&
              ta->ti == tb->ti &&
              ta->ji == tb->ji &&
              ta->jf == tb->jf &&
              ta->te == tb->te &&
              ta->tr_type == tb->tr_type);
   }
};

class Dt::Comp {
public:
   bool operator()(const Level* ta, const Level* tb) const {
      //return true if ta < tb
      if(ta->totalE == tb->totalE) return false; //ta->totalE <= tb->totalE;
      return ta->totalE < tb->totalE;
   }
};

// make thse two a template
// u. predicate for removing duplicates 
class Dt::CpRmTbLv {
public:
   bool operator()(const Level* ta) const {
      return !(ta->unph);
   }
};

class Dt::CmpRmJTR {
public:
   bool operator()(const JTR_tr* ta) const {
      return !(ta->unph);
   }
};

#endif
 
