#include "repr.hh"
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <iostream>

#ifndef _LS_H
#define _LS_H
using namespace std;

class CompLS_E;
class Comp_LS_ind;
class CompLS_err;

struct LS {
   Rep* R;
   string SEQ;
   void set_SEQ() {
       cerr << " Enter a sequence, (B, Al, Ti) : ";
       getline(cin,SEQ,'\n');
       cerr << "\n You have entered \"" << SEQ << "\", the filenames will be: LS_" << SEQ << "_Znn" << endl << endl;
       if (SEQ.size() > 2 || SEQ.size() == 0) cerr << " Most likely \"" << SEQ << "\" is not a valid sequence " << endl;
   }

   struct LS_tr{
      string UgroupL, UgroupU;
      int ind_i, ind_f;
      int Z, n;
      string Conf_L, Conf_U;
      double EL, EU;
      double LS_trE;
      double w_VAC, w_AIR;
      string TR_type;
      double SL_L, GF_L, AKI_L;
      double SL_V, GF_V, AKI_V;
      double perc_LS_error;
      LS_tr(string&,string&,string&,string&,string&,int,int,map<string,int>&);
   };

   struct LS_err{
      string UgroupL,UgroupU;
      int ZZ;
      int ind_i, ind_f;
      vector<int> vec_Z;
      vector<int> vec_N;
      vector<double> vec_GF;
      vector<float> vec_Perc_err;
      string Lconf, Uconf;
      string Lterm, Uterm;
      string Lconf_tex, Uconf_tex;
      string Lterm_tex, Uterm_tex;
      LS_err(int,int,LS_tr*,Rep*);
      void set_LS_err(int,int,float);
      typedef std::vector<LS::LS_err*>::const_iterator CI;
//      friend ostream& operator<<(ostream&, CI);
   };
   int nZ; 
   int Max_n;
   int Zdiff;
   vector<LS_err*> LS_ERR;   
   LS(const char[]);
   ~LS() {}
//   void LS_print_tex();
   void LS_mk_Err(map<int,vector<LS_tr*> >&);
   void LS_print_dat(map<int,vector<LS_tr*> >&, const char*);
   int ret_ind(string&);
};

class Find_ind_LS_tr : public unary_function<LS::LS_err*,bool> {
   int ii,fi;
public:
   explicit Find_ind_LS_tr(const int a, const int b) : ii(a), fi(b) {}
   bool operator()(const LS::LS_err* ta) const {
      return (ii == ta->ind_i && fi == ta->ind_f);
   }
};

class Comp_LS_ind {
public:
   bool operator()(const LS::LS_tr* ta,
                   const LS::LS_tr* tb) const {
      //return true if ta < tb
      if (ta->ind_i != tb->ind_i ) return ta->ind_i < tb->ind_i; 
      else if (ta->ind_f != tb->ind_f ) return ta->ind_f < tb->ind_f;
      else if (ta->n != tb->n ) return ta->n < tb->n;
      else return ta->perc_LS_error < tb->perc_LS_error; 
   }
};

class LS_trEq {
public:
   bool operator()(const LS::LS_tr* ta,
                   const LS::LS_tr* tb) const {
      //return true if ta < tb

      return  ta->Conf_L == tb->Conf_L    &&
              ta->Conf_U == tb->Conf_U    &&
              ta->EU == tb->EU            &&
              ta->LS_trE == tb->LS_trE    &&
              ta->TR_type == tb->TR_type  &&
              ta->SL_L == tb->SL_L        &&
              ta->SL_V == tb->SL_V        &&
              ta->EL == tb->EL;
   }
};


class Comp_LS_err {
public:
   bool operator()(const LS::LS_err* ta,
                   const LS::LS_err* tb) const {
      //return true if ta < tb
      if (ta->ind_i != tb->ind_i ) return ta->ind_i < tb->ind_i;
      else return ta->ind_f <= tb->ind_f;
   }
};

ostream& operator<<(ostream&, LS::LS_err*);
typedef std::vector<Rep::LSJ_Rep*>::const_iterator CI_LSJ_Rep;
typedef std::vector<Rep::LS_Rep*>::const_iterator CI_LS_Rep;
typedef std::vector<LS::LS_err*>::const_iterator CI_LSerr;


#endif
