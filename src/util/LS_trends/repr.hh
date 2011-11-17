#ifndef _REP_H
#define _REP_H
#include <string>
#include <vector>
#include<list>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

//int const sz = 256;

struct Rep {

   double RZN;
   double RZD;


   struct LS_lev{
      int ind, Z, n;
      string Conf;
      double E;
      LS_lev(string&,int,int);
   };

   struct Lev {
      string c;
      double e,ee;
      int z,n;
      Lev(LS_lev* A) : ee(0),n(A->n), z(A->Z), e(A->E), c(A->Conf) {}
      friend ostream& operator<<(ostream&, Lev*);
   };

   map<int,map<int,vector<Lev*> > > Map_LEV;
   void comp_spec();
   void print_spec();

   struct LSJ_Rep{
      string Ugroup;
      int ind;
      float JV;
      string conf_mchf;
      string term_mchf;
      string j_string_mchf;
      string parity_mchf;

      string conf_nist;
      string term_nist;
      string j_string_nist;

      string conf_tex;
      string term_tex;

      LSJ_Rep(string&);
   };

   struct LS_Rep{
      string Ugroup;
      int ind;
      string conf;
      string term;
      string parity;
      string tex_C;
      string tex_T;
      int j2;
      LS_Rep(string&);
      LS_Rep(LS_lev* A) : conf(A->Conf), ind(A->ind) {}
   };
                              
   vector<LSJ_Rep*> vec_LSJ_Rep;
   vector<LS_Rep*> vec_LS_Rep;
   vector<LS_lev*> vec_LS_lev;
   vector<LS_lev*> Vec_l;
   map<int,vector<LS_lev*> > Map_LZ;
   Rep();
   Rep(const char*);
   void print();
   double RZ(int&);
};


class Comp_LS_lev {
public:
   bool operator()(const Rep::LS_lev* ta,
                   const Rep::LS_lev* tb) const {
      //return true if ta < tb
      return  ta->E < tb->E;
   }
};

class CompLev {
public:
   bool operator()(const Rep::Lev* ta,
                   const Rep::Lev* tb) const {
      //return true if ta < tb
      return  ta->e < tb->e;
   }
};

class LevEq {
public:
   bool operator()(const Rep::Lev* ta,
                   const Rep::Lev* tb) const {
      //return true if ta < tb
      return  ta->e == tb->e &&
              ta->c == tb->c &&
              ta->n == tb->n &&
              ta->z == tb->z;
   }
};

class match_Lev : public unary_function<Rep::LS_lev*,bool> {
   int zz, nn;
   string c;  
public:
   explicit match_Lev(const Rep::LS_lev* a) : zz(a->Z),nn(a->n),
                                 c(a->Conf) {}
   bool operator()(const Rep::LS_lev* b) const {
//cout << zz << "::" << b->Z << endl;
//cout << nn << "::" << b->n << endl;
//cout << c << "::" << b->Conf << endl;
      return (zz == b->Z && c == b->Conf);
   }
};


typedef std::vector<Rep::LSJ_Rep*>::const_iterator CI_LSJ_Rep;
typedef std::vector<Rep::LS_Rep*>::const_iterator CI_LS_Rep;
typedef vector<Rep::LS_lev*>::iterator VI;
typedef map<int,vector<Rep::LS_lev*> >::iterator MI;


#endif //_REP_H
