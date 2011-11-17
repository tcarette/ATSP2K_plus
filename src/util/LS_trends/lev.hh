#ifndef _REP_H
#define _REP_H
#include <string>
#include <vector>

//int const sz = 256;

struct Lev {
   struct LS_lev{
      int ind, Z, n;
      string Conf;
      double E;
      LS_lev(string&,int,int);
   };

   struct LSJ_Lev{
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

      LSJ_Lev(string&);
   };
   struct LS_Lev{
      string Ugroup;
      int ind;
      string conf;
      string term;
      string parity;
      string tex_C;
      string tex_T;
      int j2;
      LS_Lev(string&);
   };
   vector<LSJ_Lev*> vec_LSJ_Lev;
   vector<LS_Lev*> vec_LS_Lev;
   Lev();
   void print();
};

typedef std::vector<Lev::LSJ_Lev*>::const_iterator CI_LSJ_Lev;
typedef std::vector<Lev::LS_Lev*>::const_iterator CI_LS_Lev;

#endif //_REP_H
