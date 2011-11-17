#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <list>
#include <cstring>
#include <sstream>

using namespace std;

string getOdd(string& conf);
void appendOdd(string& s);
int getL(char c);
int getNC(string& p);
bool checkOdd(vector<string>& vs);
bool getl(char c);

int getL(char c) {
   map<char,int> mapL;
   mapL['s'] = 0;
   mapL['p'] = 1;
   mapL['d'] = 2; 
   mapL['f'] = 3;
   mapL['g'] = 4;
   mapL['h'] = 5;
   mapL['i'] = 6;
   return mapL[c];
}

void appendOdd(string& s) {

   typedef vector<string*>::const_iterator LSI;

   ifstream in_lsj(s.c_str());
   string tmp = s + ".tmp.tables";
   ofstream OUT(tmp.c_str());

   if(in_lsj) cout << "...checking for odds " << s << "...\n";
   vector<string*> SL;
   const int sz = 4096;
   static char t_read[sz];

   while(in_lsj.getline(t_read,sz)) {
      string READ_LINE(t_read);
      string* A = new string(t_read);
//cout << "::" << READ_LINE << "::" << endl;
      OUT << *A << endl;
      SL.push_back(A);
   }
   OUT.close();
   in_lsj.close();

   for (int i = 0; i < SL.size(); i++) {
     string A = *SL[i];
     if (A[5] == 'S' && i > 3) {
         *SL[i-3] = getOdd(*SL[i-3]);
         *SL[i-2] = getOdd(*SL[i-2]);
//         cout << *SL[i-3] << endl;
//         cout << *SL[i-2] << endl;
     } 
   }

   ofstream OUTF(s.c_str());
   for (int ii = 0; ii < SL.size(); ii++) OUTF << *SL[ii] << endl;      

   SL.clear();
}

string getOdd(string& conf) {
   string::size_type x(0);
   if (conf.size() >= 20) x = 20;
   else x = 0;
   string::size_type y(0);
   vector<string> VS;
   while (y != string::npos) {
      y = conf.find(".",x);
      if (y != string::npos) {
         VS.push_back(conf.substr(x,y-x));
         x = y + 1;
      }
   }
   VS.push_back(conf.substr(x));
   string::size_type oa = conf.find_first_of("*",x); // check if odds are added
   if (oa != string::npos) conf[oa] = ' ';  // odds already added
   string::size_type o = conf.find_first_of(" ",x); // first the frst space
   char o_1 = conf[o-1];

   if (isdigit(o_1)) {
      o -= 1; 
   }
   if (checkOdd(VS)) {
       if (o != string::npos)
            conf[o] = '*'; // insert and '*' immediatel after the term
       else conf += "*";
   }
   return conf;
}

bool checkOdd(vector<string>& vs) {

   int nc(0);
   for (int i = 0; i < vs.size(); i++) {
      nc += getNC(vs[i]);
   }
   if (nc%2) return true;
   return false;
}

int getNC(string& p) {
   int n(1);
   int l(0);
   string::size_type bb = p.find("(");
   string::size_type eb = p.find(")");
   if (eb == string::npos || bb == string::npos) n = 1;
   else {
     string b = p.substr(bb+1,eb-bb);
     istringstream ist(b);
     ist >> n;
   } 
   for (int i = 0; i < p.size(); i++) {
      if (getl(p[i])) {
         l = getL(p[i]); 
//cout << "found match for l=" << p[i] << " n = " << l << endl;
         break;
      }
   }
   return n*l;
}

bool getl(char c) {
   if (c == 's') return true;
   if (c == 'p') return true;
   if (c == 'd') return true;
   if (c == 'f') return true;
   if (c == 'g') return true;
   if (c == 'h') return true;
   return false; 
}

//int main() {
//   string fn = "U.lsj";
//   appendOdd(fn);
//}

