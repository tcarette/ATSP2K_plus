#include "rep_lib.hh"
#include "ls.hh"
#include <string>
#include<cmath>


string& ret_str(string& s) {
   int a = s.find_first_not_of(" ");
   s.erase(0,a);
   a = s.find_last_not_of(" ");
   s.erase(a+1);
// cout << s.size() << "::" << s << "::" << endl;
   return s;
}

string& u_format(string& ss, double dr, double du) {
    if (!dr && !du) {
      ss = " & ";
      return ss;
    }

    char t1[50];
    std::ostrstream BUF1(t1,sizeof(t1));
    BUF1.precision(3);
    BUF1.setf(ios::scientific, ios::floatfield);
    BUF1 << dr << ends;
    string s1 = t1;

    char t2[50];
    std::ostrstream BUF2(t2,sizeof(t2));
    BUF2.precision(3);
    BUF2.setf(ios::scientific, ios::floatfield);
    BUF2 << du << ends;
    string s2 = t2;

    char* pe1 = &t1[6];
    char* pe2 = &t2[6];
    istrstream ie1(pe1);
    istrstream ie2(pe2);
    int e1 = 0;
    int e2 = 0;
    ie1 >> e1;
    ie2 >> e2;

    char t3[10];
    std::ostrstream BUF3(t3,sizeof(t3));
    BUF3.precision(3);
    BUF3.setf(ios::scientific, ios::floatfield);
    BUF3 << e1 << ends;
    string stmp = t3;
    string s3("");

    if (e1 < 10 && e1 > 0) s3 = " " + stmp;
    else s3 = stmp;

   if (!du) {
       s1.erase(5,8);
       ss = s1 + "(" + "-" + ")" + "      & " + s3;
       return ss;
    }

    s1.erase(5,8);
    int ediff = e1 - e2;
    if (ediff >= 4) {
       ss = s1 + "(" + "0" + ")" + "      & " + s3;
       return ss;
    }

    double maxd = (dr <= du) ? du : dr;
    double maxpd = fabs((du)/maxd)*100;

    if (maxpd >= 99.99) {
       char t10[50];
       std::ostrstream BUF1(t10,sizeof(t10));
       BUF1.precision(2);
       BUF1.setf(ios::scientific, ios::floatfield);
       BUF1 << dr << ends;
       string s10 = t10;
       s10.erase(4,7);
       ss = s10 + "(" + "100\\%" + ")" + "   & " + s3;
       return ss;
    }

    if (maxpd >= 70) {
       char t4[10];
       char t10[50];
       std::ostrstream BUF1(t10,sizeof(t10));
       BUF1.precision(2);
       BUF1.setf(ios::scientific, ios::floatfield);
       BUF1 << dr << ends;
       string s10 = t10;
       s10.erase(4,7);
       std::ostrstream BUF4(t4,sizeof(t4));
       BUF4.precision(1);
       BUF4.setf(ios::showpoint);
       BUF4.setf(ios::fixed, ios::floatfield);
       BUF4 << maxpd << ends;
       string s4 = t4;
       ss = s10 + "(" + s4 + "\\%" + ")" + "  & " + s3;
       return ss;
    }

    if (!ediff) {
       char t5[10];
       std::ostrstream BUF5(t5,sizeof(t5));
       BUF5.precision(2);
       BUF5.setf(ios::scientific, ios::floatfield);
       BUF5 << dr << ends;
       string s5 = t5;
       s5.erase(4,7);
       char t6[10];
       std::ostrstream BUF6(t6,sizeof(t6));
       BUF6.precision(2-ediff);
       BUF6.setf(ios::showpoint);
       BUF6.setf(ios::scientific, ios::floatfield);
       BUF6 << du << ends;
       string s6 = t6;
       int SZ = s6.size();
       s6.erase(SZ-4,SZ-1);
       s6.erase(1,1);
       if (!ediff)
         ss = s5 + "(" + s6 + ")" + "     & " + s3;
       else ss = s5 + "(" + s6 + ")" + "      & " + s3;
       return ss;
    }
    if (ediff == 3) {
       double DD(10);
      int I = int(du/(pow(DD,e2)));
       char t8[10];
       std::ostrstream BUF8(t8,sizeof(t8));
       BUF8 << I << ends;
       string s8 = t8;
       ss = s1 + "(" + s8 + ")" + "      & " + s3;
       return ss;
    }


    if (ediff == 2 || ediff == 1) {
       char t7[10];
       std::ostrstream BUF7(t7,sizeof(t7));
       BUF7.precision(3-ediff);
       BUF7.setf(ios::showpoint);
       BUF7.setf(ios::scientific, ios::floatfield);
       BUF7 << du << ends;
       string s7 = t7;
       int SZ = s7.size();
       s7.erase(SZ-4,SZ-1);
       s7.erase(1,1);
       if (ediff == 2) ss = s1 + "(" + s7 + ")" + "    & " + s3;
       else ss = s1 + "(" + s7 + ")" + "    & " + s3;
       return ss;
    }

//  #0123456789_123456789_123456789
//  #3.333e+05(4.343e+03)#
//  #3.333(434)  &  5 #
//  #3.33(43)    &  5 #
//  #3.33(10\\%) &  5 #

    ss = " & ";
    return ss;
}


