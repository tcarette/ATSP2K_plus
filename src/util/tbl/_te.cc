#include "_tt.hh"
#include "_te.hh"
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


Error_lsj::Error_lsj() {
   FLOG.open("tables.log", std::ios::out | std::ios::app);
   _lsj1 = "\n The correct format is:\n";
   _lsj0 = "0123456789_12345689_123456789_123456789_12345689_123456789_123456789\n";
   _lsj2 = "   0-1620.72271875  2s(2).2p(6).3s(2)_1S\n";
   _lsj3 = "   2-1619.33817124  2s(2).2p(6).3s_2S.3p_3P\n";
   _lsj4 = "  303870.44 CM-1       329.09 ANGS(VAC)       329.09 ANGS(AIR)\n";
   _lsj5 = " E1  S =  8.99519e-03   GF =  8.30277e-03   AKI =  1.70460e+08\n";
   _lsj6 = "          8.95098e-03         8.26196e-03          1.69622e+08\n";
   _lsj7 = " (where the last line (velocity values) is optional)\n\n";
   FLOG << _lsj1 << _lsj0 << _lsj2 << _lsj3 <<
           _lsj4 << _lsj5 << _lsj6 << _lsj0 << _lsj7;
}

