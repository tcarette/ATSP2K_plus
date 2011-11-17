// This program tabulates LS trends from a single file 
// concatenated for all n's and Z's.

#include "repr.hh"
#include "ls.hh"
using namespace std;

int main() {

  cout << " This program reads all .ls files  and tabulates the data"
          " for all n's and Z's\n"; 
  const char* argv[] = {"","_A_LS"};
  LS ls_tr(argv[1]);
}
