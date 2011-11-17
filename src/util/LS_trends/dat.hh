

#ifndef _DH
#define _DH


const string th1 = "\\begin{sidewaystable} \n" ;
const string th2 = "\\caption{ Accuracy indicators for allowed \n";
const string th3 = "$Singlet$$-$$Singlet$  and $Triplet$$-$$Triplet$ transitions. \n" ;
const string th4 = "The data show the differences between $S_v$ and $S_l$ in \\%.  \n" ;
const string th5 = "} \n" ;
const string th6 = "\\label{accuracyS} \n";
const string th7 = "\\begin{center} \n";
const string th8 = "\\begin{tabular}{l l l l r r r r r r r r r r r r r r r r r r r} \n";
const string th9 = "\\hline* \n";
const string th10 = "\\multicolumn{4}{c}{Transition} & \n";
const string th11= "\\multicolumn{1}{c}{7}  & \\multicolumn{1}{c}{8}  & \n";
const string th12 = "\\multicolumn{1}{c}{9} & \\multicolumn{1}{c}{10}  & \n";
const string th13 = "\\multicolumn{1}{c}{11} & \\multicolumn{1}{c}{12}  & \n";
const string th14 = "\\multicolumn{1}{c}{13} & \\multicolumn{1}{c}{14}  & \n";
const string th15 = "\\multicolumn{1}{c}{15} & \\multicolumn{1}{c}{16}  & \n";
const string th15a = "\\multicolumn{1}{c}{17} & \\multicolumn{1}{c}{}   \\\\ \n";
const string th16 = "\\hline \n";
string THeader_LS = th1 + th2 + th3 + th4 + th5 + th6 + th7 + th8 + th9 + th10 + th11 + th12
               + th13 + th14 + th15 + th15a + th16;
const string tf1 = "\\hline \n";
const string tf2 = "\\end{tabular} \n";
const string tf3 = "\\end{center} \n";
const string tf4 = "\\end{sidewaystable} \n";
string TFooter_LS = tf1 + tf2 + tf3 + tf4;

#endif
