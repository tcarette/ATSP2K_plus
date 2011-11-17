*
*     ------------------------------------------------------------------
*       B L O C K    D A T A
*     ------------------------------------------------------------------
*
      BLOCK DATA INITT
*
      COMMON/KRON/IDEL(10,10)
      COMMON/TERMS/NROWS,I(24),J(24),N(333)
*
* --- SETS QUANTUM NUMBERS OF TERMS WHICH CAN BE FORMED FROM
*     CONFIGURATIONS  L**Q . ONLY THE FIRST HALF OF THAT PART OF THE
*     TABLE, CORRESPONDING TO A GIVEN  L, IS INCLUDED, BECAUSE OF THE
*     SYMMETRY OF THE TABLE.  E.G. D**7 FORMS THE SAME TERMS AS D**3
*
*       The tables are set for a maximum value of L=9; the terms
*     for L>3 are assumed to be the same as those for L=3
*
*     S - SHELLS (ROWS 1 AND 2)
*
*     P - SHELLS (ROWS 3 TO 5)
*
*     D - SHELLS (ROWS 6 TO 10)
*
*     F - SHELLS (ROWS 11 AND 12)
*
*     G - SHELLS (ROWS 13 AND 14)
*
*     H - SHELLS (ROWS 15 AND 16)
*
*     I - SHELLS (ROWS 17 AND 18)
*
*     K - SHELLS (ROWS 19 AND 20)
*
*     L - SHELLS (ROWS 21 AND 22)
*
*     M - SHELLS (ROWS 23 AND 24)
*
      DATA NROWS/24/
*
*     THE ARRAYS I,J,N CORRESPOND TO THE ARRAYS ITAB,JTAB,NTAB
*
      DATA I/1,1, 1,3,3, 1,5,8,16,16, 1,7, 1,7, 1,7, 1,7, 1,7, 1,7, 1,7/
      DATA J/0,3, 6,9,18, 27,30,45,69,117, 165,168, 189,192,
     : 213,216, 237,240, 261,264, 285,288, 309,312/
      DATA N/1,1,2,  0,1,1,  1,3,2,  0,1,1, 2,5,1, 2,3,3, 1,3,2,
     : 3,5,2, 3,1,4,  1,5,2,  0,1,1, 2,5,1, 2,9,1, 2,3,3, 2,7,3,
     : 1,5,2,  3,3,2, 3,5,2, 3,7,2, 3,9,2, 3,11,2, 3,3,4, 3,7,4,
     : 0,1,1, 2,5,1, 2,9,1, 2,3,3, 2,7,3, 4,1,1, 4,5,1, 4,7,1, 4,9,1,
     : 4,13,1, 4,3,3, 4,5,3, 4,7,3, 4,9,3, 4,11,3, 4,5,5,
     : 1,5,2, 3,3,2, 3,5,2, 3,7,2, 3,9,2, 3,11,2, 3,3,4, 3,7,4, 5,1,2,
     : 5,5,2, 5,7,2, 5,9,2, 5,13,2, 5,5,4, 5,9,4, 5,1,6,
     : 1,7,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,9,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,11,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,13,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,15,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,17,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1,
     : 1,19,2,  2,3,3, 2,7,3, 2,11,3 ,0,1,1, 2,5,1, 2,9,1 ,2,13,1/
*
* --- READ IN OTHER INITIALIZATION DATA
*
      DATA IDEL/1,10*0,1,10*0,1,10*0,1,10*0,1,10*0,1,10*0,1,10*0,
     :          1,10*0,1,10*0,1/
*
      END
