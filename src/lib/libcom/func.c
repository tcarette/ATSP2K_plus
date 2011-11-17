#include <stdlib.h>
#include<stdio.h>

char * mmm_(size_t *n)
{
// printf("%llu \n",*n); 
return(malloc(*n));
}

