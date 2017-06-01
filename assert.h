#ifndef NDEBUG
#ifndef stderr
#include <stdio.h>
#endif
#define assert(x) if (!(x)) {fprintf(stderr,"Assertion failed: x \n"); exit(1);}
#else
#define assert(x)
#endif

