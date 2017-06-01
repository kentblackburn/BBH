/* A Bubble Sorting Routine to replace 
   to missing qsort() on the AstroVax. */

/* James Kent Blackburn
   University of Florida
   Department of Physics
   Summer 1988            */

#define TRUE  1
#define FALSE 0

/* NOTE! Arrays numbered from 0 to N-1 */

qsort(array,n,l,f)

 int n,l,(*f)();
 double array[];

{
 int j,change;
 double temp;

 if (l != sizeof(double))
   {
    printf("failed in qsort!\n");
    exit(10);
   }
 else
   {
    do{
       change = FALSE;
       for (j = 0;j < (n-1);j++)
          {
           if (array[j]>array[j+1])
             {
              temp = array[j];
              array[j] = array[j+1];
              array[j+1] = temp;
              change = TRUE;
             }
          }
      }while(change);
   }
}

