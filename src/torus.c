#include <mainheader.h>


int torus(int m, int k)
{/*torus*/
/** *******************************************************
k: it is a grid point that is determined by torus()
m: size of matrix 
***********************************************************/
  return ((m-((m-k)%m)))%m;
}/*torus*/
