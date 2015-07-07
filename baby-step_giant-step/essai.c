#include <pari/pari.h>
#include <pari/paripriv.h>
#include "bsgs.c"
#include "exhaustive_search.c"

int 
main()
{
  GEN p,a4,a6,fl;
  ulong ell,i;
  pari_init(10000000000,2);
  FILE* data = fopen("test_elkies.txt","r");
  for(i=0;i<5;i++){
  p  = gp_read_stream(data);
  a4 = gp_read_stream(data);
  a6 = gp_read_stream(data);
  ell = itou(gp_read_stream(data));
  fl = gp_read_stream(data);

  pari_timer T;
  timer_start(&T);
  ulong t = bsgs_abs(a4,a6,ell,fl,p);
  pari_printf("( bsgs    search) ell = %lu , eigenvalue = %lu (time = %ldms) \n",ell,t,timer_delay(&T)); 
  t = find_eigen_value(a4,a6,ell,fl,p);
  pari_printf("(exhaust  search) ell = %lu , eigenvalue = %lu (time = %ldms) \n",ell,t,timer_delay(&T));
  }

  pari_close();  
  return 0;
}
