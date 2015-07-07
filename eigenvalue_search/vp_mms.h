#ifndef PARI_H
#define PARI_H

#include <pari/pari.h>
#include <pari/paripriv.h>

#endif /*PARI_H*/


ulong  compute_gen(ulong );
ulong choose_q(ulong );
GEN mms_mod_compose(GEN , GEN , GEN*, GEN);
GEN compo_t(GEN , ulong , GEN*,GEN );
GEN vec_compo(GEN ,ulong , GEN*,GEN );
GEN gx(GEN ,GEN ,ulong ,GEN *,GEN );
GEN gy(GEN ,GEN ,ulong ,GEN *,GEN );
GEN T0X(GEN ,GEN,ulong ,ulong ,ulong ,ulong ,GEN*,GEN);
GEN T1X(GEN ,GEN ,ulong ,GEN ,ulong ,ulong ,GEN *,GEN );
GEN trace_calc(GEN ,GEN );
GEN mms_mod_transp(GEN ,GEN ,GEN ,GEN ,GEN );
GEN compute_trace(GEN ,GEN ,GEN ,GEN);
GEN compute_C(GEN ,GEN,GEN,GEN,GEN*,GEN);
long bsgs_t(GEN ,GEN* ,GEN ,ulong ,GEN);
long log_eigenvalue(GEN , GEN ,ulong ,ulong ,ulong ,GEN ,GEN*,GEN);
ulong find_eigenvalue(GEN ,GEN,ulong,GEN,GEN);
