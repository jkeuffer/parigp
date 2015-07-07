#ifndef PARI_H
#define PARI_H

#include <pari/pari.h>
#include <pari/paripriv.h>

#endif /*PARI_H*/


GEN compute_multiples(ulong ,GEN , GEN , GEN , GEN );
GEN mod_compose(GEN , GEN , GEN*,ulong , GEN );
GEN mod_transp(GEN ,GEN ,GEN ,GEN ,GEN);
ulong rational_equ(GEN *,GEN *, ulong ,ulong , GEN ,GEN *, GEN );
ulong bsgs_abs(GEN ,GEN ,ulong ,GEN ,GEN );

