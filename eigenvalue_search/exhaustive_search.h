#ifndef PARI_H
#define PARI_H

#include <pari/pari.h>
#include <pari/paripriv.h>

#endif /*PARI_H*/

struct eigen_ellinit
{
  GEN a4, h, T, p;
  GEN RHS, DRHS, X12, Gr, nGr,O;
  ulong pp;
};

int FqX_equal(GEN x, GEN y);
void init_eigen(struct eigen_ellinit*, GEN , GEN ,GEN, GEN , GEN );
GEN eigen_elldbl(void *, GEN );
GEN eigen_elladd(void *, GEN , GEN );
ulong find_eigen_value(GEN , GEN , ulong , GEN , GEN );

