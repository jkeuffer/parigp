#include "exhaustive_search.h"

int
FqX_equal(GEN x, GEN y) { return gequal(x,y); }


void
init_eigen(struct eigen_ellinit *Edat, GEN a4, GEN a6, GEN h,GEN T, GEN p)
{
  pari_sp ltop = avma;
  GEN RHS  = FqX_rem(mkpoln(4, gen_1, gen_0, a4, a6), h, T, p);
  GEN DRHS = FqX_rem(mkpoln(3, utoi(3), gen_0, a4), h, T, p);
  GEN lambda = FqXQ_div(DRHS, FqX_mulu(RHS, 4, T, p), h, T, p);
  GEN C = FqX_sub(FqXQ_mul(lambda, DRHS, h, T, p), monomial(gen_2,1,0), T, p);
  GEN D = FqXQ_mul(FqX_mulu(lambda, 2, T, p),FqX_sub(pol_x(0), C, T, p), h, T, p);
  GEN X12 = mkvec2(C, FqX_Fq_add(D, gen_m1, T, p));
  GEN Gr = T ? FpXQXQ_halfFrobenius(RHS, h, T, p):
               FpXQ_pow(RHS, shifti(p, -1), h, p);
  GEN nGr = FqX_neg(Gr, T, p);
  gerepileall(ltop, 5, &RHS, &DRHS, &X12, &Gr, &nGr);
  Edat->a4    = gcopy(a4);
  Edat->h     = gcopy(h);
  Edat->T     = T;
  Edat->p     = p;
  Edat->pp    = 0;
  Edat->RHS   = RHS;
  Edat->DRHS  = DRHS;
  Edat->X12   = X12;
  Edat->Gr    = Gr;
  Edat->nGr   = nGr;
  Edat->O     = mkvec2(pol_x(0), pol_1(0));
}


GEN
eigen_elldbl(void *E, GEN P)
{
  pari_sp ltop = avma;
  struct eigen_ellinit *Edat=(struct eigen_ellinit *)E;
  GEN T = Edat->T, p = Edat->p, h = Edat->h, x, y;
  if (ell_is_inf(P)) return gcopy(P);
  x = gel(P,1), y = gel(P,2);
  if (FqX_equal(x, pol_x(0)) && FqX_equal(y, pol_1(0)))
    return Edat->X12;
  else
  {
    GEN t1 = FqX_Fq_add(FqX_mulu(FqXQ_sqr(x,h,T,p),3,T, p), Edat->a4, T, p);
    GEN t2 = FqXQ_mul(FqX_mulu(y, 2, T, p), Edat->RHS, h, T, p);
    GEN lambda = FqXQ_div(t1, t2, h, T, p);
    GEN C = FqX_sub(FqXQ_mul(FqXQ_sqr(lambda, h, T, p), Edat->RHS, h, T, p),
                    FqX_mulu(x, 2, T, p), T, p);
    GEN D = FqX_sub(FqXQ_mul(lambda, FqX_sub(x, C, T, p), h, T, p), y, T, p);
    return gerepilecopy(ltop, mkvec2(C,D));
  }
}

/* Returns the addition of [P[1], P[2]*Y] and of [Q[1], Q[2]*Y]
 * Computations are done modulo Y^2 - (X^3 + a4X + a6)
 * An inversion is equivalent to 4M, so that this function requires about 7M
 * which is the same as with the method using ell-division polynomials
 * Working in mixed projective coordinates would require 11M */
 GEN
eigen_elladd(void *E, GEN P, GEN Q)
{
  pari_sp ltop = avma;
  struct eigen_ellinit *Edat=(struct eigen_ellinit *)E;
  GEN Px, Py, Qx, Qy;
  GEN T = Edat->T, p = Edat->p, h = Edat->h, lambda, C, D;
  if (ell_is_inf(P)) return gcopy(Q);
  if (ell_is_inf(Q)) return gcopy(P);
  Px = gel(P,1); Py = gel(P,2);
  Qx = gel(Q,1); Qy = gel(Q,2);
  if (FqX_equal(Px, Qx))
  {
    if (FqX_equal(Py, Qy))
      return eigen_elldbl(E, P);
    else
      return ellinf();
  }
  lambda = FqXQ_div(FqX_sub(Py, Qy, T, p), FqX_sub(Px, Qx, T, p), h, T, p);
  C = FqX_sub(FqX_sub(FqXQ_mul(FqXQ_sqr(lambda, h, T, p), Edat->RHS, h, T, p), Px, T, p), Qx, T, p);
  D = FqX_sub(FqXQ_mul(lambda, FqX_sub(Px, C, T, p), h, T, p), Py, T, p);
  return gerepilecopy(ltop, mkvec2(C,D));
}





/*Finds the eigenvalue of the Frobenius given E, ell odd prime, h factor of the
 *ell-division polynomial, p and tr the possible values for the trace
 *(useful for primes with one root)*/
 ulong
find_eigen_value(GEN a4, GEN a6, ulong ell, GEN h, GEN p)
{ 
   pari_sp ltop = avma;
   GEN BP, Dr;
   ulong t;
   struct eigen_ellinit Edat;
   init_eigen(&Edat, a4, a6, h, NULL, p);
   Dr = BP = Edat.O;
   
  /*[0,Gr], BP, Dr are not points on the curve. */
  /*To obtain the corresponding points, multiply the y-coordinates by Y */
	 pari_sp btop = avma;
	 for (t = 1; t <= (ell>>1); t++)
	   {
	     if (gequal(gel(Dr,2), Edat.Gr))  { avma = ltop; return t; }
	     if (gequal(gel(Dr,2), Edat.nGr)) { avma = ltop; return ell-t; }
	     Dr = eigen_elladd(&Edat, Dr, BP);
	     Dr = gerepileupto(btop, Dr);
	   }
	 pari_err_BUG("find_eigen_value");
	 return 0; /* NOT REACHED */
}

