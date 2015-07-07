/*incremental computation of multiples of a point*/
/*returns a vector containing n numerators and denominators of the point multiples*/

static GEN
compute_multiples(ulong n,GEN a4, GEN a6, GEN fl, GEN p)
{
  pari_sp ltop = avma;
  GEN mul;
  if(n==1){
    mul = cgetg(3,t_VEC);
    gel(mul,1) = pol_x(0); gel(mul,2) = pol_1(0);
    return mul;
  }
  if(n==2){
    mul = cgetg(5,t_VEC);
    gel(mul,1) = pol_x(0); gel(mul,2) = pol_1(0);
    gel(mul,3) = FpX_sub(FpXQ_mul(pol_x(0),FpX_mulu(mkpoln(4,gen_1,gen_0,a4,a6),4,p),fl,p),FpX_rem(Fp_elldivpol(a4,a6,3,p),fl,p),p); 
    gel(mul,4) = FpX_rem(FpX_mulu(mkpoln(4,gen_1,gen_0,a4,a6),4,p),fl,p);
    return mul;
  }
  ulong i,r;
  (n&1)?(r=(n+1)>>1):(r=n>>1);
  GEN eF = FpX_get_red(fl,p);
  GEN R = FpX_mulu(mkpoln(4,gen_1,gen_0,a4,a6),4,p);
  GEN R2 = FpXQ_sqr(R,eF,p);  
  GEN F = cgetg(n+3,t_VEC);
  GEN U = cgetg(n+1,t_VEC);
  GEN V = cgetg(n+1,t_VEC);
  gel(F,1) = pol_1(0); gel(F,2) = pol_1(0);gel(F,3) = Fp_elldivpol(a4,a6,3,p);
  gel(F,n+2) = pol_0(0);/*needed if n is even*/
  gel(V,1) = pol_1(0); gel(V,2) = pol_1(0);
  gel(V,3) = FpXQ_sqr(Fp_elldivpol(a4,a6,3,p),eF,p);  
  gel(U,1) = pol_0(0); gel(U,2) = FpX_rem(Fp_elldivpol(a4,a6,3,p),eF,p);
  gel(U,3) = Fp_elldivpol(a4,a6,4,p);
  for(i=2;i<=r;i++){
    gel(F,i<<1)=FpX_sub(FpXQ_mul(gel(V,i-1),gel(U,i+1),eF,p),FpXQ_mul(gel(V,i+1),gel(U,i-1),eF,p),p);
    if(i&1){
      gel(F,(i<<1)+1)=FpX_sub(FpXQ_mul(gel(V,i),gel(U,i+1),eF,p),FpXQ_mul(FpX_mul(gel(V,i+1),R2,p),gel(U,i),eF,p),p);
    } else{
      gel(F,(i<<1)+1)=FpX_sub(FpXQ_mul(FpX_mul(R2,gel(V,i),p),gel(U,i+1),eF,p),FpXQ_mul(gel(V,i+1),gel(U,i),eF,p),p);
    }
    gel(V,i+2) = FpXQ_sqr(gel(F,i+2),eF,p);
    gel(U,i+2) = FpXQ_mul(gel(F,i+1),gel(F,i+3),eF,p);
  }
  for(i=r+3;i<=n;i++){
    gel(V,i) = FpXQ_sqr(gel(F,i),eF,p);
    gel(U,i) = FpXQ_mul(gel(F,i-1),gel(F,i+1),eF,p);
  }
  GEN multiples = cgetg(2*n+1,t_VEC);
  for(i=1;i<=n;i++){
    if(i&1){
      gel(multiples,2*i-1) = FpX_sub(FpXQ_mul(gel(V,i),pol_x(0),eF,p),FpXQ_mul(R,gel(U,i),eF,p),p);
      gel(multiples,2*i) = gel(V,i);
    } else {
      gel(multiples,2*i-1) = FpX_sub(FpXQ_mul(R,FpX_mul(gel(V,i),pol_x(0),p),eF,p),gel(U,i),p);
      gel(multiples,2*i) = FpXQ_mul(R,gel(V,i),eF,p);
    }
  }
  return gerepilecopy(ltop,multiples);
}

/************* Modular composition :returns g(h) mod f  **************/
/* assumes that g,h,f are FpX's and  deg(f)>deg(g) deg(f)>deg(h) */

static GEN
mod_compose(GEN g, GEN h, GEN *eT,ulong n, GEN p)
{
  /*pari_sp ltop = avma;*/
  GEN tmp1,res;
  ulong i,j; 
  ulong m = usqrt(n)+1;
  tmp1 = h;
  GEN A = zeromatcopy(m,n),B = zeromatcopy(m,m);
  gcoeff(A,1,1)=gen_1;
  for(i=1;i<m;i++){
    for(j=0;j<=degree(tmp1);j++){gcoeff(A,i+1,j+1)=polcoeff0(tmp1,j,0);};
    tmp1 = FpXQ_mul(tmp1,h,*eT,p);
  }
  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      gcoeff(B,i+1,j+1) = polcoeff0(g,i*m+j,0);
    }
  }
  A = FpM_mul(B,A,p);
  res = gtopolyrev(row(A,m),0);
  for(i=1;i<m;i++){
    res = FpX_add(FpX_mul(res,tmp1,p),gtopolyrev(row(A,m-i),0),p);
  }
  /*return gerepileupto(ltop,FpX_rem(res,*eT,p));*/
  return FpX_rem(res,*eT,p);
}


/* Rational equality */

/* mod_transp computes tau o v ie the transpose of the matrix  */
/* representing the "multiplication by tau" in Fp[x]/(f) */
static GEN
mod_transp(GEN tau,GEN v,GEN h,GEN f,GEN p)
{
  /*pari_sp ltop = avma;*/
  GEN f_tilde = polrecip(f);
  long n = poldegree(f,0),degtau = degpol(tau),i;
  /*GEN a = gtopolyrev(v,0);*/
  GEN a  = RgV_to_RgX(v,0);
  GEN b_tilde;
  if(degtau==n-1){
  b_tilde = polrecip(tau);
  } else{
    b_tilde = vec_lengthen(gtovecrev(tau),n);
    for(i=degtau+2;i<=n;i++) 
      gel(b_tilde,i)=gen_0;
    b_tilde = gtopoly(b_tilde,0);
  }
  GEN tmp = gtopoly(vec_ei(n,1),0);
  GEN t1 = FpX_div(FpX_mul(f_tilde,a,p),gtopoly(vec_ei(n+1,1),0),p);
  GEN t2 = FpX_div(FpX_mul(b_tilde,a,p),tmp,p);
  GEN t3 = RgX_mullow(t1,FpX_mul(b_tilde,h,p),n-1);
  GEN ret = FpX_sub(t2,FpX_mul(pol_x(0),t3,p),p);
  /*pari_sp lbot = avma;*/
  /*
  GEN vec = gtovec(ret);
  long d = lg(vec)-1;
  if(d==n){return gerepile(ltop,lbot,vecreverse(vec));}
  else{
    return gerepile(ltop,lbot,vecreverse(shallowconcat(vec,zerovec(n-d))));
  }
  */
  GEN vec = gtovecrev(ret);
  long d = lg(vec)-1;
  if(d==n){return  vec;}
  else{
    return shallowconcat(vec,zerovec(n-d));
  }
}

/* 2m is the size of multiples */
static ulong
rational_equ(GEN *Xp,GEN *multiples, ulong m,ulong ell, GEN fl,GEN *efl, GEN p)
{
  pari_sp ltop = avma;
  ulong i, j;
  ulong n = degpol(fl);
  GEN f_tilde = polrecip(fl);
  GEN tmp =  gtopoly(vec_ei(n,1),0);
  GEN h = FpXQ_invsafe(f_tilde,tmp,p);
  GEN w = cgetg(n+1,t_VEC);
  for(i=1;i<=n;i++){
    gel(w,i)= randomi(p);
  }
  GEN Lw = cgetg(2*m+1,t_VEC);
  for(i=1;i<=m;i++){
    gel(Lw,2*i-1)=mod_transp(gel(*multiples,2*i-1),w,h,fl,p);
    gel(Lw,2*i)=mod_transp(gel(*multiples,2*i),w,h,fl,p);
  }
  GEN cj,dj;
  for(j=1;j<=m;j++){
    cj = RgX_to_RgV(mod_compose(gel(*multiples,2*j-1),*Xp,efl,n,p),n);
    dj = RgX_to_RgV(mod_compose(gel(*multiples,2*j),*Xp,efl,n,p),n);
    for(i=1;i<=m;i++){
      if(gequal(FpV_FpC_mul(gel(Lw,2*i-1),dj,p),FpV_FpC_mul(gel(Lw,2*i),cj,p))){
	avma = ltop;
	return Fl_div(i,j,ell);
      }
    }
  }
  pari_err_BUG("bsgs (rational_equ)");
  return 0;/* NOT REACHED */
}



static ulong
bsgs_abs(GEN a4,GEN a6,ulong ell,GEN fl,GEN p)
{
  pari_sp av = avma;
  ulong eig;
  ulong m = usqrt(ell)+1;
  GEN eF = FpX_get_red(fl,p);
  GEN Xp = FpXQ_pow(pol_x(0),p,eF,p);
  if(gequal(Xp,pol_x(0))){
      eig = 1;
    } else {
    GEN multiples = compute_multiples(m,a4,a6,fl,p);
    eig = rational_equ(&Xp, &multiples, m,ell,fl,&eF,p);
    }
  GEN resu = FpX_resultant(fl,mkpoln(4,gen_1,gen_0,a4,a6),p);
  long kr = kronecker(resu,p);
  if(krouu(eig,ell)==kr){avma = av;return eig;}
  else{avma = av;return ell-eig;}
}

