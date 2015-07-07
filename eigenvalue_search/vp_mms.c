#include "vp_mms.h"

/***  Computes a generator of (Z/ellZ)*  ***/
ulong 
compute_gen(ulong ell)
{
  ulong c = 2;
  while(Fl_order(c,0,ell)!=eulerphiu(ell)){c=c+1;}
  return c;
}

/***********          *************/
ulong
choose_q(ulong ell) 
{
  pari_sp av = avma;
  ulong res = 1,i=1,length;
  GEN fact = factoru(ell-1);
  GEN fact1 = gel(fact,1),fact2 = gel(fact,2);
  length=lg(fact1)-1;  
  switch(length)
    {
    case 1:
      avma = av;
      return ell-1;
    case 2: 
      res = upowuu(fact1[1],fact2[1]);
      avma = av;
      return res;
    case 3: 
      res = upowuu(fact1[1],fact2[1])*upowuu(fact1[2],fact2[2]);
      avma = av;
      return res;
    default:
      while(i<length){
	res *= upowuu(fact1[i],fact2[i]);
	i=i+2;
      }
      avma = av;
      return res;
    }
}

/************* Modular composition :returns g(h) mod f  **************/
/*Assume that g,h,f are FpX's , deg(f)>deg(g) and deg(f)>deg(h)*/

GEN
mms_mod_compose(GEN g, GEN h, GEN *eF, GEN p)
{
  pari_sp ltop = avma;
  GEN tmp1,res;
  ulong i,j;  
  ulong n =(ulong)get_FpX_degree(*eF),m = usqrt(n)+1;
  tmp1 = h;
  GEN A = zeromatcopy(m,n),B = zeromatcopy(m,m);
  gcoeff(A,1,1)=gen_1;
  for(i=1;i<m;i++){
    for(j=0;j<=degree(tmp1);j++){gcoeff(A,i+1,j+1)=polcoeff0(tmp1,j,0);};
    tmp1 = FpXQ_mul(tmp1,h,*eF,p);
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
  return gerepileupto(ltop,FpX_rem(res,*eF,p));
}

GEN
compo_t(GEN g, ulong n, GEN *eF,GEN p)
{
  GEN tmp1 = pol_x(0);
  GEN tmp2 = g;
  while(n>0){
    if(n%2==1){ tmp1 = mms_mod_compose(tmp1,tmp2,eF,p);}
    tmp2 = mms_mod_compose(tmp2,tmp2,eF,p);
    n = n>>1;
  }
  return tmp1;
}

GEN
vec_compo(GEN vec,ulong n, GEN *eF,GEN p)
{
  GEN vec1 = cgetg(3,t_VEC);
  gel(vec1,1) = FpX_red(gel(vec,1),p);
  gel(vec1,2) = FpX_red(gel(vec,2),p);
  GEN v = cgetg(3,t_VEC);
  gel(v,1) = pol_x(0);
  gel(v,2) = pol_1(0);
  while(n>0){
    if(n%2==1){
      gel(v,1) = mms_mod_compose(gel(v,1),gel(vec1,1),eF,p);
      gel(v,2) = FpXQ_mul(gel(vec1,2),mms_mod_compose(gel(v,2),gel(vec1,1),eF,p),*eF,p);
    }
    gel(vec1,2) = FpXQ_mul(gel(vec1,2),mms_mod_compose(gel(vec1,2),gel(vec1,1),eF,p),*eF,p);
    gel(vec1,1) = mms_mod_compose(gel(vec1,1),gel(vec1,1),eF,p);
    n = n>>1;
  }
  return v;
}

GEN 
gx(GEN a4,GEN a6,ulong m,GEN *eF,GEN p)
{
  pari_sp lbot = avma;
  GEN equ = FpX_rem(FpX_mulu(mkpoln(4,gen_1,gen_0,a4,a6),4,p),*eF,p);
  if(m==0) {avma = lbot; return pol_0(0);}
  if(m==1) {avma = lbot; return pol_x(0);}
  if(m%2==0) {
    GEN denom = FpXQ_invsafe(FpXQ_mul(equ,FpX_sqr(Fp_elldivpol(a4,a6,m,p),p),*eF,p),*eF,p);
    if(denom==NULL) pari_err_BUG("gx (invsafe failed)");
    GEN t1 = FpX_mul(Fp_elldivpol(a4,a6,m-1,p),Fp_elldivpol(a4,a6,m+1,p),p);
    pari_sp ltop = avma;
    return gerepile(lbot,ltop, FpX_rem(FpX_sub(pol_x(0),FpX_mul(t1,denom,p),p),*eF,p));		   
  }
  else {
    GEN denom = FpXQ_invsafe(FpX_sqr(Fp_elldivpol(a4,a6,m,p),p),*eF,p);
    GEN t1 = FpX_mul(Fp_elldivpol(a4,a6,m+1,p),Fp_elldivpol(a4,a6,m-1,p),p);
    GEN t2 = FpXQ_mul(t1,equ,*eF,p);
    pari_sp ltop = avma;
    return gerepile(lbot,ltop,FpXQ_sub(pol_x(0),FpXQ_mul(t2,denom,*eF,p),*eF,p));
  }
}

GEN
gy(GEN a4,GEN a6,ulong m,GEN *eF,GEN p)
{
  /*pari_sp av = avma;*/
  GEN equ = FpXQ_sqr(FpX_mulu(mkpoln(4,gen_1,gen_0,a4,a6),4,p),*eF,p);
  if(m==0) {/*avma = av;*/ return pol_0(0);}
  if(m==1) {/*avma = av;*/ return pol_1(0);}
  if(m==2) {
    GEN cst1 = subii(negi(powiu(a4,3)),muliu(sqri(a6),8));
    GEN cst2 = negi(muliu(mulii(a4,a6),4));
    GEN cst3 = negi(muliu(sqri(a4),5));
    GEN cst4 = muliu(a6,20);
    GEN cst5 = muliu(a4,5);
    GEN pol = FpX_red(mkpoln(7,gen_1,gen_0,cst5,cst4,cst3,cst2,cst1),p);
    GEN t1 = FpX_mulu(FpXQ_sqr(mkpoln(4,gen_1,gen_0,a4,a6),*eF,p),8,p);
    GEN denom = FpXQ_invsafe(t1,*eF,p);
    /*GEN ret = */return FpXQ_mul(pol,denom,*eF,p);
    /*return gerepilecopy(av,ret);*/
}
  if(m%2==0) {/*TO DO  traiter le cas invsafe==NULL*/
    GEN denom = FpXQ_invsafe(FpXQ_mul(equ,FpXQ_powu(Fp_elldivpol(a4,a6,m,p),3,*eF,p),*eF,p),*eF,p);
    GEN t1 = FpXQ_mul(Fp_elldivpol(a4,a6,m+2,p),FpX_sqr(Fp_elldivpol(a4,a6,m-1,p),p),*eF,p);
    GEN t2 = FpXQ_mul(Fp_elldivpol(a4,a6,m-2,p),FpX_sqr(Fp_elldivpol(a4,a6,m+1,p),p),*eF,p);
    /*GEN ret = */ return FpXQ_mul(FpXQ_sub(t1,t2,*eF,p),denom,*eF,p);
    /*return gerepileupto(av,ret);*/
  }
  else {
    GEN denom = FpXQ_invsafe(FpXQ_powu(Fp_elldivpol(a4,a6,m,p),3,*eF,p),*eF,p);
    GEN t1 = FpXQ_mul(Fp_elldivpol(a4,a6,m+2,p),FpX_sqr(Fp_elldivpol(a4,a6,m-1,p),p),*eF,p);
    GEN t2 = FpXQ_mul(Fp_elldivpol(a4,a6,m-2,p),FpX_sqr(Fp_elldivpol(a4,a6,m+1,p),p),*eF,p);
    /*GEN ret = */return FpXQ_mul(FpX_sub(t1,t2,p),denom,*eF,p);
    /*return gerepileupto(av,ret);*/
  }
}

GEN
T0X(GEN a4,GEN a6,ulong ell,ulong c,ulong q,ulong b,GEN *eF,GEN p)
{
  GEN U,V,W,tmp1,v;
  /*pari_sp av = avma;*/
  if(q%2){/*ok*/
    U = compo_t(gx(a4,a6,c,eF,p),q,eF,p);
    V = pol_x(0);
    W = pol_0(0);
    while(b>0){
      if(b%2==0)
	{
	  V = FpX_add(mms_mod_compose(V,U,eF,p),V,p);
	  U = mms_mod_compose(U,U,eF,p);
	} else {
	W = FpX_add(mms_mod_compose(W,U,eF,p),V,p);
	V = FpX_add(mms_mod_compose(V,U,eF,p),V,p);
	U = mms_mod_compose(U,U,eF,p);
      }
      b = b>>1;
    }
  } else {
    v = cgetg(3,t_VEC);gel(v,1)=gx(a4,a6,c,eF,p);gel(v,2)=gy(a4,a6,c,eF,p);  
    U = vec_compo(v,q,eF,p);
    V = cgetg(3,t_VEC);gel(V,1)=pol_x(0);gel(V,2)=pol_1(0);
    W = cgetg(3,t_VEC);gel(W,1)=pol_0(0);gel(W,2)=pol_0(0);
    while(b>0){
      if(b%2==0)
	{
	  gel(V,1) = FpX_add(gel(V,1),mms_mod_compose(gel(V,1),gel(U,1),eF,p),p);
	  tmp1 = FpXQ_mul(gel(U,2),mms_mod_compose(gel(V,2),gel(U,1),eF,p),*eF,p);
	  gel(V,2) = FpX_add(tmp1,gel(V,2),p);
	  gel(U,2) = FpXQ_mul(gel(U,2),mms_mod_compose(gel(U,2),gel(U,1),eF,p),*eF,p);
	  gel(U,1) = mms_mod_compose(gel(U,1),gel(U,1),eF,p);
	} else {
	gel(W,1) = FpX_add(gel(V,1),mms_mod_compose(gel(W,1),gel(U,1),eF,p),p);
	tmp1 = FpXQ_mul(gel(U,2),mms_mod_compose(gel(W,2),gel(U,1),eF,p),*eF,p);
	gel(W,2) = FpX_add(gel(V,2),tmp1,p);
	gel(V,1) = FpX_add(gel(V,1),mms_mod_compose(gel(V,1),gel(U,1),eF,p),p);
	gel(V,2) = FpX_add(gel(V,2),FpXQ_mul(gel(U,2),mms_mod_compose(gel(V,2),gel(U,1),eF,p),*eF,p),p);
	gel(U,2) = FpXQ_mul(gel(U,2),mms_mod_compose(gel(U,2),gel(U,1),eF,p),*eF,p);
	gel(U,1) = mms_mod_compose(gel(U,1),gel(U,1),eF,p);
      }
      b = b>>1;
    }
  }
  return W;
  /*return gerepileupto(av,W);*/
}

GEN
T1X(GEN a4,GEN a6,ulong ell,GEN e,ulong c,ulong q,GEN *eF,GEN p)
{
  /*pari_sp ltop = avma;*/
  if(q%2)
    {
      ulong qprime = (ell-1)/(2*q);
      return mms_mod_compose(e,compo_t(gx(a4,a6,c,eF,p),qprime,eF,p),eF,p);
    } else {
    ulong qprime = (ell-1)/(q);
    GEN v_tmp = cgetg(3,t_VEC);
    gel(v_tmp,1) = gx(a4,a6,c,eF,p);
    gel(v_tmp,2) = gy(a4,a6,c,eF,p);
    GEN pxy = vec_compo(v_tmp,qprime,eF,p);
    GEN ret = cgetg(3,t_VEC);
    gel(ret,1) = mms_mod_compose(gel(e,1),gel(pxy,1),eF,p);
    gel(ret,2) = FpXQ_mul(gel(pxy,2),mms_mod_compose(gel(e,2),gel(pxy,1),eF,p),*eF,p);
    return ret;
    /*return gerepileupto(ltop,ret);*/
  }
}

/******* Computing C such that C(alpha) = beta    **********/
GEN
trace_calc(GEN f,GEN p)
{
  pari_sp av = avma;
  long n  = poldegree(f,0);
  GEN lf = FpX_mul(pol_x(0),polrecip(FpX_deriv(f,p)),p);
  GEN frac = mkrfrac(lf,polrecip(f));
  GEN v = RgX_to_RgV(rfrac_to_ser(frac,n+2),n);
  return gerepileupto(av,RgV_to_FpV(v,p));
}

GEN
mms_mod_transp(GEN tau,GEN v,GEN h,GEN f,GEN p)
{
  pari_sp ltop = avma;
  GEN f_tilde = polrecip(f);
  long n = poldegree(f,0),degtau = degpol(tau),i;
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
  pari_sp lbot = avma;
 
  GEN vec = gtovecrev(ret);
  long d = lg(vec)-1;
  if(d==n){return gerepile(ltop,lbot,vec);}
  else{
    return gerepile(ltop,lbot,shallowconcat(vec,zerovec(n-d)));
  }
}

GEN
compute_trace(GEN alpha,GEN beta,GEN f,GEN p)
{
  if(signe(f)==0){err_printf("f==0!");return gen_0;}
  ulong n =(ulong) poldegree(f,0);
  GEN f_tilde = polrecip(f);
  GEN eF = FpX_get_red(f,p);
  GEN tmp =  FpXQ_powu(pol_x(0),n-1,eF,p);
  GEN h = FpXQ_invsafe(f_tilde,tmp,p);
  GEN v = mms_mod_transp(beta,trace_calc(f,p),h,f,p);/*v represents u->Tr(beta*u)*/
  ulong k1 = usqrt(n+1) , k2 = (n+1)/k1+1 ,i,j,k;
  GEN v1 = zerovec(k1+1);
  gel(v1,1) = pol_1(0);
  GEN v_tmp = zerovec(n);
  GEN v2 = zerovec(k1*k2);
  for(i=2;i<k1+2;i++)
    gel(v1,i)=FpXQ_mul(gel(v1,i-1),alpha,eF,p);
  for(i=0;i<k2;i++){
    for(j=0;j<k1;j++){
      for(k=0;k<=poldegree(gel(v1,j+1),0);k++)
	gel(v_tmp,k+1)=polcoeff0(gel(v1,j+1),k,0);
      gel(v2,i*k1+j+1)=FpV_dotproduct(v,v_tmp,p);
    }
    v_tmp = zerovec(n);
    v = mms_mod_transp(gel(v1,k1+1),v,h,f,p);
  }
  return v2;
}

GEN
compute_C(GEN alpha,GEN beta,GEN M,GEN f,GEN *eF,GEN p)
{
  if(gequalX(alpha)) 
    return beta;
  if (signe(M)==0){err_printf("M==0!");return gen_0;}
  pari_sp av = avma;  
  long q = poldegree(M,0);
  GEN v = vecslice(compute_trace(alpha,beta,f,p),1,q);
  GEN R1 = polrecip(RgX_mullow(polrecip(M),gtopolyrev(v,0),q));
  GEN R2 = FpXQ_mul(R1,FpXQ_invsafe(FpX_deriv(M,p),M,p),M,p);
  GEN r = Fp_div(pollead(beta,0),pollead(mms_mod_compose(R2,alpha,eF,p),0),p);
  return gerepileupto(av,FpX_Fp_mul(R2,r,p));
}

long 
bsgs_t(GEN C,GEN *eM,GEN Tp,ulong q,GEN p)
{
  pari_sp av = avma;
  ulong r = usqrt(q)+1,i,j;

  GEN Ci = cgetg(r+1,t_VEC);
  gel(Ci,1) = pol_x(0);
  for(i=2;i<r+1;i++) gel(Ci,i) = mms_mod_compose(gel(Ci,i-1),C,eM,p);
  
  GEN C_tilde = compo_t(C,q-r,eM,p);
  GEN C_tj = pol_x(0);
  GEN T_pj = Tp;
  for(j=0;j<r+1;j++){ 
    for(i=0;i<r;i++){
      if(gequal(gel(Ci,i+1),T_pj)){
	avma = av; return i+j*r;
      }
    }
    C_tj = mms_mod_compose(C_tj,C_tilde,eM,p);
    T_pj = mms_mod_compose(Tp,C_tj,eM,p);
  }
  pari_printf("search failed!");
  avma = av;
  return -1;
}

long
log_eigenvalue(GEN a4,GEN a6,ulong c,ulong q,ulong ell,GEN f_l,GEN *eT,GEN p)
{
  pari_sp av = avma;
  if(q%2==1){
    ulong qprime = (ell>>1)/q;
    GEN e0 = T0X(a4,a6,ell,c,q,qprime,eT,p);
    GEN e1 = T1X(a4,a6,ell,e0,c,q,eT,p);
    GEN  M =  FpXQ_minpoly(e0,*eT,p);
    GEN eM = FpX_get_red(M,p);
    long dM = degree(M);
    if (dM<0){pari_err_BUG("log_eigenvalue (minpol)");return(-1);}
    GEN C = compute_C(e0,e1,M,f_l,eT,p);
    if(!RgX_equal(mms_mod_compose(C,e0,eT,p),e1)){
      pari_err_BUG("log eigenvalue (C(e0)!=e1)");avma = av;return(-2);
    }
    GEN Tp = FpXQ_pow(pol_x(0),p,eM,p);
    long v = bsgs_t(C,&eM,Tp,dM,p);
    /* pari_printf("q odd, v= %ld\n",v);*/
    if (v>=0) {avma = av;return (qprime*v % q);}
    else {avma = av; return(-3);}
  } else {
    ulong qprime = (ell-1)/q;
    GEN e = T0X(a4,a6,ell,c,q,qprime,eT,p);
    GEN e0 = gel(e,2);
    GEN e1 = gel(T1X(a4,a6,ell,e,c,q,eT,p),2);
    GEN e0b = FpXQ_mul(mkpoln(4,gen_1,gen_0,a4,a6),FpXQ_sqr(e0,*eT,p),*eT,p);
    GEN e1b =  FpXQ_div(e1,e0,*eT,p);
    GEN N = FpXQ_minpoly(e0b,*eT,p);
    GEN M = mms_mod_compose(N,FpX_sqr(pol_x(0),p),eT,p);
    if (degree(N)<0){pari_err_BUG("log_eigenvalue (minpol)"); avma = av;return(-1);}
    GEN eN =  FpX_get_red(N,p);
    GEN eM =  FpX_get_red(M,p);
    GEN D = compute_C(e0b,e1b,N,f_l,eT,p);
    GEN C = FpX_mul(pol_x(0),gsubst(D,0,FpX_mul(pol_x(0),pol_x(0),p)),p);
    if(gequal0(C)){pari_err_BUG("log eigenvalue (compute_C)"); avma = av; return(-2);}
    GEN A = FpXQ_pow(pol_x(0),shifti(p,-1),eN,p);
    GEN Tp = FpXQ_mul(pol_x(0),mms_mod_compose(A,FpX_sqr(pol_x(0),p),&eM,p),eM,p);
    long v= bsgs_t(C,&eM,Tp,degree(M),p);
    if (v>=0){avma = av; return (qprime*v % q);}
    else {avma = av;pari_err_BUG("log eigenvalue (bsgs_t)");return(-3);}
  }
}


ulong
find_eigenvalue(GEN a4,GEN a6,ulong ell,GEN f,GEN p)
{
  pari_sp av = avma;
  GEN eF = FpX_get_red(f,p);
  ulong c = compute_gen(ell);
  ulong q1 = choose_q(ell);
  ulong q2 = (ell-1)/q1;
  long r;
  if(q2==1){
    r = log_eigenvalue(a4,a6,c,q1,ell,f,&eF,p);
    if(r<0){pari_err_BUG("find_eigenvalue");return 0;}
    return Fl_powu(c,(ulong)r,ell);
  }
  if(q1==2){/* only computing log_eigenvalue(q2) and using Dewaghe's trick*/
    r = log_eigenvalue(a4,a6,c,q2,ell,f,&eF,p);
    if(r<0){pari_err_BUG("find_eigenvalue");return r;}
    ulong r1 = Fl_powu(c,itou(Z_chinese(gen_1,stoi(r),utoi(q1),utoi(q2))),ell);
    ulong r2 = Fl_powu(c,itou(Z_chinese(gen_0,stoi(r),utoi(q1),utoi(q2))),ell);
    GEN resu = FpX_resultant(f,mkpoln(4,gen_1,gen_0,a4,a6),p);
    long kr = kronecker(resu,p);
    if(krouu(r1,ell)==kr){avma = av;return r1;}else{avma = av;return r2;}
  }
  long r1 = log_eigenvalue(a4,a6,c,q1,ell,f,&eF,p);
  long r2 = log_eigenvalue(a4,a6,c,q2,ell,f,&eF,p);
  if(r1<0 || r2<0){
    err_printf("eigenvalue search failed, r1 = %lu, r2 = %lu\n",r1,r2);
    return 0;
  }
  ulong r12 = itou(Z_chinese(stoi(r1),stoi(r2),utoi(q1),utoi(q2)));
  avma = av;
  return Fl_powu(c,r12,ell);
}
