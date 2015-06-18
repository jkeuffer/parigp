# parigp
Scripts in PARI/GP : 
mms.gp : implementation of "Computing the eigenvalue in the Schoof-Elkies-Atkin algorithm using Abelian lifts" by Mihailescu, Morain and Schost.
        usage : find_eigenvalue(E,ell,fl,flag) where:
               E is an elliptic curve
               ell is an Elkies prime
               fl is a factor of degree O(ell) of the ell-division polynomial of E
               flag = 0 or 1 : two different methods used in part of the algorithm : flag = 1 is faster
