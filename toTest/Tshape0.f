      subroutine tshape (a,b,m,w,x,t)
      implicit double precision (a-h,o-z)
c     
c     ----------------------------------------------------------------- 
c     This subroutine calculates the kinetic energy matrix, T, 
c     in a normalised Lobatto shape function basis.
c
c	  Formula for this are taken from: 
c    "QUANTUM SCATTERING VIA THE LOG DERIVATIVE VERSION
c    OF THE KOHN VARIATIONAL PRINCIPLE" 
c    D. E. Manolopoulos and R. E. Wyatt, 
c	  Chem. Phys. Lett., 1988, 152,23
c
c	  In that paper Lobatto shape (L.s.) functions are defined
c     -----------------------------------------------------------------  
c     
      dimension w(m),x(m),t(m,m)

      dimension ww(m+1),xx(m+1),tt(m+1,m+1)

c	  I think, that this is needed for the sum defined in L.S. functions to work
      n = m+1
c	  Get points and weights for n point Lobatto quadrature in (a,b)
      call lobatto (a,b,n,ww,xx)
      
c	  No idea why it's done
      do 1 i = 1,n
         ww(i) = sqrt(ww(i))
 1    continue
 
      do 4 i = 1,n
         ff = 0.0d0
         do 3 j = 1,n
c			tt(i,i) is trivially = 0, so no need for loops       
            if (j .eq. i) go to 3
c			gg will be '-' value of  derivative of i-th L.s. function at
c			j-th root, which is: i-th L.s. function evaluated at j-th
c			root divided by ( i-th root minus j-th root )
            gg = 1.0d0/(xx(i)-xx(j))
            ff = ff+gg
            
            do 2 k = 1,n
            
c			   This loop multiplies gg defined above by i-th L.s. 
c			   function evaluated at j-th root, which is itself a Lagrangian interpolation
               if (k.eq.j .or. k.eq.i) go to 2
               gg = gg*(xx(j)-xx(k))/(xx(i)-xx(k))
               
 2          continue
c			Write into tt value of derivative of j-th L.s. function
c			evaluated at i-th root. This relation is described in the paper mentioned
            tt(j,i) = ww(j)*gg/ww(i)   
            
 3       continue
c		 tt of i,i is 0 as the i-th L.s. has a maximum at the i-th root
         tt(i,i) = ff
         
 4    continue
      do 7 i = 1,m
c		 In this approach 1: roots are in decreasing order
c		 2: last root ( 0 ) doesn't get included in calculations
c		 but subroutine lobatto returns roots and weights in 
c		 increasing order so this has to be done manually
         w(i) = ww(i+1)
         x(i) = xx(i+1)

         do 6 j = 1,i
         
            hh = 0.0d0
            do 5 k = 1,n
c			   Entries in T matrix are defined as a sum over all k from 0 
c           to n+1 of: 
c           [ k-th weight ]*[ derivative 
c			   of i-th L.s. function at k-th root ] *[ derivative 
c			   of j-th L.s. function at k-th root ]

               hh = hh + tt(k,i+1)*tt(k,j+1)
 5          continue
c			t is symmetric
            t(i,j) = hh
            t(j,i) = hh
            
 6       continue
 
 7    continue
      return
      end