c	B.E. Pierzchala 08.2020

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
      
      do 4 i = 1,n
         
         do 3 j = 1,i
c			tt(i,i) is different       
            if (j .eq. i) go to 3
            
c			gg will be value of derivative of i-th Lsf at j-th root
            gg = 1.0d0/(xx(i)-xx(j))
            
            do 2 k = 1,n
            
               if (k.eq.j .or. k.eq.i) go to 2
               gg = gg*(xx(j)-xx(k))/(xx(i)-xx(k))
               
 2          continue
c			Write into tt[i,j] value of derivative of j-th L.s. function
c			evaluated at i-th root. 
            tt(i,j) = gg
c			Use relation mentioned in the paper to minimize work
            tt(j,i) = -gg*ww(j)/ww(i)
            
 3       continue
c		 tt of i,i is 0 except for first and last element
         tt(i,i) = 0.0d0
         if ( i .eq. 1) tt(i,i) = -1/(2*ww(i))
         if ( i .eq. n) tt(i,i) =  1/(2*ww(i))
         
 4    continue
      do 7 i = 1,m

c		 Why these two lines are performed we don't know
         w(i) = sqrt( ww(i+1) )
         x(i) = xx(i+1)

         do 6 j = 1,i
         
            hh = 0.0d0
            do 5 k = 1,n
c			Entries in T matrix are defined as a sum over all k from 0 
c           to n+1 of: 
c           [ k-th weight ]*[ derivative of i-th L.s. function at k-th root ] *
c			[ derivative of j-th L.s. function at k-th root ]

c			   hello is a new variable, that is incompatibile with the papers 
c			   original form, if you don't divide by henlo the result is as in theory
c			   as in hello = 1

			   hello = sqrt(ww(i+1)*ww(j+1))
               hh = hh + ww(k)*tt(i+1,k)*tt(j+1,k)/hello
 5          continue
c			t is symmetric
            t(i,j) = hh
            t(j,i) = hh
            
 6       continue
 
 7    continue
      return
      end