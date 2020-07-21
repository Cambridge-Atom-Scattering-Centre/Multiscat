!Hard limits on Multiscat parameters included everywhere at compile time

      parameter (NZFIXED_MAX = 550)        !max no of z points in fixed pot.
      parameter (NVFCFIXED_MAX = 1024)      !max no of fourier cmpts (fixed, from file)

      parameter (nmax = 1024)    !         diffraction channels
      parameter (mmax = 550)    !         z grid points
      parameter (nfcx = 1024)    !         max (?) potential fourier components

      parameter (npmax = 500)    !         NAG routine E04JAF: no of parameters
      parameter (lw=npmax*(npmax-1)/2+12)
      parameter (liw=npmax+2)   !
      parameter (liw2=npmax+1)  !         NAG routine E05CCF parameter
      parameter (ncmx=10)       !         no of experimental conditions
      parameter (npwx=100)      !         no of potential well curves
