! This script is to ask questions and upload code that we
! expect should work, but doesn't. Please add a specific 
! question as a comment to all code here

        program questions

        !(1) -Francesco
        !Found online, works with FORTRAN 77 apparently, so
        !maybe it's a difference in compilers? Something 
        !similar is used extensively in our code
        !https://www.tat.physik.uni-tuebingen.de/~kley/lehre/ftn77/tutorial/format.html

        real x
        x = 0.025
        write(*,100) 'x=', x
        100 format(A,F)

        
        
        !(2) -Francesco
        !What's up with indenting? To avoid errors here I have to ad 2 TABs 
        !before every line, but multiscat.f90 seems to get away with 2 
        !spaces?
        
        !In multiscat.f90
  character*40 inputfile,outfile,fourierfile
  
        !In my code (.f file)
       character*40 inputfile,outfile,fourierfile
        
        
        
        
        
        end 
