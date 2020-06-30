! A script for examples of code that compiles successfully.
! These are broadly speaking going to be answers to the 
! questions in the ?HowToFortran.f file, but feel free to 
! add any snippet of code that you think someone else might 
! be interested in. Please add explanations on why your code 
! works and output the result at the end.
! This file is for .f90 code.

  program answers
  
  REAL   TWO,TIME
  real x, y
  integer n

  ! How to output with print, format is a mess

  TIME = 4.3
  TWO = 2
  PRINT *, 'The time is ', TIME, ', and 2= ', TWO


  ! Basic output with write (mimics print)
  x = 0.025
  y = 1.123
  write(*,*) 'x=', x, ', and y=', y

  !(1) -Francesco 
  !Found online, works with FORTRAN 77 apparently, so
  !maybe it's a difference in compilers? Something 
  !similar is used extensively in our code
  !https://www.tat.physik.uni-tuebingen.de/~kley/lehre/ftn77/tutorial/format.html

  n = 20
  write(*,100) n ! n will be printed, 100 points to next line
100 format('N=', i2.0) 
  ! must be 2 spaces back, 'N=' is just printed but the i2.0 breaks down as: 
  ! i (integer)  2 (width in characters)  .0 (of which 0 are decimals)
  !
  ! See http://www.personal.psu.edu/jhm/f90/lectures/23.html
  ! for format specification
  ! FYI: format only works in .f90 files (at least with this syntax)

  end
