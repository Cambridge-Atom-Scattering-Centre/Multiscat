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
  character*40 var

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
  ! 'N=' is just printed but the i2.0 breaks down as: 
  ! i (integer)  2 (width in characters)  .0 (of which 0 are decimals)
  !
  ! See http://www.personal.psu.edu/jhm/f90/lectures/23.html
  ! for format specification
  ! FYI: format only works in .f90 files (at least with this syntax)
  
  !(2) -Francesco
  !Assign a string variable and format at the same time. This may be the 
  !only way to assign string varibles from o=numerical types (?)
  write(var,101)  n
101 format( 'var has been assigned this string, n= ', i2.0)
  print *, var
  ! The write command with the first argument set to a variable writes the result
  ! to the variable instead of printing to the console. The code is otherwise similar
  ! to (1). The print statement prints to console so the result of the assignment 
  ! can be seen
  end