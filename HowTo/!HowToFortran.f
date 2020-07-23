! A script for examples of code that compiles successfully.
! These are broadly speaking going to be answers to the 
! questions in the ?HowToFortran.f file, but feel free to 
! add any snippet of code that you think someone else might 
! be interested in. Please add explanations on why your code 
! works and output the result at the end.
! This file is for .f code.

        program answers

        ! How to output with print, format is a mess
        REAL   TWO,TIME
        real x, y

        TIME = 4.3
        TWO = 2
        PRINT *, 'The time is ', TIME, ', and 2= ', TWO
       
        
        ! Basic output with write (mimics print)
        x = 0.025
        y = 1.123
        write(*,*) 'x=', x, ', and y=', y
        
        end
