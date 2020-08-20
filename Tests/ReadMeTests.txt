This folder holds the files necessary to run overall tests of the fortran
component of multiscat. 
Currently 2 tests are implemented, Test1 tests a single scattering condition 
whereas Test2 tests a list of scattering conditions.

These tests do not work on Windows. They are also set up to work with 
multiscat as it was in August 2020, using exact file names. Specifically:

   - The potential file must be exactly "pot10001.in", it is shared between
     all the tests.
   - The fourier labels file must be called "FourierLabels.in"
Any changes in filenames or functionality of Multiscat, even in the low 
significant figures, will return a failed test.

Note also that tests appear to pass sometimes if the code fails to compile,
this can be easily seen by reading the console output for compiler errors.

The tests can be set up as follows:

   - Move MultiscatCaT as well as the directories holding the tests into the 
     the same directory as the fortran multiscat files.

   - Check that the "Directories" list in MultiscatCaT contains the names
     of all the directories that contain the releveant tests.

   - Change permissions of the file so it can be excecuted. This is usually
     going to need the command " chmod u+x MultiscatCaT ", but it may be 
     different depending on the situation. Look up chmod if errors come up.

   - Run " sed -i -e 's/\r$//' MultiscatCaT ". This removes the awkward 
     line endings introduced by Windows. These may be present if the file 
     has been edited on Windows at any point.

  - Run " ./MultiscatCaT " in the directory containing the fortran files.
    This runs the tests

To create a new test:

   - Make a new directory. Suppose this is called "TestN".

   - In this directory place scattering conditions file called "scatCond.in"
     and a configuration file called "Multiscat.conf".

   - Place, also in this directory, a file containing the expected results.
     This should be called "Result.out" and be formatted exactly as the 
     "diffract10001.out" file.

   - Add "TestN" to the "Directories" variable in MultiscatCaT. Include the 
     TestN directory in the test like you would Test1 ot Test2.

How tests are performed:

   - The MultiscatCat script goes, in turn, into each of the directories listed 
     in the "Directories" variable.

   - Here, it clears the directory of old files and then copies in the files it 
     needs from the parent directory, fixing the Windows line endings as necessary.

   - Then it compiles the Fortran code and runs it. The result is compared to 
     what is expected, the Reult.out file, and a " Test Passed " message is printed
     if these files are perfectly identical. Note that multiscat sends a lot of text
     to the console, so the "Test Passed" message may be buried by other text.

   - The directory is then cleane up a little and the process reperats for the next
     test folder.