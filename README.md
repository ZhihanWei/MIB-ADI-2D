# MIB-ADI-2D

# README 

## Project Description 
A package for parabolic interface problem in two-dimensional space using both matched interface method (MIB) and alternating direction implicit method (ADI) with new techniques for tangential derivative approximations.

## Implementation

Step 1: In terminal(Mac)/Command Line(Windows), go to the folder where the makefile is located, and type "make", it will compile and link to an executable file.

Step 2: Once you see ">>> compiled on (hostname of your PC) with  <<<", this mean the executable file is generated successfully.

Step 3: Type "./mib2d_e\<example number\>s\<slover number\>" to run the executable file, you'll see results on the screen.

Step 4: If you want to compile with a set of new parameters, type "make clean" to remove the previous object files and executable file. Then go to Step 1.

## Reminds

• The package is written in Fortran 90.

• Examples are should be added as file start with "data_\*\*\*.f90" including both interface informations and analytical solution informations.

