#Chombo_4/example/EBCM

# All EBCM examples should have the following.
##  A source file with code for the test.
##  The test has been run for a specfic simulation campaign.
##  A directory _test_harness
### This contains a python script used for the simulation campaign.  It is well-documented.
##  A directory _doc 
### This contains a latex document explaining the context and results from the simulation campaign.
### It contains a directory _data which  archives the the simulation campaign's data.


#Directories

# hoeb_geometry_test
##This test makes a single-box EBCM::MetaDataLevel.
## It then makes a ton of data on matrix conditioning in an hoeb/ebcm context.
## The associated document (Lawrence Berkeley National Laboratory Technical Report LBNL-2001528) technical report is quite complete.   

#hoeb_spmd

##This tests the communications infrastructure for Chombo_4/src/EBCM.
## These are the bits of Chombo that handle calculations within the context of  cut cell geometries where merger is an option.
##The associated document is complete but minimal.

#hoeb_petsc

##This test investigates the solvability of the system matrix when boundary condition equations are included in the system.
## The associated document (Lawrence Berkeley National Laboratory Technical Report LBNL-2001537)
is complete to the point of being downright
chatty.
## This test does not use petsc so the test name needs work.

#hoeb_truncation

##We intend to continue the story by investigating the truncation error of
the family of finite volume operators.
##This one will use quite a bit
of petsc and includes yet another user interface to the petsc
infrastructure.   The document and most of the code 
are vaporware as of July 7 2023.
