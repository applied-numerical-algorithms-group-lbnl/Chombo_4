

## Chombo_4/example/EBCM
*  Chombo examples have been the primary way people get started with Chombo.
*  Going forward, I will use this process to make the examples more focused and well-documented.

## EBCM example process
*  Each example will run a particular test in the source directory.
*  The test will be run for a specific simulation campaign.
*  The test will have a documents directory that explains the test and the campaign's results.
*  If there is enough data to merit a report, an LBNL technical report will result from the campaign.
*  If necessary, the test harness will be included and well-documented.

## hoeb_geometry_test
* This test makes a single-box EBCM::MetaDataLevel.
* It then makes a ton of data on matrix conditioning in an hoeb/ebcm context.
* The associated document (Lawrence Berkeley National Laboratory Technical Report LBNL-2001528) technical report is quite complete.   

## hoeb_spmd
* This tests the communications infrastructure for Chombo_4/src/EBCM.
* These are the bits of Chombo that handle calculations within the context of  cut cell geometries where merger is an option.
* The associated document is complete but minimal.

## hoeb_petsc
* This test investigates the solvability of the system matrix when boundary condition equations are included in the system.
* The associated document (Lawrence Berkeley National Laboratory Technical Report LBNL-2001537) is complete to the point of being downright chatty.
* This test does not use petsc so the test name is admittedly counter-intuitive.

## hoeb_truncation
* We intend to continue the story by investigating the truncation error of a family of finite volume elliptic operators.
* This one will use quite a bit of petsc and includes yet another user interface to the petsc infrastructure.
* The document and most of the code are vaporware as of July 7 2023.
