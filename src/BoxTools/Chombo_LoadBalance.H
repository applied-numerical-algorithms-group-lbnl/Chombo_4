#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Purpose:
//  LoadBalance() is a global function to compute an assignment
//  of boxes to processors for an AMR mesh hierarchy.  The assignment
//  is made so as to balance the computation and communication
//  workload on each processor (i.e. make it as even as possible).
//
//  ********************************************************
//  CAVEAT: this version ignores the communication workload.
//          It balances ONLY the compute workload.
//  ********************************************************
//
// Usage:
//  The meshes in the AMR hierarchy are represented using
//   a std::vector of std::vectors of Boxes.
//
//  The computational workload is a real number for each box in the
//   hierarchy, also represented as a std::vector of std::vectors.
//
//  The communication workload is a real number for each connection
//   between pairs of boxes, which specifies the cost of communication
//   between the two boxes if they are on different processors.  No
//   allowance is made for varying the communication cost depending on the
//   distance between processors (ie. all processors are equally far apart)
//   and it assumes that the effective cost of communicating between two
//   boxes on the same processor is zero.  The communication workload is
//   represented as a graph.
//
//  ********************************************************
//  CAVEAT: the communication workload argument is not
//           present in this version of the function.
//  ********************************************************
//
//  The resulting assignment is an integer for each box in the hierarchy,
//   which gives the processor number (starting at zero) on which each box
//   should reside.  Again, represented as a std::vector of std::vectors.
//
//  The other output argument is the efficiency ratio, which is a measure
//   of how close to perfect the final load balance is.  The \var{effRatio}
//   output variable is defined as the smallest load on any processor divided
//   by the largest load.  A perfect load balance has efficiency == 1.
//
//  ********************************************************
//  CAVEAT: in this version, each level in the AMR hierarchy
//          is balanced independently, so the efficiency ratio
//          that is returned is the worst over all levels.
//  ********************************************************
//
//  It is important to note that it is the sum of the computation cost
//   and communication cost that is balanced, so the values in the two
//   variables must be in the same units.  It doesn't matter what the
//   units are.
//
//  ********************************************************
//  CAVEAT: the communication workload argument is not
//           present in this version of the function.
//  ********************************************************
//
// Interface Specification:
//  int LoadBalance() arguments are:
//       std::vector<std::vector<int>>& procAssignments  //output: processor number
//                                             //   for each box
//       Real&                effRatio         //output: efficiency ratio
// const std::vector<std::vector<Box>>& Grids            //input: meshes to balance
// const std::vector<std::vector<long>>&ComputeLoads     //input: computational cost
//                                             //   of each box
//###not used###
//### const [...something...]&CommunicateLoads //input: communication cost
//                                             //   between each pair of
//                                             //   neighboring boxes
// const std::vector<int>&         RefRatios        //input: refinement ratio
//                                             //   for each level
//
// Returns: integer status_code
//    =0  means succesful completion
//  exceptions: (output variables are not defined)
//  -1011 input vectors (\var{Grids} and \var{ComputeLoads}) are not
//        the same size
//  -1012 one of the vector elements of the input vectors (\var{Grids}
//        and \var{ComputeLoads}) are not the same size
//  -1013 input vectors (\var{Grids} and \var{RefRatios}) are not
//        the same size
//
// Modification History
//  19Nov99 <dbs> initial design and coding
//

#ifndef _Chombo_LOADBALANCE_H_
#define _Chombo_LOADBALANCE_H_

#include "Chombo_SPACE.H"
#include "Chombo_REAL.H"
#include "Chombo_Box.H"
#include <vector>
#include "Chombo_BoxLayout.H"
#include "Chombo_SPMD.H"
#include "Chombo_NamespaceHeader.H"

///
/**
   procAssignments  output: processor number for each box
   effRatio         output: ratio of min load to max load
   Grids            input: meshes to balance
   ComputeLoads     input: computational cost of each box
   RefRatios        input: refinement ratio for each level
*/
int
LoadBalance( std::vector<std::vector<int> >&         a_procAssignments
             ,Real&                        a_effRatio
             ,const std::vector<std::vector<Box> >&  a_Grids
             ,const std::vector<std::vector<long> >& a_ComputeLoads
             ,const std::vector<int>&           a_RefRatios
             ,int a_nProc = numProc()
             ) ;

///
/**

Grids            in-out: input grids to balance and output proc. numbers
effRatio         output: ratio of min load to max load
ComputeLoads     input: computational cost of each box
RefRatios        input: refinement ratio for each level
*/
int
LoadBalance( std::vector<BoxLayout>&    Grids
             ,Real&                 effRatio
             ,const std::vector<std::vector<long> >& ComputeLoads
             ,const std::vector<int>&           RefRatios
             ,int nProc = numProc()
             ) ;

///
/** convenience function to load balance a std::vector of Boxes based
    on load=box.numPts() */
int LoadBalance(std::vector<int>& a_procAssignments, const std::vector<Box>& a_boxes,
                const int a_LBnumProc = numProc());

///
/** convenience function to load balance a std::vector of Loads based
    on load.  boxes also passed in for possible use with box-swapping code */
int LoadBalance(std::vector<int>&             a_procAssignments,
                const std::vector<long long>& a_computeLoads,
                const std::vector<Box>&       a_boxes,
                const int                a_LBnumProc = numProc());

///
/** Accepts "long int" computeLoads and then just calls the "long long" version. */
int LoadBalance(std::vector<int>&             a_procAssignments,
                const std::vector<long int>&  a_computeLoads,
                const std::vector<Box>&       a_boxes,
                const int                a_LBnumProc = numProc());

//
/*
   This is experimental and to keep me from doing long long <-> long conversions. (dtg)
   I think we can remove this as it simply calls the above LoadBalance() (ndk 8.4.2008)
 */
int UnLongLongLoadBalance(std::vector<int>&                      a_procAssignments,
                          const std::vector<unsigned long long>& a_computeLoads,
                          const std::vector<Box>&                a_boxes,
                          const int                         a_numProc = numProc());

///
/* Even simpler load balance scheme that just dices vector to give as close to the same number
of boxes to each rank.  
*/
int basicLoadBalance(std::vector<int>& a_procAssignments, int numBoxes, int a_numProc = numProc());

///
int
LoadBalance(std::vector<int>& procAssignments
            ,Real&                 effRatio
            ,const std::vector<Box>&  Grids
            ,const std::vector<long>& ComputeLoads);

/// convenience function to gather a distributed set of Boxes with their corresponding processor assignment
/** Assumption is that each processor has at most one valid box. This is useful when interacting with other distributed codes which might not have the entire set of distributed boxes on all processors.
 */
int LoadBalance(std::vector<int>& a_procAssignments, 
                std::vector<Box>& a_boxes,
                const Box&   a_localGridBox,
                const int    a_numProc = numProc());

#include "Chombo_NamespaceFooter.H"
#endif
