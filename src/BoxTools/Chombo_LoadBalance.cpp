#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <list>
#include <set>
using std::cout;

#include "Chombo_parstream.H"
#include "Chombo_DataIterator.H"
#include "Chombo_SPMD.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_LayoutIterator.H"
#include "Chombo_CH_Timer.H"

// Write a text file per call to LoadBalance()
//#define PRINT_EXTRA_LB_FILE
#ifndef CH_NTIMER
#include <fstream>
using std::fstream;
#endif

#include "Chombo_NamespaceHeader.H"

// local prototypes
int
min_element( const std::vector<long>& Vect );
void
min_max_elements( int& imin ,int& imax ,const std::vector<long>& Vect );

void
min_max_elements( int& imin ,int& imax ,const std::vector<long long>& Vect );

// Local class definition (needed to sort loads)
class Load
{
public:
  Load():load(0) ,grid_index(0)
  {
  }

  bool operator < (const Load& rhs) const
  {
    return load < rhs.load;
  }

  long load;      //actual load on this box
  int grid_index; //link to Grids[]
};

// Code:

///
// This version takes a std::vector<BoxLayout> and builds a matching std::vector of
// std::vector<Box> and uses it to call the full version.
///
int
LoadBalance( std::vector<BoxLayout>&    Grids            //in-out: input grids to balance
             ,Real&                 effRatio         //output: ratio of min load
             ,const std::vector<std::vector<long> >&  ComputeLoads     //input: computational cost
             ,const std::vector<int>&           RefRatios        //input: refinement ratio
             ,int nProc
             )
{
  CH_TIME("LoadBalance:std::vectorBoxLayout");
  std::vector<std::vector<Box> > boxes( Grids.size() );

  for ( int i = 0; i < Grids.size(); ++i )
    {
      for ( LayoutIterator j = Grids[i].layoutIterator() ; j.ok() ; ++j )
        {
          boxes[i].push_back( Grids[i][j()] );
        }
    }

  std::vector<std::vector<int> > procIDs( Grids.size() );
  int status = LoadBalance( procIDs ,effRatio
                            ,boxes ,ComputeLoads ,RefRatios, nProc );

  if ( status < 0 ) return status;  //LoadBalance() failed

  for ( int i = 0; i < Grids.size(); ++i )
    {
      LayoutIterator j = Grids[i].layoutIterator();
      int p;
      for ( j.reset() ,p = 0; j.ok(); ++j ,++p )
        {
          Grids[i].setProcID( j() ,procIDs[i][p] );
        }
    }

  return status;
}

///
// This version takes a single BoxLayout (i.e. one level)
// and uses the box volumes to construct the compute loads
// and calls the full version.
// a_LBnumProc specifies number of procs to use in LoadBalancing. defaults to numProc().
///
int LoadBalance(std::vector<int>& a_procAssignments, const std::vector<Box>& a_boxes,
                const int a_LBnumProc)
{
  CH_TIME("LoadBalance:std::vectorBoxEntry");
  std::vector<long long> computeLoads(a_boxes.size());

  for (int i = 0; i < a_boxes.size(); ++i)
    {
      computeLoads[i] = a_boxes[i].numPts();
    }

  int status;
  status = LoadBalance(a_procAssignments, computeLoads, a_boxes, a_LBnumProc);

  // save this older commented code
  //   std::vector<std::vector<Box> > layouts(1, boxes);
  //   std::vector<std::vector<long> > computeLoads(1, std::vector<long>(boxes.size()));
  //   std::vector<int> refRatios(1,1);
  //   std::vector<std::vector<int> > assignments(1,std::vector<int>(boxes.size(), -1));
  //   Real effRatio;

  //   for (int index = 0; index < boxes.size(); ++index)
  //     {
  //       computeLoads[0][index] = layouts[0][index].numPts();
  //     }

  //   int ret = LoadBalance(assignments, effRatio, layouts, computeLoads, refRatios);

  //   if (ret == 0)
  //     procs = assignments[0];

  return status;
}

///
// Accepts "long int" computeLoads and then just calls the "long long" version.
///
int LoadBalance(std::vector<int>&             a_procAssignments,
                const std::vector<long int>&  a_computeLoads,
                const std::vector<Box>&       a_boxes,
                const int                a_LBnumProc)
{
  std::vector<long long> cloads(a_computeLoads.size());
  for (int i=0; i<a_computeLoads.size(); i++)
    {
      cloads[i] = (long long)a_computeLoads[i];
    }
  int status = LoadBalance(a_procAssignments, cloads, a_boxes, a_LBnumProc);
  return status;
}


void split_load(std::vector<int>& procs, std::vector<long long>& loads,
                std::vector<long long>& sum, int p0, int p1, int start, int end)
{
  if(p0==(p1-1))
    {
      for(int i=start; i<end; i++){
        procs[i]=p0;
      }
      loads[p0] = sum[end-1]-sum[start];
      return;
    }
  if(p1-p0 == end-start)
    {
      int p=p0;
      for(int i=start; i<end; i++,p++)
        {
          procs[i]=p;
          if(i>0)loads[p]=sum[i]-sum[i-1];
          else loads[p]=sum[0];
        }
      return;
    }
  if(start == end-1) MayDay::Error(" we should not arrive at this terminal case in LoadBalance \n");
  long long base = 0;
  if(start!=0) base=sum[start-1]; 
  long long half = (sum[end-1]-base)/2+base;
  int scan=start;
  while(sum[scan]<half) scan++;
  scan++;
  int p = (p0+p1+1)/2;
  while(p-p0>scan-start)scan++;
  while(p1-p>end-scan) scan--;
  split_load(procs, loads, sum, p0,   p,  start,  scan);
  split_load(procs, loads, sum, p ,   p1, scan,    end);
}
  
int LoadBalance(std::vector<int>&             a_procAssignments,
                const std::vector<long long>& a_computeLoads,
                const std::vector<Box>&       a_boxes,
                const int                a_LBnumProc)
{
  CH_TIME("LoadBalance:std::vectorBoxSimple");
  // Phase 1  modified knapsack algorithm.  this one doesn't use
  // load sorting followed by round robin.  This one does bin packing
  // first by finding vector divisors, then does regular knapsack
  // optimization.

  int status = 0;
  a_procAssignments.resize(0);
  a_procAssignments.resize(a_computeLoads.size(), 0);

  const int Nboxes = a_procAssignments.size();

  if (a_LBnumProc == 1 || Nboxes == 1)
  {
    return 0;
  }

  //  first, add another simple bail out if there are fewer boxes than there are LBnumProc
  if(a_LBnumProc >= a_boxes.size())
    {
      for(int i=0; i<a_procAssignments.size(); i++) a_procAssignments[i]=i;
      return 0;
    }

 
  std::vector<long long> loads(a_LBnumProc, 0);
  /*
  //still serial version based on prefix sum   (bvs)
  std::vector<long long> sum(a_computeLoads.size());
  sum[0]=a_computeLoads[0];
  for(int i=1; i<sum.size(); i++) sum[i]=sum[i-1]+a_computeLoads[i];
  // recursive function calls for splitting
  split_load(a_procAssignments, loads, sum, 0, a_LBnumProc, 0, sum.size());

  double totalLoad = sum.back();

  */ 
  // first, compute 'total load'
  double totalLoad = 0;
  for (int ibox = 0; ibox < Nboxes; ++ibox)
    {
      totalLoad += a_computeLoads[ibox];
      a_procAssignments[ibox] = -ibox;
    }

  // figure out roughly single processor load goal
  double goal = totalLoad/a_LBnumProc;

  // debugging prints
  //pout()<<"--- LoadBalance  --- Number of Boxes: " << Nboxes << std::endl;
  //pout()<<"       a_computeLoads.size()= "<< a_computeLoads.size() << std::endl;
  //pout()<<"       a_boxes.size()=        "<< a_boxes.size() << std::endl;
  //pout() << " Box sizes: " << std::endl;
  //for (int i=0; i < Nboxes; i++) pout() <<" "<<a_computeLoads[i];
  //pout() << std::endl;

  // now, first bin assignments

  // change goal to dynamicgoal (ndk)
  // We desire equal amount of cells per processor, which
  // is the total number of cells divided by the number of processors: ie the goal.
  // But because the cells are grouped into boxes of sometimes very different sizes,
  // and because we assume that the boxes are in a locality-preferable vector,
  // using a constant goal can allow for several consecutive processors to be
  // assigned a number of cells _less_ than (or _more_ than) the goal.
  // This can lead to the very last processor receiving an abnormally large
  // (or small) workload.  Changing the algorithm to compute a new dynamic
  // goal after each processor-box assignment eliminates this problem.
  // However, there are still improvements that can be made when the
  // number of boxes to be farmed out is not much larger than the total number
  // of processors doing the work.  (ndk)

 
  int bin = 0;
  double remainingLoad = totalLoad;
  double dynamicGoal = goal;
  long long localLoad = 0;
  //pout() << " Begin loop over boxes: remainingLoad =" << remainingLoad
  //     << " dynamicGoal = " << dynamicGoal << " a_LBnumProc=" << a_LBnumProc
  //     << " Nboxes=" << Nboxes << std::endl;
  bool NboxesLeftEqualsNprocsLeft=false;
  //char ctmp[205];
  for (int ibox=0; ibox < Nboxes; ibox++)
    {
      // If the number of boxes left to sort out is <= to the number of proc bins left,
      // then force each remaining proc to have one box -- by setting NboxesLeftEqualsNprocsLeft=true
      if ( Nboxes-ibox <= a_LBnumProc-bin)
        {
          // Once set true, remains true for the rest of the boxes.
          NboxesLeftEqualsNprocsLeft = true;
        }

      if (NboxesLeftEqualsNprocsLeft)
        {
          // a_procAssignments initialized to -ibox above
          // Already something in this bin, so go to next bin.
          //  Note: Only need to set loads[bin] here for any diagnostics after loop
          if (a_procAssignments[ibox] > 0) bin++;

          if (bin > a_LBnumProc)
          {
            // Not sure this could even happen...
            MayDay::Abort("Problem in LoadBalance");
          }
          a_procAssignments[ibox] = bin;
          loads[bin] += a_computeLoads[ibox];
          bin++;
        }
      else
        {

          if (bin != a_LBnumProc-1)  // NOT last rank
            {
              localLoad += a_computeLoads[ibox];
              if (localLoad >= dynamicGoal)
                {
                  if ( ( localLoad - dynamicGoal ) > a_computeLoads[ibox]/2
                      && loads[bin]!=0)//added this last condition so no-load procs don't get skipped. mfb
                    {
                      a_procAssignments[ibox] = bin+1;
                      loads[bin+1] += a_computeLoads[ibox];
                      localLoad = a_computeLoads[ibox];
                    }
                  else
                    {
                      a_procAssignments[ibox] = bin;
                      loads[bin] += a_computeLoads[ibox];
                      localLoad = 0;
                    }
                  remainingLoad -= (double)loads[bin];
                  dynamicGoal = remainingLoad/(double)(a_LBnumProc-(bin+1));
                  bin++;
                }
              else
                {
                  loads[bin] += a_computeLoads[ibox];
                  a_procAssignments[ibox] = bin;
                }
            }
          else
            { // last rank
              a_procAssignments[ibox] = bin;
              loads[bin] += a_computeLoads[ibox];
            }
        }

      //sprintf(ctmp, " box#%04d localLoad=%10lld computeLoads[ibox]=%10lld remainingLoad=%12.1f dynamicGoal=%12.1f a_procAssignments[i]=%4d\n",
      //  ibox, localLoad, a_computeLoads[ibox], remainingLoad, dynamicGoal, a_procAssignments[ibox]);
      //pout() << ctmp;
    }
  
  return status;
}

///
// This version takes compute loads directly and does the "simple"
// load balance algorithm (modified knapsack algorithm).
// Boxes are also included, but only used if SWAP_SAME_SIZES_BOXES
// is enabled
// I think we can remove this as it simply calls the above LoadBalance() (ndk 8.4.2008)
///
int UnLongLongLoadBalance(std::vector<int>&                      a_procAssignments,
                          const std::vector<unsigned long long>& a_computeLoads,
                          const std::vector<Box>&                a_boxes,
                          const int                         a_LBnumProc)
{
  CH_TIME("UnLongLongLoadBalance:std::vectorBoxSimple");

  std::vector<long long> cloads(a_computeLoads.size());
  for (int i=0; i<a_computeLoads.size(); i++)
    {
      cloads[i] = (long long)a_computeLoads[i];
    }
  int status = LoadBalance(a_procAssignments, cloads, a_boxes, a_LBnumProc);
  return status;
}

///
// This version does the real work.
///
int
LoadBalance(std::vector<int>&          a_procAssignments
            ,Real&                a_effRatio
            ,const std::vector<Box>&   a_grids
            ,const std::vector<long>&  a_computeLoads)
{
  CH_TIME("LoadBalance:std::vectorBoxWork");
  std::vector<std::vector<Box> >   grids(1, a_grids);
  std::vector<std::vector<long> >  computeLoads(1, a_computeLoads);
  std::vector<int>            refRatios(1,2);
  std::vector<std::vector<int> > procs;
  int nproc = numProc();
  int retval = LoadBalance(procs, a_effRatio, grids, computeLoads,
                           refRatios, nproc);

  a_procAssignments = procs[0];
  return retval;
}

int
LoadBalance(std::vector<std::vector<int> >& procAssignments  //output: processor number
            ,Real&                 effRatio         //output: ratio of min load
            ,const std::vector<std::vector<Box> >&  Grids            //input: meshes to balance
            ,const std::vector<std::vector<long> >& ComputeLoads     //input: computational cost
            ,const std::vector<int>&           RefRatios        //input: refinement ratio
            ,int nProc                            // number of procs to assugn to
            )
{
  CH_TIME("LoadBalance:std::vectorBoxRealWork");
  // local variables
  Real eff_ratio; // efficiency ratio on a level
  int status = 0; // return code

  // Validate inputs
  if ( Grids.size() != ComputeLoads.size() )
    { return -1011; }
  if ( Grids.size() != RefRatios.size() )
    { return -1013; }
  for ( int lvl=0; lvl<Grids.size(); ++lvl )
    {
      if ( Grids[lvl].size() != ComputeLoads[lvl].size() )
        { return -1012; }
    }

  // set the number of elements in the output vector to the number
  // of levels and the number of elements in each element to the
  // number of boxes on each level and set the value of each
  // element to zero
  procAssignments.resize( Grids.size() );
  for ( int lvl=0; lvl<Grids.size(); ++lvl )
    {
      procAssignments[lvl].resize( Grids[lvl].size(),0 );
    }

  // check for special case of all loads on 1 processor
  if ( nProc == 1 )
    {
      for ( int lvl=0; lvl<Grids.size(); ++lvl )
        {
          for ( int i=0; i<Grids[lvl].size(); ++i )
            {
              procAssignments[lvl][i] = 0;
            }
        }
      effRatio = 1.0;
      status = 0;
    }
  else
    {
      // general case: loads on more than one processor
      effRatio = 1.0;

      // balance each level separately
      for ( int lvl=0; lvl<Grids.size(); ++lvl )
        {
          // first, build the load structure and sort by compute_load
          std::vector<Load> loads( Grids[lvl].size() );
          for ( int i=0; i<Grids[lvl].size(); ++i )
            {
              loads[i].load = ComputeLoads[lvl][i];
              loads[i].grid_index = i;
            }
          std::sort( loads.begin() ,loads.end() );
          // do the initial assignments by sequentially
          // `handing out' the loads from largest to smallest
          std::vector<long> total_loads( nProc,0 ); //total load per processor
          std::vector<std::vector<Load> > proc_loads( nProc ); //loads per processor
          int iproc_minload = 0; //processor with lowest load
          // loads are sorted in increasing order, so work backwards through the vector
          for ( int i=loads.size()-1; i>=0; --i )
            {
              // put the next load on the processor with the lowest total load
              proc_loads[iproc_minload].push_back( loads[i] );
              total_loads[iproc_minload] += loads[i].load;

              // recompute which processor has the lowest load
              //[NOTE: this would be faster if the loads were sorted]
              iproc_minload = min_element( total_loads );
            }
          // compute average load per processor, truncated to int
          long avg_load = 0;
          for ( int i=0; i<total_loads.size(); ++i ) avg_load += total_loads[i];
          avg_load /= nProc;

          // optimize the assignments by swapping a load off the
          // processor with the max load onto another processor
          // such that the load balance is improved
          int iter_count = 0, swap_count = 0;
          int iproc_maxload;
          long max_change; //largest change in load balance
          int ibmax=0,jbmax=0,ipmax=0,jpmax=0;  //box and processor indices corresponding to max_change

          while ( 1 )
            {
              max_change = 0;

              // find the processor that has the largest deviation from perfect load balance
              min_max_elements( iproc_minload ,iproc_maxload ,total_loads );
              if ( iproc_minload == iproc_maxload )
                {
                  // load balance is perfect
                  // (this won't happen except in test cases)
                  break;
                }
              CH_assert( total_loads[iproc_minload] <= avg_load &&
                      avg_load <= total_loads[iproc_maxload] );
              if ( avg_load - total_loads[iproc_minload] > total_loads[iproc_maxload] - avg_load )
                ipmax = iproc_minload;
              else
                ipmax = iproc_maxload;

              //[NOTE: dont need this here, but it may be useful for debugging.]
              eff_ratio = (Real)total_loads[iproc_minload] / (Real)total_loads[iproc_maxload];

              // deviation from perfect load balance for this proc
              long devib = total_loads[ipmax] - avg_load;

              // search all the other processors for the swap that has the maximum
              // reduction in the total deviation from the perfect load balance
              for ( int j=0; j<proc_loads.size(); ++j )
                {
                  if ( j != ipmax )
                    {
                      long devjb = total_loads[j] - avg_load;

                      // loop over all boxes on both processors
                      for ( int ibox=0; ibox<proc_loads[ipmax].size(); ++ibox )
                        {
                          for ( int jbox=0; jbox<proc_loads[j].size(); ++jbox )
                            {
                              iter_count++;
                              // how much bigger is the ibox load than the jbox load?
                              long diff = proc_loads[ipmax][ibox].load
                                - proc_loads[  j  ][jbox].load;
                              // change in total deviation from swapping boxes
                              long change = std::abs( devib ) + std::abs( devjb )
                                - std::abs( devib - diff ) - std::abs( devjb + diff );
                              // remember this pair of boxes if the change is better
                              //[NOTE: max_change starts at 0, so this is always an improvement]
                              if ( change > max_change )
                                {
                                  max_change = change;
                                  ibmax = ibox; jbmax = jbox; jpmax = j;
                                }
                            }
                        }
                    }
                }
              // if there is a swap that improves load balance, take it; else stop
              if ( max_change > 0 )
                {
                  // adjust the total loads on each processor
                  long load_diff = proc_loads[ipmax][ibmax].load
                    - proc_loads[jpmax][jbmax].load;
                  CH_assert( load_diff != 0 );
                  total_loads[ipmax] -= load_diff;
                  total_loads[jpmax] += load_diff;
                  // swap the loads
                  Load tmp = proc_loads[ipmax][ibmax];
                  proc_loads[ipmax][ibmax] = proc_loads[jpmax][jbmax];
                  proc_loads[jpmax][jbmax] = tmp;

                  swap_count++;
                }
              else
                {
                  break;
                }
            }

          // Done with this level.

          // Compute the final efficiency ratio and save it if appropriate.
          min_max_elements( iproc_minload ,iproc_maxload ,total_loads );
          eff_ratio = (Real)total_loads[iproc_minload] / (Real)total_loads[iproc_maxload];
          if ( eff_ratio < effRatio ) effRatio = eff_ratio;

          // Assign boxes to processors for this level.
          for ( int ip=0; ip<proc_loads.size(); ++ip )
            {
              for ( int jb=0; jb<proc_loads[ip].size(); ++jb )
                {
                  procAssignments[lvl][proc_loads[ip][jb].grid_index] = ip;
                }
            }

#ifndef NDEBUG
          //           if ( iter_count > 0 )
          //             {
          //               cout << "    debug: LoadBalance: level " << lvl << " used "
          //                    << iter_count << " iterations and "
          //                    << swap_count << " swaps to get efficiency ratio "
          //                    << eff_ratio << std::endl;
          //             }
#endif
        }

      // Done with all levels.
      // We could try to permute processors assignments between levels to
      // reduce communication, but it is probably not worth the effort
      // since it probably would have O(N^4) cost (N==#boxes).
    }

  return status;
}

/// convenience function to gather a distributed set of Boxes with their corresponding processor assignment
/** Assumption is that each processor has at most one valid box. This is useful when interacting with other distributed codes which might not have the entire set of distributed boxes on all processors.
 */
int LoadBalance(std::vector<int>& a_procAssignments, 
                std::vector<Box>& a_boxes,
                const Box&   a_localGridBox,
                const int    a_numProc)
{
  int status = 0;

  const int numBox = numProc();

  std::vector<Box> tempBoxes;
  std::vector<int> tempProcAssign;
  tempBoxes.resize(numBox);

  if (numBox == 1) tempBoxes[0] = a_localGridBox;

#ifdef CH_MPI
  Box* nonConstBox = const_cast<Box*>(&a_localGridBox);  
  int boxSize = sizeof(Box);
  status = MPI_Allgather(nonConstBox,  boxSize,  MPI_BYTE, &(tempBoxes[0]), 
                         boxSize , MPI_BYTE , CH4_SPMD::Chombo_MPI::comm);
#endif

  tempProcAssign.resize(tempBoxes.size());
  for (int i=0; i< tempBoxes.size(); i++)
    {
      tempProcAssign[i] = i;
    }

  // filter out any empty boxes (it's OK if there are fewer 
  // boxes than processors, but DisjointBoxLayout will 
  // choke if given an empty box, so remove them now...
  
  // this isn't the most efficient way to do this -- assumption
  // is that we're not dealing with all that many boxes, and that 
  // this is only done once anyway...
  for (int i=0; i<tempBoxes.size(); i++)
    {
      if (!tempBoxes[i].isEmpty() )
        {
          a_boxes.push_back(tempBoxes[i]);
          a_procAssignments.push_back(tempProcAssign[i]);
        }
    }
  
  return status;
}

int basicLoadBalance(std::vector<int>& a_procAssignments, int a_numBoxes, int a_numProc)
{
  CH_TIME("basicLoadBalance");
  a_procAssignments.resize(a_numBoxes);
//#pragma omp parallel
  {
    if(a_numProc == 1)
      {
//#pragma omp for
        for(int i=0; i<a_numBoxes; i++) a_procAssignments[i]=0;
      }
    else if(a_numBoxes<=a_numProc)
      {
//#pragma omp for
        for(int i=0; i<a_numBoxes; i++) a_procAssignments[i]=i;
      }
    else
      {
        float factor = (float)a_numProc/a_numBoxes;
//#pragma omp for
        for(int i=0; i<a_numBoxes; i++)
          {
            a_procAssignments[i] = i*factor;
          }
      }
  }
  return 0;
}
      
////////////////////////////////////////////////////////////////
//                utility functions                           //
////////////////////////////////////////////////////////////////

//
// Find the index of the small value in a non-empty (long) std::vector
//
int
min_element( const std::vector<long>& Vect )
{
  CH_assert( Vect.size() > 0 );
  int imin = 0;
  for ( int i=1; i<Vect.size(); ++i )
    {
      if ( Vect[i] < Vect[imin] ) imin = i;
    }
  return imin;
}

//
// Find the indices of the smallest and largest values in a non-empty (long) std::vector
//
void
min_max_elements( int& imin ,int& imax ,const std::vector<long>& Vect )
{
  CH_assert( Vect.size() > 0 );
  imin = 0; imax = 0;
  for ( int i=1; i<Vect.size(); ++i )
    {
      if ( Vect[i] < Vect[imin] ) imin = i;
      if ( Vect[i] > Vect[imax] ) imax = i;
    }
  return;
}
void
min_max_elements( int& imin ,int& imax ,const std::vector<long long>& Vect )
{
  CH_assert( Vect.size() > 0 );
  imin = 0; imax = 0;
  for ( int i=1; i<Vect.size(); ++i )
    {
      if ( Vect[i] < Vect[imin] ) imin = i;
      if ( Vect[i] > Vect[imax] ) imax = i;
    }
  return;
}

#include "Chombo_NamespaceFooter.H"
