#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_BoxPattern_H__
#define _Chombo_BoxPattern_H__

#include "Chombo_DisjointBoxLayout.H"
#include "Chombo_Pool.H"
#include <vector>
#include "Chombo_ProblemDomain.H"
#include "Chombo_CH_Timer.H"
#include "Chombo_NeighborIterator.H"
#include "Chombo_DataIterator.H"

#include <unordered_map>
#include <cstdint>

#include "Chombo_SPMD.H"
#ifdef CH_MPI  
#include "mpi.h"
#endif
using namespace Chombo4;
namespace CH4_Data_Choreography
{


  ///
  /**
     boxinfo_t
     This is the meta data for one side of a message where the data for one grid 
     is copying over a subset of the data in another grid.
     dblbox --- box of valid data.
     region --- box for copying is a subset of grow(dblbox, nghost).   Can be entirely outside dblbox.
     procid --- process id associated with this side of the message
     datind --- data index associated with this side of the message
  */
  struct boxinfo_t
  {
    Chombo4::Box       m_region;
    Chombo4::Box       m_valid;
    int                m_procID;
    Chombo4::DataIndex m_datInd;
    void 
    define(const Chombo4::Box      & a_region,
           const int               & a_procID,
           const Chombo4::DataIndex& a_datind,
           const Chombo4::Box      & a_valid)
    {
      m_datInd = a_datind;
      m_region = a_region;
      m_procID = a_procID;
      m_valid  = a_valid;
    }
  };
  
  ///
  /**
     box_interaction_t
     This holds the meta data for both sides of a communication substep
  */
  class box_interaction_t
  {
  public:
    box_interaction_t() {;}
      
    //region, procID
    boxinfo_t m_src;
    boxinfo_t m_dst;

    //for sorting
    bool operator < (const box_interaction_t& a_input) const
    {
      return m_src.m_region < a_input.m_src.m_region;
    }
    
    box_interaction_t(const Chombo4::DataIndex& a_srcInd,
                      const Chombo4::DataIndex& a_dstInd,
                      const Chombo4::Box      & a_srcRegion,
                      const Chombo4::Box      & a_dstRegion,
                      const int               & a_srcProcID,
                      const int               & a_dstProcID,
                      const Chombo4::Box      & a_srcValid,
                      const Chombo4::Box      & a_dstValid)
    {
      m_src.define(a_srcRegion, a_srcProcID, a_srcInd, a_srcValid);
      m_dst.define(a_dstRegion, a_dstProcID, a_dstInd, a_dstValid);
    }
             
  };

  ///
  /**
     BoxPattern is a class which deals with all the Box intersections associated with copyTo or exchange.
     Periodic boundary conditions are not supported.
  */
  class BoxPattern
  {
  public:
    ///
    /**
       Define for exchange--fills ghost cells around boxes if neighboring boxes are adjacent
    */
    BoxPattern(const Chombo4::DisjointBoxLayout & a_grids,
               const IntVect                    & a_ghost,
               bool a_printStuff = false)
    {
      CH_TIME("BoxPattern::exchange constructor");

      m_recvPattern.resize(CH4_SPMD::numProc());
      m_sendPattern.resize(CH4_SPMD::numProc());
      Chombo4::DataIterator dit = a_grids.dataIterator();
      Chombo4::NeighborIterator nit(a_grids);
      for (dit.begin(); dit.ok(); ++dit)
      {
        const Chombo4::Box& grid = a_grids[dit];
        int myProcID = CH4_SPMD::procID();
        Chombo4::Box gridGhost(grid);
        gridGhost.grow(a_ghost);

        for (nit.begin(dit()); nit.ok(); ++nit)
        {
          Chombo4::Box neighbor = nit.box();
          int neiProcID = a_grids.procID(nit());
          Chombo4::Box neighborGhost= neighbor;
          neighborGhost.grow(a_ghost);

          Chombo4::Box gridGhostInter(neighbor & gridGhost);
          if (!gridGhostInter.isEmpty())
          {
            //his data copies to my ghost
            Chombo4::Box srcValid = neighbor;
            Chombo4::Box dstValid = grid;
            box_interaction_t item(Chombo4::DataIndex(nit()), dit(), gridGhostInter, gridGhostInter, neiProcID, myProcID, srcValid, dstValid);
            if (myProcID == neiProcID)
            { // local move
              m_localPattern.push_back(item);
            }
            else
            {
              m_recvPattern[neiProcID].push_back(item);
            }
          }
          Chombo4::Box neighborGhostInter(grid & neighborGhost);
          if (!neighborGhostInter.isEmpty())
          {
            //my data copies to his ghost
            //local case will be taken care of when the data iterator gets to the neighbor box
            if(myProcID != neiProcID)
            {
              Chombo4::Box dstValid = neighbor;
              Chombo4::Box srcValid = grid;
              box_interaction_t item(dit(), Chombo4::DataIndex(nit()), neighborGhostInter, neighborGhostInter, myProcID, neiProcID, srcValid, dstValid);
              m_sendPattern[neiProcID].push_back(item);
            }
          }
        }
      }
      sortAllPatterns();

      
      if(a_printStuff)
      {
        printStuff(string("exchange")); 
      }
    }

    void printStuff(string a_prefix) const
    {
//      Chombo4::pout() << a_prefix << " BoxPattern constructor:"
//                      << " m_localPattern.size() = " << m_localPattern.size()
//                      << endl;
//
//      for(int iproc = 0; iproc < CH4_SPMD::numProc(); iproc++)
//      {
//        Chombo4::pout() << a_prefix << " BoxPattern constructor:"
//                        << " proc = " << iproc 
//                        << " m_recvPattern[iproc].size() = " << m_recvPattern[iproc].size()
//                        << " m_sendPattern[iproc].size() = " << m_sendPattern[iproc].size() << endl;
//                
//      }
    }
      
    ///Define for copyTo --- does NOT include ghost cells.
    BoxPattern(const Chombo4::DisjointBoxLayout & a_src,
               const Chombo4::DisjointBoxLayout & a_dst,
               bool a_printStuff = false)
    {
      CH_TIME("BoxPattern::exchange constructor");

      m_recvPattern.resize(CH4_SPMD::numProc());
      m_sendPattern.resize(CH4_SPMD::numProc());
      LayoutIterator litSrc = a_src.layoutIterator();
      LayoutIterator litDst = a_dst.layoutIterator();
      int myProcID = CH4_SPMD::procID();

      for(litSrc.begin(); litSrc.ok(); ++litSrc)
      {
        int srcProcID   = a_src.procID(litSrc());
        Chombo4::Box gridSrc   = a_src[litSrc()];
        for(litDst.begin(); litDst.ok(); ++litDst)
        {
          int dstProcID = a_dst.procID(litDst());
          Chombo4::Box gridDst = a_dst[litDst()];
          Chombo4::Box intersect(gridSrc & gridDst);
          if(!intersect.isEmpty())
          {
            Chombo4::Box dstValid = gridDst;
            Chombo4::Box srcValid = gridSrc;
            box_interaction_t item(Chombo4::DataIndex(litSrc()),
                                   Chombo4::DataIndex(litDst()),
                                   intersect, intersect,
                                   srcProcID, dstProcID, srcValid, dstValid);
            
            //if neither processor is == myprocID, nothing to do
            if((srcProcID == myProcID) || (dstProcID == myProcID))
            {
              if((srcProcID == myProcID) && (dstProcID == myProcID))
              {
                m_localPattern.push_back(item);
              }
              else if((srcProcID == myProcID) && (dstProcID != myProcID))
              {
                m_sendPattern[dstProcID].push_back(item);
              }
              else if((srcProcID != myProcID) && (dstProcID == myProcID))
              {
                m_recvPattern[srcProcID].push_back(item);
              }
              else
              {
                Chombo4::MayDay::Error("apparently I missed a case");
              }
            } //if(something is on proc)
          } //if there is an intersection
        } //loop over destination boxes
      } //loop over source boxes

      sortAllPatterns();
      if(a_printStuff)
      {
        printStuff(string("copyTo")); 
      }
    }

    void sortAllPatterns()
    {
      //sort everything in sight so that patterns on different procs show up in the same order
      std::sort(m_localPattern.begin(), m_localPattern.end());
      for(int ipat = 0; ipat < m_recvPattern.size(); ipat++)
      {
        std::sort(m_recvPattern[ipat].begin(), m_recvPattern[ipat].end() );
      }
      for(int ipat = 0; ipat < m_sendPattern.size(); ipat++)
      {
        std::sort(m_sendPattern[ipat].begin(), m_sendPattern[ipat].end() );
      }
    }

    //pattern where source == my proc, dest == my proc
    vector<box_interaction_t> m_localPattern;

    //pattern where source != my proc, dest == my proc
    /// outer loop is per proc 
    vector< vector<box_interaction_t> >m_recvPattern;

    //pattern where source == my proc, dest != my proc
    /// outer loop is per proc 
    vector< vector<box_interaction_t> >m_sendPattern;

    
  private:
    ///I like strong construction.
    BoxPattern();
    BoxPattern( const BoxPattern& a_input);
    void operator=(const BoxPattern& a_input);
  };

  
}
#endif

