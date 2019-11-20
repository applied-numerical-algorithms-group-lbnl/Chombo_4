
#ifndef _Chombo_EBPETScOperator_H_
#define _Chombo_EBPETScOperator_H_

#ifdef CH_USE_PETSC

#include "petsc.h"
#include "Proto.H"
#include "EBProto.H"
#include "Chombo_LevelData.H"
#include "Chombo_EBChombo.H"
#include "Chombo_GeometryService.H"
#include "Chombo_NamespaceHeader.H"


///
/**
   This is an interface between Chombo EB operators and Petsc matrices.
*/
template <int order>
class EBPETScOperator
{

public:
  ///
  /**
     This interface is for a single level
     a_ghost is the number of ghost cells your data will hold.
     It must be at least as big as the stencil span in each direction..
  */
  EBPETScOperator(const shared_ptr<GeometryService<order> >  & a_geoserv,
                  const DisjointBoxLayout                    & a_grids,
                  const Box                                  & a_domain,
                  const Real                                 & a_dx)
  {
    MayDay::Error("not implemented");
  }


  

  ///
  /**
     diagonal weights are for relaxation--expensive to compute so leave them out if you do not need them
     dest grow is how many cells to grow valid box for destination.
  */
  inline void registerStencil(string a_stencilName,
                              string a_domainBCName,
                              string a_ebbcName)
  {
    MayDay::Error("not implemented");
  }

  ///
  shared_ptr<Vec> getVector(const EBLevelBoxData<CELL, 1>& a_data)
  {
    shared_ptr<Vec> retval;
    MayDay::Error("not implemented");
    return retval;
  }

  ///
  shared_ptr<Mat> getMatrix(string a_stencilName,
                            string a_ebbcName)
  {
    shared_ptr<Mat> retval;
    MayDay::Error("not implemented");
    return retval;
  }

  ///after this is done, you cannot register any more stencils.  This will delete the moment information.
  void  close()
  {
    m_geoserv = shared_ptr<GeometryService<order> >();
  }

  void vectorExtract(EBLevelBoxData<CELL, 1>& a_data, 
                     const shared_ptr<Vec>  & a_vec)
  {
    MayDay::Error("not implemented");
  }
protected:


  struct stencilKey
  {
    string m_ebbcName;
    string m_stencilName;
    bool operator< (const stencilKey& a_key) const
    {
      bool retval = false;
      if(m_stencilName == a_key.m_stencilName)
      {
        retval= (m_ebbcName < a_key.m_ebbcName);
      }
      else
      {
        retval= (m_stencilName < a_key.m_stencilName);
      }
      return retval;
    }
  };
 
  /// one per box in the grid
  map< stencilKey, shared_ptr< Mat> > m_stencils;
    
  shared_ptr<GeometryService<order> > m_geoserv;
  bool                                m_isClosed;

  DisjointBoxLayout                   m_grids;
  Box                                 m_domain;
  Real                                m_dx;
  PetscInt                            m_nrows;
  
  LevelData<EBHostData<CELL, int, 1>  > m_mapToMatrixRows;

private:
  EBPETScOperator();
  EBPETScOperator(  const EBPETScOperator<order>& a_inp);
  void operator=(const EBPETScOperator<order>   & a_inp);


};

#include "Chombo_NamespaceFooter.H"
#endif
#endif
