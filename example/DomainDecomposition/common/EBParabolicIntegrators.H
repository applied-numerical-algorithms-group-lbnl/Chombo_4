#ifndef _EBParabolic_H_
#define _EBParabolic_H_
#include <cmath>
#include <memory>
#include "Proto.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBLevelBoxData.H"
#include "EBMultigrid.H"
#include "Chombo_NamespaceHeader.H"
///
/**
 */
class BaseEBParabolic
{
public:
  ///the boundary conditions for the EBMultigrid have to be correct.  Alpha and Beta do not as they are constantly reset
  BaseEBParabolic(shared_ptr<EBMultigrid>   a_diffusionSolver,
                  const DisjointBoxLayout & a_grids)
  {
    m_diffusionSolver = a_diffusionSolver;
    m_grids           = a_grids;
  }

  ///
  virtual ~BaseEBParabolic()
  {
  }

  ///this one will be different for different integrators
  virtual void advanceOneStep( EBLevelBoxData<CELL, 1>       &  a_phinew,
                               const EBLevelBoxData<CELL, 1> &  a_phiold,
                               const Real                    &  a_diffCoef,
                               const Real                    &  a_dt) = 0;

  ///this is the conservative way to advance a solution with diffusion.
  inline
  void computeDiffusion( EBLevelBoxData<CELL, 1>       &  a_diffusionTerm,
                         const EBLevelBoxData<CELL, 1> &  a_phiold,
                         const Real                    &  a_diffCoef,
                         const Real                    &  a_dt);
                        
protected:
  shared_ptr<EBMultigrid> m_diffusionSolver;  
  DisjointBoxLayout       m_grids;
private:
  
  //outlawed for all the usual reasons
  BaseEBParabolic();
  BaseEBParabolic(const  BaseEBParabolic& a_input);
  void operator= (const  BaseEBParabolic& a_input);
};


#include "Chombo_NamespaceFooter.H"
#endif


