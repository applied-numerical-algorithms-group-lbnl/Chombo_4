#ifndef __hoeb_composite_stencil__
#define __hoeb_composite_stencil__
#include "Proto.H"
#include "Chombo_EBChombo.H"
#include "Chombo_ParmParse.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_ParmParse.H"
#include "Hoeb_LAPACKMatrix.H"
#include "Hoeb_ExactSolutions.H"
#include "implem/Proto_Factorial.H"
#include "implem/Proto_AggStencil.H"

#define HOEB_MAX_ORDER 4

namespace hoeb
{

  class CompositeStencil
  {
  public:
    ///
    LocalStencil<CELL    , Real>      m_cellSten;
    LocalStencil<XFACE   , Real>      m_xfaceDiriSten;
    LocalStencil<YFACE   , Real>      m_yfaceDiriSten;
    LocalStencil<ZFACE   , Real>      m_zfaceDiriSten;
    
    LocalStencil<XFACE   , Real>      m_xfaceNeumSten;
    LocalStencil<YFACE   , Real>      m_yfaceNeumSten;
    LocalStencil<ZFACE   , Real>      m_zfaceNeumSten;


    LocalStencil<BOUNDARY, Real>      m_cutfaceDiriSten;
    LocalStencil<BOUNDARY, Real>      m_cutfaceNeumSten;

    unsigned int poutAll() const
    {
      using Chombo4::pout;
      
      pout() << "cellsten:" << endl;
      m_cellSten.print();
      pout() << "xfacediristen:" << endl;
      m_xfaceDiriSten.print();
      pout() << "yfacediristen:" << endl;
      m_yfaceDiriSten.print();
      pout() << "zfacediristen:" << endl;
      m_zfaceDiriSten.print();
    
      pout() << "xfaceneumsten:" << endl;
      m_xfaceNeumSten.print();
      pout() << "yfaceneumsten:" << endl;
      m_yfaceNeumSten.print();
      pout() << "zfaceneumsten:" << endl;
      m_zfaceNeumSten.print();


      pout() << "cutfacediristen:" << endl;
      m_cutfaceDiriSten.print();
      pout() << "cutfaceneumsten:" << endl;
      m_cutfaceNeumSten.print();
      return 0;
    }
    ///
    bool hasNeumannFaces() const
    {
      
      return ((  m_xfaceNeumSten.size() > 0) ||
              (  m_yfaceNeumSten.size() > 0) ||
              (  m_zfaceNeumSten.size() > 0) ||
              (m_cutfaceNeumSten.size() > 0));
    }
    
    ///
    bool hasDirichletFaces() const
    {
      return ((  m_xfaceDiriSten.size() > 0) ||
              (  m_yfaceDiriSten.size() > 0) ||
              (  m_zfaceDiriSten.size() > 0) ||
              (m_cutfaceDiriSten.size() > 0));
    }

    Real evaluate(const EBHostData<CELL, Real, 1>                     & a_phi,
                  const BaseExactSolution<HOEB_MAX_ORDER>             & a_exactSolution,
                  const EBGraph                                       & a_graph,
                  const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                  const Chombo4::Box                                  & a_domain,
                  const Real                                          & a_dx,
                  const Chombo4::DataIndex& a_dit) const
    {
      Real retval = 0;
      bool dividebyarea = true; //this is about values, not integrals
      //first the bits of the stencil that do not involve boundary conditions
      {
        for(int isten = 0; isten < m_cellSten.size(); isten++)
        {
          const auto& entry = m_cellSten.m_entries[isten];
          Real weight = entry.m_weight;
          EBIndex<CELL> phivof = entry.m_vof;
          Real phival = a_phi(phivof, 0);
          retval += weight*phival;
        }
      }
      //now add in contributions due to nearby boundary conditions
      //first dirichlet
      {
        const auto& xfacld = a_geoserv->getXFaceData(a_domain);
        const auto& xdat   = (*xfacld)[a_dit];
        for(int isten = 0; isten < m_xfaceDiriSten.size(); isten++)
        {
          const auto& entry  = m_xfaceDiriSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateDMO(a_graph, a_dx, xdat, face, 0, dividebyarea);
          retval += weight*phival;
        }
      }
      {
        const auto& yfacld = a_geoserv->getYFaceData(a_domain);
        const auto& ydat   = (*yfacld)[a_dit];
        for(int isten = 0; isten < m_yfaceDiriSten.size(); isten++)
        {
          const auto& entry  = m_yfaceDiriSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateDMO(a_graph, a_dx, ydat, face, 1, dividebyarea);
          retval += weight*phival;
        }
      }
#if DIM==3      
      {
        const auto& zfacld = a_geoserv->getZFaceData(a_domain);
        const auto& zdat   = (*zfacld)[a_dit];
        for(int isten = 0; isten < m_zfaceDiriSten.size(); isten++)
        {
          const auto& entry  = m_zfaceDiriSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateDMO(a_graph, a_dx, zdat, face, 2, dividebyarea);
          retval += weight*phival;
        }
      }
#endif      
      {
        const auto& ebfacld = a_geoserv->getEBFaceData(a_domain);
        const auto& ebdat   = (*ebfacld)[a_dit];
        for(int isten = 0; isten < m_cutfaceDiriSten.size(); isten++)
        {
          const auto& entry = m_cutfaceDiriSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateDIM(a_graph, a_dx, ebdat, face, dividebyarea);
          retval += weight*phival;
        }
      }
      
      //neumann
      {
        const auto& xfacld = a_geoserv->getXFaceData(a_domain);
        const auto& xdat   = (*xfacld)[a_dit];
        for(int isten = 0; isten < m_xfaceNeumSten.size(); isten++)
        {
          const auto& entry  = m_xfaceNeumSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateGradDMO(a_graph, a_dx, xdat, face, 0, dividebyarea);
          retval += weight*phival;
        }
      }
      {
        const auto& yfacld = a_geoserv->getYFaceData(a_domain);
        const auto& ydat   = (*yfacld)[a_dit];
        for(int isten = 0; isten < m_yfaceNeumSten.size(); isten++)
        {
          const auto& entry  = m_yfaceNeumSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateGradDMO(a_graph, a_dx, ydat, face, 1, dividebyarea);
          retval += weight*phival;
        }
      }
#if DIM==3      
      {
        const auto& zfacld = a_geoserv->getZFaceData(a_domain);
        const auto& zdat   = (*zfacld)[a_dit];
        for(int isten = 0; isten < m_zfaceNeumSten.size(); isten++)
        {
          const auto& entry  = m_zfaceNeumSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateGradDMO(a_graph, a_dx, zdat, face, 2, dividebyarea);
          retval += weight*phival;
        }
      }
#endif      
      {
        const auto& ebfacldx = a_geoserv->getEBNormalData_x(a_domain);
        const auto& ebfacldy = a_geoserv->getEBNormalData_y(a_domain);
        const auto& ebfacldz = a_geoserv->getEBNormalData_z(a_domain);
        const auto& ebmom    = a_geoserv->getEBFaceData(a_domain);
        for(int isten = 0; isten < m_cutfaceNeumSten.size(); isten++)
        {
          const auto& entry = m_cutfaceNeumSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          EBIndex<CELL> vof  = face.getVoF();
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateGradDIM(a_graph, a_dx,
                                                        ebfacldx, ebfacldy, ebfacldz, ebmom,
                                                        a_dit, vof, dividebyarea);
          retval += weight*phival;
        }
      }

      return retval;
    }

    //evaluate but with phi always zero
    Real evalInhomogeneousContribution(const BaseExactSolution<HOEB_MAX_ORDER>             & a_exactSolution,
                                       const EBGraph                                       & a_graph,
                                       const shared_ptr< GeometryService<HOEB_MAX_ORDER> > & a_geoserv,
                                       const Chombo4::Box                                  & a_domain,
                                       const Real                                          & a_dx,
                                       const Chombo4::DataIndex& a_dit) const
    {
      Real retval = 0;
      bool dividebyarea = true; //this is about values, not integrals

      //now add in contributions due to nearby boundary conditions
      //first dirichlet
      {
        const auto& xfacld = a_geoserv->getXFaceData(a_domain);
        const auto& xdat   = (*xfacld)[a_dit];
        for(int isten = 0; isten < m_xfaceDiriSten.size(); isten++)
        {
          const auto& entry  = m_xfaceDiriSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateDMO(a_graph, a_dx, xdat, face, 0, dividebyarea);
          retval += weight*phival;
        }
      }
      {
        const auto& yfacld = a_geoserv->getYFaceData(a_domain);
        const auto& ydat   = (*yfacld)[a_dit];
        for(int isten = 0; isten < m_yfaceDiriSten.size(); isten++)
        {
          const auto& entry  = m_yfaceDiriSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateDMO(a_graph, a_dx, ydat, face, 1, dividebyarea);
          retval += weight*phival;
        }
      }
#if DIM==3      
      {
        const auto& zfacld = a_geoserv->getZFaceData(a_domain);
        const auto& zdat   = (*zfacld)[a_dit];
        for(int isten = 0; isten < m_zfaceDiriSten.size(); isten++)
        {
          const auto& entry  = m_zfaceDiriSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateDMO(a_graph, a_dx, zdat, face, 2, dividebyarea);
          retval += weight*phival;
        }
      }
#endif      
      {
        const auto& ebfacld = a_geoserv->getEBFaceData(a_domain);
        const auto& ebdat   = (*ebfacld)[a_dit];
        for(int isten = 0; isten < m_cutfaceDiriSten.size(); isten++)
        {
          const auto& entry = m_cutfaceDiriSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateDIM(a_graph, a_dx, ebdat, face, dividebyarea);
          retval += weight*phival;
        }
      }
      
      //neumann
      {
        const auto& xfacld = a_geoserv->getXFaceData(a_domain);
        const auto& xdat   = (*xfacld)[a_dit];
        for(int isten = 0; isten < m_xfaceNeumSten.size(); isten++)
        {
          const auto& entry  = m_xfaceNeumSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateGradDMO(a_graph, a_dx, xdat, face, 0, dividebyarea);
          retval += weight*phival;
        }
      }
      {
        const auto& yfacld = a_geoserv->getYFaceData(a_domain);
        const auto& ydat   = (*yfacld)[a_dit];
        for(int isten = 0; isten < m_yfaceNeumSten.size(); isten++)
        {
          const auto& entry  = m_yfaceNeumSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateGradDMO(a_graph, a_dx, ydat, face, 1, dividebyarea);
          retval += weight*phival;
        }
      }
#if DIM==3      
      {
        const auto& zfacld = a_geoserv->getZFaceData(a_domain);
        const auto& zdat   = (*zfacld)[a_dit];
        for(int isten = 0; isten < m_zfaceNeumSten.size(); isten++)
        {
          const auto& entry  = m_zfaceNeumSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateGradDMO(a_graph, a_dx, zdat, face, 2, dividebyarea);
          retval += weight*phival;
        }
      }
#endif      
      {
        const auto& ebfacldx = a_geoserv->getEBNormalData_x(a_domain);
        const auto& ebfacldy = a_geoserv->getEBNormalData_y(a_domain);
        const auto& ebfacldz = a_geoserv->getEBNormalData_z(a_domain);
        const auto& ebmom    = a_geoserv->getEBFaceData(a_domain);
        for(int isten = 0; isten < m_cutfaceNeumSten.size(); isten++)
        {
          const auto& entry = m_cutfaceNeumSten.m_entries[isten];
          const auto& face   = entry.m_vof;
          EBIndex<CELL> vof  = face.getVoF();
          const auto& weight = entry.m_weight;
          Real phival = a_exactSolution.evaluateGradDIM(a_graph, a_dx,
                                                        ebfacldx, ebfacldy, ebfacldz, ebmom,
                                                        a_dit, vof, dividebyarea);
          retval += weight*phival;
        }
      }

      return retval;
    }
  };
}
#endif
