#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Hoeb_PetscSolver__H_
#define _Hoeb_PetscSolver__H_

#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "Chombo_EBChombo.H"
#include "EBPetscSolver.H"
#include "Hoeb_Utilities.H"
#include "Chombo_EBLevelBoxData.H"

namespace hoeb
{
/// 
/**
*/
  template <int order>
  class Hoeb_PetscSolver: public BaseEBPetscSolver<order>
  {
  public:
    typedef GraphConstructorFactory< EBBoxData<CELL, int, 1> >  devifactint_t;
    typedef GraphConstructorFactory<EBHostData<CELL, int, 1> >  host_fact_int_t;
    typedef GraphConstructorFactory<EBHostData<CELL,Real, 1> >  host_fact_real_t;
    typedef EBDictionary<order, Real, CELL, CELL>                dictionary_t;
    typedef CH4_Data_Choreography::DistributedData<EBBoxData< CELL, Real, 1> > devi_distrib_real_t;
    typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, Real, 1> > host_distrib_real_t;
    typedef CH4_Data_Choreography::DistributedData<EBBoxData< CELL, int , 1> > devi_distrib_int_t;
    typedef CH4_Data_Choreography::DistributedData<EBHostData<CELL, int , 1> > host_distrib_int_t;
    typedef CH4_Data_Choreography::DistributedData<EBGraph>                   graph_distrib_t;

    virtual ~Hoeb_PetscSolver()
    {
    }
    
    Hoeb_PetscSolver(const shared_ptr<GeometryService<order> >  & a_geoserv,
                     const shared_ptr<dictionary_t           >  & a_ebdictionary,
                     const shared_ptr<graph_distrib_t>          & a_graphs,
                     const Chombo4::DisjointBoxLayout           & a_grids,
                     const Chombo4::Box                         & a_domain,
                     string a_domainBCName[2*DIM],
                     string a_ebbcName,
                     Real a_dx,   Real a_alpha, Real a_beta, Point a_ghost,
                     int a_stenRadius, int a_maxStenRadius,
                     bool a_printStuff = false)
      :BaseEBPetscSolver<order>(a_geoserv,           
                                a_ebdictionary,      
                                a_graphs,            
                                a_grids,             
                                a_domain,            
                                a_dx,                
                                a_alpha,             
                                a_beta,              
                                a_ghost,             
                                a_printStuff)
    {
      for(int iface = 0; iface < 2*DIM; iface++)
      {
        this->m_domainBCName[iface] = a_domainBCName[iface];
      }
      this->m_ebbcName = a_ebbcName;
      m_stenRadius    = a_stenRadius;
      m_maxStenRadius = a_maxStenRadius;
      formMatrix();
    }

//    void
//    modifyRHSForInhomogeneity(EBLevelBoxData<CELL, 1>& a_rhs, bool a_printStuff = false)
//    {
//      EBLevelBoxData<CELL, 1> phi(this->m_grids, this->m_ivghost, this->m_graphs);
//      EBLevelBoxData<CELL, 1> klp(this->m_grids, this->m_ivghost, this->m_graphs);
//      phi.setVal(0.);
//
//      hoeb::getDevendranKLPInhomogeneous(klp, phi,
//                                         this->m_graphs,       this->m_grids,
//                                         this->m_domain,       this->m_dx,
//                                         this->m_ebdictionary, this->m_geoserv,
//                                         a_printStuff);
//
//      Chombo4::DataIterator dit = phi.dataIterator();
//      for(int ibox = 0; ibox < dit.size(); ibox++)
//      {
//        a_rhs[dit[ibox]] -= klp[dit[ibox]];
//      }
//    }
    
  private:
  
    virtual PetscInt formMatrix(bool a_printStuff = false)
    {
    
      CH_TIME("Hoeb_PetscSolver::formMatrix");
      {
        using Chombo4::DataIterator;
        char str[256];
        strcpy (str,"-");
        strcat (str,this->m_prestring);
#if PETSC_VERSION_GE(3,6,0)
        strcat (str,"pc_gamg_square_graph 20");
#else
        strcat (str,"pc_gamg_square_graph true");
#endif
        PetscInt  ierr;
#if PETSC_VERSION_GE(3,7,0)
        ierr = PetscOptionsInsertString(PETSC_NULL,str);CHKERRQ(ierr);
#else
        ierr = PetscOptionsInsertString(str);CHKERRQ(ierr);
#endif
        DataIterator dit = this->m_grids.dataIterator();
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          Bx srcValid = ProtoCh::getProtoBox(this->m_grids[dit[ibox]]);
          Bx dstValid = ProtoCh::getProtoBox(this->m_grids[dit[ibox]]);
          Bx srcDom = ProtoCh::getProtoBox(  this->m_domain);
          Bx dstDom = ProtoCh::getProtoBox(  this->m_domain);

          vector< EBIndex<CELL> >            dstVoFs;
          vector< Proto::LocalStencil<CELL, Real> > wstencil;
          Stencil<Real>    regStencilInterior;
          bool needDiag = false; //abandon all hope of goemetric multigrid
          hoeb::getHomogeneousDharshiStencil(this->m_stencilName, 
                                             dstVoFs, wstencil, srcValid, dstValid,
                                             srcDom, dstDom,  needDiag,
                                             this->m_geoserv,      this->m_grids,
                                             this->m_domain,       this->m_dx,   ibox,
                                             this->m_domainBCName, this->m_ebbcName,
                                             this->m_alpha, this->m_beta,
                                             m_stenRadius, m_maxStenRadius, a_printStuff);
                                             
          //now create the matrix
          for(int ivof = 0; ivof < dstVoFs.size(); ivof++)
          {
            auto vof    =   dstVoFs[ivof];
            auto& stenc =  wstencil[ivof];
            int irow = this->m_gids[dit[ibox]](vof, 0); 
            for(int ivof = 0; ivof < stenc.size(); ivof++)
            {
              auto& stenvof = stenc.m_entries[ivof].m_vof;
              auto& stenwgt = stenc.m_entries[ivof].m_weight;
              Real weight = this->m_beta*stenwgt;
              if(stenvof == vof)
              {
                weight += this->m_alpha;
              }
              int icol = this->m_gids[dit[ibox]](stenvof, 0);  
              PetscInt irowpet = irow;
              PetscInt icolpet = icol;
              ierr = MatSetValues(this->m_mat,1,&irowpet,1,&icolpet,&weight,INSERT_VALUES);
            } 
          }
        }

        ierr = MatAssemblyBegin(this->m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(  this->m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      }
      return 0;
    }

  private:
    int m_stenRadius, m_maxStenRadius;
  }; //end of class
} //namespace hoeb
#endif

#endif 
