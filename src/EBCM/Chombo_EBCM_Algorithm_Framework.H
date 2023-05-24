#ifndef _EBCM_ALGO_DEV_FRAMEWORK__
#define _EBCM_ALGO_DEV_FRAMEWORK__
#include <cmath>
#include <cstdio>
#include <iostream>

#include "EBProto.H"
#include "Chombo_EBLevelBoxData.H"
#include "Chombo_LevelData.H"
#include "Chombo_BaseFab.H"

#include "Chombo_ParmParse.H"
#include "Chombo_LoadBalance.H"
#include "Chombo_ProtoInterface.H"
#include "Chombo_BRMeshRefine.H"
#include "Chombo_GeometryService.H"
#include "Chombo_EBDictionary.H"
#include "Chombo_EBChombo.H"
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBCM_HostLevelData.H"
#include "Chombo_EigenMatrix.H"

#include "Chombo_LAPACKMatrix.H"
#include <iomanip>

namespace EBCM
{
  typedef Chombo4::Box                                   ch_box;
  typedef Chombo4::DisjointBoxLayout                     ch_dbl;
  typedef Chombo4::ProblemDomain                         ch_probdom;
  typedef Chombo4::IntVect                               ch_iv;
  typedef Chombo4::MayDay                                ch_mayday;
  typedef Proto::RealVect                                pr_rv;
  typedef Proto::IndexTM<Real, DIM>                      pr_itm_r_dim;
  typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
  typedef Proto::IndexTM<Real, DIM-1>                    pr_itm_r_dmo;
  typedef Proto::IndexTM<int , DIM-1>                    pr_itm_i_dmo;
  typedef Proto::BaseIF                                  pr_baseif;
  typedef ch_eigen::Matrix                               eigen_mat;
///
/**
   Purely functional framework for doing tests in EBCM.
   The framework has no data but lots of common  typedefs.
   All functions are static.
*/
  template <int ebcm_order>
  class Algorithm_Framework
  {
  public:
    ///test function for example/EBCM/hoeb_petsc_test
  static  void run_hoeb_petsc_tests()
    {
    CH_TIME("EBCM_Framework::runTests");
    //linear algebra machinery check
    runEigenTests();
    //the important stuff
    //see if we can make a geometry
    shared_ptr< ebcm_meta  > ebcm;
    ch_dbl                   grids;
    double                   dx;
    ch_probdom               domain;
    makeMergedGeometry(ebcm, grids, dx, domain, true, true);

    Real maxKapp, minKapp;
    //this tries out iteration and checks kappa
    checkKappa(maxKapp, minKapp, ebcm, grids);

    //makes data of inv condition number info and plots it out
    shared_ptr<ebcm_leveldata> invCondition;
    shared_ptr<condition_t>    worst;
    getConditionNumberData(invCondition, worst, ebcm);

    worst->poutLatexTable();
    Chombo4::pout() << "worst condition number at point " << worst->m_pt << endl;
    
    string plot_filename;
    ParmParse ppmain("main");  //kinda pendantic to use anything else
    ppmain.get("inv_cond_filename",   plot_filename);
    invCondition->writeToHDF5(plot_filename);
    
  }
    Algorithm_Framework()
    {}    
  

  
    typedef    EBCM::MetaDataLevel<         ebcm_order >   ebcm_meta;
    typedef    EBCM::SubVolumeVector  <     ebcm_order >   ebcm_subvol_vec;
    typedef    EBCM::EBCM_Volu<             ebcm_order >   ebcm_volu;
    typedef    EBCM::EBCM_Graph<            ebcm_order >   ebcm_graph;
    typedef    EBCM::HostLevelData<Real, 1, ebcm_order >   ebcm_leveldata;

    typedef Chombo4::GeometryService<       ebcm_order >   ch_geoserv;
    typedef   Proto::IndexedMoments<DIM  ,  ebcm_order >   pr_mom_dim;
    typedef   Proto::IndexedMoments<DIM-1,  ebcm_order >   pr_mom_dmo;
    typedef   Proto::Point pr_pt;

  
/****/

    ///
    struct  neighborhood_t
    {
  
      ch_iv                            m_startiv;
      int                              m_nghost;
      pr_rv                            m_startloc;
      shared_ptr<ebcm_subvol_vec>      m_volumes;
      //this is the actual distance stuff gets shifted
      vector<pr_rv>                    m_xbar;
      ///this is the distance used in  calculating weights
      vector<Real>                     m_distance;
      ///this is the weight of each equation (the diagonal of the weight matrix)
      vector<Real>                     m_eqweight;
  
      neighborhood_t(const ebcm_graph    & a_graph,
                     const ch_iv         & a_start,
                     int a_nghost, int a_weightPower, bool a_printStuff)
      {
        CH_TIME("Algorithm_Framework::neigborhood_t constructor");
        m_startiv = a_start;
        m_nghost  = a_nghost;
        m_startloc.setToCCLocation(a_start, a_graph.m_dx);

        ch_box region(a_start, a_start);
        region.grow(a_nghost);
        region &= a_graph.m_domain;

        m_volumes = shared_ptr<ebcm_subvol_vec>(new ebcm_subvol_vec(a_graph, region, a_printStuff));
        //begin debug
        //Chombo4::pout() << "first volmom in the neighborhood" << endl;
        //(*m_volumes)[0].m_volmom.print();
        //end debug
        int n_equations =  m_volumes->size();
        int n_unknowns  = pr_mom_dim::size();
        if(n_unknowns > n_equations)
        {
          Chombo4::pout() << a_nghost <<
            " ghost cells not getting enough volumes with start = "
                          << a_start << endl;
          ch_mayday::Error();
        }
        m_xbar.resize(    m_volumes->size());
        m_distance.resize(m_volumes->size());
        m_eqweight.resize(m_volumes->size());
    
        for(int ivec = 0; ivec < m_volumes->size(); ivec++)
        {
          const auto& curVol = (*m_volumes)[ivec];
          m_xbar[ivec] = curVol.m_centroid - m_startloc;
          m_distance[ivec] = getManhattanDistance(curVol);
          double eqnweight = 1;
          if(m_distance[ivec] > 1)
          {
            for(int iweight = 0; iweight < a_weightPower; iweight++)
            {
              eqnweight /= m_distance[ivec];
            }
          }
          m_eqweight[ivec] = eqnweight;
          //debug
          //m_eqweight[ivec] = 1;
          //end debug
        }
      }

      double getManhattanDistance(const ebcm_volu & a_volu)  const
      {
        double xbar = 1; //important that it not go to zero.
        for(int idir = 0; idir < DIM; idir++)
        {
          int idiff = a_volu.m_pt[idir] - m_startiv[idir];
          double rdiff = std::abs(double(idiff));
          xbar += rdiff;
        }
        return xbar;
      }

    };
    ///
    static inline
    shared_ptr<eigen_mat>
    getWeightMatrix(shared_ptr<neighborhood_t> a_locality, bool a_printStuff)
    {
      int nrows = a_locality->m_eqweight.size();
      shared_ptr<eigen_mat> weight_p(new eigen_mat(nrows, nrows));
      weight_p->setVal(0.);
      for(int irow = 0; irow < nrows; irow++)
      {
        (*weight_p)(irow, irow) = a_locality->m_eqweight[irow];
      }
      return weight_p;
    }


    ///
    static inline void
    shiftMomentAndFillRow(eigen_mat         & a_mat,
                          const pr_mom_dim  & a_moment,
                          const pr_rv       & a_distance,
                          const int         & a_currentRow,
                          bool a_printStuff)
    {
      CH_TIME("Algorithm_Framework::shiftMomentAndFillRow");
      if(a_printStuff)
      {
        Chombo4::pout() << "shiftMomentAndFillRow: a_currentRow = " << a_currentRow << endl;
        Chombo4::pout() << "shiftMomentAndFillRow: a_distance   = " <<   a_distance << endl;
        Chombo4::pout() << "shiftMomentAndFillRow: input mat row: " << endl;
        a_mat.poutRow(a_currentRow);
        Chombo4::pout() << "shiftMomentAndFillRow:   a_moment =" << endl;
        a_moment.print();
      }

      pr_itm_r_dim itm_dist;
      for(int idir = 0; idir < DIM; idir++)
      {
        itm_dist[idir] = a_distance[idir];
      }

      pr_mom_dim shiftedMom = a_moment;
      shiftedMom.shift(itm_dist);
      if(a_printStuff)
      {
        Chombo4::pout() << "shiftMomentAndFillRow: shiftedMom =" << endl;
        shiftedMom.print();
      }

      Real volume = shiftedMom[pr_itm_i_dim::Zero];
      bool divide = (volume > 1.0e-16);
      if(!divide)
      {
        if(a_printStuff)
        {
          Chombo4::pout() << "shiftMomentAndFillRow: smallCellRow with volume = " << volume << endl;
        }
        a_mat.setSmallCellRow(a_currentRow);
      }
      else
      {
        for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
        {
          auto momind = momit();
          unsigned int currentCol = pr_mom_dim::indexOf(momind);
          a_mat(a_currentRow, currentCol) = shiftedMom[momind]/volume;
        }
      }
      if(a_printStuff)
      {
        Chombo4::pout() << "shiftMomentAndFillRow: amat leaving = "  << endl;
        a_mat.poutRow(a_currentRow);
      }
    }

    ///
    static shared_ptr<eigen_mat>
    getMomentMatrix(shared_ptr<neighborhood_t> a_locality,
                    bool a_printStuff)
    {
      CH_TIME("Algorithm_Framework::getMomentMatrix");
      //number of equations
      int n_equations = a_locality->m_volumes->size();
      int n_unknowns  = pr_mom_dim::size();
      if(n_unknowns > n_equations)
      {
        ch_mayday::Error("getMomentMatrix::bring me more equations");
      }
      int n_rows = n_equations;
      int n_cols = n_unknowns;
    
      shared_ptr<eigen_mat> Mmat_p(new eigen_mat(n_rows, n_cols));
      for(int irow = 0; irow < n_rows; irow++)
      {
        const auto& volmom   =   (*(a_locality->m_volumes ))[irow].m_volmom;
        if(a_printStuff)
        {
          volmom.print();
        }
        pr_rv & xbar =    ( (a_locality->m_xbar))[irow];
        shiftMomentAndFillRow(*Mmat_p, volmom, xbar, irow, a_printStuff);
      }
      return Mmat_p;
    }

    ///
    static inline shared_ptr<eigen_mat>
    getAMatrix(const ebcm_volu     &  a_volu,
               const ebcm_graph    &  a_graph,
               int a_nghost, int a_weightPower,
               bool a_printStuff)
    {
  
      CH_TIME("Algorithm_Framework::AMatrix");
      shared_ptr<neighborhood_t>
        locality(new neighborhood_t(a_graph, a_volu.m_pt, a_nghost, a_weightPower, a_printStuff));
  
      shared_ptr<eigen_mat> Wmat_p = getWeightMatrix(locality, a_printStuff);
      shared_ptr<eigen_mat> Mmat_p = getMomentMatrix(locality, a_printStuff);

      eigen_mat WMmat;
      multiply(WMmat, *Wmat_p, *Mmat_p);
      eigen_mat WMTmat = WMmat;
      WMTmat.transpose();

      shared_ptr<eigen_mat> Amat_p(new eigen_mat());
      multiply(*Amat_p, WMTmat, WMmat);
    
      if(a_printStuff)
      {
        Chombo4::pout() << "Weight matrix: " << endl;
        Wmat_p->poutAll();
        Chombo4::pout() << "Moment matrix: " << endl;
        Mmat_p->poutAll();
      
        Chombo4::pout() << "WM   matrix: " << endl;
        WMmat.poutAll();
      
        Chombo4::pout() << "WM^T matrix: " << endl;
        WMTmat.poutAll();

        Chombo4::pout() << "A matrix: " << endl;
        Amat_p->poutAll();
      }

      return Amat_p;
    }

  }; //end class Algorithm_Framnework 

} //end namespace ebcm  
#endif
