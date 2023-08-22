#ifndef __Chombo_EBCM_PETSC_Framework_
#define __Chombo_EBCM_PETSC_Framework_

#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "petscmat.h"
/*
  Chombo_EBCM_Graph.H is where MetaDataLevel lives.  It made sense at the time.
  I should use this opportunity to apologize to the world for the many 
  awful choices I have made in this life in relation to file names and organization.  
  I regret all the consternation I have caused.
*/
#include "Chombo_EBCM_Graph.H"
#include "Chombo_EBLevelBoxData.H"

/// 
/**
   An EBCM interface to the magical world of PETSc, 
   wherein the evil bits of distributed linear algebra are done for us.
   I plan to use this for everything.  Geometric multigrid requires
   more infrastructure than I am willing to write so I will use PETSc
   for solves and truncation error tests and I will be hooking in
   SLEPc to do the eigenvalue analysis. --dtg
**/
namespace EBCM
{

  
  template <int ebcm_order>
  class PETSc_Framework
  {
  public:
    typedef Chombo4::Box                                   ch_box;
    typedef Chombo4::DisjointBoxLayout                     ch_dbl;
    typedef Chombo4::BoxIterator                           ch_bit;
    typedef Chombo4::IntVect                               ch_iv;
    typedef Chombo4::DataIterator                          ch_dit;
                                                           
    typedef ::Mat                                          petsc_mat;
    typedef Proto::Point                                   pr_pt;
    
    typedef  EBCM::MetaDataLevel<             ebcm_order>  ebcm_meta;
    typedef  EBCM::EBCM_Graph<                ebcm_order>  ebcm_graph;
    typedef  EBCM::HostLevelData<int,    1,   ebcm_order>  ebcm_int_sca_data;
    typedef  EBCM::HostLevelData<double, 1,   ebcm_order>  ebcm_dou_sca_data;
    typedef  EBCM::HostLevelData<int,    DIM, ebcm_order>  ebcm_int_vec_data;
    typedef  EBCM::HostLevelData<double, DIM, ebcm_order>  ebcm_dou_vec_data;
    typedef  EBCM::SubVolumeVector<           ebcm_order>  ebcm_subvol_vec;
    typedef  EBCM::EBCM_Volu<                 ebcm_order>  ebcm_volu;
    typedef  EBCM::Algorithm_Framework<       ebcm_order>  algo_framework;
    typedef Proto::IndexedMoments<DIM  ,      ebcm_order > pr_mom_dim;
    typedef Proto::IndexedMoments<DIM-1,      ebcm_order > pr_mom_dmo;
    typedef Proto::IndexTM<double, DIM>                    pr_itm_r_dim;
    typedef Proto::IndexTM<int , DIM>                      pr_itm_i_dim;
    typedef Proto::IndexTM<double, DIM-1>                  pr_itm_r_dmo;
    typedef Proto::IndexTM<int , DIM-1>                    pr_itm_i_dmo;
    
    typedef EBCM::neighborhood       <  ebcm_order>       ebcm_neighborhood;
    /// Class to replace old vofstencil.  
    /**
       This will have to be generalized for multiple vofs per cell.
     **/
    struct local_stencil_t: public std::vector< pair< pr_pt, double > >
    {
    public:
 
      ///
      /**
         adds each member of input
      */
      void
      operator+= (const local_stencil_t& a_input)
      {
        for(int ivec = 0; ivec < a_input.size(); ivec++)
        {
          add(a_input[ivec]);
        }
      }
      
      ///
      void clear()
      {
        this->resize(0);
      }

      /// 
      /**
         Append the entry to the vector if the point is not already in the vector.
         Otherwise, just increment the weight.
      */
      void add(const pair<pr_pt, double>& a_entry)
      {
        bool alreadyhere = false;
        for (int ivof = 0; ivof < this->size(); ivof++)
        {
          const pr_pt& curpt = (*this)[ivof].first;
          if(curpt == a_entry.first)
          {
            alreadyhere = true;
            (*this)[ivof].second += a_entry.second;
          }
        }
        if (!alreadyhere)
        {
          this->push_back(a_entry);
        }
      }

      ///scale the stencil by a constant
      void operator*= (const double& a_scale)
      {
        for (int ivof = 0; ivof < this->size(); ivof++)
        {
          (*this)[ivof].second *= a_scale;
        }
      }
        
    }; //end struct local_stencil_t
    
    ///
    /**
       This is a wrapper around PETSc's  ::Mat object.  The specfics of the stencil 
       are computed in the derived class.   All data is public because
       this design is cleaner than hiding behind access functions.
    **/
    class Base_PETSc_Op
    {
    public:
      shared_ptr<petsc_mat> m_mat_ptr;
      shared_ptr<ebcm_int_sca_data> m_row_map_ptr; 
      shared_ptr<ebcm_meta> m_meta_ptr;
      int                   m_ghost;
      //each box in data iterator we get total points
      std::vector<int>      m_numPtsPerBox;
      int                   m_numPtsThisProc;
      int                   m_numPtsAllProcs;
      int                   m_startPtThisProc;
      
      virtual ~Base_PETSc_Op()
      {;}
  
      Base_PETSc_Op(const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                    const int                                    & a_ghost,
                    bool                                           a_print = false)
      {
        CH_TIME("EBPetscSolver::EBPetscSolver");
        m_meta_ptr     =    a_meta_ptr;
        m_ghost        =    a_ghost;
    
        if(a_print)
        {
          Chombo4::pout() << "Base_Petsc_Op: creating map of locations in space to matrix row." << endl;
        }
        createMatrixRowMap(a_print);

        if(a_print)
        {
          Chombo4::pout() << "Base_Petsc_Op:: creating the matrix" << endl;
        }

        createOperatorPetscMatrix(a_print);
      }

      ///this gets done by the derived class
      virtual shared_ptr<local_stencil_t> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph, bool a_print = false) = 0;


    private:
      /// banning weak construction but assignment and copy constructors should work fine.
      Base_PETSc_Op();
      
    protected:
      ///makes m_row_map_ptr (operator uses this to map between volumes and petsc matrix rows)
      PetscInt createMatrixRowMap(bool a_print = false)
      {
        m_row_map_ptr = shared_ptr<ebcm_int_sca_data>(new ebcm_int_sca_data(m_meta_ptr, m_ghost));
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();

        //get all the info on points on this proc.
        int numPtsThisProc = 0;
        vector<int> numPtsPerBox(dit.size(), 0);
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          //valid is correct here.  ghost is handled via exchange.

          ebcm_subvol_vec subvollocal(graph, valid, a_print);
          numPtsThisProc     += subvollocal.size();
          numPtsPerBox[ibox]  = subvollocal.size();
        }

        Chombo4::pout() << "Base_Petsc_Op::createMatrixRowMap: numptsThisProc = " << numPtsThisProc << endl;

#ifdef CH_MPI
        std::vector<int> numPtsAllProcs(CH4_SPMD::numProc());
        MPI_Gather(&numPtsThisProc, 1, MPI_INT, &numPtsAllProcs[0], 1, MPI_INT, 0, Chombo_MPI::comm);
    
        MPI_Bcast(numPtsAllProcs.data(), numPtsAllProcs.size(), MPI_INT, 0, Chombo_MPI::comm);

        int totalNumPts = 0;
        for(int iproc = 0; iproc < numPtsAllProcs.size(); iproc++)
        {
          totalNumPts += numPtsAllProcs[iproc];
        }
        //decide which location maps to the first one of this proc
        int startgid = 0;
        for(int iproc = 0; iproc <  procID(); iproc++)
        {
          startgid += numPtsAllProcs[iproc];
        }
        Chombo4::pout()  << "Base_Petsc_Op::getMap: procID = " << procID() << ", numPtsThisProc = " << numPtsThisProc << ", startgid = " << startgid <<", totalNumPts = " << totalNumPts <<  endl;
#else
        int startgid = 0;
        PetscInt totalNumPts = numPtsThisProc;
#endif
        m_startPtThisProc      = startgid;
        m_numPtsThisProc       = numPtsThisProc;
        m_numPtsAllProcs       = totalNumPts;
        m_numPtsPerBox         = numPtsPerBox;
        auto& mapdata = *(m_row_map_ptr->m_data);
        int curgid = startgid;

        //get all the info on points on this proc into the map
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          ebcm_subvol_vec volvec(graph, valid, a_print);
          auto      & mapfab =      mapdata[dit[ibox]];
          for(int ivec = 0; ivec < volvec.size(); ivec++)
          {
            const auto& volume = volvec[ivec];
            mapfab(volume.m_pt, 0) = curgid;
            for(int icell = 0; icell < volume.m_cells.size(); icell++)
            {
              const auto& cell = volume.m_cells[icell];
              mapfab(cell, 0) = curgid;
            }
            curgid++;
          }
        }

        Chombo4::pout() << "Base_Petsc_Op::getMap calling exchange" << endl;

        m_row_map_ptr->exchange(a_print);
        
        return 0;
      }//end function createMatrixRowMap

      ///makes m_mat_ptr (operator uses this to map between volumes and petsc matrix rows)
      PetscInt createOperatorPetscMatrix(bool a_print = false)
      {
        CH_TIME("Base_Petsc_Op::createOperatorPetscMatrix");
        // create matrix
        PetscInt nnzrow = 0;
#ifdef CH_MPI
        MPI_Comm wcomm = Chombo_MPI::comm;
#else
        MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
        PetscInt ierr;
        ierr = MatCreate(wcomm, & (*m_mat_ptr));CHKERRQ(ierr);
        ierr = MatSetOptionsPrefix(*m_mat_ptr,"");CHKERRQ(ierr);

        ierr = MatSetSizes(*m_mat_ptr,m_numPtsThisProc,m_numPtsThisProc,m_numPtsAllProcs,m_numPtsAllProcs);CHKERRQ(ierr);
        ierr = MatSetType( *m_mat_ptr,MATDENSE);CHKERRQ(ierr);
#ifdef CH_MPI
        ierr = MatMPIDenseSetPreallocation(*m_mat_ptr, NULL);CHKERRQ(ierr);
#else
        PetscInt matsize = m_numPtsThisProc*m_numPtsThisProc;
        PetscScalar* matspace = new PetscScalar[matsize];
        ierr = MatSeqDenseSetPreallocation(*m_mat_ptr, matspace);CHKERRQ(ierr);
#endif    

        ierr = MatSetFromOptions( *m_mat_ptr ); CHKERRQ(ierr);

        //sets the actual values in *m_mat_ptr
        ierr = formMatrix(a_print); CHKERRQ(ierr);
        return 0;
      }

      ///sets the actual values in *m_mat_ptr
      PetscInt formMatrix(bool a_print = false)
      {
    
        CH_TIME("Base_Petsc_Op::formMatrix");
        
        const auto& graphs = *(m_meta_ptr->m_graphs);
        const auto& grids  =  (m_meta_ptr->m_grids);
        auto dit = grids.dataIterator();

        //get all the info on points on this proc.
        int numPtsThisProc = 0;
        vector<int> numPtsPerBox(dit.size(), 0);
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =  grids[dit[ibox]];
          const auto& graph = graphs[dit[ibox]];
          ebcm_subvol_vec volvec(graph, valid, a_print);
          for(int ivol = 0; ivol < volvec.size(); ivol++)
          {
            const auto& volu = volvec[ivol];
            shared_ptr<local_stencil_t>  local = getLocalStencil(volu.m_pt, graph);
            const auto    & maplev = *(m_row_map_ptr->m_data);
            const auto    & mapfab = maplev[dit[ibox]];
            PetscInt irow = mapfab(volu.m_pt, 0);
            for(int isten = 0; isten < local->size(); isten++)
            {
              const auto& stenpair = (*local)[isten];
              PetscInt icol = mapfab(stenpair.first, 0);
              double weight = stenpair.second;
              /// This is PETSc for
              /// (*m_mat)(irow, icol) = weight;
              MatSetValues(*(this->m_mat_ptr),1,&irow,1,&icol,&weight,INSERT_VALUES);
            }
          }
        }
        return 0;
      }

    }; //end class Base_Petsc_Op

    /// 
    /**
       Class to force compilation of all this lovely code. It is useful but unusable.  
       If a class can be a zen koan, this may qualify.
    **/
    class DummyDerivedOp: public Base_PETSc_Op
    {
    public:
      DummyDerivedOp(const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                     const int                                    & a_ghost,
                     bool                                           a_print = false)
        :Base_PETSc_Op(a_meta_ptr, a_ghost, a_print)
      {
      }

      virtual ~DummyDerivedOp()
      {
      }
      
      virtual shared_ptr<local_stencil_t> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph, bool a_print = false)
      {
        shared_ptr<local_stencil_t> retval(new local_stencil_t());
        return retval;
      }
    };
    

    ///
    template <CENTERING cent>
    static shared_ptr<eigen_mat>
    getHelmholtzQMatrix( const EBCM_Face<cent, ebcm_order>  & a_face,
                         shared_ptr<ebcm_neighborhood>        a_locality,
                         const int                          & a_facedir,
                         const double                       & a_dx,
                         bool a_print)
    {
      pr_mom_dmo facmom =  a_face.m_facmom;
      pr_mom_dim volmom = ebcm_basics::getVolMomFromFacMom(facmom, a_facedir);
      double momzero    = facmom[pr_itm_i_dmo::Zero];

      pr_itm_r_dim xbar = a_locality->getXBar(a_face, a_dx);
      volmom.shift(xbar);

      unsigned int Pcols = pr_mom_dim::size();
      shared_ptr<eigen_mat> retval(new eigen_mat(1, Pcols));


      for(MomentIterator<DIM, ebcm_order> momit; momit.ok(); ++momit)
      {
        auto momind = momit();
        double matval = 0;
        if(momind[a_facedir] > 0)
        {
          int pd = momind[a_facedir];
          pr_itm_i_dim derind = momind - pr_itm_i_dim::BASISV(a_facedir);
          Real derval = volmom[derind];
          if(std::abs(momzero) > 1.0e-10)
          {
            matval = pd*derval/momzero;
          }
        }
        unsigned int currentCol = pr_mom_dim::indexOf(momind);
        (*retval)(0,currentCol) = matval;
      }
    
      return retval;
    }
    ///Returns a local stencil for \int (\grad \phi) \cdot (\ehat^facedir)
    /**
       Get stencil for integrated gradient in a  coordinate direction. 
       This is the same as the average gradient * face area.
       \int (\grad \phi) \cdot (\ehat^d) dA = A_f < \frac{\partial \phi}{\partial x_d} >
    */
    static local_stencil_t
    getStencilFromMat(const eigen_mat                                & a_matrix,
                      shared_ptr< EBCM::neighborhood< ebcm_order > > & a_locality)
    {
      int nrows_volumes   = a_locality->m_weig_volumes.size();
      local_stencil_t retval;
      const auto& volvec = *a_locality->m_volumes;
      for( int irow = 0; irow < nrows_volumes; irow++)
      {
        pair<pr_pt, double> entry;
        const ebcm_volu& volume = volvec[irow];
        entry.first = volume.m_pt;
        entry.second = a_matrix(irow, 0);
        retval.add(entry);
      }
      return retval;
    }

    ///
    static shared_ptr<eigen_mat>
    getWeightMatrix(shared_ptr< EBCM::neighborhood< ebcm_order > >  a_locality,
                    bool a_print = false)
    {
      shared_ptr<eigen_mat> retval;
      ch_mayday::Error("not implemented");
      return retval;
    }
    ///
    template<CENTERING cent>
    static local_stencil_t
    getHelmholtzFDotEhatStencil(const EBCM_Face<cent, ebcm_order>  & a_face,
                                const int                          & a_stenrad,
                                const int                          & a_weightpower,
                                const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                const string                       & a_dom_bc,
                                const string                       & a_eb_bc,
                                const int                          & a_facedir,
                                const double                       & a_dx,
                                bool a_print = false)
    {
      /**
         From the document (Lawrence Berkeley National Laboratory Technical Report LBNL-2001537):  
         Q : Q matrix (Unfortunately, I did not give it a fancy mathy name in the document.
               It is the one that converts  Taylor coefficients to fluxes).   
         S : stencil matrix.  
         P : phi vector.  
         W : weighting matrix.
         M : moment matrix.
         A: system matrix =  ((WM)^T(WM))
         (S P) = <\nabla \psi> = Q C= (Q ((WM)^T(WM))^{-1} (WM)^T  P) Q A^{-1} WM^T P.
         S = Q C = (Q ((WM)^T(WM))^{-1} (WM)^T  ) = Q A^{-1} (WM)^T.
      **/
      shared_ptr< EBCM::neighborhood< ebcm_order > > locality;
      shared_ptr< eigen_mat > Amat_p =
        EBCM::Algorithm_Framework<ebcm_order>::
        getSystemMatrix(locality, a_face, a_graph,
                        a_stenrad, a_weightpower,
                        a_dom_bc, a_eb_bc, a_print);
      shared_ptr<eigen_mat> Wmat_p = algo_framework::getWeightMatrix(locality, a_print);
      shared_ptr<eigen_mat> Mmat_p = algo_framework::getMomentMatrix(locality, a_dom_bc, a_eb_bc,a_print);

      eigen_mat WMmat;
      multiply(WMmat, *Wmat_p, *Mmat_p);
      eigen_mat WMTmat = WMmat;
      WMTmat.transpose();

      eigen_mat Ainvmat = *Amat_p;
      Ainvmat.invert();
      
      shared_ptr<eigen_mat> Qmat_p = getHelmholtzQMatrix(a_face, locality, a_facedir, a_dx, a_print);
      eigen_mat AinvWMTmat;
      multiply( AinvWMTmat, Ainvmat, WMTmat);
      eigen_mat Smat;
      multiply( Smat, *Qmat_p, AinvWMTmat);
      local_stencil_t retval = getStencilFromMat(Smat, locality);
      return retval;
    }

    ///
    static shared_ptr<eigen_mat>
    getHelmholtzEBQMatrix( const IrregFace< ebcm_order >      & a_face,
                           shared_ptr<ebcm_neighborhood>        a_locality,
                           const double                       & a_dx,
                           bool a_print)
    {
      
      shared_ptr<eigen_mat> retval;

      ch_mayday::Error("not implemented");
      return retval;
    }
    ///Returns a local stencil for \int (\grad \phi) \cdot (\ehat^facedir)
    /**
       Get stencil for integrated gradient in a  coordinate direction. 
       This is the same as the average gradient * face area.
       \int (\grad \phi) \cdot (\ehat^d) dA = A_f < \frac{\partial \phi}{\partial x_d} >
    */
    static local_stencil_t
    getIrregHelmholtzFDotNStencil(const IrregFace<ebcm_order> & a_face,
                                  const int                          & a_stenrad,
                                  const int                          & a_weightpower,
                                  const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                  const string                       & a_dom_bc,
                                  const string                       & a_eb_bc,
                                  const double                       & a_dx,
                                  bool a_print = false)
    {

      shared_ptr< EBCM::neighborhood< ebcm_order > > locality;
      shared_ptr< eigen_mat > Amat_p =
        EBCM::Algorithm_Framework<ebcm_order>::
        getEBSystemMatrix(locality, a_face, a_graph,
                          a_stenrad, a_weightpower,
                          a_dom_bc, a_eb_bc, a_print);

      shared_ptr<eigen_mat> Wmat_p = algo_framework::getWeightMatrix(locality, a_print);
      shared_ptr<eigen_mat> Mmat_p = algo_framework::getMomentMatrix(locality, a_dom_bc, a_eb_bc, a_print);

      eigen_mat WMmat;
      multiply(WMmat, *Wmat_p, *Mmat_p);
      eigen_mat WMTmat = WMmat;
      WMTmat.transpose();

      eigen_mat Ainvmat = *Amat_p;
      Ainvmat.invert();
      
      shared_ptr<eigen_mat> Qmat_p = getHelmholtzEBQMatrix(a_face, locality,  a_dx, a_print);
      eigen_mat AinvWMTmat;
      multiply( AinvWMTmat, Ainvmat, WMTmat);
      eigen_mat Smat;
      multiply( Smat, *Qmat_p, AinvWMTmat);
      local_stencil_t retval = getStencilFromMat(Smat, locality);
      return retval;
    }
    
    ///
    /**
       Class which forms a helmholtz operator.
       This creates an operator of the form
       L(phi) = alpha phi + beta * div(F(phi))
    **/
    class HelmholtzOp: public Base_PETSc_Op
    {
    public:

      typedef Proto::Side         pr_sd;
      typedef Proto::Box          pr_box;
      
      double                                 m_alpha;
      double                                 m_beta;
      int                                    m_weight_power;
      int                                    m_sten_rad;
      string m_dom_bc;
      string m_eb_bc;
      /// 
      /**
         alpha        = identity coeff
         beta         = coeff of laplcian
         ghost        = neighbor ghost size (grow the face box by this much)
         weight_power = goes to neighborhood to weight vofs.
         meta         = ebcm's very data intensive meta data
      **/
      HelmholtzOp(const double                                 & a_alpha,
                  const double                                 & a_beta,
                  const int                                    & a_sten_rad,
                  const int                                    & a_weight_power,
                  const shared_ptr< ebcm_meta  >               & a_meta_ptr,
                  string a_dom_bc, string a_eb_bc, 
                  bool                                           a_print = false)
      :Base_PETSc_Op(a_meta_ptr, a_sten_rad, a_print)
      {
        m_alpha        = a_alpha;       
        m_beta         = a_beta;
        m_sten_rad     = a_sten_rad;
        m_weight_power = a_weight_power;
        m_dom_bc       = a_dom_bc;
        m_eb_bc        = a_eb_bc;
      }

      virtual ~HelmholtzOp()
      {
      }

      ///this is the function that the driver cod calls
      virtual shared_ptr<local_stencil_t> getLocalStencil(const pr_pt& a_pt, const ebcm_graph& a_graph, bool a_print = false)
      {
        shared_ptr<local_stencil_t> retval =
          privateStencil(a_pt,          
                         m_sten_rad,     
                         m_weight_power, 
                         a_graph,       
                         m_dom_bc,      
                         m_eb_bc,       
                         a_print);        
        return retval;
      }

      void applyOp(shared_ptr<ebcm_dou_sca_data> a_lph,
                   shared_ptr<ebcm_dou_sca_data> a_phi) const
      {
        ch_mayday::Error("not implemented");
      }
    private:
      ///this is where the canal locks are dug in the basement
      virtual shared_ptr<local_stencil_t> privateStencil(const pr_pt                         & a_pt, 
                                                         const int                          & a_stenrad,
                                                         const int                          & a_weightpower,
                                                         const EBCM::EBCM_Graph<ebcm_order> & a_graph,
                                                         const string                       & a_dom_bc,
                                                         const string                       & a_eb_bc,
                                                         bool                                 a_print        = false)
      {
        typedef Proto::Box pr_bx;
        
        shared_ptr<local_stencil_t> retval(new local_stencil_t());
        pr_bx pointbx(a_pt, a_pt);
        pr_bx grbx = pointbx.grow(pr_pt::Ones(this->m_ghost));
        
        auto volume = a_graph.getVolumeCoveringPoint(a_pt);
        auto dx = volume.m_dx;
        {
          for(int iface = 0; iface < volume.m_xfaces.size(); iface++)
          {
            const auto& face = volume.m_xfaces[iface];
            double fluxsign = faceSign(volume, face);
            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,       
                                                               dx, 0, a_print);        
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_yfaces.size(); iface++)
          {
            const auto& face = volume.m_yfaces[iface];
            double fluxsign = faceSign(volume, face);
            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,       
                                                               dx, 1, a_print);        

                                                               
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_zfaces.size(); iface++)
          {
            const auto& face = volume.m_zfaces[iface];
            double fluxsign = faceSign(volume, face);
            local_stencil_t flux = getHelmholtzFDotEhatStencil(face,
                                                               a_stenrad,     
                                                               a_weightpower, 
                                                               a_graph,       
                                                               a_dom_bc,      
                                                               a_eb_bc,       
                                                               dx, 2, a_print);        
                                                               
            flux *= fluxsign;
            (*retval) += flux;
          }
        }
        {
          for(int iface = 0; iface < volume.m_ifaces.size(); iface++)
          {
            const auto& face = volume.m_ifaces[iface];
            local_stencil_t flux = getIrregHelmholtzFDotNStencil(face,
                                                                 a_stenrad,     
                                                                 a_weightpower, 
                                                                 a_graph,       
                                                                 a_dom_bc,      
                                                                 a_eb_bc,
                                                                 dx, a_print);
                                                                 
            (*retval) += flux;
          }
        }
        
        return retval;
      }

    };
    ///
    static void
    setPhiExact(shared_ptr<ebcm_dou_sca_data>       a_phi_exac_p,
                const ch_dbl                 &      a_grids,
                const ch_probdom             &      a_domain,
                const double                 &      a_dx,
                shared_ptr< ebcm_meta  >     &      a_meta_data)
    {
      ch_mayday::Error("not implemented");
    }
    ///
    static void
    setLphExact(shared_ptr<ebcm_dou_sca_data>       a_lph_exac_p,
                const ch_dbl                 &      a_grids,
                const ch_probdom             &      a_domain,
                const double                 &      a_dx,
                shared_ptr< ebcm_meta  >     &      a_meta_data)
    {
      ch_mayday::Error("not implemented");
    }


    static void
    getHelmholtzError(int a_coarseningFactor)
    {
      
      shared_ptr< ebcm_meta  > meta_data;
      ch_dbl                   grids;
      double                   dx;
      ch_probdom               domain;
      ///stuff to parmparse.
      double alpha     = 4586;
      double beta      = 4586;
      int weight_power = 4586;
      int nghost       = 4586;
      string dom_bc;
      string eb_bc;
      ParmParse pp("trunc_test");
      pp.get("alpha", alpha);
      pp.get("beta" , beta);
      pp.get("nghost" , beta);
      pp.get("weight_power", weight_power);
      pp.get("dom_bc", dom_bc);
      pp.get("eb_bc",   eb_bc);
      algo_framework::makeMergedGeometry(meta_data, grids, dx, domain, a_coarseningFactor, true, true);

      
      HelmholtzOp helmholtzOperator(alpha, beta, nghost, weight_power, meta_data, dom_bc, eb_bc, false);
      shared_ptr<ebcm_dou_sca_data> phiExac(new ebcm_dou_sca_data(meta_data, nghost));
      shared_ptr<ebcm_dou_sca_data> lphExac(new ebcm_dou_sca_data(meta_data, nghost));
      shared_ptr<ebcm_dou_sca_data> lphCalc(new ebcm_dou_sca_data(meta_data, nghost));
      shared_ptr<ebcm_dou_sca_data> truncat(new ebcm_dou_sca_data(meta_data, nghost));
      setPhiExact(phiExac, grids, domain, dx, meta_data);
      setLphExact(lphExac, grids, domain, dx, meta_data);

      helmholtzOperator.applyOp(lphCalc, phiExac);
      truncat->setVal(0,);
      (*truncat) += (*lphExac);
      (*truncat) -= (*lphCalc);
      double retval = truncat->maxNorm();
      return retval;
    }
    //this is the function called by Chombo_4/example/EBCM/hoeb_truncation
    static  void run_hoeb_truncation_tests()
    {

      
      ///apply the helmholtz operator and compute Error
      double fineErrorNorm = getHelmholtzError(1);
      double coarErrorNorm = getHelmholtzError(2);
      double order = (coarErrorNorm/fineErrorNorm)/std::log(2);
      Chombo4::pout() << "Helmholtz truncation error order = " << order << endl;
    }

  }; //end petsc_framwork class
}
#endif

#endif 
