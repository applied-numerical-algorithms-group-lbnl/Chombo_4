
#ifndef _Proto_EBStencilImplem_H_
#define _Proto_EBStencilImplem_H_
#include "Proto_EBHostData.H"
namespace Proto
{
  ///
  /**
     regular stencil gets applied in the interior of the domain.
     The aggstencil overwrites it near the EB.
     The domainBC stencils fill ghost cells outside the domain.  
     The vector length needs to be the same as the one-sided span of the stencil.
  */
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  EBStencil<order, data_t, srcCenter, dstCenter>::
  EBStencil(string              a_stencilName,
            string              a_domainBCName[2*DIM],
            string              a_ebbcName,
            const Box         & a_srcValid,
            const Box         & a_dstValid,
            const Box         & a_srcDomain,
            const Box         & a_dstDomain,
            const Point       & a_srcGhost,
            const Point       & a_dstGhost,
            const EBGraph     & a_srcGraph,
            const EBGraph     & a_dstGraph,
            const VoluData    & a_voludata,
            const EBFaData    & a_ebfadata,
            const XFacData    & a_xfacdata,
            const YFacData    & a_yfacdata,
            const ZFacData    & a_zfacdata,
            double              a_dx,
            bool                a_periodiconly,
            bool                a_needDiagonalWeights,
            Point               a_dstGrown)
  {
    m_periodic_domain_bcs_only = a_periodiconly;
    vector< EBIndex<dstCenter> >              dstVoFs;
    vector< LocalStencil<srcCenter, data_t> > stencil;

    m_hasDomainBCs = 
      archive_t::getStencil(dstVoFs, stencil, m_regStencilInterior, 
                            m_regStencilBCS, m_BCApplyBoxes, m_regApplyBox,
                            m_bcsOnly, m_irregOnly,
                            a_stencilName,  a_domainBCName, a_ebbcName,
                            a_srcValid,  a_dstValid,
                            a_srcDomain, a_dstDomain,     
                            a_srcGhost,  a_dstGhost,    
                            a_srcGraph,  a_dstGraph,
                            a_voludata,  a_ebfadata,    
                            a_xfacdata,  a_yfacdata, a_zfacdata,    
                            a_dx, a_periodiconly, a_dstGrown);

    if(m_hasDomainBCs) defineApplyITZ();

    if(!m_bcsOnly)
    {
      m_irregStencil = shared_ptr<AggStencil<srcCenter, dstCenter, data_t> >
        (new AggStencil<srcCenter, dstCenter, data_t>(dstVoFs, stencil, a_srcGraph, a_dstGraph,
                                                      a_srcValid, a_dstValid, a_srcGhost, a_dstGhost));
      if(a_needDiagonalWeights)
      {
        computeDiagonalWeights(dstVoFs, stencil, a_dstGraph, a_dstValid, a_dstGhost);
      }
    }
  }

  //Constructor specifically for Weighted Least Squares stencils
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  EBStencil<order, data_t, srcCenter, dstCenter>::
  EBStencil(string              a_stencilName,
            string              a_domainBCName[2*DIM],
            string              a_ebbcName,
            const Box         & a_srcValid,
            const Box         & a_dstValid,
            const Box         & a_srcDomain,
            const Box         & a_dstDomain,
            const Point       & a_srcGhost,
            const Point       & a_dstGhost,
            const EBGraph     & a_srcGraph,
            const EBGraph     & a_dstGraph,
            const VoluData    & a_voludata,
            const EBFaData    & a_ebfadata,
            const XFacData    & a_xfacdata,
            const YFacData    & a_yfacdata,
            const ZFacData    & a_zfacdata,
	    const EBNorData   & a_ebnorxdata,
            const EBNorData   & a_ebnorydata,
            const EBNorData   & a_ebnorzdata,
            double              a_dx,
            bool                a_periodiconly,
            bool                a_needDiagonalWeights,
            Point               a_dstGrown)
  {
    m_periodic_domain_bcs_only = a_periodiconly;
    vector< EBIndex<dstCenter> >              dstVoFs;
    vector< LocalStencil<srcCenter, data_t> > stencil;

    m_hasDomainBCs = 
      archive_t::getStencil(dstVoFs, stencil, m_regStencilInterior, 
                            m_regStencilBCS, m_BCApplyBoxes, m_regApplyBox, m_bcsOnly,
                            a_stencilName,  a_domainBCName, a_ebbcName,
                            a_srcValid,  a_dstValid,
                            a_srcDomain, a_dstDomain,     
                            a_srcGhost,  a_dstGhost,    
                            a_srcGraph,  a_dstGraph,
                            a_voludata,  a_ebfadata,    
                            a_xfacdata,  a_yfacdata, a_zfacdata,
			    a_ebnorxdata, a_ebnorydata,
                            a_ebnorzdata,
                            a_dx, a_periodiconly, a_dstGrown);

    if(m_hasDomainBCs) defineApplyITZ();

    if(!m_bcsOnly)
    {
      m_irregStencil = shared_ptr<AggStencil<srcCenter, dstCenter, data_t> >
        (new AggStencil<srcCenter, dstCenter, data_t>(dstVoFs, stencil, a_srcGraph, a_dstGraph,
                                                      a_srcValid, a_dstValid, a_srcGhost, a_dstGhost));
      if(a_needDiagonalWeights)
      {
        computeDiagonalWeights(dstVoFs, stencil, a_dstGraph, a_dstValid, a_dstGhost);
      }
    }
  }
  ///
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  void 
  EBStencil<order, data_t, srcCenter, dstCenter>::
  computeDiagonalWeights(const vector< EBIndex<dstCenter> >              & a_dstVoFs,
                         const vector< LocalStencil<srcCenter, data_t> > & a_stencil,
                         const EBGraph                                   & a_graph,
                         const Box                                       & a_validBox,
                         const Point                                     & a_dstGrown)
  {
    //not sure this makes any sense if the centerings are different
    if(srcCenter == dstCenter)
    {
      Box grownBox = a_validBox.grow(a_dstGrown);

      m_diagonalWeights = shared_ptr< EBBoxData<dstCenter, data_t, 1> >
        (new EBBoxData<dstCenter, data_t, 1>(grownBox, a_graph));
    
      //have to compute stuff on the host first
      EBHostData<dstCenter, data_t, 1> hostDiagonal(grownBox, a_graph);
      //first set the regular data to the value of the regular stencil's diagonal.
      data_t regDiagonal = m_regStencilInterior.diagonalValue();
      hostDiagonal.setVal(regDiagonal);
    
      //now fix up values on irregular cells
      for(int idst = 0; idst < a_dstVoFs.size(); idst++)
      {
        const EBIndex<dstCenter>& dstvof               = a_dstVoFs[idst];
        const LocalStencil<srcCenter, data_t>& irrsten = a_stencil[idst];
        if(!irrsten.m_betaWeightSet)
        {
          PROTO_ASSERT(false, "EBStencil asking for diagonal weight that has not been set");
        }
        data_t diagval = irrsten.m_betaWeight;
        hostDiagonal(dstvof, 0) = diagval;
      }

      //now copy the diagonal value to the device
      EBHostData<dstCenter, data_t, 1>& hostfab = hostDiagonal;
      EBBoxData< dstCenter, data_t, 1>& devifab = *m_diagonalWeights;
      for(unsigned int iwhichfab = 0; iwhichfab < 2; iwhichfab++)
      {
        size_t   fabsizeHost, fabsizeDevi; 
        data_t*  dataptrHost, *dataptrDevi;
        hostfab.dataInfo(dataptrHost, fabsizeHost, iwhichfab);
        devifab.dataInfo(dataptrDevi, fabsizeDevi, iwhichfab);

        if(fabsizeDevi != fabsizeHost)
        {
          PROTO_ASSERT(false, "host and device sizes need to match here");
        }
        size_t datSize = sizeof(data_t)*fabsizeDevi;
        protoMemcpy(dataptrDevi, dataptrHost, datSize, protoMemcpyHostToDevice);
      }
    }
  }

  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  void
  EBStencil<order, data_t, srcCenter, dstCenter>::
  defineApplyITZ() 
  {
#ifdef PROTO_CUDA
/*      auto& vecStencil = sFused.getStencils();  
      for(int ibc = 0; ibc < 2*DIM; ibc++)
        for(int ivec = 0; ivec < m_regStencilBCS[ibc].size(); ivec++)
          vecStencil.push_back(&m_regStencilBCS[ibc][ivec]);
   */   
      sFused.define(m_regStencilBCS, m_BCApplyBoxes, 2*DIM);
#endif
  }

  ///
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  void 
  EBStencil<order, data_t, srcCenter, dstCenter>::
  applyITZ(EBBoxData<       dstCenter, data_t, 1>  &   a_dst,
           const EBBoxData< srcCenter, data_t, 1>  &   a_src,
           double                                      a_scale)  const
  {
    PR_TIME("EBStencil::applyITZ");
    
    //first set the boundary conditions at the domain box
    //this fills the ghost cells when the box abuts the domain.
    if(m_hasDomainBCs && (!m_irregOnly))
    {
      PR_TIME("bcs_apply");
#ifdef PROTO_CUDA
      const BoxData<data_t, 1>& castSrc = const_cast<BoxData<data_t, 1> &>(a_src.getRegData());
      sFused.cudaApplyFused(castSrc, (BoxData<data_t, 1> &) (a_src.getRegData()), true, 1.0);
#else
      BoxData<data_t, 1>& castSrc = const_cast<BoxData<data_t, 1> &>(a_src.getRegData());
      for(int ibc = 0; ibc < 2*DIM; ibc++)
      {
        for(int ivec = 0; ivec < m_regStencilBCS[ibc].size(); ivec++)
        {
          Box bcsrcbox = m_BCApplyBoxes[ibc][ivec];
          m_regStencilBCS[ibc][ivec].apply(castSrc, castSrc,
                                           bcsrcbox, true, 1.0);
        }
      }
#endif
    }
    if(!m_bcsOnly)
    {
      PR_TIME("actual_apply");
      bool incrementOnly = false;
      //already set stuff to zero if necessary
      bool initToZero    = true;
      if(!m_irregOnly)
      {
        m_regStencilInterior.apply(a_src.getRegData(), 
                                   a_dst.getRegData(),
                                   m_regApplyBox, initToZero, a_scale);
      }

      m_irregStencil->apply(a_dst, a_src, a_scale, incrementOnly);
    }
    
  }


  ///
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  void 
  EBStencil<order, data_t, srcCenter, dstCenter>::
  apply(EBBoxData<       dstCenter, data_t, 1>  &  a_dst,
        const EBBoxData< srcCenter, data_t, 1>  &  a_src,
        bool                                       a_initToZero,
        double                                     a_scale)  const
  {
    PR_TIME("EBStencil::apply");
    if(a_initToZero)
    {
      applyITZ(a_dst, a_src, a_scale);
    }
    else
    {
      EBBoxData< dstCenter, data_t, 1>   incr(a_dst.inputBox(), a_dst.ebgraph());
      incr.setVal(0.);
      applyITZ(incr, a_src, a_scale);
      a_dst += incr;
    }
  }
              
}

#endif
