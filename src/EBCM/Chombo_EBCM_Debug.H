#ifndef __EBCM_DBGFunctionss__
#define __EBCM_DBGFunctionss__

#include "Chombo_EBCM_HostLevelData.H"
#include "Chombo_EBCM_FluxLevelData.H"
///debugging use only.   


namespace EBCM_Debug
{
  ///
  /**
     class to print a window of data in index xpace to screen 
     Only static data members.
   **/
  template<typename data_t, int  ncomp, int order>
  class DataWindow
  {
  public:
    typedef  Proto::Point                                          pr_pt;
    typedef  Proto::Box                                            pr_box;
    typedef  Proto::HostBoxData< data_t, ncomp>                    pr_fab;
    typedef  EBCM::HostLevelData<data_t, ncomp, order>             ebcm_level_data;
    typedef  EBCM::FluxLevelData<data_t, ncomp, order>             ebcm_level_flux;
    typedef  EBCM::EBCM_Graph<order>                               ebcm_graph;
    typedef  Chombo4::DisjointBoxLayout                            ch_dbl;
    typedef  Chombo4::Box                                          ch_box;
    static void writeWindowLevel(const ebcm_level_data* a_data_ptr, pr_pt a_lo, pr_pt a_hi, string a_prefix)
    {
      using std::cout;
      using std::endl;
      using Chombo4::pout;
      pout() << "writeWindow: prefix = " << a_prefix << endl;
      if(a_data_ptr != NULL)
      {
        pout() << "writeWindowLevel:  window from " << a_lo << " to " << a_hi << endl;
        pout()    << setprecision(8)
                << setiosflags(ios::showpoint)
                << setiosflags(ios::scientific);
        pr_box window(a_lo, a_hi);
        const auto& graphs = *(a_data_ptr->m_meta->m_graphs);
        const auto& grids  =  (a_data_ptr->m_meta->m_grids);
        auto dit = grids.dataIterator();
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const auto& valid =           grids[dit[ibox]];
          pr_box intersect = window & valid;
          if(!intersect.empty())
          {
            const auto& graph =        graphs[dit[ibox]];
            const auto&  data = (*a_data_ptr)[dit[ibox]];
            pr_pt lo = intersect.low();
            pr_pt hi = intersect.high();
            
            pout() << "writeWindowLevel:  intersection for ibox  "<< ibox  << " goes from " << a_lo << " to " << a_hi << endl;
            writeWindowBox(&data, intersect, graph);
          } //end (has intersection on this box)
        } //end loop over boxes
      } //end if data_ptr real
    } //end funtion writeWindowLevel

    static void writeFluxWindowLevel(const ebcm_level_flux* a_data_ptr, pr_pt a_lo, pr_pt a_hi, int a_idir, string a_prefix)
    {
      using std::cout;
      using std::endl;
      using Chombo4::pout;
      
      pout() << "writeWindow: prefix = " << a_prefix << endl;
      if(a_data_ptr != NULL)
      {
        pout() << "writeWindowLevel:  window from " << a_lo << " to " << a_hi << endl;
        pout()    << setprecision(8)
                << setiosflags(ios::showpoint)
                << setiosflags(ios::scientific);
        pr_box window(a_lo, a_hi);
        pr_box window_face = window;
        window_face.m_high[a_idir] += 1;
        window_face.grow(a_idir, 1);
        const auto& graphs   = *(a_data_ptr->m_meta->m_graphs);
        const auto& grids    =  (a_data_ptr->m_meta->m_grids);
        const auto& flux_lev =  (a_data_ptr->m_data)[a_idir];
        auto dit = grids.dataIterator();
        for(int ibox =0; ibox < dit.size(); ibox++)
        {
          const ch_box& valid = grids[dit[ibox]];
          pr_box valid_face = getProtoBox(valid);
          valid_face.m_high[a_idir] += 1;

          pr_box intersect = window_face & valid_face;
          if(!intersect.empty())
          {
            const auto& graph =        graphs[dit[ibox]];
            const auto&  data =    (*flux_lev)[dit[ibox]];
            pr_pt lo = intersect.low();
            pr_pt hi = intersect.high();
            
            pout() << "writeWindowLevel: FLUX_DIR = "<< a_idir<< "\n intersection for ibox  "<< ibox  << " goes from " << a_lo << " to " << a_hi << endl;
            writeWindowBox(&data, intersect, graph);
          } //end (has intersection on this box)
        } //end loop over boxes
      } //end if data_ptr real
    } //end funtion writeWindowLevel
    
    static void writeWindowBox(const pr_fab    * a_data_ptr,
                               const pr_box    & a_window,
                               const ebcm_graph& a_graph)
    {
      using std::cout;
      using std::endl;
      using Chombo4::pout;
      if(a_data_ptr != NULL)
      {
        pr_box intersect = a_data_ptr->box();
        intersect &= a_window;
        const auto&  data =   (*a_data_ptr);
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          pout()<< "writeWindowBox: intersection = " << intersect << ", comp = " << icomp << endl;
          pr_pt lo = intersect.low();
          pr_pt hi = intersect.high();
              
#if DIM==3
          for(int z_co = hi[2]; co >= lo[2]; z_co--)
          {
#endif
            for(int y_co = hi[1]; y_co >= lo[1]; y_co--)
            {
              for(int x_co = lo[0]; x_co <= hi[0]; x_co++)
              {
                pr_pt cur_pt(D_DECL(x_co, y_co, z_co));
                pout() << "(";
                for(int idir = 0; idir <DIM; idir++)
                {
                  int ipt = cur_pt[idir];
                  //points should be positive and <= 100 for now
                  if(ipt < 10)
                    pout() << "  ";
                  else if (ipt < 100)
                    pout() << " ";
                  pout() << cur_pt[idir];
                  if(idir == (DIM-1))
                    pout() << "):";
                  else
                    pout() << ",";
                }
                if(a_graph.isCovered(cur_pt))
                {
                  pout()  <<   "---------------" << "  ";
                }
                else 
                {
                  const double& data_val = data(cur_pt, icomp);
                  if(std::isnan(data_val))
                  {
                    pout()  << "=-=-=NaN=-=-=-=" << "  ";
                  }
                  else
                  {
                    if(data_val >= 0)
                    {
                      //keeps positive and negative lined up 
                      pout() << " ";
                    } 
                    pout()  << data_val << "  ";
                  }
                } //end if data is not covered
              } //end loop over x
              pout()  << endl;
            } //end loop over y
#if DIM==3
          } //end loop over z
#endif

        } //end loop over components
      }//end if data is real
    }//end function writeWindowBox
  }; //end class DataWindow
} //end namespace EBCM_Debug

//include guard
#endif   
