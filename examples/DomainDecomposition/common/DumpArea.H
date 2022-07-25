#ifndef __EBDumpARea__
#define __EBDumpARea__
#define PI 3.141592653589793
#include "Chombo_EBLevelFluxData.H"
#include "Chombo_NamespaceHeader.H"

///will not work on device--debugging use only
class DumpArea
{
public:
  template <int order>
  static void dumpAsOneBox( const EBLevelBoxData<CELL, 1>*    a_input,
                            shared_ptr<GeometryService<order> >   a_geoserv)
  {
    shared_ptr<EBLevelBoxData<CELL, 1> > oneboxCopy=
      EBLevelBoxData<CELL, 1>::getOneBoxCopyOfData(*a_input, a_geoserv);
    DumpArea::dumpLDCell1Area(&(*oneboxCopy));
  }
  
  static void 
  genDumpCell1(const EBBoxData<CELL, Real, 1>* dataPtr, Point lodeb, Point hideb)
  {
    if(dataPtr != NULL)
    {
      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      EBGraph graph = dataPtr->ebgraph();
      const EBBoxData<CELL, Real, 1> & data = *dataPtr;
      Bx databox = dataPtr->box();
      databox &= graph.getDomain();
      if((databox.contains(lodeb)) || (databox.contains(hideb)))
      {
        cout << "data region contains:" << endl;
        for(int j = hideb[1]; j >= lodeb[1]; j--)
        {
          for(int i = lodeb[0]; i <= hideb[0]; i++)
          {
            Point pt(i,j);
            if(databox.contains(pt))
            {
              vector<EBIndex<CELL> > vofs = graph.getVoFs(pt);
              cout << pt << ":";
              if(vofs.size() >  0)
              {
                if(data(vofs[0], 0) >= 0)
                {
                  cout << " ";
                }
                cout  << data(vofs[0], 0) << "  ";
              }
              else
              {
                cout  << "---------------" << "  ";
              }
            }
          }
          cout << endl;
        }
      }
    }
  }
  static void 
  dumpLDCell1Area(const EBLevelBoxData<CELL, 1>* dataPtr)
  {
    if(dataPtr != NULL)
    {
      const EBLevelBoxData<CELL, 1>& data = *dataPtr;
      DataIterator dit = data.dataIterator();
      int nx = dataPtr->disjointBoxLayout().physDomain().domainBox().size(0);
      int ny = dataPtr->disjointBoxLayout().physDomain().domainBox().size(1);
      int nxfrac = (3*nx)/32;
      int nyfrac = (8*ny)/32;
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        IntVect ivlo(D_DECL(nxfrac-3, nyfrac-5, 0));
        IntVect ivhi(D_DECL(nxfrac+3, nyfrac+5, 0));
//        ivlo = IntVect::Zero;
//        ivhi = 9*IntVect::Unit;
//        ivhi[0] = 6;
        Box area(ivlo, ivhi);

        
        Box valid = data.disjointBoxLayout()[dit[ibox]];
        valid &= area;
        if(!valid.isEmpty())
        {
          Point lovalid = ProtoCh::getPoint(valid.smallEnd());
          Point hivalid = ProtoCh::getPoint(valid.bigEnd());
          cout << "ibox = " << ibox << ":" << endl;
          genDumpCell1(&data[dit[ibox]], lovalid, hivalid);
        }
      }
    }
  }

};

#include "Chombo_NamespaceFooter.H"

#endif
