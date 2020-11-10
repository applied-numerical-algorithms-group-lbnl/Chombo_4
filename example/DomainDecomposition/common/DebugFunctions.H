#ifndef __DBGFunctionss__
#define __DBGFunctionss__
#define PI 3.141592653589793
#include "Chombo_EBLevelFluxData.H"
#include "Chombo_NamespaceHeader.H"
///will not work on device--debugging use only


using Proto::Var;
void 
dumpVecRV(Vector<RealVect>* a_dataPtr)
{
  if(a_dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    for(int ivec = 0; ivec < a_dataPtr->size(); ivec++)
    {
      
      cout << ivec << ":" << (*a_dataPtr)[ivec];
      cout << endl;
    }
  }
}

void
dumpPS(PointSet* dataPtr)
{
  if(dataPtr != NULL)
  {
    PointSet irregset = *dataPtr;
    int irow = 0;
    for(PointSetIterator ivsit(irregset); ivsit.ok(); ++ivsit)
    {
      cout << ivsit() << " ";
      irow++;
      if(irow > 10) 
      {
        cout << endl;
        irow = 0;
      }
    }
    cout << endl;
  }
}

#if DIM==2
void 
dumpBDDIM(BoxData<Real, DIM>* dataPtr)
{
  if(dataPtr != NULL)
  {
    Point lodeb = dataPtr->box().low();
    Point hideb = dataPtr->box().high();
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    BoxData<Real, DIM> & data = *dataPtr;
    for(int ivar = 0; ivar <DIM; ivar++)
    {
      
      cout << "data region contains for variable " << ivar << endl;
      for(int j = hideb[1]; j >= lodeb[1]; j--)
      {
        for(int i = lodeb[0]; i <= hideb[0]; i++)
        {
          Point pt(i,j);
          cout << pt << ":";
          if(data(pt, ivar) >= 0)
          {
            cout << " ";
          }
          cout  << data(pt, ivar) << "  ";
        }
        cout << endl;
      }
    }
  }
}

void 
dumpBD1(BoxData< Real, 1>* dataPtr)
  
{
if(dataPtr != NULL)
  {
    Point lodeb = dataPtr->box().low();
    Point hideb = dataPtr->box().high();
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    BoxData< Real, 1> & data = *dataPtr;

    cout << "data region contains:" << endl;
    for(int j = hideb[1]; j >= lodeb[1]; j--)
    {
      for(int i = lodeb[0]; i <= hideb[0]; i++)
      {
        Point pt(i,j);
        cout << pt << ":";
        if(data(pt, 0) >= 0)
        {
          cout << " ";
        }
        cout  << data(pt, 0) << "  ";
      }
      cout << endl;
    }
  }
}

void 
genDumpCellDIM(EBBoxData<CELL, Real, DIM>* dataPtr, Point lodeb, Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBGraph graph = dataPtr->ebgraph();
    EBBoxData<CELL, Real, DIM> & data = *dataPtr;
    Bx databox = dataPtr->box();
    databox &= graph.getDomain();
    if((databox.contains(lodeb)) || (databox.contains(hideb)))
    {
      for(int ivar = 0; ivar <DIM; ivar++)
      {
      
        cout << "data region contains for variable " << ivar << endl;
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
                if(data(vofs[0], ivar) >= 0)
                {
                  cout << " ";
                }
                cout  << data(vofs[0], ivar) << "  ";
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
}

void 
genDumpCell1(EBBoxData<CELL, Real, 1>* dataPtr, Point lodeb, Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBGraph graph = dataPtr->ebgraph();
    EBBoxData<CELL, Real, 1> & data = *dataPtr;
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

void 
genDumpHost1(EBHostData<CELL, Real, 1>* dataPtr, Point lodeb, Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBGraph graph = dataPtr->ebgraph();
    EBHostData<CELL, Real, 1> & data = *dataPtr;
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


void 
fullDumpCell1(EBBoxData<CELL, Real, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBGraph graph = dataPtr->ebgraph();

    Bx databox = dataPtr->box();
    databox &= graph.getDomain();
    Point lodeb = databox.low();
    Point hideb = databox.high();
    genDumpCell1(dataPtr, lodeb, hideb);
  }
}


void 
fullDumpCellDIM(EBBoxData<CELL, Real, DIM>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBGraph graph = dataPtr->ebgraph();
    Bx databox = dataPtr->box();
    databox &= graph.getDomain();
    Point lodeb = databox.low();
    Point hideb = databox.high();
    genDumpCellDIM(dataPtr, lodeb, hideb);
  }
}
void 
dumpLD1(EBLevelBoxData<CELL, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    EBLevelBoxData<CELL, 1>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      fullDumpCell1(&data[dit[ibox]]);
    }
  }
}
void 
dumpLDDIM(EBLevelBoxData<CELL, DIM>* dataPtr)
{
  if(dataPtr != NULL)
  {
    EBLevelBoxData<CELL, DIM>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      fullDumpCellDIM(&data[dit[ibox]]);
    }
  }
}



void 
dumpEB1(EBBoxData<CELL, Real, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);



    EBGraph graph = dataPtr->ebgraph();

    Bx databox = dataPtr->box();
    databox &= graph.getDomain();
    Point lodeb = databox.low();
    Point hideb = databox.high();
    genDumpCell1(dataPtr, lodeb, hideb);
  }
}

void 
dumpHost1(EBHostData<CELL, Real, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);



    EBGraph graph = dataPtr->ebgraph();

    Bx databox = dataPtr->box();
    databox &= graph.getDomain();
    Point lodeb = databox.low();
    Point hideb = databox.high();
    genDumpHost1(dataPtr, lodeb, hideb);
  }
}




void 
dumpEBDIM(EBBoxData<CELL, Real, DIM>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Point lodeb = dataPtr->box().low();
    Point hideb = dataPtr->box().high();
    genDumpCellDIM(dataPtr, lodeb, hideb);
  }
}


void 
genDumpXFace1(shared_ptr<EBBoxData<XFACE, Real, 1> >* dataPtr,
              Point lodeb,
              Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBBoxData<XFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(0, -1);


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
            bool printLines = false;
            if(vofs.size() >  0)
            {
              vector< EBIndex<XFACE> >faces = graph.getXFaces(vofs[0], Side::Lo);
              if(faces.size() == 0)
              {
                printLines = true;
              }
              else
              {
                if (data(faces[0], 0) >= 0) cout << " ";
                cout  << data(faces[0], 0) << "  ";
              }
            }
            else
            {
              printLines = true;
            }
            if(printLines)
            {
              cout << "---------------" << "  ";
            }
          }
        }
        cout << endl;
      }
    }
  }
}
void 
dumpXFace(shared_ptr<EBBoxData<XFACE, Real, 1> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBBoxData<XFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(0, -1);

    Bx area = databox & graph.getDomain();
    Point lodeb = area.low();
    Point hideb = area.high();
    genDumpXFace1(dataPtr, lodeb, hideb);
  }
}


void 
genDumpXFaceDIM(shared_ptr<EBBoxData<XFACE, Real, DIM> >* dataPtr, 
                Point lodeb,
                Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBBoxData<XFACE, Real, DIM> & data = *(*dataPtr);
    for(int idir = 0 ; idir < DIM; idir++)
    {
      Bx databox = data.box();
      if((databox.contains(lodeb)) || (databox.contains(hideb)))
      {
 
        pout() << "variable " << idir << endl;
        EBGraph graph = data.ebgraph();
        Bx databox = data.box();
        databox.grow(0, -1);

//        Bx area = databox & graph.getDomain();

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
              bool printLines = false;
              if(vofs.size() >  0)
              {
                vector< EBIndex<XFACE> >faces = graph.getXFaces(vofs[0], Side::Lo);
                if(faces.size() == 0)
                {
                  printLines = true;
                }
                else
                {
                  if (data(faces[0], idir) >= 0) cout << " ";
                  cout  << data(faces[0], idir) << "  ";
                }
              }
              else
              {
                printLines = true;
              }
              if(printLines)
              {
                cout << "---------------" << "  ";
              }
            }
          }
          cout << endl;
        }
      }
    }
  }
}

void 
dumpXFaceDIM(shared_ptr<EBBoxData<XFACE, Real, DIM> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);
    Bx databox = (*dataPtr)->box();
    databox.grow(0, -1);
    EBGraph graph = (*dataPtr)->ebgraph();

    Bx area = databox & graph.getDomain();

    Point lodeb = area.low();
    Point hideb = area.high();
    genDumpXFaceDIM(dataPtr, lodeb, hideb);
  }
}


void 
genDumpYFace1(shared_ptr<EBBoxData<YFACE, Real, 1> >* dataPtr,
              Point lodeb,
              Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBBoxData<YFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(0, -1);

    if((databox.contains(lodeb)) || (databox.contains(hideb)))
    {
//    Bx area = databox & graph.getDomain();

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
          bool printLines = false;
          if(vofs.size() >  0)
          {
            vector< EBIndex<YFACE> >faces = graph.getYFaces(vofs[0], Side::Lo);
            if(faces.size() == 0)
            {
              printLines = true;
            }
            else
            {
              if (data(faces[0], 0) >= 0) cout << " ";
              cout  << data(faces[0], 0) << "  ";
            }
          }
          else
          {
            printLines = true;
          }
          if(printLines)
          {
            cout << "---------------" << "  ";
          }
        }
      }
      cout << endl;
    }
  }
  }
}


void 
dumpYFace(shared_ptr<EBBoxData<YFACE, Real, 1> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBBoxData<YFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(0, -1);

    Bx area = databox & graph.getDomain();             // 
    Point lodeb = area.low();
    Point hideb = area.high();
    genDumpYFace1(dataPtr, lodeb, hideb);
  }
}


void 
genDumpYFaceDIM(shared_ptr<EBBoxData<YFACE, Real, DIM> >* dataPtr,
                Point lodeb, Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBBoxData<YFACE, Real, DIM> & data = *(*dataPtr);
    for(int idir = 0; idir < DIM; idir++)
    {
      pout() << "variable  = "  << idir << endl;
      EBGraph graph = data.ebgraph();
      Bx databox = data.box();
      databox.grow(0, -1);
      if((databox.contains(lodeb)) || (databox.contains(hideb)))
      {

//        Bx area = databox & graph.getDomain();

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
              bool printLines = false;
              if(vofs.size() >  0)
              {
                vector< EBIndex<YFACE> >faces = graph.getYFaces(vofs[0], Side::Lo);
                if(faces.size() == 0)
                {
                  printLines = true;
                }
                else
                {
                  if (data(faces[0], idir) >= 0) cout << " ";
                  cout  << data(faces[0], idir) << "  ";
                }
              }
              else
              {
                printLines = true;
              }
              if(printLines)
              {
                cout << "---------------" << "  ";
              }
            }
          }
          cout << endl;
        }
      }
    }
  }
}


void 
dumpYFaceDIM(shared_ptr<EBBoxData<YFACE, Real, DIM> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBBoxData<YFACE, Real, DIM> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(1, -1);

    Bx area = databox & graph.getDomain();
    Point lodeb = area.low();
    Point hideb = area.high();
    genDumpYFaceDIM(dataPtr, lodeb, hideb);
  }
}

void 
dumpFlux1(EBFluxData<Real, 1>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    pout() << "x face data:" << endl;
    dumpXFace(&(dataPtr->m_xflux));
    pout() << "y face data:" << endl;
    dumpYFace(&(dataPtr->m_yflux));
  }
}

void 
dumpFluxDIM(EBFluxData<Real, DIM>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    pout() << "x face data:" << endl;
    dumpXFaceDIM(&(dataPtr->m_xflux));
    pout() << "y face data:" << endl;
    dumpYFaceDIM(&(dataPtr->m_yflux));
  }
}

void 
dumpLDFlux1(EBLevelFluxData<1>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    EBLevelFluxData<1>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      
      pout() << "level flux data for ibox = " << ibox << ":" << endl;
      dumpFlux1(&(data[dit[ibox]]));
    }
  }
}
void 
dumpLDFluxDIM(EBLevelFluxData<DIM>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    EBLevelFluxData<DIM>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      
      pout() << "level flux data for ibox = " << ibox << ":" << endl;
      dumpFluxDIM(&(data[dit[ibox]]));
    }
  }
}

void dummyfunc()
{
  dumpLDDIM(NULL);
  dumpLDFluxDIM(NULL);
  dumpLD1(NULL);
  dumpLDFlux1(NULL);
  dumpXFace(NULL);
  dumpYFace(NULL);
  dumpFlux1(NULL);
  dumpEB1(NULL);
  dumpHost1(NULL);
  EBLevelBoxData<CELL,1>* phi1=nullptr;
  EBLevelBoxData<CELL,2>* phi2=nullptr;
  phi1->dumpFile(0,0.1);
  phi2->dumpFile(0,0.1);
}


Point lodebarea(3 , 7);
Point hidebarea(5,  11);                                    

void dumpEB1Area(EBBoxData<CELL, Real, 1>* dataPtr)
{
  genDumpCell1(dataPtr, lodebarea, hidebarea);
}
void dumpEBDIMArea(EBBoxData<CELL, Real, DIM>* dataPtr)
{
  genDumpCellDIM(dataPtr, lodebarea, hidebarea);
}

void 
dumpLD1Area(EBLevelBoxData<CELL, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    EBLevelBoxData<CELL, 1>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      genDumpCell1(&data[dit[ibox]], lodebarea, hidebarea);
    }
  }
}
void 
dumpLDDIMArea(EBLevelBoxData<CELL, DIM>* dataPtr)
{
  if(dataPtr != NULL)
  {
    EBLevelBoxData<CELL, DIM>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      genDumpCellDIM(&data[dit[ibox]], lodebarea, hidebarea);
    }
  }
}
void 
dumpXFaceDIMArea(shared_ptr<EBBoxData<XFACE, Real, DIM> >* dataPtr)
{
  genDumpXFaceDIM(dataPtr, lodebarea, hidebarea);
}
void 
dumpXFace1Area(shared_ptr<EBBoxData<XFACE, Real, 1> >* dataPtr)
{
  genDumpXFace1(dataPtr, lodebarea, hidebarea);
}
void 
dumpYFaceDIMArea(shared_ptr<EBBoxData<YFACE, Real, DIM> >* dataPtr)
{
  genDumpYFaceDIM(dataPtr, lodebarea, hidebarea);
}
void 
dumpYFace1Area(shared_ptr<EBBoxData<YFACE, Real, 1> >* dataPtr)
{
  genDumpYFace1(dataPtr, lodebarea, hidebarea);
}

void 
dumpFlux1Area(EBFluxData<Real, 1>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    pout() << "x face data:" << endl;
    dumpXFace1Area(&(dataPtr->m_xflux));
    pout() << "y face data:" << endl;
    dumpYFace1Area(&(dataPtr->m_yflux));
  }
}
void 
dumpFluxDIMArea(EBFluxData<Real, DIM>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    pout() << "x face data:" << endl;
    dumpXFaceDIMArea(&(dataPtr->m_xflux));
    pout() << "y face data:" << endl;
    dumpYFaceDIMArea(&(dataPtr->m_yflux));
  }
}

void 
dumpLDFluxDIMArea(EBLevelFluxData<DIM>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    EBLevelFluxData<DIM>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      pout() << "level flux data for ibox = " << ibox << ":" << endl;
      dumpFluxDIMArea(&(data[dit[ibox]]));
    }
  }
}
void 
dumpLDFlux1Area(EBLevelFluxData<1>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    EBLevelFluxData<1>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      pout() << "level flux data for ibox = " << ibox << ":" << endl;
      dumpFlux1Area(&(data[dit[ibox]]));
    }
  }
}

void dummyfunc_too()
{
  dumpLDDIMArea(NULL);
  dumpLDFluxDIMArea(NULL);
  dumpLD1Area(NULL);
  dumpLDFlux1Area(NULL);
  dumpXFace1Area(NULL);
  dumpYFace1Area(NULL);
  dumpXFaceDIMArea(NULL);
  dumpYFaceDIMArea(NULL);
  dumpFlux1Area(NULL);
  dumpEB1Area(NULL);
  dumpEBDIMArea(NULL);
  dumpBD1(NULL);
  dumpBDDIM(NULL);
}
#endif
#include "Chombo_NamespaceFooter.H"

#endif
