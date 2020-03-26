#ifndef __HostDBGFunctionss__
#define __HostDBGFunctionss__
#define PI 3.141592653589793
#include "Proto_EBHostData.H"
#include "Proto_EBHostFluxData.H"
#include "Chombo_NamespaceHeader.H"
using Proto::EBHostData;
using Proto::EBHostFluxData;
#if DIM==2
void 
genDumpHostCellDIM(EBHostData<CELL, Real, DIM>* dataPtr, Point lodeb, Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBGraph graph = dataPtr->ebgraph();
    EBHostData<CELL, Real, DIM> & data = *dataPtr;
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
genDumpHostCell1(EBHostData<CELL, Real, 1>* dataPtr, Point lodeb, Point hideb)
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
fullDumpHostCell1(EBHostData<CELL, Real, 1>* dataPtr)
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
    genDumpHostCell1(dataPtr, lodeb, hideb);
  }
}


void 
fullDumpHostCellDIM(EBHostData<CELL, Real, DIM>* dataPtr)
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
    genDumpHostCellDIM(dataPtr, lodeb, hideb);
  }
}
void 
dumpHostLD1(LevelData<EBHostData<CELL, Real, 1> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBHostData<CELL, Real, 1> >& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      fullDumpHostCell1(&data[dit[ibox]]);
    }
  }
}
void 
dumpHostLDDIM(LevelData<EBHostData<CELL, Real, DIM> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBHostData<CELL, Real, DIM> >& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      fullDumpHostCellDIM(&data[dit[ibox]]);
    }
  }
}



void 
dumpHostEB1(EBHostData<CELL, Real, 1>* dataPtr)
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
    genDumpHostCell1(dataPtr, lodeb, hideb);
  }
}




void 
dumpHostEBDIM(EBHostData<CELL, Real, DIM>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Point lodeb = dataPtr->box().low();
    Point hideb = dataPtr->box().high();
    genDumpHostCellDIM(dataPtr, lodeb, hideb);
  }
}


void 
genDumpHostXFace1(shared_ptr<EBHostData<XFACE, Real, 1> >* dataPtr,
              Point lodeb,
              Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBHostData<XFACE, Real, 1> & data = *(*dataPtr);
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
dumpHostXFace(shared_ptr<EBHostData<XFACE, Real, 1> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBHostData<XFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(0, -1);

    Bx area = databox & graph.getDomain();
    Point lodeb = area.low();
    Point hideb = area.high();
    genDumpHostXFace1(dataPtr, lodeb, hideb);
  }
}


void 
genDumpHostXFaceDIM(shared_ptr<EBHostData<XFACE, Real, DIM> >* dataPtr, 
                Point lodeb,
                Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBHostData<XFACE, Real, DIM> & data = *(*dataPtr);
    for(int idir = 0 ; idir < DIM; idir++)
    {
      Bx databox = data.box();
      if((databox.contains(lodeb)) || (databox.contains(hideb)))
      {
 
        pout() << "variable " << idir << endl;
        EBGraph graph = data.ebgraph();
        Bx databox = data.box();
        databox.grow(0, -1);

        Bx area = databox & graph.getDomain();

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
dumpHostXFaceDIM(shared_ptr<EBHostData<XFACE, Real, DIM> >* dataPtr)
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
    genDumpHostXFaceDIM(dataPtr, lodeb, hideb);
  }
}


void 
genDumpHostYFace1(shared_ptr<EBHostData<YFACE, Real, 1> >* dataPtr,
              Point lodeb,
              Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBHostData<YFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(0, -1);

    if((databox.contains(lodeb)) || (databox.contains(hideb)))
    {
    Bx area = databox & graph.getDomain();

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
dumpHostYFace(shared_ptr<EBHostData<YFACE, Real, 1> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBHostData<YFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(0, -1);

    Bx area = databox & graph.getDomain();             // 
    Point lodeb = area.low();
    Point hideb = area.high();
    genDumpHostYFace1(dataPtr, lodeb, hideb);
  }
}


void 
genDumpHostYFaceDIM(shared_ptr<EBHostData<YFACE, Real, DIM> >* dataPtr,
                Point lodeb, Point hideb)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBHostData<YFACE, Real, DIM> & data = *(*dataPtr);
    for(int idir = 0; idir < DIM; idir++)
    {
      pout() << "variable  = "  << idir << endl;
      EBGraph graph = data.ebgraph();
      Bx databox = data.box();
      databox.grow(0, -1);
      if((databox.contains(lodeb)) || (databox.contains(hideb)))
      {

      Bx area = databox & graph.getDomain();

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
dumpHostYFaceDIM(shared_ptr<EBHostData<YFACE, Real, DIM> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBHostData<YFACE, Real, DIM> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
    databox.grow(1, -1);

    Bx area = databox & graph.getDomain();
    Point lodeb = area.low();
    Point hideb = area.high();
    genDumpHostYFaceDIM(dataPtr, lodeb, hideb);
  }
}

void 
dumpHostFlux1(EBHostFluxData<Real, 1>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    pout() << "x face data:" << endl;
    dumpHostXFace(&(dataPtr->m_xflux));
    pout() << "y face data:" << endl;
    dumpHostYFace(&(dataPtr->m_yflux));
  }
}

void 
dumpHostFluxDIM(EBHostFluxData<Real, DIM>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    pout() << "x face data:" << endl;
    dumpHostXFaceDIM(&(dataPtr->m_xflux));
    pout() << "y face data:" << endl;
    dumpHostYFaceDIM(&(dataPtr->m_yflux));
  }
}

void 
dumpHostLDFlux1(LevelData<EBHostFluxData<Real, 1> >*  dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBHostFluxData<Real, 1> > & data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      
      pout() << "level flux data for ibox = " << ibox << ":" << endl;
      dumpHostFlux1(&(data[dit[ibox]]));
    }
  }
}
void 
dumpHostLDFluxDIM(LevelData<EBHostFluxData<Real, DIM> >*  dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBHostFluxData<Real, DIM> >& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      
      pout() << "level flux data for ibox = " << ibox << ":" << endl;
      dumpHostFluxDIM(&(data[dit[ibox]]));
    }
  }
}



Point lodebareahost(1, 10);
Point hidebareahost(4, 20);                                    

void dumpHostEB1Area(EBHostData<CELL, Real, 1>* dataPtr)
{
  genDumpHostCell1(dataPtr, lodebareahost, hidebareahost);
}
void dumpHostEBDIMArea(EBHostData<CELL, Real, DIM>* dataPtr)
{
  genDumpHostCellDIM(dataPtr, lodebareahost, hidebareahost);
}

void 
dumpHostLD1Area(LevelData<EBHostData<CELL, Real, 1> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBHostData<CELL, Real, 1> >& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      genDumpHostCell1(&data[dit[ibox]], lodebareahost, hidebareahost);
    }
  }
}
void 
dumpHostLDDIMArea(LevelData<EBHostData<CELL, Real, DIM> >* dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBHostData<CELL, Real, DIM> >& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      genDumpHostCellDIM(&data[dit[ibox]], lodebareahost, hidebareahost);
    }
  }
}
void 
dumpHostXFaceDIMArea(shared_ptr<EBHostData<XFACE, Real, DIM> >* dataPtr)
{
  genDumpHostXFaceDIM(dataPtr, lodebareahost, hidebareahost);
}
void 
dumpHostXFace1Area(shared_ptr<EBHostData<XFACE, Real, 1> >* dataPtr)
{
  genDumpHostXFace1(dataPtr, lodebareahost, hidebareahost);
}
void 
dumpHostYFaceDIMArea(shared_ptr<EBHostData<YFACE, Real, DIM> >* dataPtr)
{
  genDumpHostYFaceDIM(dataPtr, lodebareahost, hidebareahost);
}
void 
dumpHostYFace1Area(shared_ptr<EBHostData<YFACE, Real, 1> >* dataPtr)
{
  genDumpHostYFace1(dataPtr, lodebareahost, hidebareahost);
}

void 
dumpHostFlux1Area(EBHostFluxData<Real, 1>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    pout() << "x face data:" << endl;
    dumpHostXFace1Area(&(dataPtr->m_xflux));
    pout() << "y face data:" << endl;
    dumpHostYFace1Area(&(dataPtr->m_yflux));
  }
}
void 
dumpHostFluxDIMArea(EBHostFluxData<Real, DIM>*  dataPtr)
{
  if(dataPtr != NULL)
  {
    pout() << "x face data:" << endl;
    dumpHostXFaceDIMArea(&(dataPtr->m_xflux));
    pout() << "y face data:" << endl;
    dumpHostYFaceDIMArea(&(dataPtr->m_yflux));
  }
}

void 
dumpHostLDFluxDIMArea(LevelData<EBHostFluxData<Real, DIM>  >*  dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBHostFluxData<Real, DIM> > & data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      
      pout() << "level flux data for ibox = " << ibox << ":" << endl;
      dumpHostFluxDIMArea(&(data[dit[ibox]]));
    }
  }
}
void 
dumpHostLDFlux1Area(LevelData<EBHostFluxData<Real, 1> > *  dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBHostFluxData<Real, 1> > & data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      
      pout() << "level flux data for ibox = " << ibox << ":" << endl;
      dumpHostFlux1Area(&(data[dit[ibox]]));
    }
  }
}

void dummyfunc_too_bad_for_me()
{
  dumpHostLDDIM(NULL);
  dumpHostLDFluxDIM(NULL);
  dumpHostLD1(NULL);
  dumpHostLDFlux1(NULL);
  dumpHostXFace(NULL);
  dumpHostYFace(NULL);
  dumpHostFlux1(NULL);
  dumpHostEB1(NULL);
  dumpHostLDDIMArea(NULL);
  dumpHostLDFluxDIMArea(NULL);
  dumpHostLD1Area(NULL);
  dumpHostLDFlux1Area(NULL);
  dumpHostXFace1Area(NULL);
  dumpHostYFace1Area(NULL);
  dumpHostXFaceDIMArea(NULL);
  dumpHostYFaceDIMArea(NULL);
  dumpHostFlux1Area(NULL);
  dumpHostEB1Area(NULL);
  dumpHostEBDIMArea(NULL);
}
#endif
#include "Chombo_NamespaceFooter.H"

#endif
