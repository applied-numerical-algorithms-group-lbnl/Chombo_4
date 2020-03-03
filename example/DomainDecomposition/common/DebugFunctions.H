#ifndef __DBGFunctionss__
#define __DBGFunctionss__
#define PI 3.141592653589793
#include "Chombo_NamespaceHeader.H"
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
fullDumpCell1(EBBoxData<CELL, Real, 1>* dataPtr)
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
    Point lodeb = databox.low();
    Point hideb = databox.high();
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
            cout  << data(vofs[0], 0) << "  ";
            if(data(vofs[0], 0) >= 0)
            {
              cout << " ";
            }
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


void 
fullDumpCellDIM(EBBoxData<CELL, Real, DIM>* dataPtr)
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
    Point lodeb = databox.low();
    Point hideb = databox.high();
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
              cout  << data(vofs[0], ivar) << "  ";
              if(data(vofs[0], ivar) >= 0)
              {
                cout << " ";
              }
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
    EBBoxData<CELL, Real, 1> & data = *dataPtr;
    Bx databox = dataPtr->box();
    databox &= graph.getDomain();
    Point lodeb = databox.low();
    Point hideb = databox.high();
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
            if(data(vofs[0], 0) >= 0.0) cout << " ";
            cout  << data(vofs[0], 0) << "  ";
          }
          else
          {
            cout  << "----------------" << "  ";
          }
        }
      }
      cout << endl;
    }
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
    Bx area(lodeb, hideb);

    EBGraph graph = dataPtr->ebgraph();
    EBBoxData<CELL, Real, DIM> & data = *dataPtr;
    Bx databox = dataPtr->box();
    for(int  icomp = 0; icomp < DIM; icomp++)
    {
      cout << "comp = :" << icomp << endl;
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
              if(data(vofs[0], icomp) >= 0.0) cout << " ";
              cout  << data(vofs[0], icomp) << "  ";
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

    Bx area = databox & graph.getDomain();
    Point lodeb = area.low();
    Point hideb = area.high();


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

void dummyfunc()
{
  dumpLDDIM(NULL);
  dumpLD1(NULL);
  dumpXFace(NULL);
  dumpYFace(NULL);
  dumpEB1(NULL);
}
#endif
#include "Chombo_NamespaceFooter.H"

#endif
