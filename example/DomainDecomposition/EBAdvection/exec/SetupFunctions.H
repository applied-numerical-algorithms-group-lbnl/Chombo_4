#ifndef __SetupFunctionss__
#define __SetupFunctionss__
#define PI 3.141592653589793
#include "Chombo_NamespaceHeader.H"
using Proto::Var;

#if DIM==2
Point lodeb(15, 15);
Point hideb(21, 21);
void 
dumpBlob(BoxData<Real, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Bx area(lodeb, hideb);

    BoxData<Real, 1> & data = *dataPtr;
    Bx databox = dataPtr->box();
    cout << "data region contains:" << endl;
    for(int j = hideb[1]; j >= lodeb[1]; j--)
    {
      for(int i = lodeb[0]; i <= hideb[0]; i++)
      {
        Point pt(i,j);
        if(databox.contains(pt))
        {
          cout << pt << ":" << data(pt, 0) << "  ";
        }
      }
      cout << endl;
    }
  }
}
void 
dumpEB1(EBBoxData<CELL, Real, 1>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Bx area(lodeb, hideb);

    EBGraph graph = dataPtr->ebgraph();
    EBBoxData<CELL, Real, 1> & data = *dataPtr;
    Bx databox = dataPtr->box();
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
          }
          else
          {
            cout  << "------------" << "  ";
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
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

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
              cout  << data(vofs[0], icomp) << "  ";
            }
            else
            {
              cout  << "------------" << "  ";
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
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Bx area(lodeb, hideb);

    EBBoxData<XFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();
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
              cout  << data(faces[0], 0) << "  ";
            }
          }
          else
          {
            printLines = true;
          }
          if(printLines)
          {
            cout  << "------------" << "  ";
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
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Bx area(lodeb, hideb);

    EBBoxData<YFACE, Real, 1> & data = *(*dataPtr);
    EBGraph graph = data.ebgraph();
    Bx databox = data.box();

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
              cout  << data(faces[0], 0) << "  ";
            }
          }
          else
          {
            printLines = true;
          }
          if(printLines)
          {
            cout  << "------------" << "  ";
          }
        }
      }
      cout << endl;
    }
  }
}
#endif
//=================================================
PROTO_KERNEL_START 
void InitializeSpotF(int       a_p[DIM],
                     Var<Real, 1>       a_phi,
                     Real      a_X0,
                     Real      a_rad,
                     Real      a_dx)
{
  Real val = 0;
  Real rlocsq = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    Real xrel = (a_p[idir] + 0.5)*a_dx - a_X0;
    rlocsq += xrel*xrel;
  }
  Real rloc = sqrt((Real) rlocsq);
  
  if(rloc < a_rad)
  {
    Real carg = 0.5*PI*(rloc)/a_rad;
    Real cosval = cos(carg);
    val = cosval*cosval;
  }
  a_phi(0) = val;
}
PROTO_KERNEL_END(InitializeSpotF, InitializeSpot)


//=================================================
PROTO_KERNEL_START 
void InitializeLineF(int          a_p[DIM],
                     Var<Real, 1> a_phi,
                     Real      a_px,
                     Real      a_rad,
                     Real      a_dx)
{
  Real val = 0;
  Real xloc[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    xloc[idir] = (a_p[idir] + 0.5)*a_dx;
    
  }
  
  Real funcval = std::abs(a_px - xloc[0]);
  Real rloc = funcval;
  
  if(rloc < a_rad && (xloc[1] < 0.75) )
  {
    Real carg = 0.5*PI*(rloc)/a_rad;
    Real cosval = cos(carg);
    val = std::abs(cosval);
//    val = 1.;
  }
  a_phi(0) = val;
}
PROTO_KERNEL_END(InitializeLineF, InitializeLine)

//=================================================
PROTO_KERNEL_START 
void InitializeVCellF(int       a_p[DIM],
                      Var<Real, DIM>       a_vel,
                      Real      a_cen,  //geom center
                      Real      a_rad,  //geom rad
                      Real      a_mag,  //max vel
                      Real      a_maxr, //radius for max vel
                      Real      a_dx)
{
  Real rlocsq = 0;
  Real xrel[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    xrel[idir] = (a_p[idir] + 0.5)*a_dx - a_cen;
    rlocsq += xrel[idir]*xrel[idir];
  }

  Real rloc = sqrt((Real) rlocsq);
  Real carg = 0.5*PI*(rloc - a_maxr)/a_maxr;
  Real cosval = cos(carg);
  Real velmag = a_mag*cosval*cosval;
//debug
  velmag = 1;
#if DIM==2
  a_vel(0) =  velmag*xrel[1];
  a_vel(1) = -velmag*xrel[0];
#else
  a_vel(0) =  velmag*( xrel[1] + xrel[2]);
  a_vel(1) =  velmag*(-xrel[0] - xrel[2]);
  a_vel(2) = -velmag*( xrel[0] - xrel[1]);
#endif
#if 0
  a_vel(0) =  0;
  a_vel(1) =  1;
#endif
//end debug
}
PROTO_KERNEL_END(InitializeVCellF, InitializeVCell)


//=================================================
PROTO_KERNEL_START 
void InitializeVCellConstF(int       a_p[DIM],
                           Var<Real, DIM>       a_vel,
                           Real      a_xvel,
                           Real      a_yvel)
{
  a_vel(0) =  a_xvel;
  a_vel(1) =  a_yvel;
#if DIM==3  
  a_vel(2) =  0;
#endif
//end debug
}
PROTO_KERNEL_END(InitializeVCellConstF, InitializeVCellConst)
#include "Chombo_NamespaceFooter.H"

#endif
