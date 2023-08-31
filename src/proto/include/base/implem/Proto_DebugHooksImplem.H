#ifndef __DEBUGHOOKSIMPLEM__
#define __DEBUGHOOKSIMPLEM__
void writeData1(BoxData<double, 1>* a_dataPtr)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
    WriteData<1>(*a_dataPtr, -1, 1.0, string("var"), string("debugdat"));
  }
}


void printData1(BoxData<double, 1>* a_dataPtr)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
    a_dataPtr->print();
   Box dom = a_dataPtr->box();
   for(BoxIterator bit = dom.begin(); bit != dom.end(); ++bit)
   {
     std::cout << *bit << " " << (*a_dataPtr)(*bit, 0)  << std::endl;
   }
  }
}
void printColumn1(BoxData<double, 1>* a_dataPtr, int column)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
    a_dataPtr->print();
   Box dom = a_dataPtr->box();
   for(BoxIterator bit = dom.begin(); bit != dom.end(); ++bit)
   {
     Point pt = *bit;
     if(pt[0] == column)
     {
       std::cout << *bit << " " << (*a_dataPtr)(*bit, 0)  << std::endl;
     }
   }
  }
}
void printColumnDim(BoxData<double, DIM>* a_dataPtr, int column)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
   a_dataPtr->print();
   Box dom = a_dataPtr->box();
   for(BoxIterator bit = dom.begin(); bit != dom.end(); ++bit)
   {
     Point pt = *bit;
     if(pt[0] == column)
     {
       std::cout << *bit << " ";
       for(int idir = 0; idir < DIM; idir++)
       {
         std::cout << (*a_dataPtr)(*bit, idir) << " ";
       }
       std::cout << std::endl;
     }
   }
  }
}

void printDataDim(BoxData<double, DIM>* a_dataPtr)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
   a_dataPtr->print();
   Box dom = a_dataPtr->box();
   for(BoxIterator bit = dom.begin(); bit != dom.end(); ++bit)
   {
     std::cout << *bit << " ";
     for(int idir = 0; idir < DIM; idir++)
     {
       std::cout << (*a_dataPtr)(*bit, idir) << " ";
     }
     std::cout << std::endl;
   }
  }
}


void  viewData1(BoxData<double, 1>* a_dataPtr)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
    WriteData<1>(*a_dataPtr, -1, 1.0, string("var"), string("debugdat"));
    int sysret = system("visit -o debugdat.vtk");
    std::cout << "system call returned " << sysret << std::endl;;
  }
}


void  maxMinData1(BoxData<double, 1>* a_dataPtr)
{
  if(a_dataPtr != NULL)
  {
    std::cout << " data max = " << a_dataPtr->max() << ", min = " << a_dataPtr->min() << ", absmax = "   << a_dataPtr->absMax() << std::endl;
  }
}
void  maxMinDataDim(BoxData<double, DIM>* a_dataPtr)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
    for(int idir = 0; idir < DIM; idir++)
    {
      std::cout << idir << " data max = " << a_dataPtr->max(idir) << ", min = " << a_dataPtr->min(idir) << ", absmax = "   << a_dataPtr->absMax(idir) << std::endl;
    }
  }
}


void writeDataDim(BoxData<double, DIM>* a_dataPtr)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
    WriteData<DIM>(*a_dataPtr, -1, 1.0, string("var"), string("debugdat"));
  }
}

void  viewDataDim(BoxData<double, DIM>* a_dataPtr)
{
using Proto::WriteData;
using Proto::Box;
using Proto::BoxIterator;
using Proto::Point;
  if(a_dataPtr != NULL)
  {
    WriteData<DIM>(*a_dataPtr, -1, 1.0, string("var"), string("debugdat"));
    int sysret = system("visit -o debugdat.vtk");
    std::cout << "system call returned " << sysret << std::endl;
  }
}

void ProtoFunctionToGetStuffIntoTheSymbolTable()
{
  writeData1(NULL);
  viewData1(NULL);
  maxMinData1(NULL);
  writeDataDim(NULL);
  viewDataDim(NULL);
  maxMinDataDim(NULL);
}

#endif
