#include "Chombo_EBCM_GDB_Utilities.H"
#include "Chombo_ParmParse.H"
///
/**
   These functions exist to get templated functions into the symbol table.
   It is okay to call them but nothing should happen
**/
int ebcm_getPolynomialOrder()
{
  int order = 4586;
  Chombo4::ParmParse pp("main");
  pp.get("polynomial_order", order);
  Chombo4::pout() <<  "polynomial order = " << order << endl;
  return order;
}
  
void ebcm_gdb_dummyfunc()
{
  int order = ebcm_getPolynomialOrder();
  EBCM::dumpEigen(NULL);
  if(order == 1)
  {
    EBCM::GDB_Debug_Framework<1>::allTheFuncsNULL();
  }
  else if(order == 2)
  {
    EBCM::GDB_Debug_Framework<2>::allTheFuncsNULL();
  }
  else if(order == 3)
  {
    EBCM::GDB_Debug_Framework<3>::allTheFuncsNULL();
  }
  else if(order == 4)
  {
    EBCM::GDB_Debug_Framework<4>::allTheFuncsNULL();
  }
  else
  {
    Chombo4::pout() << "cannot find option for order = " << order << endl;
  }
}

