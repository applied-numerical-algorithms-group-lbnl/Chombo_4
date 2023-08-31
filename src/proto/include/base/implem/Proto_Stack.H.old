#ifndef _STACK_ALLOC_H_
#define _STACK_ALLOC_H_

#include <cstddef>
#include <stack>
#include <iostream>
#include <memory>
#include "Proto_MemType.H"
using std::shared_ptr;
///  class Stack implements a memory allocation as chunks of a big array.
/** Stack is meant to be the default memory allocation for automatic variables that have runtime storage requirements, 
    like FArrayBox's dataPtr.  It relies on the automatic variable semantics of LIFO.  This is literally a static storage array and
    a stack-based LIFO control system.  When code is built with OpenMP, the static arrays are threadprivate.
*/

template <Proto::MemType MEMTYPE>
class StackAlloc
{
    /// base class constructor
public:
  inline StackAlloc()
  {
      protoMalloc(MEMTYPE, m_stack, PR_STACK_ALLOC);
      m_top = m_stack;
      m_maxusage= 0;
  }


  inline ~StackAlloc()
  {
      protoFree(MEMTYPE, m_stack);
//    std::cout << "max stack memory usage at destruction = " << m_maxusage << std::endl;
  }

  /**
     Allocate a dynamic memory arena of size a_sz.
     A pointer to this memory should be returned.
  */
  inline void* alloc(size_t a_sz)
  {
    PROTO_ASSERT(false, "Stack::alloc | Error: This function should not be called.");
      PR_assert(a_sz>0);
    void* rtn = (void*)(m_top);
    //compute the safe increment into the stack and reset the top value
    //size_t increment = (2 + (a_sz-1)/sizeof(size_t))*sizeof(size_t); // safe round up divide and move ahead to next

#ifdef PROTO_CUDA
     const size_t align = 256;
     size_t increment = align*((a_sz+align-1)/align);
#else
     size_t increment = a_sz;
#endif
//    std::cout << "alloc increment = " << increment << std::endl;
    m_top += increment;
    m_increments.push(increment);

    size_t curusage = (m_top - m_stack);
    if(curusage > m_maxusage)
    {
      m_maxusage = curusage;
    }
    
    //make certain we have not gone off the end of the world
    //if((uintptr_t)m_top > ((uintptr_t)m_stack+PR_STACK_ALLOC)) throw std::bad_alloc();
    if(curusage > PR_STACK_ALLOC) throw std::bad_alloc();

    return rtn;
  }

  /**
     virtual function for deleting the arena pointed to by a_pt.
     rewinds top by size of last allocaiton.
  */
  inline void free(void* a_pt)
  {
    //remember how the increment was stored in the extra bit at the end of the allocation?
    //this is only technically correct if objects are deallocated in the reverse order of allocation
    if(!m_increments.empty())
    {
//    size_t  decrement = m_increments.top();
//    std::cout << "free decrement = " << decrement << std::endl;
      m_top  -= m_increments.top();
      m_increments.pop();
    }
  }

  inline void reset()
  {
    while (!m_increments.empty())
    {
      m_increments.pop();
    }
    m_top = m_stack;
  }

  /**
     accesor to m_top, debugging tool
  */
  inline char* getTop()
  {
    return m_top;
  } 

  /**
     accesor to m_stack, debugging tool
  */
  inline char* getStack()
  {
    return m_stack;
  }


  /**
     accessor to m_increments, using for checking
  */
  std::stack<size_t>& getIncrements()
  {
    return m_increments;
  }

private:
  char *            m_stack;
  char *            m_top;
  size_t             m_maxusage;
  std::stack<size_t> m_increments;

//#pragma omp threadprivate(m_stack, m_top, m_maxusage, m_increments)

};

// Defined in DisjointBoxLayout. If you use this function, you create a second stack
/*
template<Proto::MemType MEMTYPE>
StackAlloc<MEMTYPE>& getStackAllocator()
{
  static StackAlloc<MEMTYPE> stack;
  return stack;
}*/


#endif //end include guard
