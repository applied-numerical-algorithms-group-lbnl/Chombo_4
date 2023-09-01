#ifndef _PROTO_STACK_H_
#define _PROTO_STACK_H_

#include <memory>
#include <stack>
#include "Proto_PAssert.H"
#include "Proto_MemType.H"

#define PR_STACK_ALLOC    4294967296 //4GB
//#define PR_STACK_OVERFLOW 268435456  //256MB
#define PR_STACK_OVERFLOW 536870912  //512MB
namespace Proto {

    template<MemType MEM = MEMTYPE_DEFAULT>
    class Stack
    {
        public:
        
        /// Get Stack Singleton
        static Stack<MEM>& getStack()
        {
            static Stack<MEM> s_stack(PR_STACK_ALLOC);
            return s_stack;
        }
        
        /// Constructor
        inline Stack(size_t a_capacity)
        {
            protoMalloc(MEM, m_stack, a_capacity);
            m_top = m_stack;
            m_capacity = a_capacity;
            m_hasOverflow = false;
            m_useStack = false;
        }
        
        /// Destructor
        inline ~Stack()
        {
            protoFree(MEM, m_stack);
        }
        
        /// Public Allocate
        inline void* alloc(size_t a_sz)
        {
            PROTO_ASSERT(a_sz > 0,
                "Stack::alloc | Error: Tried to allocate an array of size 0 to Stack.");
            PROTO_ASSERT(m_useStack,
                "Stack::alloc | Error: Tried to use the stack when it is not turned on. Did you forget to call PR_STACK_ON?");
#ifdef PROTO_CUDA
            const size_t align = 256;
            size_t increment = align*((a_sz+align-1)/align);
#else
            size_t increment = a_sz;
#endif
            if (m_useStack)
            {
                return doAlloc(increment);
            } else {
                return NULL;
            }
        }

        /// Toggle Stack
        inline void toggleStack(bool a_flag)
        {
            PROTO_ASSERT(a_flag != m_useStack,
                "Stack::toggleStack | Error: Toggling stack usage to %d when it is already %d. Did you forget to call PR_STACK_OFF?", a_flag, m_useStack);
            m_useStack = a_flag;
            if (!a_flag) { freeStack(); }
        }

        /// Query Stack Active
        inline bool enabled()
        {
            return m_useStack;
        }
        
        /// Query Unused
        inline size_t unused()
        {
            size_t u = m_top - m_stack;
            if (m_hasOverflow)
            {
                return (m_capacity - u) + m_overflow->unused();
            } else {
                return (m_capacity - u);
            }
        }
        
        /// Query Unused
        inline size_t used()
        {
            size_t u = m_top - m_stack;
            if (m_hasOverflow)
            {
                return m_overflow->used() + u; 
            } else {
                return u;
            }
        }

        private:
        
        /// Private Allocate
        inline void* doAlloc(size_t a_increment)
        {
            if (m_top + a_increment - m_stack < m_capacity)
            {
                void* rtn = (void*)(m_top);
                m_top += a_increment;
                m_increments.push(a_increment);
                return rtn;
            } else {
                if (m_hasOverflow)
                {
                    return m_overflow->doAlloc(a_increment);
                } else {
                    m_overflow = std::make_shared<Stack<MEM>>(PR_STACK_OVERFLOW);
                    m_hasOverflow = true;
                    return m_overflow->doAlloc(a_increment);
                }
            }
        }
        
        /// Free Stack
        inline void freeStack()
        {
            if (m_hasOverflow)
            {
                m_overflow->freeStack();
            }
            while(!m_increments.empty())
            {
                m_increments.pop();
            }
            m_top = m_stack;
        }


        bool               m_useStack;
        char*              m_stack;
        char*              m_top;
        size_t             m_capacity;
        std::stack<size_t> m_increments;
        int                m_bufferIndex;
        bool               m_hasOverflow;
        std::shared_ptr<Stack<MEM>>  m_overflow;

    }; // end class Stack

#define PR_STACK_ON Stack<MEMTYPE_DEFAULT>::getStack().toggleStack(true);
#define PR_STACK_OFF Stack<MEMTYPE_DEFAULT>::getStack().toggleStack(false);

} // end namespace Proto
#endif // end include guard
