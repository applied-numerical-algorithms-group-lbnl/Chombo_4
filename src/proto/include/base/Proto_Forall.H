#pragma once
#ifndef _PROTO_FORALL_H_
#define _PROTO_FORALL_H_

#include "Proto_BoxData.H"

namespace Proto {

/** @name Pointwise Operators

   The suite of forall functions facilitate writing functions that operate pointwise on BoxData.
   To this end, the user must write a function with one of the following structures:
   @code
   PROTO_KERNEL_START
   void F_temp(Var<T,C,MEMTYPE,D,E>&, Args...)
   { ... }
   PROTO_KERNEL_END(F_temp, F)
   // OR
   PROTO_KERNEL_START
   void F_p_temp(Point&, Var<T,C,MEMTYPE,D,E>&, Args...)
   { ... }
   PROTO_KERNEL_END(F_p_temp, F_p)
   @endcode
   - PROTO_KERNEL_START and PROTO_KERNEL_END are required for cross-platform code.
   - The "#_temp" symbols are temporaries; the actual function symbol is the one without "_temp"
   - The template arguments of the first Var argument must match the output BoxData
   - The Point argument in the second signature corresponds to the Point of function application
   - Args... may include any number of Var& or read-only scalars (including Point or Box).
   - The elements of Args... may have arbitrary tensor structure and const-ness
   - non-const objects in Args... have input-output semantics
   - The order and template arguments of the Vars in Args... must match the BoxData inputs of forall
   - If F is a member function of a class F MUST BE DECLARED STATIC
   - F or F_p may be an anonymous (lambda) function defined using the PROTO_LAMBDA macro

   Refer to the following code snippet for some sample valid forall input functions:
   \snippet Snippets.cpp proto_forall_func
*/
///@{
    /// Pointwise Operator
    /**
      Computes the function <code> a_F </code> at each Point of this BoxData.
      This version of forall returns a new BoxData corresponding to the first argument of 
      <code> a_F </code> which must be a <code> Var </code>.
      The domain of the created output is equal to the intersection of all other BoxData inputs 
      (the function will fail if there are no other inputs).
      This function MUST specify the template parameters of the output BoxData. 
      See the example code snippet below.

      Example Usage:
      Input Function:
      \snippet Snippets.cpp proto_forall_constoprim
      Calling forall:
      \snippet Snippets.cpp proto_forall_1   

      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
      \tparam T       (Mandatory!) Data type of return BoxData
      \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
      \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
      \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1

    */  
    template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEMTYPE=MEMTYPE_DEFAULT, typename Func, typename... Srcs>
    inline BoxData<T,C,MEMTYPE,D,E> forall(const Func& a_F, Srcs&&... a_srcs);

    /// Pointwise Operator (Instrumented)
    /**
      Overload of the pointwise operator with additional instrumentation inputs
     
      \param a_num_flops_point  Number of flops used to compoute a_F at each Point
      \param a_timername        The name of a timer
      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
      \tparam T       (Mandatory!) Data type of return BoxData
      \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
      \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
      \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
    */
    template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEMTYPE=MEMTYPE_DEFAULT, typename Func, typename... Srcs>
    inline BoxData<T,C,MEMTYPE,D,E> forallOp(
            unsigned long long int a_num_flops_point,   
            const char* a_timername, 
            const Func& a_F, Srcs&&... a_srcs);

    /// Pointwise Operator (Explicit Range Box)
    /**
      Computes the function <code> a_F </code> at each Point of this BoxData.
      This version of forall returns a new BoxData corresponding to the first argument of 
      <code> a_F </code> which must be a <code> Var </code>.
      This function MUST specify the template parameters of the output BoxData. 
      See the example code snippet below.

      In general, this function should not be used unless absolutely necessary. 
      Some valid use cases are:
      - Creating and initializing a BoxData without any BoxData inputs
      - Evaluating a pointwise function on a Box that is a proper subset of the intersection
        of all input BoxData
      <code>a_box</code> must be a subset of the union of all input BoxData domains.

      Example Usage:
      Input Function:
      \snippet Snippets.cpp proto_forall_constoprim
      Calling forall:
      \snippet Snippets.cpp proto_forall_2   

      \param a_F      Pointwise function
      \param a_box    Range on which a_F will be evaluated
      \param a_srcs   Inputs (BoxData and primitive data)
      \tparam T       (Mandatory!) Data type of return BoxData
      \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
      \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
      \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
      */
    template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEMTYPE=MEMTYPE_DEFAULT, typename Func, typename... Srcs>
    inline BoxData<T,C,MEMTYPE,D,E> forall(const Func& a_F, Box  a_box, Srcs&&... a_srcs); 

    /// Pointwise Operator (Explicit Range Box, Instrumented)
    /**
      Overload of the pointwise operator with additional instrumentation inputs
     
      \param a_num_flops_point  Number of flops used to compoute a_F at each Point
      \param a_timername        The name of a timer
      \param a_F      Pointwise function
      \param a_box    Range on which a_F will be evaluated
      \param a_srcs   Inputs (BoxData and primitive data)
      \tparam T       (Mandatory!) Data type of return BoxData
      \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
      \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
      \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
    */
    template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEMTYPE=MEMTYPE_DEFAULT, typename Func, typename... Srcs>
    inline BoxData<T,C,MEMTYPE,D,E> forallOp(
            unsigned long long int  a_num_flops_point,   
            const char*             a_timername, 
            const Func& a_F, Box  a_box, Srcs&&... a_srcs); 

    /// Pointwise Operator (Pointwise Dependence)
    /**
      Computes the function <code> a_F </code> at each Point of this BoxData.
      This version of forall allows the input function to be dependent on the Point 
      at which it is applied. Hence, the first argument of a_F is a Point&, followed by 
      the Var corresponding to the output BoxData.
      This function MUST specify the template parameters of the output BoxData.
      See the example code snippet below.

      Example Usage:
      Input Function:
      \snippet Snippets.cpp proto_forall_pointwise
      Calling forall:
      \snippet Snippets.cpp proto_forall_3   

      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
      \tparam T       (Mandatory!) Data type of return BoxData
      \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
      \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
      \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
      */
    template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEMTYPE=MEMTYPE_DEFAULT, typename Func, typename... Srcs>
    inline BoxData<T,C,MEMTYPE,D,E> forall_p(const Func& a_F, Srcs&&... a_srcs);

    /// Pointwise Operator (Pointwise Dependence, Instrumented)
    /**
      Overload of the pointwise operator with additional instrumentation inputs
     
      \param a_num_flops_point  Number of flops used to compoute a_F at each Point
      \param a_timername        The name of a timer
      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
      \tparam T       (Mandatory!) Data type of return BoxData
      \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
      \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
      \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
    */
    template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEMTYPE=MEMTYPE_DEFAULT, typename Func, typename... Srcs>
    inline BoxData<T,C,MEMTYPE,D,E> forallOp_p(
            unsigned long long int a_num_flops_point,   
            const char* a_timername, 
            const Func& a_F, Srcs&&... a_srcs);

    /// Pointwise Operator (Pointwise Dependence, Explicit Range Box)
    /**
      Computes the function <code> a_F </code> at each Point of this BoxData.
      This version of forall allows the input function to be dependent on the 
      Point at which it is applied. Hence, the first argument of a_F is a Point&, 
      followed by the Var corresponding to the output BoxData.
      This function MUST specify the template parameters of the output BoxData.
      See the example code snippet below.

      In general, this function should not be used unless absolutely necessary.
      <code>a_box</code> must be a subset of the union of all input BoxData domains.

      Example Usage:
      Input Function:
      \snippet Snippets.cpp proto_forall_pointwise
      Calling forall:
      \snippet Snippets.cpp proto_forall_4   

      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
      \param a_box    Range on which a_F will be evaluated
      \tparam T       (Mandatory!) Data type of return BoxData
      \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
      \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
      \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
      */
    template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEMTYPE=MEMTYPE_DEFAULT, typename Func, typename... Srcs>
    inline BoxData<T,C,MEMTYPE,D,E> forall_p(const Func& a_F, Box a_box, Srcs&&... a_srcs); 

    /// Pointwise Operator (Pointwise Dependence, Explicit Range Box, Instrumented)
    /**
      Overload of the pointwise operator with additional instrumentation inputs
     
      \param a_num_flops_point  Number of flops used to compoute a_F at each Point
      \param a_timername        The name of a timer
      \param a_F      Pointwise function
      \param a_box    Range on which a_F will be evaluated
      \param a_srcs   Inputs (BoxData and primitive data)
      \tparam T       (Mandatory!) Data type of return BoxData
      \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
      \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
      \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
    */
    template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
        MemType MEMTYPE=MEMTYPE_DEFAULT, typename Func, typename... Srcs>
    inline BoxData<T,C,MEMTYPE,D,E> forallOp_p(
            unsigned long long int a_num_flops_point,   
            const char* a_timername, 
            const Func& a_F, Box a_box, Srcs&&... a_srcs); 

    /// In-Place Pointwise Operator
    /**
      Computes the function <code> a_F </code> at each Point of this BoxData.
      The "InPlace" versions of forall execute on existing BoxData and do not produce a new array.
      The range of the function evaluation is equal to the intersection of all BoxData inputs.

      Example Usage:
      Input Function:
      \snippet Snippets.cpp proto_forall_constoprim
      Calling forall:
      \snippet Snippets.cpp proto_forall_5   

      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
      */
    template<typename Func, typename... Srcs>
        inline void forallInPlace(const Func& a_F, Srcs&&... a_srcs);
    
    /// In-Place Pointwise Operator (Instrumented)
    /**
      Overload of the pointwise operator with additional instrumentation inputs
     
      \param a_num_flops_point  Number of flops used to compoute a_F at each Point
      \param a_timername        The name of a timer
      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
    */
    template<typename Func, typename... Srcs>
    inline void forallInPlaceOp(
            unsigned long long int a_num_flops_point,   
            const char* a_timername, 
            const Func& a_F, Srcs&&... a_srcs);

    /// In-Place Pointwise Operator (Explicit Range Box)
    /**
      Computes the function <code> a_F </code> at each Point of this BoxData.
      The "InPlace" versions of forall execute on existing BoxData and do not produce a new array.
      <code> a_F </code> will be applied at all points of the input <code> a_box </code>.

      In general, this function should not be used unless you want to restrict the range of 
      a_F's application to be something smaller than the intersection of all the inputs.   
      <code>a_box</code> must be a subset of the union of all input BoxData domains.

      Example Usage:
      Input Function:
      \snippet Snippets.cpp proto_forall_constoprim
      Calling forall:
      \snippet Snippets.cpp proto_forall_6   

      \param a_F      Pointwise function
      \param a_box    Range of computation.
      \param a_srcs   Inputs (BoxData and primitive data)
      */
    template<typename Func, typename... Srcs>
    inline void forallInPlace(const Func& a_F, Box a_box, Srcs&&... a_srcs); 

    /// In-Place Pointwise Operator (Explicit Range Box, Instrumented)
    /**
      Overload of the pointwise operator with additional instrumentation inputs
     
      \param a_num_flops_point  Number of flops used to compoute a_F at each Point
      \param a_timername        The name of a timer
      \param a_F      Pointwise function
      \param a_box    Range of computation.
      \param a_srcs   Inputs (BoxData and primitive data)
    */
    template<typename Func, typename... Srcs>
    inline void forallInPlaceOp(
            unsigned long long int a_num_flops_point,   
            const char* a_timername, 
            const Func& a_F, Box a_box, Srcs&&... a_srcs); 

    /// In-Place Pointwise Operator (Pointwise Dependence)
    /**   
      Computes the function <code> a_F </code> at each Point of this BoxData.
      This version of forall allows the input function to be dependent on the 
      Point at which it is applied. Hence, the first argument of <code> a_F </code> 
      is a <code> Point& </code>, followed by the normal <code> Var </code> inputs

      Example Usage:
      Input Function:
      \snippet Snippets.cpp proto_forall_pointwise
      Calling forall:
      \snippet Snippets.cpp proto_forall_7 

      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
      */
    template<typename Func, typename... Srcs>
    inline void forallInPlace_p(const Func& a_F, Srcs&&... a_srcs);

    /// In-Place Pointwise Operator (Pointwise Dependence, Instrumented)
    /**
      Overload of the pointwise operator with additional instrumentation inputs
     
      \param a_num_flops_point  Number of flops used to compoute a_F at each Point
      \param a_timername        The name of a timer
      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
    */
    template<typename Func, typename... Srcs>
    inline void forallInPlaceOp_p(
            unsigned long long int a_num_flops_point,   
            const char* a_timername, 
            const Func& a_F, Srcs&&... a_srcs);

    /// In-Place Pointwise Operator (Pointwise Dependence, Explicit Box Range)
    /**
      Computes the function <code> a_F </code> at each Point of this BoxData.
      This version of forall allows the input function to be dependent on the 
      Point at which it is applied. Hence, the first argument of <code> a_F </code> 
      is a <code> Point& </code>, followed by the normal <code> Var </code> inputs

      In general, this function should not be used unless you want to restrict the 
      domain of a_F's application to be something smaller than the intersection of all the inputs.   
      <code>a_box</code> must be a subset of the union of all input BoxData domains.

      Example Usage:
      Input Function:
      \snippet Snippets.cpp proto_forall_pointwise
      Calling forall:
      \snippet Snippets.cpp proto_forall_8 

      \param a_F      Pointwise function
      \param a_srcs   Inputs (BoxData and primitive data)
      \param a_box    Range of computation
      */
    template<typename Func, typename... Srcs>
    inline void forallInPlace_p(const Func& a_F,  Box a_box, Srcs&&... a_srcs); 

    /// In-Place Pointwise Operator (Pointwise Dependence, Explicit Box Range, Instrumented)
    /**
      Overload of the pointwise operator with additional instrumentation inputs
     
      \param a_num_flops_point  Number of flops used to compoute a_F at each Point
      \param a_timername        The name of a timer
      \param a_F      Pointwise function
      \param a_box    Range of computation.
      \param a_srcs   Inputs (BoxData and primitive data)
    */
    template<typename Func, typename... Srcs>
    inline void forallInPlaceOp_p(
            unsigned long long int a_num_flops_point,   
            const char* a_timername, 
            const Func& a_F,  Box a_box, Srcs&&... a_srcs); 

    // End of BoxData Doxygen Module
    /*@}*/

    //========================================================================
    // CUDA API ||
    //==========++
#ifdef PROTO_CUDA
    template<typename Func, typename... Srcs>
    inline void protoForall(const Func& a_F, Box a_box, Srcs&&... a_srcs); 
    
    template<typename Func, typename... Srcs>
    inline void protoForallStream(
            protoStream_t& a_stream,
            const Func& a_F, Box a_box, Srcs&&... a_srcs); 

    template<typename Func, typename... Srcs>
    inline void protoForall_p(const Func& a_F, Box a_box, Srcs&&... a_srcs); 
#endif
#include "implem/Proto_ForallImplem.H"
} // end namespace Proto
#endif // end include guard
