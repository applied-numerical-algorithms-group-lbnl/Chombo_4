#pragma once

#include "Proto_BoxData.H"
//#include "Proto_LevelBoxData.H"

namespace Operator {

    using namespace Proto;    

    /// 4th Order Convolution
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
        Distinguishes between 4th and 2nd order inputs.
         
        \param a_ctr    4th  order cell centered quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _convolve(
            const BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_2nd);
    
    /// 4th Order Convolution
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
        Equal to _convolve(ctr, ctr).

        \param a_ctr    4th  order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> convolve(
            const BoxData<T,C,MEM,D,E>& a_ctr);
    
    /// 4th Order Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
        Distinguishes between 4th and 2nd order inputs.
       
        \param a_avg    4th  order cell averaged quantity (output)
        \param a_ctr    4th  order cell centered quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _convolve(
            BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_2nd);
    
    /// 4th Order Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
        Equal to _convolve(avg, ctr, ctr).

        \param a_avg    4th  order cell averaged quantity (output)
        \param a_ctr    4th  order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void convolve(
            BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_ctr);
    
    /// 4th Order Convolution (Level, In Place)
    /**
        Compute the 4th order accurate average of a quantity over a level. 
        Distinguishes between 4th and 2nd order inputs.
        TODO: implement appropriate behavior for non-cell centered LevelData
         
        \param a_avg    4th  order cell averaged quantity (output)
        \param a_ctr    4th  order cell centered quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    /*
    template<class T, unsigned int C, MemType MEM, Centering CTR>
    inline void _convolve(
            LevelBoxData<T,C,MEM,CTR>& a_avg,
            const LevelBoxData<T,C,MEM,CTR>& a_ctr,
            const LevelBoxData<T,C,MEM,CTR>& a_2nd);
    */
    /// 4th Order Convolution (Level, In Place)
    /**
        Compute the 4th order accurate average of a quantity over a level. 
        Equal to _convolve(avg, ctr, ctr).
        TODO: implement appropriate behavior for non-cell centered LevelData
        
        \param a_avg    4th  order cell averaged quantity (output)
        \param a_ctr    4th  order cell centered quantity
    */
    /*
    template<class T, unsigned int C, MemType MEM, Centering CTR>
    inline void convolve(
            LevelBoxData<T,C,MEM,CTR>& a_avg,
            const LevelBoxData<T,C,MEM,CTR>& a_ctr);
    */
    /// 4th Order Deconvolution
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
        Distinguishes between 4th and 2nd order inputs

        \param a_avg    4th  order cell averaged quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _deconvolve(
            const BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_2nd);
    
    /// 4th Order Deconvolution
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
        Equal to _deconvolve(avg, avg).
            
        \param a_avg    4th  order cell averaged quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> deconvolve(
            const BoxData<T,C,MEM,D,E>& a_avg);
    
    /// 4th Order Deconvolution (In Place)
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
        Distinguishes between 4th and 2nd order inputs
        
        \param a_ctr    4th  order cell centered quantity (output)
        \param a_avg    4th  order cell averaged quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _deconvolve(
            BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_2nd);

    /// 4th Order Deconvolution (In Place)
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
        Equal to _deconvolve(ctr, avg, avg)
        
        \param a_ctr    4th  order cell centered quantity (output)
        \param a_avg    4th  order cell averaged quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void deconvolve(
            BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_avg);
    
    /// 4th Order Deconvolution (Level, In Place)
    /**
        Compute the 4th order accurate cell centered quantity over a level. 
        Distinguishes between 2nd and 4th order inputs.
        TODO: implement appropriate behavior for non-cell centered LevelData
         
        \param a_ctr    4th  order cell centered quantity (output)
        \param a_avg    4th  order cell averaged quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    /*
    template<class T, unsigned int C, MemType MEM, Centering CTR>
    inline void _deconvolve(
            LevelBoxData<T,C,MEM,CTR>& a_ctr,
            const LevelBoxData<T,C,MEM,CTR>& a_avg,
            const LevelBoxData<T,C,MEM,CTR>& a_2nd);
    */
    /// 4th Order Deconvolution (Level, In Place)
    /**
        Compute the 4th order accurate cell centered quantity over a level. 
        Equal to _deconvolve(ctr, avg, avg);
        TODO: implement appropriate behavior for non-cell centered LevelData
         
        \param a_ctr    4th  order cell centered quantity (output)
        \param a_avg    4th  order cell averaged quantity
    */
    /*
    template<class T, unsigned int C, MemType MEM, Centering CTR>
    inline void deconvolve(
            LevelBoxData<T,C,MEM,CTR>& a_ctr,
            const LevelBoxData<T,C,MEM,CTR>& a_avg);
    */
    /// 4th Order Face Convolution
    /**
        Compute the 4th order accurate average of a quantity over a face.
        Distinguishes between 2nd and 4th order inputs.

        \param a_ctr    4th  order face centered quantity
        \param a_2nd    2nd+ order face centered quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _convolveFace(
            const BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Convolution
    /**
        Compute the 4th order accurate average of a quantity over a face.
        Equal to _convolveFace(ctr, ctr, dir);

        \param a_ctr    4th  order face centered quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> convolveFace(
            const BoxData<T,C,MEM,D,E>& a_ctr,
            int a_dir);
    
    /// 4th Order Face Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a face. 
        Distinguishes between 2nd and 4th order inputs.

        \param a_avg    4th  order face averaged quantity (output)
        \param a_ctr    4th  order face centered quantity
        \param a_2nd    2nd+ order face centered quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _convolveFace(
            BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a face. 
        Equal to _convolveFace(avg, ctr, ctr, dir)
        
        \param a_avg    4th  order face averaged quantity (output)
        \param a_ctr    4th  order face centered quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void convolveFace(
            BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_ctr,
            int a_dir);
    
    /// 4th Order Face Deconvolution
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
        Distinguishes between 2nd and 4th order inputs.

        \param a_avg    4th  order face averaged quantity
        \param a_2nd    2nd+ order face averaged quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _deconvolveFace(
            const BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Deconvolution
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
        Equal to _deconvolveFace(avg, avg, dir)
        
        \param a_avg    4th  order face averaged quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> deconvolveFace(
            const BoxData<T,C,MEM,D,E>& a_avg,
            int a_dir);
   
    /// 4th Order Face Deconvolution (In Place)
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
        Distinguishes between 2nd and 4th order inputs.

        \param a_ctr    4th  order face centered quantity (output)
        \param a_avg    4th  order face averaged quantity
        \param a_2nd    2nd+ order face averaged quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _deconvolveFace(
            BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Deconvolution (In Place)
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
        Equal to _deconvolveFace(ctr, avg, avg, dir)
         
        \param a_ctr    4th  order face centered quantity (output)
        \param a_avg    4th  order face averaged quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void deconvolveFace(
            BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_avg,
            int a_dir);
   
    /// 4th Order Cell Averaged Product
    /**
        Compute the 4th order accurate cell averaged product rule <P>*<Q> 
        Distinguishes between 2nd and 4th order inputs.
    
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
        \param a_p2nd   2nd+ order cell centered quantity P
        \param a_q2nd   2nd+ order cell centered quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _cellProduct(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd);
    
    /// 4th Order Cell Averaged Product
    /**
        Compute the 4th order accurate cell averaged product rule <P>*<Q> 
        Equal to _cellProduct(pavg, qavg, pavg, qavg)
    
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> cellProduct(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg);
    
    /// 4th Order Cell Averaged Product (In Place)
    /**
        Compute the 4th order accurate cell averaged product rule <P>*<Q> 
        Distinguishes between 2nd and 4th order inputs.
        
        \param a_prod   4th  order cell averaged product P*Q (output)
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
        \param a_p2nd   2nd+ order cell centered quantity P
        \param a_q2nd   2nd+ order cell centered quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _cellProduct(
            BoxData<T,C,MEM,D,E>& a_prod,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd);
    
    /// 4th Order Cell Averaged Product (In Place)
    /**
        Compute the 4th order accurate cell averaged product rule <P>*<Q> 
        Equal to _cellProduct(prod, pavg, qavg, pavg, qavg)

        \param a_prod   4th  order cell averaged product P*Q (output)
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void cellProduct(
            BoxData<T,C,MEM,D,E>& a_product,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg);
    
    /// 4th Order Face Averaged Product
    /**
        Compute the 4th order accurate face averaged product rule <P>_d * <Q>_d 
        Distinguishes between 2nd and 4th order inputs.
        
        \param a_prod   4th  order face averaged product P*Q (output)
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_p2nd   2nd+ order face centered quantity P
        \param a_q2nd   2nd+ order face centered quantity Q
        \param a_dir    face normal 
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _faceProduct(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Face Averaged Product
    /**
        Compute the 4th order accurate face averaged product rule <P>_d * <Q>_d 
        Equal to _faceProduct(pavg, qavg, pavg, qavg, dir)
        
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> faceProduct(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            int a_dir);
    
    /// 4th Order Face Averaged Product (In Place)
    /**
        Compute the 4th order accurate face averaged product rule <P>_d * <Q>_d 
        Distinguishes between 2nd and 4th order inputs.

        \param a_prod   4th  order face averaged product P*Q (output)
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_p2nd   2nd+ order face centered quantity P
        \param a_q2nd   2nd+ order face centered quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _faceProduct(
            BoxData<T,C,MEM,D,E>& a_prod,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Face Averaged Product (In Place)
    /**
        Compute the 4th order accurate face averaged product rule <P>_d * <Q>_d 
        Equal to _faceProduct(prod, pavg, qavg, pavg, qavg, dir)
        
        \param a_prod   4th  order face averaged product P*Q (output)
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void faceProduct(
            BoxData<T,C,MEM,D,E>&       a_prod,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            int a_dir);
    
    /// 4th Order Cell Averaged Quotient
    /**
        Compute the 4th order accurate cell averaged quotient rule <P> / <Q> 
        Distinguishes between 2nd and 4th order inputs.

        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
        \param a_p2nd   2nd+ order cell centered quantity P
        \param a_q2nd   2nd+ order cell centered quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _cellQuotient(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd);
    
    /// 4th Order Cell Averaged Quotient
    /**
        Compute the 4th order accurate cell averaged quotient rule <P> / <Q> 
        Equal to _cellQuotient(pavg, qavg, pavg, qavg)

        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> cellQuotient(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg);
    
    /// 4th Order Cell Averaged Quotient (In Place)
    /**
        Compute the 4th order accurate cell averaged quotient rule <P> / <Q> 
        Distinguishes between 2nd and 4th order inputs.

        \param a_quot   4th  order cell averaged quotient P / Q (output)
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
        \param a_p2nd   2nd+ order cell centered quantity P
        \param a_q2nd   2nd+ order cell centered quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _cellQuotient(
            BoxData<T,C,MEM,D,E>&       a_quot,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd);
    
    /// 4th Order Cell Averaged Quotient (In Place)
    /**
        Compute the 4th order accurate cell averaged quotient rule <P> / <Q> 
        Equal to _cellQuotient(qout, pavg, qavg, pavg, qavg)

        \param a_quot   4th  order cell averaged quotient P / Q (output)
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void cellQuotient(
            BoxData<T,C,MEM,D,E>&       a_quot,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg);
    
    /// 4th Order Face Averaged Quotient
    /**
        Compute the 4th order accurate face averaged quotient rule <P>_d / <Q>_d 
        Distinguishes between 2nd and 4th order inputs.

        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_p2nd   2nd+ order face centered quantity P
        \param a_q2nd   2nd+ order face centered quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _faceQuotient(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Face Averaged Quotient
    /**
        Compute the 4th order accurate face averaged quotient rule <P>_d / <Q>_d 
        Equal to _faceQuotient(pavg, qavg, pavg, qavg, dir)
            
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> faceQuotient(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            int a_dir);

    /// 4th Order Face Averaged Quotient (In Place)
    /**
        Compute the 4th order accurate face averaged quotient rule <P>_d / <Q>_d 
        Distinguishes between 2nd and 4th order inputs.
                
        \param a_quot   4th  order face averaged quotient P / Q
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_p2nd   2nd+ order face centered quantity P
        \param a_q2nd   2nd+ order face centered quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _faceQuotient(
            BoxData<T,C,MEM,D,E>& a_quotient,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Face Averaged Quotient (In Place)
    /**
        Compute the 4th order accurate face averaged quotient rule <P>_d / <Q>_d 
        Equal to _faceQuotient(quot, pavg, qavg, pavg, qavg, dir)
        
        \param a_quot   4th  order face averaged quotient P / Q
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void faceQuotient(
            BoxData<T,C,MEM,D,E>& a_quotient,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            int a_dir);
}
#include "implem/Proto_OperatorImplem.H"
