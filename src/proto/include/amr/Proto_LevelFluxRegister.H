#pragma once
#ifndef _PROTO_FLUX_REGISTER_
#define _PROTO_FLUX_REGISTER_

#include "Proto_Register.H"
#include "Proto_BoxData.H"
#include "Proto_Copier.H"
namespace Proto
{
    // Forward Declarations
    template<typename T, unsigned int C, MemType MEM>
    class LevelFluxRegister;
    template<typename T, unsigned int C, MemType MEM>
    class FluxRegisterCopier;

    // =======================================================================
    // FLUX REGISTER COPIER OP

    /// Flux Register Copier Op
    /**
     */
    template<typename T, unsigned int C=1, MemType MEM = MEMTYPE_DEFAULT> 
    class FluxRegisterCopierOp
    {
        friend class FluxRegisterCopier<T,C,MEM>;
        private:
        LevelFluxRegister<T, C, MEM>* m_register;

        public:
        inline FluxRegisterCopierOp(){};
        inline FluxRegisterCopierOp(LevelFluxRegister<T, C, MEM>& a_register);

        /// Get Serialized Size
        inline int linearSize(const Box& a_bx,
                const DataIndex<BoxPartition>& a_index) const;

        /// Serial Write To Buffer
        inline void linearOut(void* a_buf, 
                const Box& a_bx,
                const DataIndex<BoxPartition>& a_index) const;

        /// Serial Read From Buffer
        inline void linearIn(void* a_buf, 
                const Box& a_bx,
                const DataIndex<BoxPartition>& a_index) const;

        /// Local Copy
        inline void localCopy(
                const Box& a_domain, const DataIndex<BoxPartition>& a_domainIndex,
                const Box& a_range,  const DataIndex<BoxPartition>& a_rangeIndex) const;

    }; // end class FluxRegisterCopierOp

    // =======================================================================
    // FLUX REGISTER COPIER

    /// Flux Register Copier
    /**
     */
    template<typename T, unsigned int C=1, MemType MEM=MEMTYPE_DEFAULT>
    class FluxRegisterCopier : public Copier<FluxRegisterCopierOp<T, C, MEM>, BoxPartition, MEM, MEM>
    {
        public:

        inline FluxRegisterCopier(){};

        /// Build Copier Motion Plan
        inline void buildMotionPlans(FluxRegisterCopierOp<T, C, MEM>& a_op);
    }; //end class FluxRegisterCopier

    // =======================================================================
    // FLUX REGISTER

    /// Flux Register
    /**
     */
    template<typename T, unsigned int C, MemType MEM>
    class LevelFluxRegister
    {
        friend class FluxRegisterCopierOp<T,C,MEM>;
        friend class FluxRegisterCopier<T,C,MEM>;

        public:
        /// Default constructor
        LevelFluxRegister(){};
        /// Constructor
        LevelFluxRegister(const DisjointBoxLayout& a_dblCoarse,
                const DisjointBoxLayout& a_dblFine,
                const Point& a_refRatio,
                std::array<T, DIM> a_dxCoarse);

        /// Define function
        void define(const DisjointBoxLayout& a_dblCoarse,
                const DisjointBoxLayout& a_dblFine,
                const Point& a_refRatio,
                std::array<T, DIM> a_dxCoarse);

        /// Increment Coarse
        /**
          Increments all the registers contained in the coarse patch a_di, direction a_dir. 
          Assumes that a_flux is defined on all the a_dir faces in the patch.
         */
        void incrementCoarse(
                const BoxData<T,C,MEM>& a_flux,
                const DataIndex<BoxPartition>& a_di,
                const T& a_weight,
                unsigned int a_dir);

        /// Increment Fine
        /**
          Increments all the registers contained in the fine patch a_di, direction a_dir. 
          Assumes that a_flux is defined on all the a_dir faces in the patch.
         */
        void incrementFine(const BoxData<T,C,MEM>& a_flux,
                const DataIndex<BoxPartition>& a_di,
                const T& a_weight,
                unsigned int a_dir); 

        /// Reflux
        /**
          Performs the reflux correction on input coarse data. 
         */
        void reflux(LevelBoxData<T,C,MEM>& a_coarseData, const T& a_weight);
        
        /// Reset
        void reset(T a_val = 0);
        
        /// compute the hash key for a coarse register.
        int key(const Box& a_bx,const DataIndex<BoxPartition>& a_di);

        void print() const;
      
        /// Access to layouts.
        DisjointBoxLayout crseLayout() const {return m_crseLayout;};
        DisjointBoxLayout fineLayout() const {return m_fineLayout;};
        protected:
        
        BoxData<T,C,MEM>& sourceData(const Box& a_bx,const DataIndex<BoxPartition>& di);
        BoxData<T,C,MEM>& destData(const Box& a_bx,const DataIndex<BoxPartition>& di);

        DisjointBoxLayout  m_crseLayout;
        DisjointBoxLayout  m_fineLayout;
        std::array<T, DIM> m_dxCrse;
        FluxRegisterCopier<T,C,MEM> m_copier;
        Point              m_refRatio;

        // Registers for coarse data
        vector<vector<Register<T,C,MEM> > >              m_crseRegisters;
        // Registers for fine data on coarse layout
        vector<vector<shared_ptr<BoxData<T,C,MEM> > > >  m_crseIncrement;
        // Registers for fine data
        vector<vector<Register<T,C,MEM> > >              m_fineRegisters; 
        // crseIndices is a map because there are in general multiple registers for each (dir, side)
        vector<unordered_map<unsigned int,unsigned int> >m_crseIndices;
        // fineIndices are not a map because there is at most one register for each (dir, side)
        vector<array<array<int,2>, DIM> >                m_fineIndices;
        array<array<Stencil<T>,2>,DIM>                   m_avgFlux;
    }; //end class LevelFluxRegister
#include "implem/Proto_LevelFluxRegisterImplem.H"
} // end namespace Proto
#endif // end include guard
