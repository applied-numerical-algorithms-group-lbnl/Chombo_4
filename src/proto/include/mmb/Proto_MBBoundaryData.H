#pragma once
#ifndef _PROTO_MB_BOUNDARY_DATA_
#define _PROTO_MB_BOUNDARY_DATA__

#include "Proto.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM>
    struct MBBoundaryData
    {
        MBBoundaryData();

        MBBoundaryData(
            unsigned int a_localBlock,
            unsigned int a_adjBlock,
            Box a_localBox,
            Box a_adjBox);
        
        inline void define(
            unsigned int a_localBlock,
            unsigned int a_adjBlock,
            Box a_localBox,
            Box a_adjBox);

        inline void print();

        unsigned int localBlock;
        unsigned int adjBlock;
        std::shared_ptr<BoxData<T, C, MEM>> localData;
        std::shared_ptr<BoxData<T, C, MEM>> adjData;
    };
    
    template<typename T, unsigned int C, MemType MEM>
    MBBoundaryData<T, C, MEM>::MBBoundaryData()
    {
        Box empty;
        define(0, 0, empty, empty);
    }

    template<typename T, unsigned int C, MemType MEM>
    MBBoundaryData<T, C, MEM>::MBBoundaryData(
            unsigned int a_localBlock,
            unsigned int a_adjBlock,
            Box a_localBox,
            Box a_adjBox)
    {
        define(a_localBlock, a_adjBlock, a_localBox, a_adjBox);
    }
    
    template<typename T, unsigned int C, MemType MEM>
    void MBBoundaryData<T, C, MEM>::define(
            unsigned int a_localBlock,
            unsigned int a_adjBlock,
            Box a_localBox,
            Box a_adjBox)
    {
        localBlock = a_localBlock;
        adjBlock = a_adjBlock;
        localData = std::make_shared<BoxData<T, C, MEM>>(a_localBox);
        adjData = std::make_shared<BoxData<T, C, MEM>>(a_adjBox);
    }
    
    template<typename T, unsigned int C, MemType MEM>
    void MBBoundaryData<T, C, MEM>::print()
    {
        std::cout << "MBBoundaryData | localBlock: " << localBlock;
        std::cout << " | localBoundary: " << localData->box();
        std::cout << " | adjBlock: " << adjBlock;
        std::cout << " | adjBoundary: " << adjData->box() << std::endl;
    }
} // end namespace Proto
#endif //end include guard
