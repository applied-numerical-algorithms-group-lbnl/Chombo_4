#pragma once
#ifndef _PROTO_BoxDataOp_H_
#define _PROTO_BoxDataOp_H_


using std::shared_ptr;
namespace Proto 
{

    enum BoxDataOp
    {
        Invalid = -1,
        Copy    = 0,
        Add     = 1,
        Assign,
        Subtract, 
        Multiply, 
        Divide, 
        NUM_OP_TYPES
    };
} //end Proto namespace
#endif //end include guard
