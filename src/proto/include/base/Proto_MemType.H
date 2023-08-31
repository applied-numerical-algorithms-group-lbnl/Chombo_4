#pragma once
#ifndef _PROTO_MEMTYPE_H_
#define _PROTO_MEMTYPE_H_

namespace Proto
{
  enum MemType {INVALID=0,HOST=1, DEVICE=2,BOTH=3};
  // for better error reporting
  inline std::string parseMemType(MemType a_mem)
  {
      switch (a_mem)
      {
          case INVALID: return "INVALID"; break;
          case HOST: return "HOST"; break;
          case DEVICE: return "DEVICE"; break;
          case BOTH: return "BOTH"; break;
      }
      return "INVALID";
  }
#ifndef MEMTYPE_DEFAULT
#ifdef PROTO_CUDA
#define MEMTYPE_DEFAULT MemType::DEVICE
#else
#define MEMTYPE_DEFAULT MemType::HOST
#endif
#endif
}

#endif
