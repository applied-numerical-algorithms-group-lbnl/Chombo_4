#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_CH_TIMER_H_
#define _Chombo_CH_TIMER_H_

#include "base/Proto_Timer.H"
#include "Chombo_SPMD.H"



#define CH_TIMER(name, tpointer)   PR_TIMER(name, tpointer) 
#define CH_TIME(name)              PR_TIME(name)   
#define CH_FLOPS(flops)            PR_FLOPS(flops)
#define CH_TIMELEAF(name)          PR_TIMELEAF(name)                                                   
#define CH_TIMERS(name)            PR_TIMERS(name)  
#define CH_START(tpointer)         PR_START(tpointer)
#define CH_STOP(tpointer)          PR_STOP(tpointer)  
#define CH_TIMER_REPORT()          PR_TIMER_REPORT()
#define CH_TIMER_RESET()           PR_TIMER_RESET() 
#define CH_TIMER_PRUNE(threshold)  PR_TIMER_PRUNE(threshold)

inline void CH_TIMER_SETFILE(const char* a_prefix_string = "chombo.time.table")
{
#ifndef CH_MPI
  std::string name(a_prefix_string);
#else
  char char_str[1024];

  int outInterv = 1;
  char* charInterv = getenv("CH_OUTPUT_INTERVAL");
  if (charInterv != NULL)
  {
    outInterv =  atoi(charInterv);
  }
  int thisProc = CH4_SPMD::procID();
  if ((thisProc % outInterv) != 0)
  {
    sprintf(char_str,"/dev/null");
  }
  else
  {
    sprintf(char_str, "%s.%d",a_prefix_string, CH4_SPMD::procID());
  }

  std::string name(char_str);
#endif
  PR_TIMER_SETFILE(name.c_str());
}





#endif // CH_TIMER_H
