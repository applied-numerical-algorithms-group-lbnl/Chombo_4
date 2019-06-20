#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_TIMER_H_
#define _CH_TIMER_H_



#include "Proto_Timer.H"

namespace Proto
{
#ifdef CH_TURN_OFF_TIMERS

#define PR_TIMER(name, tpointer) 
#define PR_TIME(name)   
#define PR_FLOPS(flops)
#define PR_TIMELEAF(name)                                                   
#define PR_TIMERS(name)  
#define PR_START(tpointer)
#define PR_STOP(tpointer)  
#define PR_TIMER_REPORT()
#define PR_TIMER_RESET() 
#define PR_TIMER_PRUNE(threshold)
#define PR_TIMER_SETFILE(filename)

#else

#define CH_TIMER(name, tpointer)                                        \
  const char* TimerTag_##tpointer = name ;                              \
  Proto::TraceTimer* tpointer = NULL ;                                  \
  if(Proto::TraceTimer::getTID()==0)                                    \
  {                                                                     \
    tpointer = Proto::TraceTimer::staticGetTimer(TimerTag_##tpointer) ; \
  }

#define CH_TIME(name)                                           \
  const char* TimerTagA = name ;                                \
  char CH_TimermutexA = 0;                                      \
  Proto::TraceTimer* tpointer = NULL;                        \
  if(Proto::TraceTimer::getTID()==0)                            \
  {                                                             \
    CH_tpointer = Proto::TraceTimer::staticGetTimer(TimerTagA); \
  }                                                             \
  Proto::AutoStart autostart(CH_tpointer, &PR_TimermutexA)

#define CH_FLOPS(flops)                         \
  if(tpointer)tpointer->addFlops(flops);

#define CH_TIMELEAF(name)                                       \
  const char* TimerTagA = name ;                                \
  Proto::TraceTimer* tpointer = NULL;                        \
  if(Proto::TraceTimer::getTID()==0)                            \
  {                                                             \
    tpointer = Proto::TraceTimer::staticGetTimer(TimerTagA); \
  }                                                             \
  Proto::AutoStartLeaf autostart(tpointer)

#define CH_TIMERS(name)                                                 \
  const char* TimerTagA = name ;                                        \
  char PR_TimermutexA = 0;                                              \
  char PR_Timermutex = 0;                                               \
  Proto::TraceTimer* tpointer = NULL;                                \
  if(Proto::TraceTimer::getTID()==0)                                    \
  {                                                                     \
    tpointer = Proto::TraceTimer::staticGetTimer(TimerTagA);         \
  }                                                                     \
  Proto::AutoStart autostart(tpointer, &PR_TimermutexA, &PR_Timermutex)


#define CH_START(tpointer)                      \
  if(Proto::TraceTimer::getTID()==0)            \
  {                                             \
    tpointer->start(&PR_Timermutex);            \
  }

#define CH_STOP(tpointer)                       \
    if(Proto::TraceTimer::getTID()==0)          \
    {                                           \
      tpointer->stop(&PR_Timermutex);           \
    }
//#define CH_STOPV(tpointer, val ) val = tpointer->stop(&PR_Timermutex)

#define CH_TIMER_REPORT() Proto::TraceTimer::staticReport()

#define CH_TIMER_RESET() Proto::TraceTimer::staticReset()

#define CH_TIMER_PRUNE(threshold) Proto::TraceTimer::staticPruneTimersParentChildPercent(threshold)

#define CH_TIMER_SETFILE(filename) Proto::TraceTimer::staticSetTimerFileName(filename);
#endif

}//namespace proto


#endif // CH_TIMER_H
