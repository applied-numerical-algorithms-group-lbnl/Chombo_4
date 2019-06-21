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
#define CH_TIMER_SETFILE(filename) PR_TIMER_SETFILE(filename)






#endif // CH_TIMER_H
