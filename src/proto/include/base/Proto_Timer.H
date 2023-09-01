
#ifndef _PR_TIMER_H_
#define _PR_TIMER_H_

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include "implem/Proto_ClockTicks.H"
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <cstdio>
#include <cstring>
#include <sys/time.h>

using std::string;

/** TraceTimer class is a self-tracing code instrumentation system

    TraceTimer class is a self-tracing code instrumentation system
    for Chombo (or any other package really).  The user interface is specified
    by a small set of macros.  The usage model is that you just leave these
    timers in the code, for good.  Initially, your application will have 'main'
    and a few hewavy functions instrumented, and the lower level Chombo library
    instrumentation.  As your tool or application matures, it
    will garner a larger set of instrumentation giving clear views of your
    code performance.  After a routine has been
    cleverly and lovingly optimized, you leave in the timers, to spot when
    some later bug fix or *improvement* undoes your previous labors.

    \note
    You should never need to use or interact with the the classes TraceTimer or
    AutoStart.  Use the macros. They call the right functions and classes for you.

    The first macro is what people will use the most:
    \code
    PR_TIME("label");
    \endcode

    This is the simplest interface for timers.  you place this macro call in a function
    you wish to be timed.  It handles making the timer, calling 'start' when you
    enter the function, and calling 'stop' when you leave the function.  A good
    idea is to use a 'label' specific enough to be unambiguous without being
    overwhelming.  for instance:

    \code
    void AMRLevelPolytropicGas::define(AMRLevel*            a_coarserLevelPtr,
    const ProblemDomain& a_problemDomain,
    int                  a_level,
    int                  a_refRatio)
    {
    PR_TIME("AMRLevelPolytropicGas::define");
    .
    .
    }
    \endcode

    In this case, we have a class with many constructors and define functions that
    all funnel into a single general function.  We can just call this 'define' and
    not worry about naming/instrumenting all the different overloaded instances. If
    you slip up and use the same label twice, that is not a real problem, the two
    locations will be timed and tracked properly (even if one is a sibling or parent
    of the other). The only place it will make things a little harder is in the output
    where you might have the same name show up and look confusing.
    <br><br>
    In serial, you will see a file called <em>time.table</em> (in parallel, you will get a
    <em>time.table.n</em> (where n is the rank number) files).  If
    you want fewer files, you can do 
    setenv PR_OUTPUT_INTERVAL nproc 
    and it will only output every nproc processors time.table.n files
    (where n%nproc == 0).     I won't go into this file
    format.  It is kind of gprof-ish, with what I consider improvements.  The real
    benefit here is profiling that understands our Chombo context, a smaller information
    set to observe, and the fact that, so far in my testing, the timers have negligible
    impact on the run time or memory use of the code.
    <br><br>
    By default, Chombo compiles in the instructions for the timers wherever the macros
    appear.  If the compiler macro <b>PR_NTIMER</b> is defined, then all the PR_TIME* macros
    evaluate to empty expressions at compile time.

    \par So, you put some PR_TIME calls in your code and ran it, and nothing happened:
    Chombo looks for the environment variable <b>PR_TIMER</b>. If it is set to anything (even
    if it is set to 'false' or 'no' or whatever) then the timers will be active and
    reporting will happen.  If this environment variable is not set, then all the timers
    check a bool and return after doing nothing.
    \par
    One point of interest with using the environment variable: In parallel jobs using
    mpich, only processor 0 inherits the environment variables from the shell where
    you invoke 'mpirun', the rest read your .cshrc (.bashrc, etc.) file to get their
    environment.  To time all your processes, you need to make sure the <b>PR_TIMER</b>
    environment variable gets to all your processes.



    \par Auto hierarchy:
    The timers automatically figure out their parent/child relationships.  They
    also can be placed in template code.  This has some consequences.  First,
    if you have a low level function instrumented that has no timers near it in
    the code call stack, you will see it show up as a child of a high level timer.
    the root timer "main" will catch all orphaned timers.  So, even though you
    might make no call to, say, 'exchange' in your 'main' function, you might
    very well call a function, that calls a function, that calls 'exchange'. Since
    no code in between was instrumented, this exchange is accounted for at 'main'.
    This might look strange, but it should prove very powerful. An expensive orphan
    is exactly where you should consider some more timers, or reconsidering code
    design.

    \par
    For performance reasons, child timers have only one parent.  As a consequence
    each PR_TIME("label") label can show up at multiple places in your output. Each
    instance has it's own timer.  So, each path through the call graph that arrives
    at a low-level function has a unique lineage, with it's own counter and time.
    Thus, I can instrument LevelData::copyTo once, but copyTo can appear in many
    places in the time.table file.


    The next level up in complexity is the set of *four* macros for when you want
    sub-function resolution in your timers. For instance, in a really huge function
    that you have not figured out how to re-factor, or built with lots of bad cut n paste
    code 're-use'.
    \code
    PR_TIMERS("parent");
    PR_TIMER("child1", t1);
    PR_TIMER("child2", t2);
    PR_START(t1);
    PR_STOP(t1);
    PR_START(t2);
    PR_STOP(t2);
    PR_START(t1);
    PR_STOP(t1);
    \endcode

    PR_TIMERS has the same semantic as PR_TIME, except that you can declare an
    arbitrary number of children after it in the same function scope.  The
    children here do not autostart and autostop, you have to tell them where to
    start and stop timing.  The children can themselves be parents for timers
    in called functions, of course. The children obey a set of mutual exclusions. The
    following generate run time errors:
    - double start called
    - double stop called
    - start called when another child is also started
    - you leave the function with a child not stopped

    the following will generate compile time errors:
    - more than one PR_TIME macro in a function
    - invoking PR_TIMER("child", t) without having first invoked PR_TIMERS
    - re-using the timer handle ie. PR_TIMER("bobby", t1); PR_TIMER("sally", t1)
    - mixing PR_TIME macro with PR_TIMER
    - mixing PR_TIME macro with PR_TIMERS

    You do not have to put any calls in your main routine to activate the clocks
    or generate a report at completion, this is handled with static iniitalization
    and an atexit function.
    <br><br>
    There is a larger argument of manual instrumentation being counter to good development.
    Profiling the code is supposed to tell you where to expend your optimization effort.
    Manual instrumentation opens the door to people wasting time *assuming* what parts of the
    code are going to take up lots of time and instrumenting them, before seeing any real
    performance data.  Good judgement is needed.  We have a body of knowledge about Chombo
    that will inform us about a good minimal first set of functions to instrument.
*/
namespace Proto 
{
#ifndef PR_TURN_OFF_TIMERS
  inline double TimerGetTimeStampWC()
  {

    struct timeval tv;   //  Values from call to gettimeofday
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return((double)tv.tv_sec + 0.000001 * (double)tv.tv_usec);

  }

  class TraceTimer
  {
  public:
    virtual ~TraceTimer()
    {
      for(unsigned int i=0; i<m_children.size(); ++i)
      {
        delete m_children[i];
      }
    }

    TraceTimer()
    {
      m_filename = string("/dev/null");
    }
    inline void start(char* mutex);
    inline unsigned long long int stop(char* mutex);
    inline void report(bool a_closeAfter=false);
    inline void reset();


    inline void leafStart();
    inline void leafStop();

    unsigned long long int time() const
    {
      return m_accumulated_WCtime;
    }

    int rank() const
    {
      return m_rank;
    }
    long long int count() const
    {
      return m_count;
    }

    void prune();
    bool isPruned() const {return m_pruned;}


    ///need to set the filename so it is not /dev/null
    void setTimerFileName(string a_filename)
    {
      
      std::vector<TraceTimer*>* troots = getRootTimerPtr();
      (*troots)[0]->m_filename = a_filename;
    }

    inline TraceTimer* getTimer(const char* name); // don't use
    inline const std::vector<TraceTimer*>& children() const ;//don't use.

    inline void PruneTimersParentChildPercent(double percent);
    inline void sampleMemUsage() ;

    inline void addFlops(long long int flops) {m_flops+=flops;}


    const char*        m_name;
    mutable int        m_rank;
    int tid = 0;
    int                m_thread_id;
    long long int      m_flops;
    long long int      m_count;
//#pragma omp atomic

    static std::vector<TraceTimer*>* getRootTimerPtr()
    {
      static std::vector<TraceTimer*>* retval = new std::vector<TraceTimer*>();
      static bool initialized = false;
      if(!initialized)
      {
        std::vector<TraceTimer*>* current = getCurrentTimerPtr();
        //this is the code from initialize
        const char* rootName = "root";
        TraceTimer* rootTimer = new TraceTimer(rootName, NULL, 0);
        rootTimer->m_thread_id = 0;
        char mutex = 0;

        retval->resize(1);
        (*retval)[0] = rootTimer;

        (*current).resize(1);
        (*current)[0]=rootTimer;

        (*retval)[0]->tid = 0;
        rootTimer->start(&mutex);
        rootTimer->zeroTime = TimerGetTimeStampWC();
        rootTimer->zeroTicks = PR_ticks();
      }
      initialized = true;
      return retval;
    }


    static std::vector<TraceTimer*>* getCurrentTimerPtr()
    {
      static std::vector<TraceTimer*>* retval = new std::vector<TraceTimer*>(1);
      return retval;
    }


    //oh the evil crap we have to do to avoid static initialization
    static void staticReport()
    {
      std::vector<TraceTimer*>* vecptr = getRootTimerPtr();
      if(vecptr->size() > 0)
      {
        (*vecptr)[0]->report();
      }
    }

    //oh the evil crap we have to do to avoid static initialization
    static void staticReset()
    {
      std::vector<TraceTimer*>* vecptr = getRootTimerPtr();
      if(vecptr->size() > 0)
      {
        (*vecptr)[0]->reset();
      }
    }

    //oh the evil crap we have to do to avoid static initialization
    static void staticPruneTimersParentChildPercent(double percent)
    {
      std::vector<TraceTimer*>* vecptr = getRootTimerPtr();
      if(vecptr->size() > 0)
      {
        (*vecptr)[0]->PruneTimersParentChildPercent(percent);
      }
    }

    //oh the evil crap we have to do to avoid static initialization
    static void staticSetTimerFileName(string a_filename)
    {
      std::vector<TraceTimer*>* vecptr = getRootTimerPtr();
      if(vecptr->size() > 0)
      {
        (*vecptr)[0]->setTimerFileName(a_filename);
      }
    }

    //more evil crap
    static int getTID()
    {
      int retval = -1;
      std::vector<TraceTimer*>* vecptr = getRootTimerPtr();
      if(vecptr->size() > 0)
      {
        retval =  (*vecptr)[0]->tid;
      }
      return retval;
    }
    
    //and all for the lack of a perfect reasonable language feature 
    static TraceTimer* staticGetTimer(const char* name)
    {
      TraceTimer* retval = NULL;
      std::vector<TraceTimer*>* vecptr = getRootTimerPtr();
      if(vecptr->size() > 0)
      {
        retval =  (*vecptr)[0]->getTimer(name);
      }
      return retval;
    }

  private:
    string m_filename; //only meaningful at root

    TraceTimer(const char* a_name, TraceTimer* parent, int thread_id)
    //:m_pruned(false), m_parent(parent), m_name(a_name), m_count(0),
    // m_accumulated_WCtime(0),m_last_WCtime_stamp(0), m_thread_id(thread_id),
    // m_flops(0)
    {
        m_pruned = false;
        m_parent = parent;
        m_name = a_name;
        m_count = 0;
        m_accumulated_WCtime = 0;
        m_last_WCtime_stamp = 0;
        m_thread_id = thread_id;
        m_flops = 0;
    }



    bool               m_pruned;
    TraceTimer*        m_parent;
    std::vector<TraceTimer*> m_children;

    unsigned long long int      m_accumulated_WCtime;
    unsigned long long int      m_last_WCtime_stamp;
    
    double zeroTime = 0;
    unsigned long long int zeroTicks = 0;



    void reportTree(FILE* out, const TraceTimer& node, int depth);
    const TraceTimer* activeChild() const;

    int m_depth;
    inline void macroTest();
    inline void macroTest2();

    inline void currentize() const ;

    inline void reportFullTree(FILE* out, const TraceTimer& timer,
                        unsigned long long int totalTime, int depth, double secondspertick);
    inline void reportOneTree(FILE* out, const TraceTimer& timer, double secondspertick);


    inline void reset(TraceTimer& timer);
    inline void PruneTimersParentChildPercent(double threshold, TraceTimer* parent);
    inline void sumFlops(TraceTimer& timer);
  };


  class AutoStartLeaf
  {
  public:
    AutoStartLeaf(TraceTimer* a_timer):m_timer(a_timer) {if(a_timer)a_timer->leafStart();}
    ~AutoStartLeaf(){if(m_timer)m_timer->leafStop();}
    bool active();
  private:
    TraceTimer* m_timer;
  };

  class AutoStart
  {
  public:
    AutoStart(TraceTimer* a_timer, char* mutex, char* childMutex)
      :m_mutex(mutex), m_childMutex(childMutex),m_timer(a_timer)
    {if(a_timer)a_timer->start(mutex);}
    AutoStart(TraceTimer* a_timer, char* mutex)
      :m_mutex(mutex), m_childMutex((char*)(getOKPtr()) ),m_timer(a_timer)
    {
      if(a_timer)a_timer->start(mutex);
    }
    inline ~AutoStart();
    bool active();
  private:
    static char* getOKPtr()
    {
      static bool init = false;
      static char* retval = new char[1024];
      if(!init)
      {
        sprintf(retval, "0");
      }
      init = true;
      return retval;
    }
    char* m_mutex;
    char* m_childMutex;
    TraceTimer* m_timer;
    static const char ok = 0;

  };

  inline
  AutoStart::~AutoStart()
  {
#  ifndef NDEBUG
    if(*m_childMutex == 1)
    {
      std::cerr<<" stop called in timer unexpectedly";
      abort();
    }
#  endif
    if(m_timer)m_timer->stop(m_mutex);
  }


  inline void TraceTimer::leafStart()
  {
    if(m_pruned) return;
    ++m_count;
    m_last_WCtime_stamp = PR_ticks();
  }

  inline void TraceTimer::leafStop()
  {
    if(m_pruned) return;
    m_accumulated_WCtime +=  PR_ticks() - m_last_WCtime_stamp;
    m_last_WCtime_stamp=0;
  }
#endif

#ifdef PR_TURN_OFF_TIMERS

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

#define PR_TIMER(name, tpointer)                                        \
  const char* TimerTag_##tpointer = name ;                              \
  ::Proto::TraceTimer* tpointer = NULL ;                                  \
  if(::Proto::TraceTimer::getTID()==0)                                    \
  {                                                                     \
    tpointer = ::Proto::TraceTimer::staticGetTimer(TimerTag_##tpointer) ;       \
  }

#define PR_TIME(name)                                           \
  const char* TimerTagA = name ;                                \
  char PR_TimermutexA = 0;                                      \
  ::Proto::TraceTimer* PR_tpointer = NULL;                      \
  if(::Proto::TraceTimer::getTID()==0)                          \
  {                                                             \
    PR_tpointer = ::Proto::TraceTimer::staticGetTimer(TimerTagA);       \
  }                                                             \
  ::Proto::AutoStart autostart(PR_tpointer, &PR_TimermutexA)

#define PR_FLOPS(flops)                         \
  if(PR_tpointer)PR_tpointer->addFlops(flops);

#define PR_TIMELEAF(name)                                       \
  const char* TimerTagA = name ;                                \
  ::Proto::TraceTimer* PR_tpointer = NULL;                        \
  if(::Proto::TraceTimer::getTID()==0)                            \
  {                                                             \
    PR_tpointer = ::Proto::TraceTimer::staticGetTimer(TimerTagA);       \
  }                                                             \
  ::Proto::AutoStartLeaf autostart(PR_tpointer)

#define PR_TIMERS(name)                                                 \
  const char* TimerTagA = name ;                                        \
  char PR_TimermutexA = 0;                                              \
  char PR_Timermutex = 0;                                               \
  ::Proto::TraceTimer* PR_tpointer = NULL;                                \
  if(::Proto::TraceTimer::getTID()==0)                                    \
  {                                                                     \
    PR_tpointer = ::Proto::TraceTimer::staticGetTimer(TimerTagA);               \
  }                                                                     \
  ::Proto::AutoStart autostart(PR_tpointer, &PR_TimermutexA, &PR_Timermutex)


#define PR_START(tpointer)                      \
  if(::Proto::TraceTimer::getTID()==0)            \
  {                                             \
  tpointer->start(&PR_Timermutex);              \
}

#define PR_STOP(tpointer)                       \
  if(::Proto::TraceTimer::getTID()==0)            \
  {                                             \
    tpointer->stop(&PR_Timermutex);             \
  }
  //#define PR_STOPV(tpointer, val ) val = tpointer->stop(&PR_Timermutex)

#define PR_TIMER_REPORT() ::Proto::TraceTimer::staticReport()

#define PR_TIMER_RESET() ::Proto::TraceTimer::staticReset()

#define PR_TIMER_PRUNE(threshold) ::Proto::TraceTimer::staticPruneTimersParentChildPercent(threshold)

#define PR_TIMER_SETFILE(filename) ::Proto::TraceTimer::staticSetTimerFileName(filename);
#endif
}//namespace proto

#include "implem/Proto_TimerImplem.H"
#endif // PR_TIMER_H
