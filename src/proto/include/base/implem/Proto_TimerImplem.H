
#ifndef __PROTO_TIMER_IMPLEM_H__
#define __PROTO_TIMER_IMPLEM_H__
//using namespace std;
#ifndef PR_TURN_OFF_TIMERS

#ifdef PROTO_CUDA
#ifndef PROTO_HIP
#include <nvToolsExt.h>
#include <nvToolsExtCuda.h>
#else
#include <roctx.h>
#endif
#endif

namespace Proto
{



  struct elem
  {
    elem()
    {
      val=NULL;
      time=0;
    }

    const TraceTimer* val;
    unsigned long long int      time;
    elem(const TraceTimer* p, unsigned long long int t):val(p), time(t){;}
    inline bool operator < (const elem& rhs) const
    {
      if(val->isPruned()) return false;
      return time > rhs.time;
    }
    inline static void buildList(std::list<elem>& tlist, const TraceTimer& timer)
    {
      tlist.push_back(elem(&timer, timer.time()));
      const std::vector<TraceTimer*>& children = timer.children();
      for (unsigned int i=0; i<children.size(); i++)
      {
        buildList(tlist, *(children[i]));
      }
    }
  };




//  const char* currentTimer()
//  {
//    return TraceTimer::currentTimer();
//  }

  inline void writeOnExit()
  {
    std::vector<TraceTimer*>* vecrootptr = TraceTimer::getRootTimerPtr();
    if(vecrootptr->size() > 0)
    {
      (*vecrootptr)[0]->report(true);
    }
  }


  inline void TraceTimer::currentize() const
  {
    if(m_pruned) return;

    if(m_last_WCtime_stamp != 0){
      for(unsigned int i=0; i<m_children.size(); i++)
      {
        m_children[i]->currentize();
      }
      unsigned long long int current = PR_ticks();
      (unsigned long long int&)m_accumulated_WCtime += current - m_last_WCtime_stamp;
      (unsigned long long int&)m_last_WCtime_stamp = current;

    }
  }

  inline int computeRank(std::list<elem>& tracerlist, const TraceTimer& a_this) 
  {
    tracerlist.clear();
    elem::buildList(tracerlist, a_this);
    tracerlist.sort();
    int r=0;
    std::list<elem>::iterator it;
    for(it=tracerlist.begin(); it!=tracerlist.end(); ++it)
    {
      const elem& e = *it;
      //printf("%s %e %d\n",e.val->m_name, e.time, r);
      e.val->m_rank = r;
      ++r;
    }
    return r;
  }

  inline const TraceTimer* TraceTimer::activeChild() const
  {
    TraceTimer* child;
    for(unsigned int i=0; i<m_children.size(); i++)
    {
      child = m_children[i];
      if(child->m_last_WCtime_stamp != 0) return child;
    }
    return NULL;

  }

  inline const std::vector<TraceTimer*>& TraceTimer::children() const
  {
    return m_children;
  }

  inline
  void subReport(FILE* out, 
                 const char* header, 
                 unsigned long long int totalTime, 
                 std::list<elem>& tracerlist,
                 double         & secondspertick)
  {
    size_t length = strlen(header);
    fprintf(out, "=======================================================\n");
    unsigned long long int subTime = 0;
    std::list<elem>::iterator it;
    for(it=tracerlist.begin(); it!= tracerlist.end(); ++it)
    {
      const char* name = (*it).val->m_name;
      if(strncmp(header, name, length) == 0){
        if((*it).val->isPruned()){
          //fprintf(out, "             pruned  %s  \n", name);

        } else {
          unsigned long long int t = (*it).val->time();
          unsigned long long int f = (*it).val->m_flops;
          int rank = (*it).val->rank();
          subTime += t;
          fprintf(out, "  %8.5f %8lld  %s [%d]  %lld \n", t*secondspertick, (*it).val->m_count, name, rank, f);
        }
      }
    }
    if(subTime > 0)
      fprintf(out, "  %8.5f   %4.1f%%    Total\n", subTime*secondspertick, (double)subTime/totalTime*100.0);
  }

  inline void TraceTimer::report(bool a_closeAfter)
  {

    static std::list<elem> tracerlist(false);
    static double secondspertick = 0;



    static std::vector<TraceTimer*>* roots = getRootTimerPtr();
    TraceTimer& root = *((*roots)[0]); // in MThread code, loop over roots
    root.currentize();
    int numCounters = computeRank(tracerlist, root);

    double elapsedTime = TimerGetTimeStampWC() - root.zeroTime;
    unsigned long long int elapsedTicks = PR_ticks() - root.zeroTicks;
    secondspertick = elapsedTime/(double)elapsedTicks;


    static FILE* out = fopen(root.m_filename.c_str(), "w");
    static int reportCount = 0;
    fprintf(out, "-----------\nTimer report %d (%d timers)\n--------------\n",
            reportCount, numCounters);
    reportCount++;
    sumFlops(root);
    reportFullTree(out, root, root.m_accumulated_WCtime, 0, secondspertick); //uses recursion
    std::list<elem>::iterator it;
    for(it=tracerlist.begin(); it!=tracerlist.end(); ++it)
      reportOneTree(out, *((*it).val), secondspertick);
    subReport(out, "FORT_", root.m_accumulated_WCtime , tracerlist, secondspertick);
    subReport(out, "MPI_" , root.m_accumulated_WCtime , tracerlist, secondspertick);
    fflush(out);
    if(a_closeAfter) fclose(out);

  }
  inline void TraceTimer::reset()
  {
    static std::vector<TraceTimer*>* roots = getRootTimerPtr();
    TraceTimer& root = *((*roots)[0]);
    root.currentize();
    reset(root);
  }

  inline void TraceTimer::reset(TraceTimer& node)
  {
    node.m_count = 0;
    node.m_accumulated_WCtime = 0;
    for(unsigned int i=0; i<node.m_children.size(); i++)
    {
      reset(*(node.m_children[i]));
    }
  }


  inline void sorterHelper(const std::vector<TraceTimer*>& children, std::vector<int>& order)
  {
    int n = children.size();
    order.resize(n);
    for(int i=0; i<n; ++i) order[i]=i;
    bool swaps = true;
    while(swaps)
    {
      swaps = false;
      for(int i=0; i<n-1; ++i){
        if(children[order[i]]->time()  < children[order[i+1]]->time())
        {
          int tmp = order[i];
          order[i] = order[i+1];
          order[i+1] = tmp;
          swaps = true;
          break;
        }
      }
    }
  }


  inline void TraceTimer::sumFlops(TraceTimer& timer)
  {
    if(timer.m_pruned) return;

    for(unsigned int i=0; i<timer.m_children.size(); ++i){
      TraceTimer& child = *(timer.m_children[i]);
      sumFlops(child);
      timer.m_flops+=child.m_flops;
    }
  }

  inline
  void TraceTimer::reportFullTree(FILE* out, const TraceTimer& timer,
                                  unsigned long long int totalTime, int depth,
                                  double secondspertick)
  {
    if(timer.m_pruned) return;
    unsigned long long int time = timer.m_accumulated_WCtime;

    if(depth < 20){
      for(int i=0; i<depth; ++i) fprintf(out,"   ");
      double percent = ((double)time)/totalTime * 100.0;
      fprintf(out, "[%d] %s %.5f %4.1f%% %lld %lld %6.1f MFlops \n", timer.m_rank, timer.m_name, time*secondspertick, 
              percent, timer.m_count, timer.m_flops, timer.m_flops/(time*secondspertick*1000000));
    }
    std::vector<int> ordering;
    sorterHelper(timer.m_children, ordering);
    for(unsigned int i=0; i<timer.m_children.size(); ++i){
      reportFullTree(out, *(timer.m_children[ordering[i]]), totalTime, depth+1, secondspertick);
    }

  }
  inline
  void TraceTimer::reportOneTree(FILE* out, const TraceTimer& timer, double secondspertick)
  {
    if(timer.m_pruned) return;
    unsigned long long int time = timer.m_accumulated_WCtime;
    unsigned long long int subTime = 0;

    fprintf(out,"---------------------------------------------------------\n");

    fprintf(out,"[%d]%s %.5f %lld\n", timer.m_rank, timer.m_name, time*secondspertick, timer.m_count);
    const std::vector<TraceTimer*>& children = timer.m_children;
    std::vector<int> ordering;
    sorterHelper(children, ordering);
    for(unsigned int i=0; i<children.size(); ++i)
    {
      const TraceTimer& child = *(children[ordering[i]]);
      if(!child.m_pruned)
      {
        unsigned long long int childtime = child.m_accumulated_WCtime;
        if(childtime > 0){
          subTime += childtime;
          double percent = ((double)childtime) / time * 100.0;
          fprintf(out,"    %4.1f%% %7.4f %8lld %s [%d]\n",
                  percent, childtime*secondspertick, child.m_count, child.m_name, child.m_rank);
        }
      } else {
        fprintf(out,"           pruned           \n");
        i=children.size();
      }
    }
    if(time > 0 && children.size() > 0){
      double totalPercent = ((double)subTime)/ time * 100.0;
      fprintf(out, "    %4.1f%%                  Total \n", totalPercent);
    }

  }



// some compilers complain if there isn't at least 1 non-inlined
// function for every class.  These two are mostly to make those compilers happy



  inline bool AutoStart::active() { return true;}

  inline bool AutoStartLeaf::active() { return true;}


  inline TraceTimer* TraceTimer::getTimer(const char* name)
  {
    std::vector<TraceTimer*>* current = getCurrentTimerPtr();
    int thread_id = 0; // this line will change in MThread-aware code.
    TraceTimer* parent = (*current)[thread_id];
    if(parent->m_pruned) return parent;
    std::vector<TraceTimer*>& children = parent->m_children;
    unsigned int i=0;
    for(; i<children.size(); ++i){
      TraceTimer* timer =  children[i];
      if(timer->m_name == name) return timer;
    }
    TraceTimer* newTimer = new TraceTimer(name, parent, thread_id);
    children.push_back(newTimer);
    return newTimer;
  }


  inline void TraceTimer::prune()
  {
    unsigned int i=0;
    for(; i<m_children.size(); ++i){
      TraceTimer* timer =  m_children[i];
      timer->prune();
    }
    m_pruned = true;
  }

  inline void TraceTimer::start(char* mutex)
  {
    if(m_pruned) return;
# ifndef NDEBUG
    if(*mutex == 1) {
      char buf[1024];
      sprintf(buf, "double TraceTimer::start called: %s ",m_name);
      std::cerr<<buf;
      abort();
    }
# endif


    ++m_count;
    *mutex = 1;
    std::vector<TraceTimer*>* current = getCurrentTimerPtr();
#ifdef PROTO_CUDA
#ifndef PROTO_HIP
    nvtxRangePushA(m_name);
#else
   roctxRangePushA(m_name);
#endif
#endif
    (*current)[m_thread_id] = this;
    m_last_WCtime_stamp = PR_ticks();
#ifdef  PR_USE_MEMORY_TRACKING
    if(s_traceMemory)
    {
      long long int peak;
      overallMemoryUsage(m_last_Memory_Stamp, peak);
      if(peak > m_peak) m_peak = 0;
    }
#endif
  }
  inline unsigned long long int overflowLong()
  {
    unsigned long long int retval = (unsigned long long int)1<<50;
    return retval;
  }
  inline unsigned long long int TraceTimer::stop(char* mutex)
  {
    if(m_pruned) return 0;
    std::vector<TraceTimer*>* current = getCurrentTimerPtr();
#ifdef PROTO_CUDA
#ifndef PROTO_HIP
    nvtxRangePop();
#else
    roctxRangePop();
#endif
#endif
    
#ifndef NDEBUG
    if((*current)[0] != this)
    {
      char buf[1024];
      sprintf(buf, "TraceTimer::stop called while not parent: %s ",m_name);
      std::cerr<<buf;
      abort();
    }
#endif
    unsigned long long int diff = PR_ticks();
    diff -= m_last_WCtime_stamp;
    if(diff > overflowLong()) diff = 0;
    m_accumulated_WCtime += diff;

    m_last_WCtime_stamp=0;

    (*current)[m_thread_id] = m_parent;
    *mutex=0;
    return diff;
  }


  inline void TraceTimer::macroTest2()
  {
    PR_TIME("macroTest2");
  }
  inline void TraceTimer::macroTest()
  {
    PR_TIMERS("testy");
    PR_TIMER("billy", t1);
    PR_TIMER("sally", t2);
    PR_TIMER("renaldo", billy);
    PR_START(t1);
    PR_STOP(t1);
    PR_START(t2);
    PR_STOP(t2);
    PR_START(billy);
    PR_STOP(billy);
    PR_TIMER_REPORT();
    PR_TIMER_RESET();
    PR_TIMER_PRUNE(0.01);
  }

  inline void TraceTimer::PruneTimersParentChildPercent(double threshold, TraceTimer* parent)
  {
    if(parent->isPruned()) return;
    unsigned long long int time = parent->time();
    const std::vector<TraceTimer*>& children = parent->children();

    for(unsigned int i=0; i<children.size(); ++i)
    {
      TraceTimer* child = children[i];
      if(!child->isPruned())
      {
        unsigned long long int childtime = child->time();
        if(((double)childtime)/time < threshold) child->prune();
        else PruneTimersParentChildPercent(threshold, child);
      }

    }
  }


  inline void TraceTimer::PruneTimersParentChildPercent(double percent)
  {
    static std::vector<TraceTimer*>* roots = getRootTimerPtr();
    TraceTimer* root = (*roots)[0]; // in MThread code, loop over roots
    root->currentize();
    PruneTimersParentChildPercent(percent, root);
  }

} //namespace proto
#endif
#endif
