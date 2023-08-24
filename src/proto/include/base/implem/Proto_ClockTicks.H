/*                       __      
       ____  _________  / /_____ 
      / __ \/ ___/ __ \/ __/ __ \
     / /_/ / /  / /_/ / /_/ /_/ /
    / .___/_/   \____/\__/\____/ 
   /_/                           
*/
#ifndef _PR_CLOCKTICKS_H_
#define _PR_CLOCKTICKS_H_

// (dfm 4/29/08)note that NamespaceHeader.H gets included 3x, once for
// each "if"...

#include <chrono>

inline unsigned long long int PR_ticks()
{

  using namespace std::chrono;
  return steady_clock::now().time_since_epoch().count();

}
#define PR_TICKS



#endif
