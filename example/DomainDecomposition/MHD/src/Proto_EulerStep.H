#ifndef _PROTO_EULERSTEP_H_
#define _PROTO_EULERSTEP_H_

namespace Proto {

///  Generic Euler Algorithm
/**  class X is the type of your state data. class F provides an evaluation
     of the X' for a specified state of X and the time.  If your X' is 
     determined entirely from the current state values then F can ignore.
     F must implement the function:
     @code
     void operator(X& k, double a_time, double a_dt, const X& a_state, const X& a_shift);
     @endcode
     F interprets the real state as a_state+a_shift

*/
template <class X, class F, class dX>
class EulerStep
{
public:
    void advance(double& a_time, double& a_dt, X& a_state);
protected:
    dX m_k; 
    F m_f;
};

template <class X, class F, class dX>
void EulerStep<X, F, dX>::advance(double& a_time, double& a_dt, X& a_state)
{
    m_k.init(a_state);                  // init must allocate storage, and initialize it to zero.
    m_f(m_k, a_time, a_dt, a_state);    // compute k1
	a_state.increment(m_k);	
}

} //end Proto namespace
  
#endif //end include guard
