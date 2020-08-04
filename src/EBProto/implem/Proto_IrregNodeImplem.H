namespace Proto
{

  ///
  template<int order>
  void 
  IrregNode<order>::
  setMomentsToZero()
  {
    m_volumeMoments.setToZero();
    m_EBMoments.setToZero();
    for(int idir  = 0; idir < DIM; idir++)
    {
      m_normalPartialDeriv[idir].setToZero();
    }
    for(int iface  = 0; iface < 2*DIM; iface++)
    {
      m_faceMoments[iface].setToZero();
    }
  }

  ///
  template <int order>
  void 
  IrregNode<order>::
  setMomentsToRegular(const Real& a_dx)
  {
    //no eb for regular
    m_EBMoments.setToZero();
    for(int idir  = 0; idir < DIM; idir++)
    {
      m_normalPartialDeriv[idir].setToZero();
    }
    m_volumeMoments.setRegular(a_dx);
    for(int iface  = 0; iface < 2*DIM; iface++)
    {
      m_faceMoments[iface].setRegular(a_dx);
    }
  }
}
