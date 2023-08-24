
template<typename T>
InterpStencilTP<T>::InterpStencilTP(int a_order, int a_refRatio)
{
    define(a_order, a_refRatio);
}

template<typename T>
void InterpStencilTP<T>::define(int a_order, int a_refRatio)
{
    if (a_refRatio  == 2)
    {
      array<array<T,5> , 2 > coefs;
      if (a_order == 4)
        { 
          coefs =
            {{
                {-3.0/128,22.0/128,128.0/128,-22.0/128,3.0/128},
                {3.0/128,-22.0/128,128.0/128,22.0/128,-3.0/128}
              }};
        }
      else if (a_order== 5)
        {
          coefs =
            {{
                {-3.0/128,22.0/128,128.0/128,-22.0/128,3.0/128},
                {3.0/128,-22.0/128,128.0/128,22.0/128,-3.0/128}
              }};
          // PC:this is the same as the 4th order stencil! 
        }
      else
        {
          MayDay<void>::Abort(
                              "InterpStencilTP::Build | Error: Invalid order. Must be 4 or 5.");  
        }
      for (int dir = 0; dir < DIM; dir++)
        {
            Point refRatio = Point::Ones();
            refRatio[dir] = a_refRatio;
            m_interp[dir].define(refRatio);
            for (int k = 0; k < 2; k++)
            {
                m_interp[dir](Point::Basis(dir, k)) =
                    build1D5pt(coefs[k], dir);
            }
            m_interp[dir].close();
        }
    } else if (a_refRatio == 4)
    {
      array<array<T,5> , 4 > coefs;
      if (a_order == 4)
        {
          coefs =
            {{
                {-35.0/1024,294.0/1024,960.0/1024,-230.0/1024,35.0/1024},
                {-13.0/1024,58.0/1024,1088.0/1024,-122.0/1024,13.0/1024},
                {13.0/1024,-122.0/1024,1088.0/1024,58.0/1024,-13.0/1024},
                {35.0/1024,-230.0/1024,960.0/1024,294.0/1024,-35.0/1024}
              }};
        }
      else if (a_order == 5)
        {
          coefs =
            {{
                {-77.0/2048,616.0/2048,1878.0/2048,-432.0/2048,63.0/2048},
                {-19.0/2048,88.0/2048,2218.0/2048,-272.0/2048,33.0/2048},
                {33.0/2048,-272.0/2048,2218.0/2048,88.0/2048,-19.0/2048},
                {63.0/2048,-432.0/2048,1878.0/2048,616.0/2048,-77.0/2048}                
              }};
        }
        else
          {
            MayDay<void>::Abort(
                                "InterpStencilTP::Build | Error: Invalid order. Must be 4 or 5.");  
          }
        for (int dir = 0; dir < DIM; dir++)
        {
            Point refRatio = Point::Ones();
            refRatio[dir] = a_refRatio;
            m_interp[dir].define(refRatio);
            for (int k = 0; k < 4; k++)
            {
                m_interp[dir](Point::Basis(dir, k)) =
                    build1D5pt(coefs[k], dir);
            }
            m_interp[dir].close();
        }
    } else {
        MayDay<void>::Abort(
            "InterpStencilTP::Build | Error: Invalid refinement ratio. Must be 2 or 4.");
    }
}

template<typename T>
Stencil<T> InterpStencilTP<T>::build1D5pt(
    array<double,5>& a_coefs, int a_dir)
{
    return
          a_coefs[0]*Shift::Basis(a_dir,-2) +
          a_coefs[1]*Shift::Basis(a_dir,-1) +
          a_coefs[2]*Shift::Zeros() +
          a_coefs[3]*Shift::Basis(a_dir,1) +
          a_coefs[4]*Shift::Basis(a_dir,2);
}

template<typename T>
Point InterpStencilTP<T>::ghost() const
{
    Point span_0 = m_interp[0].ghost();
    for (int jj = 1; jj < DIM; jj++)
    {
        Point span = m_interp[jj].ghost();
        for (int dir = 0; dir < DIM; dir++)
        {
            span_0[dir] = max(span[dir], span_0[dir]);
        }
    }
    return span_0;
}

#if 0
template<typename T>
template<unsigned int C, MemType MEM>
void InterpStencilTP<T>::apply(
    BoxData<T, C, MEM>& a_output,
    BoxData<T, C, MEM>& a_input,
    T                   a_scale)
{
    std::array<BoxData<T, C, MEM>, DIM+1> stage;
    stage[0].define(a_input.box());
    a_input.copyTo(stage[0]);
    for (int dir = 0; dir < DIM; dir++)
    {
        std::cout << "stage " << dir << " box: " << stage[dir].box() << std::endl;
        stage[dir+1] = m_interp[dir](stage[dir]);
    }
    std::cout << "stage " << DIM << " box: " << stage[DIM].box() << std::endl;
    stage[DIM].copyTo(a_output);
}
#endif
template<typename T>
template<unsigned int C, MemType MEM>
void InterpStencilTP<T>::apply(
    BoxData<T, C, MEM>& a_output,
    BoxData<T, C, MEM>& a_input,
    T                   a_scale)
{
  PR_TIMERS("InterpStencilTP::apply");
    std::array<BoxData<T, C, MEM>, DIM+1> stage;
    // cout << "before first TP stencil" << endl;
    stage[1] = m_interp[0](a_input);
    for (int dir = 1; dir < DIM; dir++)
    {
      //std::cout << "stage " << dir << " box: " << stage[dir].box() << std::endl;
        stage[dir+1] = m_interp[dir](stage[dir]);
    }
    //std::cout << "stage " << DIM << " box: " << stage[DIM].box() << std::endl;
    Stencil<T> idOp = 1.0*Shift::Zeros();
    a_output |= idOp(stage[DIM]);
}
