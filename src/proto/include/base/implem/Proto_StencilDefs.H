#ifndef PRroto_stencildefs
#define PRroto_stencildefs

//========================================================================
// Stencil Definitions
//======================
template<typename T>
Stencil<T> Stencil<T>::Laplacian()
{
    Stencil<T> ret = ((T)(-2*DIM))*Shift::Zeros();
    for (int dir = 0; dir < DIM; dir++)
    {
      ret += ((T)(1))*Shift::Basis(dir, 1);
      ret += ((T)(1))*Shift::Basis(dir,-1);
    }
    return ret;
}
#if DIM == 2
template<typename T>
Stencil<T> Stencil<T>::Laplacian_9()
{
  T factor = -20./6.;
    Stencil<T> ret = factor*Shift::Zeros();
    Box K = Box::Cube(3).shift(Point::Ones(-1));
    T a = 1.0/6.0;
    T b = 2.0/3.0;
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        Point p = *iter;
        int pointSum = abs(p[0]) + abs(p[1]);
        if (pointSum == 1)
        {
            ret += b*Shift(p); 
        } else if (pointSum == 2)
        {
            ret += a*Shift(p); 
        }
    }
    return ret;
}
#elif DIM == 3
template<typename T>
Stencil<T> Stencil<T>::Laplacian_19()
{
    Stencil<T> ret = ((T)(-24.0))*Shift::Zeros();
    Box K = Box::Cube(3).shift(Point::Ones(-1));
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        int sum = (*iter).sumAbs();
        if (sum == 1) {
            ret += ((T)(2.0))*Shift(*iter);
        } else if (sum == 2) {
            ret += ((T)(1.0))*Shift(*iter);
        }
    }
    ret *= (1.0/6.0);
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::Laplacian_27()
{
    Stencil<T> ret = ((T)(-64.0/15.0))*Shift::Zeros();
    Box K = Box::Cube(3).shift(Point::Ones(-1));
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        Point p = *iter;
        int pointSum = p.sumAbs();
        if (pointSum == 1)
        {
            ret += ((T)(7.0/15.0))*Shift(p);
        } else if (pointSum == 2) {
            ret += ((T)(1.0/10.0))*Shift(p);
        } else if (pointSum == 3) {
            ret += ((T)(1.0/30.0))*Shift(p);
        }
    }
    return ret;
}
#endif
template<typename T>
Stencil<T> Stencil<T>::LaplacianFace(int a_dir, int a_order)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::Laplacian(int a_order, int a_dir) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    Stencil<T> ret;
    switch (a_order)
    {
        case 2:
        {
            ret = ((T)(-2*(DIM-1)))*Shift(Point::Zeros());
            for (int dir = 0; dir < DIM; dir++)
            {
                if (a_dir != dir)
                {
                    ret += ((T)(1))*Shift(Point::Basis(dir));
                    ret += ((T)(1))*Shift(Point::Basis(dir)*(-1));
                }
            }
            break;
        }
        default:
        {
            std::cout << "Stencil::Laplacian(int a_order, int a_dir) invalid for a_order = ";
            std::cout << a_order << ". This configuration is not yet implemented as a built-in Stencil." << std::endl;
            abort();
            break;
        }
    }
    return ret;
}

template<typename T>
Stencil<T>
Stencil<T>::CellToFace(int a_dir, Side::LoHiSide a_side, int a_order)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
            "Stencil::CellToFace(int a_dir, int a_order) invalid for a_dir = %i.\
            a_dir should be in [0,DIM=%i).",a_dir, DIM);
    Stencil<T> ret;
    switch (a_order)
    {
        case 4:
        {
            T coefs[4] = {-1,7,7,-1};
            ret = 
                coefs[0]*Shift::Basis(a_dir, -2) + coefs[3]*Shift::Basis(a_dir, 1) +
                coefs[1]*Shift::Basis(a_dir, -1)  + coefs[2]*Shift::Zeros();
            ret *= (T)(1.0/12.0);
            break;
        }
        case 5:
        {
          T coefs[5] = {2,-13,47,27,-3};
            ret =  
              coefs[0]*Shift::Basis(a_dir, -3) + coefs[4]*Shift::Basis(a_dir, 1) +
              coefs[1]*Shift::Basis(a_dir, - 2) + coefs[3]*Shift::Zeros() +
              coefs[2]*Shift::Basis(a_dir, -1);
            ret *= (T)(1/60.0);
            break;
        }
        default:
        {
            std::abort();
            break;
        }
    }
    if (a_side == Side::Hi)
    {
        ret.invert(a_dir);
    }
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::DiffCellToFace(int a_dir, Side::LoHiSide a_side, int a_order)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::DiffCellToFace(int a_dir, int a_order) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    Stencil<T> ret;
    switch (a_order)
    {
        case 4:
        {
            T coefs[4] = {1, -15, 15, -1};
            ret = coefs[0]*Shift::Basis(a_dir, -2) +
                  coefs[1]*Shift::Basis(a_dir, -1) +
                  coefs[2]*Shift::Zeros() +
                  coefs[3]*Shift::Basis(a_dir, 1);
            ret *= (T)(1.0/12.0);
            break;
        }
        default:
        {
            std::cout << "Stencil::DiffCellToFace(int a_dir, int a_order) invalid for a_order = ";
            std::cout << a_order << ". This configuration is not yet implemented as a built-in Stencil." << std::endl;
            std::abort();
            break;
        }
    }
    if (a_side == Side::Hi)
    {
        ret.invert(a_dir);
    }
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::CellToFaceL(int a_dir, Side::LoHiSide a_side, int a_order)
{
    auto S = CellToFace(a_dir, a_side, a_order);
    if (a_order % 2 == 0)
    {
        return S;
    } else {
        if (a_side == Side::Lo)
        {
            return CellToFace(a_dir, Side::Lo, a_order);
        } else {
            auto S = CellToFace(a_dir, Side::Lo, a_order);
            return S * (1.0*Shift::Basis(a_dir, +1));
        }
    }
}

template<typename T>
Stencil<T> Stencil<T>::CellToFaceH(int a_dir, Side::LoHiSide a_side, int a_order)
{
    auto S = CellToFace(a_dir, a_side, a_order);
    if (a_order % 2 == 0)
    {
        return S;
    } else {
        if (a_side == Side::Hi)
        {
            return CellToFace(a_dir, Side::Hi, a_order);
        } else {
            auto S = CellToFace(a_dir, Side::Hi, a_order);
            return S * (1.0*Shift(Point::Basis(a_dir, -1)));
        }
    }
}
template<typename T>
Stencil<T> Stencil<T>::AvgDown(int a_refRatio)
{
    return Stencil<T>::AvgDown(Point::Ones(a_refRatio));
}

template<typename T>
Stencil<T> Stencil<T>::AvgDown(Point a_refRatio)
{
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_refRatio[ii] > 0,
                "Stencil::AvgDown(Ppoint a_refRatio) invalid for a_refRatio = %i. a_refRatio should be positive.",
                a_refRatio[ii]);
    }
    
    Stencil<T> ret;
    Box kernel(a_refRatio);
    T coef = 1.0 / kernel.size();
    for (auto iter = kernel.begin(); iter != kernel.end(); ++iter)
    {
        ret += coef*Shift(*iter);
    }
    ret.srcRatio() = a_refRatio;
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::Sum(int a_refRatio)
{
    return Stencil<T>::Sum(Point::Ones(a_refRatio));
}

template<typename T>
Stencil<T> Stencil<T>::Sum(Point a_refRatio)
{
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_refRatio[ii] > 0,
                "Stencil::AvgDown(Ppoint a_refRatio) invalid for a_refRatio = %i. a_refRatio should be positive.",
                a_refRatio[ii]);
    }
    
    Stencil<T> ret;
    Box kernel(a_refRatio);
    T coef = 1.0;
    for (auto iter = kernel.begin(); iter != kernel.end(); ++iter)
    {
        ret += coef*Shift(*iter);
    }
    ret.srcRatio() = a_refRatio;
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::AvgDownFace(int a_dir, Side::LoHiSide a_side, int a_refRatio)
{
    return Stencil<T>::AvgDownFace(a_dir, a_side, Point::Ones(a_refRatio));
}

template<typename T>
Stencil<T> Stencil<T>::AvgDownFace(int a_dir, Side::LoHiSide a_side, Point a_refRatio)
{
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_refRatio[ii] > 0,
                "Stencil::AvgDown(Ppoint a_refRatio) invalid for a_refRatio = %i. a_refRatio should be positive.",
                a_refRatio[ii]);
    }
    
    Stencil<T> ret;
    Box kernel(a_refRatio);
    kernel = kernel.face(a_dir, a_side);
    T coef = 1.0 / kernel.size();
    for (auto iter = kernel.begin(); iter != kernel.end(); ++iter)
    {
        ret += coef*Shift(*iter);
    }
    ret.srcRatio() = a_refRatio;
    return ret;
}
template<typename T>
Stencil<T> Stencil<T>::FluxDivergence(int a_dir)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::FluxDivergence(int a_dir) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    return ((T)(-1))*Shift(Point::Zeros()) + ((T)(1.0))*Shift(Point::Basis(a_dir));
}

template<typename T>
Stencil<T> Stencil<T>::Derivative(int a_n, int a_dir, int a_order)
{
  static const int MAXORDH = 6;
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::Derivative(int a_n, int a_dir, int a_order) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    PROTO_ASSERT((a_n > 0) && (a_n <= 2*MAXORDH),
    "Stencil::Derivative(int a_n, int a_dir, int a_order) invalid for a_n = %i.\
    a_n must be in [1,%i].",a_n,2*MAXORDH);
    PROTO_ASSERT((a_order > 0) && (a_order <= 2*MAXORDH),
    "Stencil::Derivative(int a_n, int a_dir, int a_order) invalid for a_order = %i.\
    a_order must be in [1,%i].",a_order,2*MAXORDH);
    PROTO_ASSERT(2*MAXORDH +2 - a_n - a_order >= 0,
    "Stencil::Derivative(int a_n, int a_dir, int a_order) invalid for a_n = %i and a_order = %i.\
    Implementation of requested accuracy not available.\n The maximum accuracy available for a_n = %i is\
    %i.",a_n, a_order, a_n, 2*MAXORDH +2 - (a_n/2)*2);

    if (a_order % 2 != 0){a_order ++;}

    Point edir = Point::Basis(a_dir);
    Stencil<T> factor = ((T)(1.0))*Shift(edir) + ((T)(-2.0))*Shift(Point::Zeros()) + ((T)(1.0))*Shift(edir*(-1.0));
    Stencil<T> deriv;
    Stencil<T> temps[MAXORDH];
    if (a_n % 2 == 0) //even derivative
    {
        T evenMatrixCoefs[MAXORDH][MAXORDH] =
          {{1.,-1./12,1./90,-1./560,1./3150,-1./16632},
           {0,1.,-1./6,7./240,-41./7560,479./453600},
           {0,0,1.,-1./4,13./240,-139./12096},
           {0,0,0,1.,-1./3,31./360},{0,0,0,0,1.,-5./12},
           {0,0,0,0,0,1.}};
        
        int p = (a_n - 2)/2;
        temps[0] = factor;
        for (int ord = 1; ord < MAXORDH; ord++)
        {
            temps[ord] = temps[ord-1]*factor;
        }
        for (int ii = 0; ii < (a_n + a_order - 2)/2; ii++)
        {
            deriv += evenMatrixCoefs[p][ii]*temps[ii];
        }
    } else {//odd derivative
        T oddMatrixCoefs[MAXORDH][MAXORDH] =
          {{1.,-1./6,1./30,-1./140,1./630,-1./2772},
           {0,1.,-1./4,7./120,-41./3024,479./151200},
           {0,0,1.,-1./3,13./144,-139./6048},
           {0,0,0,1.,-5./12,31./240},
           {0,0,0,0,1.,-1./2},
           {0,0,0,0,0,1.}};
        
        int p = (a_n - 1)/2;
        temps[0] = ((T)(0.5))*Shift(edir) + ((T)(-0.5))*Shift(edir*(-1));
        for (int ord = 1; ord < MAXORDH; ord++)
        {
            temps[ord] = temps[ord-1]*factor;
        }
        for (int ii = 0; ii < (a_n + 1 + a_order - 2)/2; ii++)
        {
            deriv += oddMatrixCoefs[p][ii]*temps[ii];
        }
    }
    return deriv;
    
}


#endif
