template<class T>
InterpStencil<T> InterpStencil<T>::PiecewiseConstant(Point a_ratio)
{
    InterpStencil<T> ret(a_ratio);
    Box K = Box(a_ratio);
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        Stencil<T>& S = ret(*iter);
        S = 1.0*Shift::Zeros();
        S.destRatio() = a_ratio;
        S.destShift() = (*iter);
    }
    return ret;
}

template<class T>
InterpStencil<T> InterpStencil<T>::PiecewiseLinear(Point a_ratio)
{
    InterpStencil<T> ret(a_ratio);
    Box K = Box(a_ratio);
    Box U = Box::Cube(2);
    for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
    {
        Point p = *kiter;
        Point q = a_ratio - p;
        Stencil<double>& S = ret(*kiter);
        for (auto uiter = U.begin(); uiter != U.end(); ++ uiter)
        {
            Point e = *uiter;
            double coef = 1.0;
            for (int ii = 0; ii < DIM; ii++)
            {
                if (e[ii] == 0)
                {
                    coef *= q[ii]/((double)a_ratio[ii]);
                } else {
                    coef *= p[ii]/((double)a_ratio[ii]);
                }
            }
            if (coef != 0)
            {
                S += coef*Shift(e);
            }
        }
        S.destRatio() = a_ratio;
        S.destShift() = p;
    }
    return ret;    
}

// helper function for Quadratic InterpStencil
inline double qcoef(int a_ni, int a_r, int a_side)
{
    PROTO_ASSERT(std::abs(a_side) <= 1,
    "qcoef(ni, r, side) invalide for side = %i. Value must be in [-1,1]", a_side);
    double n = a_ni; double r = a_r;
    double coef = 0;
    if (a_side == -1)
    {
        coef = 1.0 + 3.0*n + 3.0*n*n - 3.0*r - 6.0*n*r + 2.0*r*r;
    } 
    else if (a_side == 0)
    {
        coef = -2.0 - 6.0*n - 6.0*n*n + 3.0*r + 6.0*n*r + 5.0*r*r;
    }
    else if (a_side == 1)
    {
        coef = -(-1.0 - 3.0*n - 3.0*n*n + r*r);
    }
    coef /= (6.0*r*r);
    return coef;
}

template<class T>
InterpStencil<T> InterpStencil<T>::Quadratic(int a_ratio)
{
    //can make it anisotropic later...
    InterpStencil<T> ret(Point::Ones(a_ratio));
    Box K = Box(Point::Ones(a_ratio));
    Box U = Box::Cube(3).shift(Point::Ones(-1));

    double r = a_ratio;
     
    for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
    {
        Point n = *kiter;

        Stencil<double>& S = ret(*kiter);
        for (auto iter = U.begin(); iter != U.end(); ++iter)
        {
            T coef = 1;
            for (int ii = 0; ii < DIM; ii++)
            {
                coef *= (T)qcoef(n[ii],r,(*iter)[ii]);
            }
            S += coef*Shift(*iter);
        }
        S.destRatio() = Point::Ones(a_ratio);
        S.destShift() = n;
    }
    return ret;    
}

// helper function for Build(...)
/*
    Returns <x^p> over the coarse cell shifted by s
*/
inline double coarseCoef(Point p, Point s, double h)
{
    double coef = 1.0;
    for (int d = 0; d < DIM; d++)
    {
        double high = h*(2.0*s[d] + 1)/2.0;
        double low = h*(2.0*s[d] - 1)/2.0;
        coef *= 1.0/(p[d]+1)
            *(std::pow(high, p[d]+1) - std::pow(low, p[d]+1));
    }
    coef /= std::pow(h,DIM);
    return coef;
}

// helper function for Build(...)
/**
    Returns <x^p> over the fine cell s with coarse spacing h
*/
inline double fineCoef(Point p, Point s, double h, Point refRatio)
{
    std::array<double, DIM> hf;
    double denominator = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        hf[dir] = h / refRatio[dir];
        denominator *= hf[dir];
    }
    //double hf = h/refRatio;
    auto low = s - (refRatio / 2);
    //auto low = Proto::Point::Ones(-refRatio/2) + s;
    auto high = low + Proto::Point::Ones();
    
    double coef = 1.0;
    for (int d = 0; d < DIM; d++)
    {
        coef *= 1.0/(p[d]+1)
            *(std::pow(high[d]*hf[d], p[d]+1) - std::pow(low[d]*hf[d], p[d]+1));
    }
    //coef /= std::pow(hf,DIM);
    coef /= denominator;
    return coef;
}

inline double fineCoef(Point p, Point s, double h, int refRatio)
{
    return fineCoef(p, s, h, Point::Ones(refRatio));
}


// LAPACK functions for Build

template<class T>
InterpStencil<T> InterpStencil<T>::Build(int a_smax, Box a_skernel,
        int a_expMax, int a_refRatio)
{
    return Build(a_smax, a_skernel, a_expMax, Point::Ones(a_refRatio));
}

template<class T>
InterpStencil<T> InterpStencil<T>::Build(int a_smax, Box a_skernel,
        int a_expMax, Point a_refRatio)
{
    std::vector<Point> shifts;
    for (auto iter = a_skernel.begin(); iter != a_skernel.end(); ++iter)
    {
        int dist = (*iter).sumAbs();
        if (dist <= a_smax && dist > 0)
        {
            shifts.push_back(*iter);
        }
    }
    std::vector<Point> orders;
    Box pkernel = Box::Cube(a_expMax+1);
    for (auto iter = pkernel.begin(); iter != pkernel.end(); ++iter)
    {
        int sum = (*iter).sum();
        if (sum <= a_expMax){orders.push_back(*iter);}
    }
    return InterpStencil<T>::Build(shifts,orders,a_refRatio);
}

template<class T>
InterpStencil<T> InterpStencil<T>::Build(
        std::vector<Point>& a_shifts, 
        int a_maxOrder,
        Point a_refRatio)
{
    std::vector<Point> orders;
    Box K = Box(Point::Ones(a_maxOrder+1));
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        if ((*iter).sum() <= a_maxOrder)
        {
            orders.push_back(*iter);
        }
    }

    return InterpStencil<T>::Build(a_shifts, orders, a_refRatio);
}

template<class T>
InterpStencil<T> InterpStencil<T>::Build(
        std::vector<Point>& a_shifts, 
        const std::vector<Point>& a_orders,
        Point a_refRatio)
{
    int ref0 = a_refRatio[0];
    for (int dir = 1; dir < DIM; dir++)
    {
        if (a_refRatio[dir] != ref0)
        {
            MayDay<void>::Abort(
                "InterpStencil::Build | Error: Anisotropic refinement ratios are not currently compatible with Build.");
        }
    }

    // ==============================================================
    // COMPUTE COEFFICIENTS
    // --------------------
    
    // Prevent the zero shift or duplicate shifts from being used
  using std::set;
    std::set<Point> tmp;
    for (auto iter = a_shifts.begin(); iter != a_shifts.end(); ++iter)
    {
        if ((*iter) != Proto::Point::Zeros())
        {
            tmp.insert(*iter);
        }
    }
    a_shifts.clear();
    for (auto iter = tmp.begin(); iter != tmp.end(); ++iter)
    {
        a_shifts.push_back(*iter);
    }
    int M = a_shifts.size();
    int N = a_orders.size();
    /*
    std::cout << "Shifts: " << std::endl;
    for (auto iter = a_shifts.begin(); iter != a_shifts.end(); ++iter)
    {
        std::cout << "\t" << *iter << std::endl;
    }
    std::cout << "Orders: " << std::endl;
    for (auto iter = a_orders.begin(); iter != a_orders.end(); ++iter)
    {
        std::cout << "\t" << *iter << std::endl;
    }
    */
    double* X = (double*)malloc(sizeof(double)*M*N);
    for (int ii = 0; ii < M; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            X[jj*M+ii] = coarseCoef(a_orders[jj], a_shifts[ii], 1.0);
        }
    }
    
    double* X0 = (double*)malloc(sizeof(double)*N);
    for (int jj = 0; jj < N; jj++)
    {
        X0[jj] = coarseCoef(a_orders[jj],Point::Zeros(),1.0);
    }
    
    double* XTX2 = (double*)malloc(sizeof(double)*N*N);
    char _trans_a, _trans_b;
    int _m, _n, _k, _lda, _ldb, _ldc;
    double _alpha, _beta;

    _trans_a = 'T'; _trans_b = 'N';
    _m = N; _n = N; _k = M;
    _lda = M; _ldb = M; _ldc = N;
    _alpha = 2.0; _beta = 0.0;
    dgemm_(&_trans_a, &_trans_b, &_m, &_n, &_k, &_alpha,
            X, &_lda, X, &_ldb, &_beta, XTX2, &_ldc);
    
    double* Y = (double*)malloc(sizeof(double)*(N+1)*(N+1));
    for (int ii = 0; ii <= N; ii++)
    {
        for (int jj = 0; jj <= N; jj++)
        {
            if (ii < N && jj < N) {
                Y[jj*(N+1)+ii] = XTX2[jj*N+ii];
            } else if (ii >= N && jj < N){
                Y[jj*(N+1)+ii] = X0[jj];
            } else if (jj >= N && ii < N){
                Y[jj*(N+1)+ii] = X0[ii];
            } else {
                Y[jj*(N+1)+ii] = 0;
            }
        }
    }

    double* Z = (double*)malloc(sizeof(double)*(N+1)*(M+1));
    for (int ii = 0; ii <= N; ii++)
    {
        for (int jj = 0; jj <= M; jj++)
        {
            if (ii < N && jj < M) {
                Z[jj*(N+1)+ii] = 2.0*X[ii*M+jj];//X transpose
            } else if (ii >= N && jj < M){
                Z[jj*(N+1)+ii] = 0;
            } else if (jj >= M && ii < N){
                Z[jj*(N+1)+ii] = 0;
            } else {
                Z[jj*(N+1)+ii] = 1;
            }
        }
    }

    int* pivots = (int*)malloc(sizeof(int)*(N+1));
    double* _work = (double*)malloc(sizeof(double)*(N+1));
    int _info;
    _n = N+1,
    dgetrf_(&_n, &_n, Y, &_n, pivots, &_info);
    PROTO_ASSERT(_info == 0, "Y is singular.");
    dgetri_(&_n, Y, &_n, pivots, _work, &_n, &_info);
    PROTO_ASSERT(_info == 0, "Failed to invert Y.");
    
    double* YinvZ = (double*)malloc(sizeof(double)*(N+1)*(M+1));
    _trans_a = 'N'; _trans_b = 'N';
    _n = M+1; _k = N+1; _m = N+1;
    _alpha = 1.0;
    _lda = _m; _ldb = _m; _ldc = _m;
    dgemm_(&_trans_a, &_trans_b, &_m, &_n, &_k,
            &_alpha, Y, &_lda, Z, &_ldb, &_beta, YinvZ, &_ldc);

    double* Q = (double*)malloc(sizeof(double)*N*(M+1));
    for (int ii = 0; ii < N; ii++)
    {
        for (int jj = 0; jj < M+1; jj++)
        {
            Q[N*jj+ii] = YinvZ[(N+1)*jj+ii];
        }
    }
   
    //int Mf = ipow<DIM>(a_refRatio);
    
    int Mf = 1;
    for (int dir = 0; dir < DIM; dir++)
    {
        Mf *= a_refRatio[dir];
    }
    Box fkernel(a_refRatio);
    
    //auto fkernel = Box::Cube(a_refRatio);
    double* Xf = (double*)malloc(sizeof(double)*Mf*N);
    int ii = 0;
    for (auto fiter = fkernel.begin(); fiter != fkernel.end();
             ++fiter, ++ii)
    {
        for (int jj = 0; jj < N; jj++)
        {
            Xf[jj*Mf+ii] = fineCoef(a_orders[jj], *fiter, 1.0, a_refRatio);
        }
    }
   
    double* R = (double*)malloc(sizeof(double)*Mf*(M+1));
    _trans_a = 'N'; _trans_b = 'N';
    _n = M+1; 
    _k = N; _m = Mf;
    _alpha = 1.0; _beta = 0.0;
    _lda = _m; _ldb = _k; _ldc = _m;
    dgemm_(&_trans_a, &_trans_b, &_m, &_n, &_k, &_alpha,
            Xf, &_lda, Q, &_ldb, &_beta, R, &_ldc);

    free(X);
    free(X0);
    free(XTX2);
    free(Y);
    free(Z);
    free(YinvZ);
    free(_work);
    free(pivots);
    free(Q);
    free(Xf);
    
    // ==============================================================
    // BUILD STENCILS
    // --------------

    //InterpStencil<T> ret(Point::Ones(a_refRatio));
    InterpStencil<T> ret(a_refRatio);
    ii = 0;
    for (auto fiter = fkernel.begin(); fiter != fkernel.end();
             ++fiter, ++ii)
    {
        Stencil<T>& S = ret(*fiter);
        for (int jj = 0; jj < a_shifts.size(); jj++)
        {
            S += ((T)R[Mf*jj+ii])*Shift(a_shifts[jj]);
        }
        S += ((T)R[Mf*M+ii])*Shift::Zeros();
        //S.destRatio() = Point::Ones(a_refRatio);
        S.destRatio() = a_refRatio;
        S.destShift() = (*fiter);
    }
    
    free(R);
    
    return ret;
}

template<class T>
InterpStencil<T> InterpStencil<T>::Build(int a_order, int a_refRatio)
{
    return Build(a_order, Point::Ones(a_refRatio));
}

template<class T>
InterpStencil<T> InterpStencil<T>::Build(int a_order, Point a_refRatio)
{
    InterpStencil<T> ret;
    if (a_order <= 1)
    {
        return InterpStencil<T>::PiecewiseConstant(a_refRatio);
    } else if (a_order <= 2)
    {
        return InterpStencil<T>::Build(1, Box::Kernel(1), 1, a_refRatio);
    } else if (a_order <= 3)
    {
        return InterpStencil<T>::Build(DIM, Box::Kernel(1), 2, a_refRatio);
    } else if (a_order <= 4)
    {
        return InterpStencil<T>::Build(DIM, Box::Kernel(2), 3, a_refRatio);
    } else if (a_order <= 5)
    {
        return InterpStencil<T>::Build(4, Box::Kernel(2), 4, a_refRatio);
    } else {
        std::cout << "InterpStencil::Build | Error: No output defined for order = " << a_order << std::endl; 
    }
    return ret; //should never get here, this is just to stifle a compiler warning
}
