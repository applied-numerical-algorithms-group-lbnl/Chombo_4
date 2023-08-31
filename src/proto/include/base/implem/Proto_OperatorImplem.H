using namespace Proto;

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_convolve(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    a_avg |= Stencil<T>::Laplacian()(a_2nd, 1.0/24.0);
    a_avg += a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::convolve(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr)
{
    _convolve(a_avg, a_ctr, a_ctr);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_convolve(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    BoxData<T,C,MEM,D,E> a_avg(a_ctr.box() & a_2nd.box().grow(-1));
    _convolve(a_avg, a_ctr, a_2nd);
    return a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::convolve(const BoxData<T,C,MEM,D,E>& a_ctr)
{
    return _convolve(a_ctr, a_ctr);
}
/*
template<class T, unsigned int C, MemType MEM, Centering CTR>
void
Operator::_convolve(
        LevelBoxData<T,C,MEM,CTR>& a_avg,
        const LevelBoxData<T,C,MEM,CTR>& a_ctr,
        const LevelBoxData<T,C,MEM,CTR>& a_2nd)
{
    //TODO: switch functionality based on value of CTR
    for (auto iter = a_avg.begin(); iter.ok(); ++iter)
    {
        auto& avg =  a_avg[*iter];
        auto& ctr =  a_ctr[*iter];
        auto& ctr2 = a_2nd[*iter];
        _convolve(avg, ctr, ctr2);
    }
}

template<class T, unsigned int C, MemType MEM, Centering CTR>
void
Operator::convolve(
        LevelBoxData<T,C,MEM,CTR>& a_avg,
        const LevelBoxData<T,C,MEM,CTR>& a_ctr)
{
    _convolve(a_avg, a_ctr, a_ctr);
}
*/
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_deconvolve(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    a_ctr |= Stencil<T>::Laplacian()(a_2nd, -1.0/24);
    a_ctr += a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::deconvolve(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg)
{
    _deconvolve(a_ctr, a_avg, a_avg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_deconvolve(
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    BoxData<T,C,MEM,D,E> a_ctr(a_avg.box() & a_2nd.box().grow(-1));
    _deconvolve(a_ctr, a_avg, a_2nd);
    return a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::deconvolve(const BoxData<T,C,MEM,D,E>& a_avg)
{
    return _deconvolve(a_avg, a_avg);
}
/*
template<class T, unsigned int C, MemType MEM, Centering CTR>
void
Operator::_deconvolve(
        LevelBoxData<T,C,MEM,CTR>& a_ctr,
        const LevelBoxData<T,C,MEM,CTR>& a_avg,
        const LevelBoxData<T,C,MEM,CTR>& a_2nd)
{
    //TODO: switch functionality based on value of CTR
    for (auto iter = a_avg.begin(); iter.ok(); ++iter)
    {
        auto& ctr =  a_ctr[*iter];
        auto& avg =  a_avg[*iter];
        auto& ctr2 = a_2nd[*iter];
        deconvolve(ctr, avg, ctr2);
    }
}

template<class T, unsigned int C, MemType MEM, Centering CTR>
void
Operator::deconvolve(
        LevelBoxData<T,C,MEM,CTR>& a_ctr,
        const LevelBoxData<T,C,MEM,CTR>& a_avg)
{
    _deconvolve(a_ctr, a_avg, a_avg);
}
*/
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_convolveFace(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1).grow(a_dir,1)),
    "Error in Operator::convolveFace | Insufficient source data.");
    a_avg |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, 1.0/24.0);
    a_avg += a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::convolveFace(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        int a_dir)
{
    _convolveFace(a_avg, a_ctr, a_ctr, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_convolveFace(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    auto b = a_2nd.box().grow(-1).grow(a_dir, 1);
    b &= a_ctr.box();
    BoxData<T,C,MEM,D,E> avg(b);
    _convolveFace(avg, a_ctr, a_2nd, a_dir);
    return avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::convolveFace(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        int a_dir)
{
    return _convolveFace(a_ctr, a_ctr, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_deconvolveFace(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1).grow(a_dir,1)),
    "Error in Operator::convolveFace | Insufficient source data.");
    a_ctr |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, -1.0/24.0);
    a_ctr += a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::deconvolveFace(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        int a_dir)
{
    _deconvolveFace(a_ctr, a_avg, a_avg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_deconvolveFace(
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    auto b = a_2nd.box().grow(-1).grow(a_dir, 1);
    b &= a_avg.box();
    BoxData<T,C,MEM,D,E> ctr(b);
    _deconvolveFace(ctr, a_avg, a_2nd, a_dir);
    return ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::deconvolveFace(
        const BoxData<T,C,MEM,D,E>& a_avg,
        int a_dir)
{
    return _deconvolveFace(a_avg, a_avg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_cellProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    auto b = (a_p2nd.box().grow(-1) & a_q2nd.box().grow(-1));
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,MEM,D,E> pq(b);
    _cellProduct(pq, a_pavg, a_qavg, a_p2nd, a_q2nd);
    return pq;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::cellProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    return _cellProduct(a_pavg, a_qavg, a_pavg, a_qavg);
}
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_cellProduct(BoxData<T,C,MEM,D,E>& a_product,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellProduct | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellProduct | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellProduct | incompatible patch sizes.");
    for (int ii = 0; ii < DIM; ii++)
    {
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);

        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, int v_dir)
            {
                for (int cc = 0; cc < C; cc++)
                for (int dd = 0; dd < D; dd++)
                for (int ee = 0; ee < E; ee++)
                {
                    if (v_dir == 0)
                    {
                        v_out(cc,dd,ee) = v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                    } else {
                        v_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                    }
                }
            }, a_product, Dp, Dq, ii);
    }

    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += v_p(cc,dd,ee)*v_q(cc,dd,ee);
            }
        }, a_product, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::cellProduct(BoxData<T,C,MEM,D,E>& a_product,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    _cellProduct(a_product, a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_faceProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    Proto::Box b = (a_p2nd.box() & a_q2nd.box()).grow(-1).grow(a_dir, 1);
    
    BoxData<T,C,MEM,D,E> out(b);
    _faceProduct(out, a_pavg, a_qavg, a_p2nd, a_q2nd, a_dir);
    return out;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::faceProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    return _faceProduct(a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_faceProduct(BoxData<T,C,MEM,D,E>& a_dest,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    bool do_replace = true;
    for (int ii = 0; ii < DIM; ii++)
    {
        if (ii == a_dir){continue;}
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, bool v_replace)
            {
                for (int cc = 0; cc < C; cc++)
                {
                for (int dd = 0; dd < D; dd++)
                {
                for (int ee = 0; ee < E; ee++)
                {
                    if (v_replace)
                    {
                        v_out(cc,dd,ee) = v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                    } else {
                        v_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                    }
                }
                }
                }
            }, a_dest, Dp, Dq, do_replace);
        do_replace = false;
    }
    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += v_p(cc,dd,ee)*v_q(cc,dd,ee);
            }
        }, a_dest, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::faceProduct(BoxData<T,C,MEM,D,E>& a_dest,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    _faceProduct(a_dest, a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}


template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_cellQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    auto b = (a_p2nd.box().grow(-1) & a_q2nd.box().grow(-1));
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,MEM,D,E> quotient(b);
    _cellQuotient(quotient, a_pavg, a_qavg, a_p2nd, a_q2nd);
    return quotient;

}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::cellQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    return _cellQuotient(a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_cellQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellQuotient | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    //BoxData<T,C,MEM,D,E> out(a_p2nd.box().grow(-1));
    //out.setVal(0);
    for (int ii = 0; ii < DIM; ii++)
    {
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q, int v_dir)
            {
                for (int cc = 0; cc < C; cc++)
                for (int dd = 0; dd < D; dd++)
                for (int ee = 0; ee < E; ee++)
                {
                    if (v_dir == 0)
                    {
                        v_out(cc,dd,ee) = v_dp(cc,dd,ee)*v_dq(cc,dd,ee);
                    } else {
                        v_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee);
                    }
                    v_out(cc,dd,ee) -= v_p(cc,dd,ee)/v_q(cc,dd,ee)*v_dq(cc,dd,ee)*v_dq(cc,dd,ee);
                    v_out(cc,dd,ee) *= -1.0/(12.0*v_q(cc,dd,ee)*v_q(cc,dd,ee));
                }
            }, a_quotient, Dp, Dq, a_p2nd, a_q2nd, ii);
    }
    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += (v_p(cc,dd,ee) / v_q(cc,dd,ee));
            }
        }, a_quotient, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::cellQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    _cellQuotient(a_quotient, a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_faceQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    Proto::Box b = (a_p2nd.box() & a_q2nd.box()).grow(-1).grow(a_dir, 1);
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,MEM,D,E> out(b);
    _faceQuotient(out, a_pavg, a_qavg, a_p2nd, a_q2nd, a_dir);
    return out;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::faceQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    return _faceQuotient(a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_faceQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellQuotient | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    
    for (int ii = 0; ii < DIM; ii++)
    {
        if (ii == a_dir){continue;}
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q, int v_dir)
            {
                for (int cc = 0; cc < C; cc++)
                for (int dd = 0; dd < D; dd++)
                for (int ee = 0; ee < E; ee++)
                {
                    if (v_dir == 0)
                    {
                        v_out(cc,dd,ee) = v_dp(cc,dd,ee)*v_dq(cc,dd,ee);
                    } else {
                        v_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee);
                    }
                    v_out(cc,dd,ee) -= v_p(cc,dd,ee)/v_q(cc,dd,ee)*v_dq(cc,dd,ee)*v_dq(cc,dd,ee);
                    v_out(cc,dd,ee) *= -1.0/(12.0*v_q(cc,dd,ee)*v_q(cc,dd,ee));
                }
            }, a_quotient, Dp, Dq, a_p2nd, a_q2nd, ii);
    }
    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += (v_p(cc,dd,ee) / v_q(cc,dd,ee));
            }
        }, a_quotient, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::faceQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    _faceQuotient(a_quotient, a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

