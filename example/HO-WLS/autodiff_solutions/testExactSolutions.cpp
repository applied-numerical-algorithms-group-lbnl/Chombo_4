// main() provided by linkage to Catch2WithMain
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include "ExactSolutions.hpp"
#include "EBGeomInfo.hpp"
#include "EBHashDataWriteVTK.hpp"

#include "../geometry/exampleEBGeomInfo.hpp"

#include "catch2ParallelMain.hpp"

TEST_CASE("Exact solutions for points and regular cells are correct", "[ebTools]")
{
  using namespace hoeb;
  using namespace exactSolution;
  // LogLevel verbosity = g_defaultLogLevel;
  // LogLevel verbosity = LogLevel::Detail;

  constexpr std::array<Real, 10> testEvalPts{ 0., 0.5,  1.,  -1.,           1.3,
                                              2., -5.1, 1e4, hbr::Real_tol, hoeb::Exp() };

  SECTION("One dimension, quadratic monomial")
  {
    // x^2 function
    exactSolution::Monomial<1> fnV(stc::IVec<1>{ 2 });
    // direct point evaluations are possible, no autodiff invocation
    CHECK(fnV(stc::RVec<1>{ 0. }) == 0.);
    CHECK(fnV(stc::RVec<1>{ 1. }) == 1.);
    CHECK(fnV(stc::RVec<1>{ 2. }) == 4.);

    // unit size moments in 1D
    auto mom = hoeb::make_moments<1, 2>::regular();
    stc::RVec<1> dx = stc::make_RVec<1>::unit();
    // test autodiff invocations the assorted test points
    for (Real xi : testEvalPts)
    {
      // points evaluations
      CHECK_THAT(exactSolution::evalPoint<1>(fnV, stc::RVec<1>{ xi }),
                 Catch::Matchers::WithinRel(std::pow(xi, 2), hbr::Real_eps));
      // point derivatives
      CHECK_THAT(exactSolution::evalPoint<1>(fnV, stc::RVec<1>{ xi }, stc::IVec<1>{ 1 }),
                 Catch::Matchers::WithinRel(2. * std::pow(xi, 1), hbr::Real_tol));
      CHECK_THAT(exactSolution::evalPoint<2>(fnV, stc::RVec<1>{ xi }, stc::IVec<1>{ 2 }),
                 Catch::Matchers::WithinRel(2., hbr::Real_tol));
      CHECK_THAT(exactSolution::evalPoint<4>(fnV, stc::RVec<1>{ xi }, stc::IVec<1>{ 3 }),
                 Catch::Matchers::WithinRel(0., hbr::Real_tol));

      // averages for unit length cells
      CHECK_THAT(exactSolution::evalAverage<2>(fnV, mom, stc::RVec<1>{ xi }),
                 Catch::Matchers::WithinRel(std::pow(xi + 0.5, 3) / 3. - std::pow(xi - 0.5, 3) / 3.,
                                            hbr::Real_tol));
      // averages of derivatives
      Real exactDAvg = std::pow(xi + 0.5, 2) - std::pow(xi - 0.5, 2);
      CHECK_THAT(exactSolution::evalAverage<3>(fnV, mom, stc::RVec<1>{ xi }, dx, stc::IVec<1>{ 1 }),
                 Catch::Matchers::WithinRel(exactDAvg, hbr::Real_tol)
                     || Catch::Matchers::WithinAbs(exactDAvg, hbr::Real_eps));

      Real exactD2Avg = 2 * std::pow(xi + 0.5, 1) - 2 * std::pow(xi - 0.5, 1);
      CHECK_THAT(exactSolution::evalAverage<4>(fnV, mom, stc::RVec<1>{ xi }, dx, stc::IVec<1>{ 2 }),
                 Catch::Matchers::WithinRel(exactD2Avg, hbr::Real_tol)
                     || Catch::Matchers::WithinAbs(exactD2Avg, hbr::Real_eps));
    }
  }

  SECTION("One dimension, quartic polynomial")
  {
    // 3*x^4 - 2x^3 + 1
    exactSolution::Polynomial<1, 3> fn(
        { { stc::IVec<1>{ 4 }, 3. }, { stc::IVec<1>{ 3 }, -2. }, { stc::IVec<1>{ 0 }, 1. } });
    // unit size moments in 1D
    auto mom = hoeb::make_moments<1, 4>::regular();
    // exact solutions
    auto exactFn = [](Real xi) { return 3. * std::pow(xi, 4) - 2. * std::pow(xi, 3) + 1; };
    auto exactFnD2
        = [](Real xi) { return 3. * 4. * 3. * std::pow(xi, 2) - 2. * 3. * 2. * std::pow(xi, 1); };

    auto exactFnI
        = [](Real xi) { return 3. * std::pow(xi, 5) / 5. - 2. * std::pow(xi, 4) / 4. + xi; };
    auto exactCellAvg = [&](Real xi) { return exactFnI(xi + 0.5) - exactFnI(xi - 0.5); };
    auto exactFnID2 = [](Real xi) { return 3. * 4. * std::pow(xi, 3) - 2. * 3. * std::pow(xi, 2); };
    auto exactCellAvgD2 = [&](Real xi) { return exactFnID2(xi + 0.5) - exactFnID2(xi - 0.5); };

    for (Real xi : testEvalPts)
    {
      // points
      CHECK_THAT(exactSolution::evalPoint<4>(fn, stc::RVec<1>{ xi }),
                 Catch::Matchers::WithinRel(exactFn(xi), hbr::Real_tol));
      // point second derivatives
      CHECK_THAT(exactSolution::evalPoint<4>(fn, stc::RVec<1>{ xi }, stc::IVec<1>{ 2 }),
                 Catch::Matchers::WithinRel(exactFnD2(xi), hbr::Real_tol));
      // averages for dx=1
      CHECK_THAT(exactSolution::evalAverage<4>(fn, mom, stc::RVec<1>{ xi }),
                 Catch::Matchers::WithinRel(exactCellAvg(xi), hbr::Real_tol));
      // averages of derivative
      // CHECK_THAT(exactSolution::evalAverage<6>(fn, mom, stc::RVec<1>{xi}, stc::IVec<1>{2}),
      //            Catch::Matchers::WithinRel(exactCellAvgD2(xi), hbr::Real_tol));
      // averages for non-unit dx
      stc::RVec<1> dx{0.5};
      Real exactAvg = (exactFnI(xi + .5*dx[0]) - exactFnI(xi - .5*dx[0])) / dx[0];
      CHECK_THAT(exactSolution::evalAverage<4>(fn, mom, stc::RVec<1>{ xi }, dx),
                 Catch::Matchers::WithinRel(exactAvg, hbr::Real_tol));
    }
  }

  SECTION("Three dimensions, cell averages of monmials")
  {
    auto mom = hoeb::make_moments<3, 4>::regular();
    auto momit = mom.getMomentIterator();
    for (momit.begin(); momit.ok(); ++momit)
    {
      auto p = momit.momentIndex();
      exactSolution::Monomial fn(p);
      stc::RVec<3> pt{ 0.5, 0.5, 0.5 };
      //
      CHECK_THAT(exactSolution::evalAverage<4>(fn, mom, pt),
                 Catch::Matchers::WithinRel(1. / (stc::product(p + 1)), hbr::Real_tol));
    }
  }

  SECTION("Three dimensions, quartic polynomial")
  {
    // 3*x^2*y^2 + y^4 - x*y*z^2 - 2z^3
    exactSolution::Polynomial<3, 4> fn({ { stc::IVec<3>{ 2, 2, 0 }, 3. },
                                         { stc::IVec<3>{ 0, 4, 0 }, 1. },
                                         { stc::IVec<3>{ 1, 1, 2 }, -1. },
                                         { stc::IVec<3>{ 0, 0, 3 }, -2. } });
    // unit size moments
    auto mom = hoeb::make_moments<3, 4>::regular();
    stc::RVec<3> dx = stc::make_RVec<3>::unit();
    // about zero
    stc::RVec<3> pt{ 0, 0, 0 };
    CHECK_THAT(exactSolution::evalAverage<4>(fn, mom, pt),
               Catch::Matchers::WithinRel(std::pow(0.5, 7) * 64. / 15., hbr::Real_tol));

    // generate random points to test
    stc::RVec<3> testPt{ GENERATE(take(4, random(-1e3, 1e3))), GENERATE(take(4, random(-1e3, 1e3))),
                         GENERATE(take(4, random(-1e3, 1e3))) };

    //  test point values
    auto exactFn = [](stc::RVec<3> x)
    {
      return 3 * std::pow(x[0], 2) * std::pow(x[1], 2) + std::pow(x[1], 4)
             - x[0] * x[1] * std::pow(x[2], 2) - 2 * std::pow(x[2], 3);
    };
    CHECK_THAT(exactSolution::evalPoint<1>(fn, testPt),
               Catch::Matchers::WithinRel(exactFn(testPt), hbr::Real_tol));

    auto exactFnDx = [](stc::RVec<3> x)
    { return 3 * 2 * std::pow(x[0], 1) * std::pow(x[1], 2) - x[1] * std::pow(x[2], 2); };
    CHECK_THAT(exactSolution::evalPoint<1>(fn, testPt, stc::IVec<3>{ 1, 0, 0 }),
               Catch::Matchers::WithinRel(exactFnDx(testPt), hbr::Real_tol));

    auto exactFnDxx = [](stc::RVec<3> x) { return 3 * 2 * std::pow(x[1], 2); };
    CHECK_THAT(exactSolution::evalPoint<2>(fn, testPt, stc::IVec<3>{ 2, 0, 0 }),
               Catch::Matchers::WithinRel(exactFnDxx(testPt), hbr::Real_tol));

    auto exactFnDy = [](stc::RVec<3> x)
    { return 3 * std::pow(x[0], 2) * 2 * x[1] + 4 * std::pow(x[1], 3) - x[0] * std::pow(x[2], 2); };
    CHECK_THAT(exactSolution::evalPoint<1>(fn, testPt, stc::IVec<3>{ 0, 1, 0 }),
               Catch::Matchers::WithinRel(exactFnDy(testPt), hbr::Real_tol));

    auto exactFnDyy
        = [](stc::RVec<3> x) { return 3 * std::pow(x[0], 2) * 2 + 4 * 3 * std::pow(x[1], 2); };
    CHECK_THAT(exactSolution::evalPoint<2>(fn, testPt, stc::IVec<3>{ 0, 2, 0 }),
               Catch::Matchers::WithinRel(exactFnDyy(testPt), hbr::Real_tol));

    auto exactFnDz = [](stc::RVec<3> x)
    { return -x[0] * x[1] * 2 * std::pow(x[2], 1) - 2 * 3 * std::pow(x[2], 2); };
    CHECK_THAT(exactSolution::evalPoint<1>(fn, testPt, stc::IVec<3>{ 0, 0, 1 }),
               Catch::Matchers::WithinRel(exactFnDz(testPt), hbr::Real_tol));

    auto exactFnDzz
        = [](stc::RVec<3> x) { return -x[0] * x[1] * 2 - 2 * 3 * 2 * std::pow(x[2], 1); };
    CHECK_THAT(exactSolution::evalPoint<2>(fn, testPt, stc::IVec<3>{ 0, 0, 2 }),
               Catch::Matchers::WithinRel(exactFnDzz(testPt), hbr::Real_tol));

    auto exactFnDxyz = [](stc::RVec<3> x) { return -2 * std::pow(x[2], 1); };
    CHECK_THAT(exactSolution::evalPoint<3>(fn, testPt, stc::IVec<3>{ 1, 1, 1 }),
               Catch::Matchers::WithinRel(exactFnDxyz(testPt), hbr::Real_tol));
    // test point operators
    auto grad = exactSolution::evalGradient(fn, testPt);
    stc::RVec<3> exactGrad = { exactFnDx(testPt), exactFnDy(testPt), exactFnDz(testPt) };
    for (int d = 0; d != 3; d++)
    {
      CHECK_THAT(grad[d], Catch::Matchers::WithinRel(exactGrad[d], hbr::Real_tol));
    }
    // 6 x^2 y + 6 x y^2 - 2 x y z - x z^2 + 4 y^3 - y z^2 - 6 z^2
    CHECK_THAT(exactSolution::evalDivergence(fn, testPt),
               Catch::Matchers::WithinRel(exactFnDx(testPt) + exactFnDy(testPt) + exactFnDz(testPt),
                                          hbr::Real_tol));
    // 6 x^2 - 2 x y + 18 y^2 - 12 z
    CHECK_THAT(exactSolution::evalLaplacian(fn, testPt),
               Catch::Matchers::WithinRel(
                   exactFnDxx(testPt) + exactFnDyy(testPt) + exactFnDzz(testPt), hbr::Real_tol));

    // test cell averages
    auto exactFnI = [](stc::RVec<3> x)
    {
      return 3 * (std::pow(x[0], 2) + 1. / 12.) * (std::pow(x[1], 2) + 1. / 12.) + std::pow(x[1], 4)
             + 0.5 * std::pow(x[1], 2) + 0.0125 - (std::pow(x[2], 2) + 1. / 12.) * x[0] * x[1]
             - 2 * std::pow(x[2], 3) - 0.5 * x[2];
    };
    CHECK_THAT(exactSolution::evalAverage<5>(fn, mom, testPt, dx),
               Catch::Matchers::WithinRel(exactFnI(testPt), hbr::Real_tol));

    // test cell average operators
    auto exactFnIDiv = [](stc::RVec<3> x)
    {
      return 0.5 * x[0] + 6. * std::pow(x[1], 2) * x[0] + x[1] * (0.5 + 6. * std::pow(x[0], 2))
             + (-1. / 12. - std::pow(x[2], 2) - 2. * x[2] * x[1]) * x[0] + x[1]
             + 4. * std::pow(x[1], 3) + (-1. / 12. - std::pow(x[2], 2)) * x[1] - 0.5
             - 6. * std::pow(x[2], 2);
    };

    CHECK_THAT(exactSolution::evalAvgDivergence(fn, mom, testPt, dx),
               Catch::Matchers::WithinRel(exactFnIDiv(testPt), hbr::Real_tol));

    auto exactFnILap = [](stc::RVec<3> x) {
      return 2. - 12. * x[2] + 18. * std::pow(x[1], 2) - 2. * x[1] * x[0] + 6. * std::pow(x[0], 2);
    };
    CHECK_THAT(exactSolution::evalAvgLaplacian(fn, mom, testPt, dx),
               Catch::Matchers::WithinRel(exactFnILap(testPt), hbr::Real_tol));
  }

  SECTION("Three dimensions, one-dimensional product of sin")
  {
    // sin(x)*sin(y)*sin(z)
    exactSolution::SinProduct fn(stc::RVec<3>{ 0.1, -0.2, 0. }, 1., stc::RVec<3>{ 1, 1.5, 2 });
    // unit size moments
    auto mom = hoeb::make_moments<3, 4>::regular();
    // about the center
    stc::RVec<3> pt{ 0.1, -0.2, 0. };
    CHECK_THAT(exactSolution::evalPoint<1>(fn, pt), Catch::Matchers::WithinRel(0., hbr::Real_tol));

    CHECK_THAT(exactSolution::evalAverage<4>(fn, mom, pt),
               Catch::Matchers::WithinRel(0., hbr::Real_tol));

    // test random points
    stc::RVec<3> testPt{ GENERATE(take(5, random(-hoeb::Pi(), hoeb::Pi()))), GENERATE(take(5, random(-hoeb::Pi(), hoeb::Pi()))),
                         GENERATE(take(5, random(-hoeb::Pi(), hoeb::Pi()))) };
    // test point values
    auto exactFn = [](stc::RVec<3> x)
    { return sin((x[0] - 0.1)) * sin(1.5 * (x[1] + 0.2)) * sin(2 * x[2]); };
    CHECK_THAT(exactSolution::evalPoint<4>(fn, testPt),
               Catch::Matchers::WithinRel(exactFn(testPt), hbr::Real_tol));
  }
}

TEST_CASE("Exact solution over a domain", "[ebTools]")
{
  // not really a test so much as a visual demo currently
  using namespace hoeb;
  // make the geometry on a unit domain
  EBGeomInfo<g_Order> geom;
  setExampleRegularGeometry<g_Order>(geom, 16);
  auto graph = geom.getGraph();

  // volume data
  EBHashData<Real> soln(graph->domainLayout(GridType::Vol));

  // generate a solution
  {
    exactSolution::Polynomial<g_SpaceDim, 1> fn;
    fn.setTerm(0, stc::IVec<g_SpaceDim>{ 1, 2, 0 }, 3);
    geom.setCellAverages(soln, fn, soln.layout());
    writeToFileVTK(soln, geom, "polynomialFunction", 0);
  }

  // generate a solution
  {
    exactSolution::SinProduct fn(stc::RVec<g_SpaceDim>{ 0, 0.25, 0 }, 0.5,
                                 stc::RVec<g_SpaceDim>{ hoeb::Pi(), 2 * hoeb::Pi(), 0 });
    geom.setCellAverages(soln, fn, soln.layout());
    writeToFileVTK(soln, geom, "sinFunction", 0);
  }

  // generate a solution
  {
    exactSolution::CosProduct fn(stc::RVec<g_SpaceDim>{ 0, 0.25, 0 }, 0.5,
                                 stc::RVec<g_SpaceDim>{ hoeb::Pi(), 2 * hoeb::Pi(), 0 });
    geom.setCellAverages(soln, fn, soln.layout());
    writeToFileVTK(soln, geom, "cosFunction", 0);
  }

  {
    exactSolution::Gaussian fn(stc::RVec<g_SpaceDim>{ 0.5, 0.4, 0.7 }, 1,
                               stc::RVec<g_SpaceDim>{ 0.4, 0.3, 1 });
    geom.setCellAverages(soln, fn, soln.layout());
    writeToFileVTK(soln, geom, "gaussianFunction", 0);
  }

  {
    exactSolution::PseudoGaussian fn(stc::RVec<g_SpaceDim>{ 0.5, 0.4, 0.7 }, 1,
                                     stc::RVec<g_SpaceDim>{ 0.4, 0.3, 1 });
    geom.setCellAverages(soln, fn, soln.layout());
    writeToFileVTK(soln, geom, "pseudogaussianFunction", 0);
  }
}
