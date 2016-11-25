#include <valarray>
#include <catch.hpp>
#include "opt.h"

using Coss::VecD;
using Coss::Vec_func;

TEST_CASE("lambda_minmax")
{
    std::vector<double> val{1.0, 4.0, 5.4, 45.0};
    std::valarray<double> val2{1.0, 4.0, 5.4, 45.0};
    VecD A{-1.0, -1.0, -1.0};
    VecD B = -A;
    VecD P{0.5, 0.5, 0.5};
    VecD n{1.0, 1.0, 1.0};

    double lambda_min, lambda_max;
    Coss::Opt::lambda_minmax(A, B, P, n, lambda_min, lambda_max);
    REQUIRE(lambda_max == 0.5);
    REQUIRE(lambda_min == -1.5);
}

TEST_CASE("xj_minmax")
{
    Vec_func F = [](VecD const& x) { return -x.norm(); };
    auto A = VecD(3, -1.0);
    auto B = VecD(3, 1.0);
    auto X = VecD(3, 0.0);
    double c = -0.5;
    size_t k = 0, j = 0;
    double xj_min, xj_max;
//    xj_minmax(F, A, B, k, j, c, X, xj_min, xj_max);
//
//    REQUIRE(xj_min == -0.5);
//    REQUIRE(xj_max == 0.5);
}