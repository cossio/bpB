#include <catch.hpp>
#include <Coss/gsl/function_wrapper.h>
#include <Coss/gsl/integration_workspace_wrapper.h>
#include <gsl/gsl_fft_real.h>

TEST_CASE("playground", "[play]")
{
	double z[] = {1, 2, 3, 4, 5, 6, 7, 8};

	gsl_fft_real_radix2_transform(z, 1, 8);

	CHECK(true);

}


TEST_CASE("gsl_integrate", "[play]")
{
	/*
	 * Just checking that the eps_err argument to gsl_integration
	 * functions doesn't make that much difference.
	 */

	auto integrand = [&](double x) { return exp(x) * sin(x); };
	Coss::gsl::gsl_function_pp<decltype(integrand)> Fp(integrand);
	gsl_function* F = static_cast<gsl_function*>(&Fp);

	double pts[2] = {1, 10};
	Coss::gsl::gsl_integration_workspace_wrapper W;

	double integral_1, abs_error_1;
	int gsl_err_code_1 = gsl_integration_qagp(F, pts, 2, 0, 1e-04, W.limit(), W.w(), &integral_1, &abs_error_1);

	double integral_2, abs_error_2;
	int gsl_err_code_2 = gsl_integration_qagp(F, pts, 2, 0, 1e-10, W.limit(), W.w(), &integral_2, &abs_error_2);

	double integral_3, abs_error_3;
	int gsl_err_code_3 = gsl_integration_qagp(F, pts, 2, 0, 1e-12, W.limit(), W.w(), &integral_3, &abs_error_3);

	CHECK(true);
}