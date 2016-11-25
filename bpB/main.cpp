#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include <gsl/gsl_errno.h>
#include <Coss/py/CPyInstance.h>

int main(int argc, char* const argv[])
{
	// global setup...

	// so that GSL doesn't kill the program
	gsl_set_error_handler_off();

	// run Catch tests
	int result = Catch::Session().run(argc, argv);

	// global clean-up...

	return result;
}