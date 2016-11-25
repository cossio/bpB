#include <catch.hpp>
#include "MsgN.h"

TEST_CASE("msgN", "[msgN]")
{
	Coss::Peak p(0, -1, 1, 1, 0);

	Coss::bpB::MsgN msg(0.0, 1.0);
	CHECK(msg.alpha_trunc() == Approx(1.0).epsilon(1e-10));
	CHECK(msg.beta_trunc() == Approx(1.0).epsilon(1e-10));

	msg = Coss::bpB::MsgN(-20, 40, 5.71428571428571428571428571428571428571, 12.41736131072147168726645587812589824662, 0.0, p);
	CHECK(msg.alpha_trunc() == Approx(30).epsilon(1e-10));
	CHECK(msg.beta_trunc() == Approx(40).epsilon(1e-10));
}