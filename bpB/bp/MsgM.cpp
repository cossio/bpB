#include <catch.hpp>
#include "MsgM.h"

namespace Coss { namespace bpB {

TEST_CASE("msgM 1", "[bp]")
{
	MsgM m(-1, 2);
	CHECK(m.alpha() == 2.0);
	CHECK(m.beta() == 2.0);
	CHECK(m.A() == -1.0);
	CHECK(m.B() == 2.0);
	CHECK(m.mean() == 0.5);
	CHECK(m.var() == 0.45);
	CHECK(m.mean_trunc() == Approx(m.mean()).epsilon(1e-10));
	CHECK(m.var_trunc() == Approx(m.var()).epsilon(1e-10));
}

TEST_CASE("msgM 2", "[bp]")
{
	MsgM m(-1, 1, -1, 2, 4, 5);

	CHECK(m.left() == -1);
	CHECK(m.right() == 1);

	CHECK(m.left_01_scaled() == 0);
	CHECK(m.right_01_scaled() == Approx(0.6666666666666666).epsilon(1e-10));

	CHECK(m.mean() == Approx(0.3333333333333333).epsilon(1e-10));
	CHECK(m.var() == Approx(0.22222222222222357).epsilon(1e-10));
	CHECK(m.mean_trunc() == Approx(0.2501485442661913).epsilon(1e-10));
	CHECK(m.var_trunc() == Approx(0.1626188133475536).epsilon(1e-10));

	CHECK(m.log_norm() == Approx(3.154108706175628).epsilon(1e-10));
	CHECK(m.log_norm_trunc() == Approx(3.0620549164850637).epsilon(1e-10));
}

TEST_CASE("msgM_regularize", "[bp]")
{
	MsgM m(-1, 1, -1, 2, 4, 5);
	double mean_tr = m.mean_trunc(), var_tr = m.var_trunc();
	m.regularize();
	CHECK(m.A() == m.left());
	CHECK(m.B() == m.right());
	CHECK(m.mean() == Approx(mean_tr).epsilon(1e-7));
	CHECK(m.var() == Approx(var_tr).epsilon(1e-7));
}

}}	// namespace Coss::bpB