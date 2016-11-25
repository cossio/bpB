#include <catch.hpp>
#include "Net.h"

namespace Coss { namespace bpB {

TEST_CASE("metabolic network")
{
	Net net;

	net.add_eqn(0);
	net.add_eqn(1);

	net.add_var(0.0, 1.0, 0);
	net.add_var(0.0, 1.0, 1);
	net.add_var(0.0, 1.0, 2);
	net.add_var(0.0, 1.0, 3);
	net.add_var(0.0, 1.0, 4);

	net.connect(net.eqns(0), net.vars(0), 1.0);
	net.connect(net.eqns(0), net.vars(1), -1.0);
	net.connect(net.eqns(0), net.vars(4), 1.0);

	net.connect(net.eqns(1), net.vars(1), -1.0);
	net.connect(net.eqns(1), net.vars(2), 1.0);
	net.connect(net.eqns(1), net.vars(3), -1.0);

	CHECK(!net.connected(net.eqns(0), net.vars(3)));
	CHECK(net.connected(net.eqns(0), net.vars(0)));
	CHECK(net.connected(net.eqns(1), net.vars(2)));
	CHECK(net.find_link(net.eqns(1), net.vars(4)) == nullptr);
	CHECK(net.find_link(net.eqns(0), net.vars(1))->s == -1.0);
	CHECK(net.M() == 2);
	CHECK(net.N() == 5);
	CHECK(net.E() == 6);

	net.disconnect(net.eqns(1), net.vars(2));
	CHECK(net.E() == 5);
	CHECK(!net.connected(net.eqns(1), net.vars(2)));
}

}}	// namespace Coss::bpB