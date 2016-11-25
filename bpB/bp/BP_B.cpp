#include <catch.hpp>
#include "BP_B.h"
#include "Entropies.h"

namespace Coss { namespace bpB {

TEST_CASE("beta_characteristic", "[bp]")
{
	CHECK(BP_B::beta_characteristic_function(500, 300, -40, 20, 300).real() == Approx(-2.141068745361687097131248403093826635525690e-436).epsilon(1e-7));
	CHECK(BP_B::beta_characteristic_function(500, 300, -40, 20, 300).imag() == Approx(9.8815339224478225949115503857210156333358465e-436).epsilon(1e-7));
	CHECK(BP_B::beta_characteristic_function(500, 300, -40, 20, 5).real() == Approx(9.682682724656784859630505430828461767365544186483209e-7).epsilon(1e-7));
	CHECK(BP_B::beta_characteristic_function(500, 300, -40, 20, 5).imag() == Approx(1.310091237980244460854584978564166729791167515435863e-6).epsilon(1e-7));
	CHECK(BP_B::beta_characteristic_function(100, 300, -40, 20, 5).real() == Approx(4.679833228128399905484669050375702374970962614730019e-10).epsilon(1e-7));
	CHECK(BP_B::beta_characteristic_function(100, 300, -40, 20, 5).imag() == Approx(1.436069514363859144778660683555687051225986409761646e-9).epsilon(1e-7));
}

TEST_CASE("update_m", "[bp]")
{
	Net net;
	net.add_eqn(0);
	net.add_var(0.0, 1.0, 0);
	net.add_var(0.0, 1.0, 1);
	net.add_var(0.0, 1.0, 2);
	net.add_var(0.0, 1.0, 3);
	Link* link0 = net.connect(net.eqns(0), net.vars(0), -1.2);
	Link* link1 = net.connect(net.eqns(0), net.vars(1), -2.0);
	Link* link2 = net.connect(net.eqns(0), net.vars(2), 1.0);
	Link* link3 = net.connect(net.eqns(0), net.vars(3), 4.2);
	CHECK(net.M() == 1);
	CHECK(net.N() == 4);
	CHECK(net.E() == 4);

	BP_B bp(&net);
	bp.damp(0).seed(0).p(1);

	Peak p(0.5, 0.0, 1.0, 10.2, .001);

	bp.msg_n(link1) = MsgN(0.2, 1, 0.44444444444444453, 0.024691358024691353, -5.634789603169249, p);
	bp.msg_n(link2) = MsgN(0.2, 1, 0.2986223321664844, 0.05324095444000649, 10.533040446903026, p);
	bp.msg_n(link3) = MsgN(0.2, 1, 0.6444404558915912, 0.04395343141366635, 2.420861833918419, p);

	double A = -(link1->s * bp.msg_n(link1).right() + link2->s * bp.msg_n(link2).left() + link3->s * bp.msg_n(link3).left()) / link0->s;
	double B = -(link1->s * bp.msg_n(link1).left() + link2->s * bp.msg_n(link2).right() + link3->s * bp.msg_n(link3).right()) / link0->s;

	SECTION("regularize=false")
	{
		bp.update_m(link0);
		MsgM const& m = bp.msg_m(link0);
		CHECK(m.A() == Approx(A).epsilon(1e-10));
		CHECK(m.B() == Approx(B).epsilon(1e-10));
		CHECK(m.mean() == Approx(1.7636527983518988).epsilon(1e-10));
		CHECK(m.var() == Approx(0.6439895254693376).epsilon(1e-10));
		CHECK(m.left() == std::max(m.A(), m.lb()));
		CHECK(m.right() == std::min(m.B(), m.ub()));
		CHECK(m.lb() == net.vars(0)->lb);
		CHECK(m.ub() == net.vars(0)->ub);
	}

	bp.update_m(link0);
	MsgM tmp = bp.msg_m(link0);

	SECTION("regularize=true")
	{
		bp.update_m(link0, true);
		MsgM const& m = bp.msg_m(link0);
		CHECK(m.mean() == Approx(tmp.mean_trunc()).epsilon(1e-10));
		CHECK(m.var() == Approx(tmp.var_trunc()).epsilon(1e-10));
		CHECK(m.A() == Approx(tmp.left()).epsilon(1e-10));
		CHECK(m.B() == Approx(tmp.right()).epsilon(1e-10));
		CHECK(m.left() == m.A());
		CHECK(m.right() == m.B());
		CHECK(m.lb() == net.vars(0)->lb);
		CHECK(m.ub() == net.vars(0)->ub);
	}
}

TEST_CASE("update_n", "[bp]")
{
	Net net;
	net.add_eqn(0);
	net.add_eqn(1);
	net.add_eqn(2);
	net.add_eqn(3);
	net.add_var(0.0, 1.0, 0);
	Link* link0 = net.connect(net.eqns(0), net.vars(0), -1.0);
	Link* link1 = net.connect(net.eqns(1), net.vars(0), -1.0);
	Link* link2 = net.connect(net.eqns(2), net.vars(0), 1.0);
	Link* link3 = net.connect(net.eqns(3), net.vars(0), 1.0);
	CHECK(net.M() == 4);
	CHECK(net.N() == 1);
	CHECK(net.E() == 4);

	BP_B bp(&net);
	bp.seed(0).damp(0).p(1);

	bp.msg_m(link1) = MsgM(0, 1, -1, 1, 2, 2);
	bp.msg_m(link2) = MsgM(0, 1, 0.2, 1, 4, 5);
	bp.msg_m(link3) = MsgM(0, 1, 0, 1, 1.5, 3);

	bp.update_n(link0);
	
	MsgN const& n0 = bp.msg_n(link0);

	CHECK(n0.left() == Approx(0.2).epsilon(1e-15));
	CHECK(n0.right() == Approx(1.0).epsilon(1e-15));

	CHECK(n0.mean_trunc() == Approx(0.4859492982772886).epsilon(1e-10));
	CHECK(n0.var_trunc() == Approx(0.011334229149905858).epsilon(1e-10));
	CHECK(n0.log_norm_trunc() == Approx(-9.636365226156267).epsilon(1e-10));
}

TEST_CASE("1_eqn", "[bp]")
{
	Net net;
	net.add_eqn(0);
	net.add_var(0.0, 1.0, 0);
	net.add_var(0.0, 1.0, 1);
	net.connect(net.eqns(0), net.vars(0), 1.0);
	net.connect(net.eqns(0), net.vars(1), -1.0);
	CHECK(net.M() == 1);
	CHECK(net.N() == 2);
	CHECK(net.E() == 2);

	BP_B bp(&net);
	bp.seed(0).damp(0).p(1);
	double total_change = 1;
	size_t iterations = 0;
	for (; total_change > 1e-6 && iterations < 100; ++iterations)
		total_change = bp.update_all();
	CHECK(iterations < 100);

	MsgM const& m0 = bp.msg_m(net.find_link(net.eqns(0), net.vars(0)));
	MsgM const& m1 = bp.msg_m(net.find_link(net.eqns(0), net.vars(1)));
	CHECK(m0.alpha() == Approx(1.0).epsilon(1e-2));
	CHECK(m0.beta() == Approx(1.0).epsilon(1e-2));

	Entropies entropies = bp.entropies();
	CHECK(entropies.S == Approx(std::log(std::sqrt(2.0))).epsilon(1e-5));
	CHECK(entropies.logZ == Approx(0.0).epsilon(1e-5));

//	Entropies entropies_noFT = bp.entropies(false);
//	CHECK(entropies_noFT.S == Approx(std::log(std::sqrt(2.0))).epsilon(1e-5));
//	CHECK(entropies_noFT.logZ == Approx(0.0).epsilon(1e-5));
}

TEST_CASE("small tree", "[bp]")
{
	Net net;
	net.add_eqn(0);
	net.add_eqn(1);
	net.add_var(0, 1, 0);
	net.add_var(0, 1, 1);
	net.add_var(0, 1, 2);
	net.add_var(0, 1, 3);
	net.connect(net.eqns(0), net.vars(0), 1.0);
	net.connect(net.eqns(0), net.vars(1), 1.0);
	net.connect(net.eqns(0), net.vars(2), -1.0);
	net.connect(net.eqns(1), net.vars(2), 1.0);
	net.connect(net.eqns(1), net.vars(3), -1.0);
	CHECK(net.M() == 2);
	CHECK(net.N() == 4);
	CHECK(net.E() == 5);

	BP_B bp(&net);
	bp.seed(0).damp(0).p(1);
	double total_change = 1;
	size_t iterations = 0;
	for (; total_change > 1e-6 && iterations < 100; ++iterations)
		total_change = bp.update_all();
	CHECK(iterations < 100);

	MsgM const& m0 = bp.msg_m(net.find_link(net.eqns(0), net.vars(0)));
	MsgM const& m1 = bp.msg_m(net.find_link(net.eqns(0), net.vars(1)));
	MsgM const& m2 = bp.msg_m(net.find_link(net.eqns(0), net.vars(2)));
	MsgM const& m3 = bp.msg_m(net.find_link(net.eqns(1), net.vars(2)));
	MsgM const& m4 = bp.msg_m(net.find_link(net.eqns(1), net.vars(3)));

	Entropies entropies = bp.entropies();
	CHECK(entropies.S == Approx(std::log(std::sqrt(5.0) / 2.0)).epsilon(1e-5));
	CHECK(entropies.logZ == Approx(log(0.5)).epsilon(1e-5));

//	Entropies entropies_noFT = bp.entropies(false);
//	CHECK(entropies_noFT.S == Approx(std::log(std::sqrt(5.0) / 2.0)).epsilon(1e-45));
//	CHECK(entropies_noFT.logZ == Approx(log(0.5)).epsilon(1e-5));
}

TEST_CASE("1 eqn regularize", "[bp]")
{
	Net net;
	net.add_eqn(0);
	net.add_var(0.0, 1.0, 0);
	net.add_var(0.0, 1.0, 1);
	net.connect(net.eqns(0), net.vars(0), 1.0);
	net.connect(net.eqns(0), net.vars(1), -1.0);
	CHECK(net.M() == 1);
	CHECK(net.N() == 2);
	CHECK(net.E() == 2);

	BP_B bp(&net);
	bp.seed(0).damp(0).p(1);
	double total_change = 1;
	size_t iterations = 0;
	for (; total_change > 1e-6 && iterations < 100; ++iterations)
		total_change = bp.update_all(true);
	CHECK(iterations < 100);

	MsgM const& m0 = bp.msg_m(net.find_link(net.eqns(0), net.vars(0)));
	MsgM const& m1 = bp.msg_m(net.find_link(net.eqns(0), net.vars(1)));
	CHECK(m0.alpha() == Approx(1.0).epsilon(1e-2));
	CHECK(m0.beta() == Approx(1.0).epsilon(1e-2));

	Entropies entropies = bp.entropies();
	CHECK(entropies.S == Approx(std::log(std::sqrt(2.0))).epsilon(1e-5));
	CHECK(entropies.logZ == Approx(0.0).epsilon(1e-5));

//	Entropies entropies_noFT = bp.entropies(false);
//	CHECK(entropies_noFT.S == Approx(std::log(std::sqrt(2.0))).epsilon(1e-5));
//	CHECK(entropies_noFT.logZ == Approx(0.0).epsilon(1e-5));
}

TEST_CASE("small tree, regularize", "[bp]")
{
	Net net;
	net.add_eqn(0);
	net.add_eqn(1);
	net.add_var(0, 1, 0);
	net.add_var(0, 1, 1);
	net.add_var(0, 1, 2);
	net.add_var(0, 1, 3);
	net.connect(net.eqns(0), net.vars(0), 1.0);
	net.connect(net.eqns(0), net.vars(1), 1.0);
	net.connect(net.eqns(0), net.vars(2), -1.0);
	net.connect(net.eqns(1), net.vars(2), 1.0);
	net.connect(net.eqns(1), net.vars(3), -1.0);
	CHECK(net.M() == 2);
	CHECK(net.N() == 4);
	CHECK(net.E() == 5);

	BP_B bp(&net);
	bp.seed(0).damp(0).p(1);
	double total_change = 1;
	size_t iterations = 0;
	for (; total_change > 1e-6 && iterations < 200; ++iterations)
		total_change = bp.update_all(true);
	CHECK(iterations < 200);

	Entropies entropies = bp.entropies();
	CHECK(entropies.S == Approx(std::log(std::sqrt(5.0) / 2.0)).epsilon(1e-5));
	CHECK(entropies.logZ == Approx(log(0.5)).epsilon(1e-5));

//	Entropies entropies_noFT = bp.entropies(false);
//	CHECK(entropies_noFT.S == Approx(std::log(std::sqrt(5.0) / 2.0)).epsilon(1e-5));
//	CHECK(entropies_noFT.logZ == Approx(log(0.5)).epsilon(1e-5));
}

}}	// namespace Coss::bpB