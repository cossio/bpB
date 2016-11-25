#pragma once
#include <map>
#include "../net/Net.h"

namespace Coss { namespace bpB {

struct Entropies
{
	std::map<Eqn const*, double> Sa;
	std::map<Var const*, double> Si;
	double S = 0, logZ = 0;
};

}}	// namespace Coss::bpB