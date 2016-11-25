#pragma once
#include <cassert>

namespace Coss { namespace bpB
{

class Net;
class Var;
class Eqn;

class Link
{
	friend class Net;

public:
	Net* const owner;
	Eqn* const eqn;
	Var* const var;
	double const s;

	// uncopyable and unmovable
	Link(Link const&) = delete;
	Link(Link&&) = delete;
	Link& operator=(Link const&) = delete;
	Link& operator=(Link&&) = delete;

private:
	Link(Net* owner, Eqn* eqn, Var* var, double s)
			: owner(owner), eqn(eqn), var(var), s(s)
	{ assert(owner != nullptr && eqn != nullptr && var != nullptr && s != 0.0); }
};

}}	// namespace Coss::bpB