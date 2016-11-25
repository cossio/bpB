#pragma once
#include <vector>
#include <string>
#include <cassert>

namespace Coss { namespace bpB {

class Net;
class Link;

class Var
{
	friend class Net;

public:
	Var(Var const&) = delete;
	Var(Var&&) = delete;
	Var& operator=(Var const&) = delete;
	Var& operator=(Var&&) = delete;

	Net* const owner;
	double const lb, ub;
	size_t const id;
	std::string const name;

	std::vector<Link*> const& links() const
	{ return _links; }

	size_t deg() const
	{ return _links.size(); };

private:
	std::vector<Link*> _links;

	Var(Net* owner, double lb, double ub, size_t id = 0, std::string name = "")
			: owner(owner), lb(lb), ub(ub), id(id), name(name)
	{ assert(owner != nullptr && lb <= ub); }
};

}}	// namespace Coss::bpB