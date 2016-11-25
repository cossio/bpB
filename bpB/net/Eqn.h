#pragma once
#include <vector>
#include <string>

namespace Coss { namespace bpB {

class Net;
class Link;

class Eqn
{
	friend class Net;

private:
	std::vector<Link*> _links;

	explicit Eqn(Net* owner, size_t id = 0, std::string name = "")
			: owner(owner), id(id), name(name)
	{ }

public:
	Eqn(Eqn const&) = delete;
	Eqn(Eqn&&) = delete;
	Eqn& operator=(Eqn const&) = delete;
	Eqn& operator=(Eqn&&) = delete;

	Net* const owner;
	size_t const id;
	std::string const name;

	inline std::vector<Link*> const& links() const
	{ return _links; }

	inline size_t deg() const
	{ return _links.size(); }
};

}}	// namespace Coss::bpB