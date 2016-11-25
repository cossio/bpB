#pragma once
#include <algorithm>
#include <vector>
#include <cassert>
#include "Link.h"
#include "Eqn.h"
#include "Var.h"

namespace Coss { namespace bpB {

class Net
{
private:
	bool _frozen = false;
public:
	Net() = default;

	~Net()
	{
		for (auto eqn : _eqns)
			delete eqn;
		for (auto var : _vars)
			delete var;
	}

	// uncopyable and unmovable
	Net(Net const&) = delete;
	Net(Net&&) = delete;
	Net& operator=(Net const&) = delete;
	Net& operator=(Net&&) = delete;

	std::vector<Eqn*> const& eqns() const
	{ return _eqns; }

	std::vector<Var*> const& vars() const
	{ return _vars; }

	Eqn* eqns(size_t a) const
	{ return _eqns.at(a); }

	Var* vars(size_t i) const
	{ return _vars.at(i); }

	bool frozen() const
	{ return _frozen; }

	Eqn* add_eqn(size_t id = 0, std::string name = "")
	{
		assert(!_frozen);
		auto eqn = new Eqn(this, id, name);
		_eqns.push_back(eqn);
		return eqn;
	}

	Var* add_var(double lb = 0, double ub = 1, size_t id = 0, std::string name = "")
	{
		auto var = new Var(this, lb, ub, id, name);
		assert(var->owner == this);
		_vars.push_back(var);
		return var;
	}

	void remove_eqn(Eqn* eqn)
	{
		assert(eqn != nullptr && eqn->owner == this);
		while (!eqn->links().empty())
			disconnect(eqn, eqn->links().back()->var);
		auto itr = std::find(_eqns.begin(), _eqns.end(), eqn);
		assert(itr != _eqns.end());
		_eqns.erase(itr);
		delete eqn;
	}

	void remove_var(Var* var)
	{
		assert(var != nullptr && var->owner == this);
		while (!var->links().empty())
			disconnect(var->links().back()->eqn, var);
		auto itr = std::find(_vars.begin(), _vars.end(), var);
		assert(itr != _vars.end());
		_vars.erase(itr);
		delete var;
	}

	Link* connect(Eqn* eqn, Var* var, double s)
	{
		assert(eqn != nullptr && var != nullptr && eqn->owner == this && var->owner == this);
		disconnect(eqn, var);
		Link* link = new Link(this, eqn, var, s);
		eqn->_links.push_back(link);
		var->_links.push_back(link);
		return link;
	}

	void disconnect(Eqn* eqn, Var* var)
	{
		assert(eqn != nullptr && var != nullptr && eqn->owner == this && var->owner == this);
		Link const* link = find_link(eqn, var);
		if (link == nullptr)
			return;
		eqn->_links.erase(std::remove(eqn->_links.begin(), eqn->_links.end(), link), eqn->_links.end());
		var->_links.erase(std::remove(var->_links.begin(), var->_links.end(), link), var->_links.end());
		delete link;
	}

	Link const* find_link(Eqn const* eqn, Var const* var) const
	{
		assert(eqn != nullptr && var != nullptr && eqn->owner == this && var->owner == this);
		auto itr = std::find_if(eqn->_links.begin(), eqn->_links.end(), [var](Link const* e) { return e->var == var; });
		return itr == eqn->_links.end() ? nullptr : *itr;
	}

	Link* find_link(Eqn const* eqn, Var const* var)
	{
		assert(eqn != nullptr && var != nullptr && eqn->owner == this && var->owner == this);
		auto itr = find_link_itr(eqn, var);
		return itr == eqn->_links.end() ? nullptr : *itr;
	}

	bool connected(Eqn const* eqn, Var const* var) const
	{
		assert(eqn != nullptr && var != nullptr);
		assert(eqn->owner == this && var->owner == this);
		return find_link(eqn, var) != nullptr;
	}

	bool empty() const
	{ return eqns().empty() && vars().empty(); }

	size_t M() const
	{ return eqns().size(); }

	size_t N() const
	{ return vars().size(); }

	size_t E() const
	{
		size_t e = 0;
		for (auto eqn : eqns())
			e += eqn->deg();
		return e;
	}

	void clear()
	{
		while (!eqns().empty())
			remove_eqn(_eqns.back());
		while (!vars().empty())
			remove_var(_vars.back());
	}

private:
	std::vector<Eqn*> _eqns;
	std::vector<Var*> _vars;

	std::vector<Link*>::const_iterator find_link_itr(Eqn const* eqn, Var const* var) const
	{
		assert(eqn != nullptr && var != nullptr);
		assert(eqn->owner == this && var->owner == this);
		return std::find_if(eqn->_links.begin(), eqn->_links.end(), [var](Link const* e) { return e->var == var; });
	}
};

}}	// namespace Coss::bpB
