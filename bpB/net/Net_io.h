#pragma once
#include <fstream>
#include <string>
#include <Coss/Vec.h>
#include <Coss/VecD_io.h>
#include "Net.h"

namespace Coss { namespace bpB
{

struct Net_io
{
	// converts a stoichiometric matrix to a network of metabolites and reactions
	static void sto_to_net(Vec <VecD> const& sto, VecD const& lb, VecD const& ub, Net& net)
	{
		net.clear();

		size_t m = sto.size();
		size_t n = sto.back().size();
		assert(n == lb.size() && n == ub.size());
		assert(m < n);

		for (size_t i = 0; i < n; ++i)
			net.add_var(lb[i], ub[i], i);

		for (size_t a = 0; a < m; ++a)
			net.add_eqn(a);

		for (size_t a = 0; a < m; ++a)
			for (size_t i = 0; i < n; ++i)
				if (!almost_zero(sto[a][i]))
					net.connect(net.eqns(a), net.vars(i), sto[a][i]);
	}

	// reads a network of metabolites and reactions directly from the .sto, .lb and .ub files
	static void read_net(std::string const& path, Net& net)
	{
		VecD lb = Vec_io::read_vec_D(path + ".lb");
		VecD ub = Vec_io::read_vec_D(path + ".ub");
		assert(lb.size() == ub.size());
		size_t n = lb.size();
		Vec<VecD> mat = Vec_io::read_mat_D(path + ".sto");
		size_t m = mat.size();
		for (VecD const& row : mat)
			assert(row.size() == n);

		sto_to_net(mat, lb, ub, net);

		assert(net.N() > 0 && net.M() > 0);
	}

	// reads network from a sparse stochiometric matrix file
	static void read_sparse_net(std::string const& path, Net& net)
	{
		auto dat = Vec_io::read_mat_D(path + ".sto");
		size_t m = 0, n = 0;

		// determine dimensions of stochiometric matrix (m,n)
		for (auto line : dat)
		{
			assert(line[0] > 0 && line[1] > 0);
			if (line[0] > m)
				m = line[0];
			if (line[1] > n)
				n = line[1];
		}

		std::cout << "Sparse matrix has dimensions " << m << " x " << n << "." << std::endl;

		// make sure net is empty
		net.clear();

		// add equations
		for (size_t a = 0; a < m; ++a)
			net.add_eqn(a);

		// add variables
		auto lb = Vec_io::read_vec_D(path + ".lb");
		auto ub = Vec_io::read_vec_D(path + ".ub");
		assert(n == lb.size() && n == ub.size());
		for (size_t i = 0; i < n; ++i)
			net.add_var(lb[i], ub[i], i);

		for (auto line : dat)
		{
			size_t a = round(line[0]) - 1;
			size_t i = round(line[1]) - 1;
			assert(a < m && i < n);
			double s = line[2];
			assert(s != 0.0);
			net.connect(net.eqns()[a], net.vars()[i], s);
		}
	}
};

}}	// namespace Coss::bpB