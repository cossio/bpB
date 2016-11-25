#pragma once
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include "../net/Net.h"
#include "../net/Net_io.h"
#include "Entropies.h"

namespace Coss { namespace bpB {

class BP_B_io
{
public:

	static void read_M_msgs(std::string const& path, BP_B& bp_beta)
	{
		std::cout << "reading input messages..." << std::endl;

		std::ifstream file(path);
		std::string line;
		while (getline(file, line))
		{
			std::istringstream iss(line);
			size_t a, i;
			iss >> a >> i;
			// find a->i link
			auto eqn = bp_beta.net->eqns(a);
			auto var = bp_beta.net->vars(i);
			assert(bp_beta.net->connected(eqn, var));
			auto link = bp_beta.net->find_link(eqn, var);
			MsgM& m = bp_beta.msg_m(link);
			double A, B, alpha, beta;
			iss >> A >> B >> alpha >> beta;
			m = MsgM(m.lb(), m.ub(), A, B, alpha, beta);

			for (Link* e : var->links())
				if (e->eqn != eqn)
					bp_beta.msg_n(e).needs_update = true;
		}
		file.close();
	}

	static void export_msgs(BP_B const& bp_beta, std::string const& path)
	{
		std::ofstream file(path);
		auto net = bp_beta.net;
		for (size_t a = 0; a < net->M(); ++a)
		{
			for (size_t i = 0; i < net->N(); ++i)
			{
				// find link (if it exists) between met(a) and rxn(i)

				Link const* link = net->find_link(net->eqns(a), net->vars(i));
				if (link == nullptr)   // there is no link between eqn(a) and var(i)
					continue;
				MsgM const& m = bp_beta.msg_m(link);
				file << a << "\t" << i << "\t" <<
						m.A() << "\t" << m.B() << "\t" <<
						m.alpha() << "\t" << m.beta() << std::endl;
			}
		}
		file.close();
	}

	static void export_Ks(BP_B const& bp_beta, std::string const& path)
	{
		std::ofstream file(path + ".Ki");
		auto net = bp_beta.net;
		for (Var const* var : net->vars())
			file << bp_beta.logK(var) << std::endl;
		file.close();

//		file.open(path + ".Kia");
//
//		for (Met const* met : net->mets())
//			for (Link const* e : met->links())
//				file << met->id<< "\t" << e->rxn->id << "\t" << bp_beta.logKia(e) << std::endl;
//		file.close();

//		file.open(path + ".Kai");
//		for (Met const* met : net->mets())
//			for (Link const* e : met->links())
//				file << met->id << "\t" << e->rxn->id << "\t" << bp_beta.logKai(e) << std::endl;
//		file.close();
	}

	static void export_Ss(BP_B const& bp, Entropies const& v, std::string const& path)
	{
		std::ofstream file(path + ".Si");
		for (Var const* var : bp.net->vars())
			file << v.Si.at(var) << std::endl;
		file.close();

		file.open(path + ".Sa");
		for (Eqn const* eqn : bp.net->eqns())
			file << v.Sa.at(eqn) << std::endl;
		file.close();
	}

};

}}	// namespace Coss::bpB