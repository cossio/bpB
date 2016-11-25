#include <string>
#include <iostream>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <boost/program_options.hpp>
#include <Coss/VecD_io.h>
#include <Coss/py/CPyInstance.h>
#include "bpB/net/Net.h"
#include "bpB/net/Net_io.h"
#include "bpB/bp/BP_B.h"
#include "bpB/bp/BP_B_io.h"

// uncomment to disable assert()
// #define NDEBUG

int main(int const argc, char const* const argv[])
{
	// Turn off GSL error handling !! TODO: remove this later
	gsl_set_error_handler_off();

	std::string path;
	int iterations, seed;
	double tol, damp, p;
	bool regularize = false;

	boost::program_options::variables_map vm;
	boost::program_options::options_description desc("Program usage", 1024, 512);
	desc.add_options()
			("help,h", "produces this help message")
			("path", boost::program_options::value<std::string>(&path)->required(), "network path")
			("iterations,i", boost::program_options::value<int>(&iterations)->default_value(10000), "iterations")
			("tolerance,t", boost::program_options::value<double>(&tol)->default_value(1e-4), "used as stopping criteria")
			("damping,d", boost::program_options::value<double>(&damp)->default_value(0.3), "damping used in update_m")
			("msg,m", "if set, reads initial messages from path .msg")
			("sparse,s", "if set, reads stochiometric matrix in sparse representation")
			("prob,p", boost::program_options::value<double>(&p)->default_value(0.1), "0 < p <= 1 controls randomness of message passing schedule")
			("regularize,r", "if set, regularizes the bounds of m messages")
			("writemsg,w", "if set, writes messages at each iteration")
			("vols,v", "if set, outputs the volume computed at each iteration")
			("seed,e", boost::program_options::value<int>(&seed)->default_value(1), "seed for random generators");
	boost::program_options::positional_options_description p_desc;
	p_desc.add("path", 1);

	try
	{
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(desc).positional(p_desc).run(), vm);
		if (vm.count("help"))
		{
			std::cout << desc << std::endl;
			return 1;
		}
		boost::program_options::notify(vm);
	}
	catch (std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
	catch (...)
	{
		std::cerr << "Unknown error!" << std::endl;
		return 1;
	}

	if (iterations < 0)
	{
		std::cerr << "Error: iterations must be non-negative." << std::endl;
		return 1;
	}

	if (p <= 0 || p > 1)
	{
		std::cerr << "Error: p must be 0 <= p < 1." << std::endl;
		return 1;
	}

	Coss::bpB::Net net;

	if (vm.count("sparse"))
	{
		std::cout << "--sparse set. Reading sparse stoichiometric matrix representation." << std::endl;
		Coss::bpB::Net_io::read_sparse_net(path, net);
	}
	else
		Coss::bpB::Net_io::read_net(path, net);

	std::cout << "Stoichiometric matrix read with " << net.M() << " metabolites and " << net.N() << " reactions." << std::endl;

	Coss::bpB::BP_B bp(&net);
	bp.damp(damp).seed(seed).p(p);

	if (vm.count("regularize"))
		regularize = true;

	if (vm.count("msg"))
	{
		std::cout << "reading initialization messages" << std::endl;
		Coss::bpB::BP_B_io::read_M_msgs(path + ".msg", bp);
	}

	// initialize n messages from initial m messages
	bp.update_all_n();
	std::cout << "m and n messages initialized" << std::endl;

	if (vm.count("writemsg"))
		Coss::bpB::BP_B_io::export_msgs(bp, path + "_0.msg");

	if (vm.count("vols"))
	{
		auto entropies = bp.entropies();
		std::cout << "starting logvol=" << entropies.S << std::endl;
	}

	double max_change = 0;

	std::cout << "begin message passing. Using damp=" << damp << "." << std::endl;
	int itr;
	for (itr = 1; itr <= iterations; ++itr)
	{
		max_change = bp.update_all(regularize);
		if (vm.count("vols"))
		{
			auto entropies = bp.entropies();
			std::cout << "iteration " << itr << ", max_change=" << max_change
			<< ", logvol=" << entropies.S << std::endl;
		}
		else
			std::cout << "iteration " << itr << ", max_change=" << max_change
				 << std::endl;

		if (vm.count("writemsg"))
			Coss::bpB::BP_B_io::export_msgs(bp, path + "_" + std::to_string(itr) + ".msg");

		if (max_change <= tol)
			break;
	}

	std::cout << "exporting final messages to " << path + ".msg" << std::endl;
	Coss::bpB::BP_B_io::export_msgs(bp, path + ".msg");

	// compute volume
	std::cout << "computing entropies ..." << std::endl;
	auto entropies = bp.entropies();
	std::cout << "final S = " << entropies.S << std::endl;
	Coss::bpB::BP_B_io::export_Ss(bp, entropies, path);
	Coss::bpB::BP_B_io::export_Ks(bp, path);

	std::ofstream bpvol_file(path + ".bpvol");
	bpvol_file << entropies.S << std::endl;
	bpvol_file.close();

	std::ofstream logz_file(path + ".bplogz");
	logz_file << entropies.logZ << std::endl;
	logz_file.close();

	std::cout << "logvol written to " << path + ".bpvol" << std::endl;
	std::cout << "logZ written to " << path + ".bplogz" << std::endl;
}
