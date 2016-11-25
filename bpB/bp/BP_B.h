#pragma once
#include <map>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <Coss/Beta.h>
#include <Coss/gsl/Integrate_QAGP.h>
#include <Coss/Peak.h>
#include <Coss/FT_conv.h>
#include <Coss/log_conv.h>
#include <Coss/py/Hyp1F1.h>
#include "../net/Net.h"
#include "MsgN.h"
#include "MsgM.h"
#include "Entropies.h"

namespace Coss { namespace bpB {

class BP_B
{
private:
	double _damp = 0.4;
	size_t _seed = 0;

	boost::random::mt19937 _rand_eng;
	boost::random::bernoulli_distribution<double> _dist;

	gsl::Integrate_QAGP const _I;
	Coss::FT_conv ftconv;
	Coss::gsl::FT_QAWO ft;

	std::map<Link const* const, MsgM> _msgs_m;
	std::map<Link const* const, MsgN> _msgs_n;

public:

	constexpr static double min_ab = 1.0 + 1e-4;
	constexpr static double max_lambda = 2e4;

	Net const* const net;

	explicit BP_B(Net const* net, size_t limit = 1000000, size_t levels = 50)
			: ftconv(limit, levels), net(net), ft(limit, levels)
	{
		assert(net != nullptr);

		seed(_seed);
		p(1.0);

		for (auto rxn : net->vars())
		{
			for (auto link : rxn->links())
			{
				_msgs_m.insert({link, MsgM(rxn->lb, rxn->ub)});
				_msgs_n.insert({link, MsgN(rxn->lb, rxn->ub)});
			}
		}
	}

	double update_m(Link const* link, bool regularize = false)
	{
		// because we aply damping here, we cannot rely on a needs_update_m variable.
		// even if other messages connected to this node don't change, the damping means
		// that calling this function again can make a change. Therefore we always update_m.
		// and we always set _needs_update_m to true.
		// We could set _needs_update_m = false for damp = 0, but we don't do this for consistency
		// and easier debugging.

		if (link->eqn->links().size() < 2)
			return 0;

		double A, B, mean, var;
		A = B = mean = var = 0.0;

		for (Link const* e : link->eqn->links())
		{
			if (e == link)
				continue;

			double r = e->s / link->s;
			assert(Coss::is_finite(r) && std::abs(r) > 0);
			mean -= r * msg_n(e).mean_trunc();
			var += r * r * msg_n(e).var_trunc();
			assert(msg_n(e).left() < msg_n(e).right());
			A -= r > 0 ? r * msg_n(e).right() : r * msg_n(e).left();
			B -= r > 0 ? r * msg_n(e).left() : r * msg_n(e).right();
			assert(A < B);
			assert(A < mean && mean < B);
			assert(var < (mean - A) * (B - mean));
			assert(var > 0);
		}

		assert(std::max(msg_m(link).lb(), A) < std::min(msg_m(link).ub(), B));

		// apply damping
		mean = (1 - damp()) * mean + damp() * msg_m(link).mean();
		var = (1 - damp()) * var + damp() * msg_m(link).var();
		A = (1 - damp()) * A + damp() * msg_m(link).A();
		B = (1 - damp()) * B + damp() * msg_m(link).B();

		double lambda = std::min((mean - A) * (B - mean) / var - 1, max_lambda);
		assert(Coss::is_finite(lambda));
		assert(0 < lambda && lambda <= max_lambda);

		// force a,b >= min_ab. Note that this may change mean, var
		double a = std::max(lambda * (mean - A) / (B - A), min_ab);
		double b = std::max(lambda * (B - mean) / (B - A), min_ab);
		assert(a > 1 && b > 1);

		//update m
		MsgM m_old = msg_m(link);
		assert(m_old.alpha() > 1 || m_old.beta() > 1);

		msg_m(link) = MsgM(msg_m(link).lb(), msg_m(link).ub(), A, B, a, b);

		if (regularize)
		{
			// sets A >= lb, B <= ub
			msg_m(link).regularize();

//			if (m.alpha() < 1.5)
//				a = 0.9 * m_old.alpha() + 0.1 * m.alpha();
//			if (m.beta() < 1.5)
//				b = 0.9 * m_old.beta() + 0.1 * m.beta();
			assert(msg_m(link).A() >= msg_m(link).lb() && msg_m(link).B() <= msg_m(link).ub());
		}
		else
			assert(msg_m(link).alpha() >= 1 && msg_m(link).beta() >= 1);

		//m = MsgM(m.lb(), m.ub(), m.A(), m.B(), a, b);

		double change = std::abs(m_old.A() - msg_m(link).A()) + std::abs(m_old.B() - msg_m(link).B()) +
						std::abs(m_old.mean() - msg_m(link).mean()) + std::abs(sqrt(m_old.var()) - std::sqrt(msg_m(link).var()));

		for (Link* e : link->var->links())
			if (e != link)
				msg_n(e).needs_update = true;

		return change;
	}

	void update_n(Link const* link)
	{
		assert(link != nullptr);

		if (link->var->links().size() < 2)
			return;

		double left = link->var->lb;
		double right = link->var->ub;
		double Amax = -INF;
		double Bmin = INF;

		for (Link* e : link->var->links())
		{
			if (e == link)
				continue;
			left  = std::max(left,  msg_m(e).left());
			right = std::min(right, msg_m(e).right());
			assert(left < right);
			assert(link->var->lb <= left && right <= link->var->ub);
			Amax = std::max(Amax, msg_m(e).A());
			Bmin = std::min(Bmin, msg_m(e).B());
			assert(Amax < Bmin);
		}

		std::function<double(double)> log_f = [this, link](double x) { return log_n(link, x); };
		assert(log_f(Amax) == log(0) && log_f(Bmin) == log(0));
		Peak peak = Peak::find_log_peak(log_f, Amax, Bmin);

		double log_norm_trunc = _I.log_integrate(log_f, peak, left, right);
		double logM1_trunc = _I.log_integrate([log_f, Amax](double x) { return log_f(x) + log(x - Amax); }, peak, left, right) - log_norm_trunc;
		double logM2_trunc = _I.log_integrate([log_f, Amax](double x) { return log_f(x) + 2 * log(x - Amax); }, peak, left, right) - log_norm_trunc;

		double mean_trunc = exp(logM1_trunc) + Amax;
		double var_trunc = exp(log_sub(logM2_trunc, 2 * logM1_trunc));

		assert(left < mean_trunc && mean_trunc < right);
		assert(var_trunc > 0);
		assert(var_trunc < (right - mean_trunc) * (mean_trunc - left));
		assert(is_finite(log_norm_trunc));

		msg_n(link) = MsgN(left, right, mean_trunc, var_trunc, log_norm_trunc, peak);
		msg_n(link).needs_update = false;
	}

	void update_all_n()
	{
		for (Eqn* eqn : net->eqns())
			for (Link* e : eqn->links())
				if (msg_n(e).needs_update)
					update_n(e);
	}

	double update_all_m(bool regularize = false)
	{
		double max_change = 0;
		for (Eqn* eqn : net->eqns())
			for (Link* e : eqn->links())
				if (_dist(_rand_eng))
					max_change = std::max(update_m(e, regularize), max_change);
		return max_change;
	}

	double update_all(bool regularize = false)
	{
		double max_change = 0;
		for (std::size_t iter = 0; iter < updates_per_iteration(); ++iter)
		{
			max_change = std::max(max_change, update_all_m(regularize));
			update_all_n();
		}
		return max_change;
	}

	MsgM& msg_m(Link const* link)
	{
		assert(link != nullptr && link->owner == net);
		return _msgs_m.at(link);
	}

	MsgN& msg_n(Link const* link)
	{
		assert(link != nullptr && link->owner == net);
		return _msgs_n.at(link);
	}

	MsgM const& msg_m(Link const* link) const
	{
		assert(link != nullptr && link->owner == net);
		return _msgs_m.at(link);
	}

	MsgN const& msg_n(Link const* link) const
	{
		assert(link != nullptr && link->owner == net);
		return _msgs_n.at(link);
	}

	double damp() const
	{ return _damp; }

	double p() const
	{ return _dist.p(); }

	std::size_t seed() const
	{ return _seed; }

	std::size_t updates_per_iteration() const
	{ return static_cast<std::size_t>(std::ceil(1.0 / p())); }

	BP_B& damp(double value)
	{
		assert(0 <= value && value < 1);
		_damp = value;
		return *this;
	}

	BP_B& p(double value)
	{
		assert(0 < value && value <= 1);
		_dist.param(boost::random::bernoulli_distribution<double>::param_type(value));
		return *this;
	}

	BP_B& seed(size_t value)
	{
		_rand_eng = boost::random::mt19937(value);
		_seed = value;
		return *this;
	}

	double logZ_formulita() const
	{
		double sum = 0;

		for (Eqn const* eqn : net->eqns())
		{
			for (Link const* e : eqn->links())
			{
				sum += (1.0 / eqn->deg()) * logKai(e);
				sum += (1.0 - 1.0 / eqn->deg()) * logKia(e);
			}
		}

		for (Var const* var : net->vars())
		{
			double invA = 0;    // sum 1/da
			for (Eqn const* met : net->eqns())
				invA += 1.0 / met->deg();

			sum -= (var->deg() - invA - 1.0) * logK(var);
		}

		double logvol1 = sum;

		sum = 0;

		for (Eqn const* eqns : net->eqns())
		{
			sum += logK(eqns);
			assert(Coss::is_finite(sum));

			for (Link* e : eqns->links())
			{
				sum += logKia(e);
				assert(Coss::is_finite(sum));
			}
		}

		for (Var const* var : net->vars())
		{
			// 1.0 -> convert size_t (which is unsigned) to signed value
			sum += (1.0 - var->deg()) * logK(var);
			assert(Coss::is_finite(sum));

		}

		assert(Coss::is_finite(sum));

		return sum;
	}

	Entropies entropies(bool use_FT = true, bool use_1F1 = true) const
	{
		Entropies R;

		for (Var* var : net->vars())
			R.Si[var] = S(var);

		for (Eqn* eqn : net->eqns())
		{
			R.Sa[eqn] = S(eqn, use_FT, use_1F1);
			std::cout << "Sa(" << eqn->id << ") = " << R.Sa[eqn] << std::endl;
		}

		for (Eqn* eqn : net->eqns())
			R.S += R.Sa[eqn];
		for (Var* var : net->vars())
			R.S -= (var->links().size() - 1) * R.Si[var];

		R.logZ = R.S;
		for (Eqn* eqn : net->eqns())
			R.logZ -= std::log(pi(eqn));

		return R;
	}

	double Amax(Var const* var) const
	{
		double Amax = -INF;
		for (Link* e : var->links())
			Amax = std::max(Amax, msg_m(e).A());
		return Amax;
	}

	double Bmin(Var const* var) const
	{
		double Bmin = INF;
		for (Link* e : var->links())
			Bmin = std::min(Bmin, msg_m(e).B());
		return Bmin;
	}

	double left(Var const* var) const
	{ return std::max(Amax(var), var->lb); }

	double right(Var const* var) const
	{ return std::min(Bmin(var), var->ub); }

	double S(Var const* var) const
	{
		Peak peak = Peak::find_log_peak([this, var](double x) { return log_b_unnorm(var, x); }, Amax(var), Bmin(var));
		double logKi = logK(var);
		Re_func integrand = [this, var, logKi](double x) {
			return exp(log_b_unnorm(var, x) - logKi) * (log_b_unnorm(var, x) - logKi);
		};
		double result = -_I.integrate(integrand, peak, left(var), right(var));
		assert(is_finite(result));
		return result;
	}

    double S(Eqn const* eqn, bool use_FT = true, bool use_1F1 = true) const
    {
		double result = 0;

		for (Link const* e : eqn->links())
		{
			double logKi = logK(e->var);
			Peak peak = Peak::find_log_peak([this, e](double x) { return log_b_unnorm(e->var, x); },
											Amax(e->var), Bmin(e->var));
			result -= _I.integrate([this, e, logKi](double x) { return exp(log_b_unnorm(e->var, x) - logKi) * log_n(e, x); },
								   peak, left(e->var), right(e->var));
		}

		return result + logK(eqn, use_FT, use_1F1);
    }

	double logKia(Link const* link) const
	{
		double result = msg_n(link).log_norm_trunc();
		for (Link* e: link->var->links())
			if (e != link)
				result -= msg_m(e).log_norm_trunc();
		assert(is_finite(result));
		return result;
	}

	double logKai(Link const* link) const
	{
		MsgM const& m = msg_m(link);
		return Coss::Beta::log_beta(m.alpha(), m.beta(), m.left_01_scaled(), m.right_01_scaled()) -
				Coss::Beta::log_beta(m.alpha(), m.beta()) - std::log(std::abs(link->s));
	}

	double logK(Var const* var) const
	{
		assert(var != nullptr);
		assert(var->links().size() > 0);
		auto log_f = [this, var](double x) { return log_b_unnorm(var, x); };
		Peak peak = Peak::find_log_peak(log_f, Amax(var), Bmin(var));
		double result = _I.log_integrate(log_f, peak, left(var), right(var));
		assert(is_finite(result));
		return result;
	}

	double logK(Eqn const* eqn, bool use_FT = true, bool use_1F1 = true) const
	{ return use_FT ? logK_FT(eqn, use_1F1) : logK_logconv(eqn); }

	double logK_FT(Eqn const* eqn, bool use_1F1 = true) const
	{ return use_1F1 ? logK_FT_1F1(eqn) : logK_FT_conv(eqn); }

	double logK_FT_1F1(Eqn const* eqn) const
	{
		assert(eqn->links().size() > 0);
		std::function<std::complex<double>(double)> FT = [this, eqn](double w) {
			std::complex<double> result = 1.0;
			for (Link const* link : eqn->links())
				result *= std::abs(link->s) * msg_n_FT(link, link->s * w);
			assert(is_finite(result.real()) && is_finite(result.imag()));
			return result;
		};

		std::function<Log_Num(double)> FT_re = [FT] (double w) { return Log_Num(FT(w).real()); };
		std::function<Log_Num(double)> FT_im = [FT] (double w) { return Log_Num(FT(w).imag()); };

		Log_Num result = ft.FT_inv(FT_re, FT_im, 0.0, 1e-7);
		//assert(result.sign() == PositiveSign);

		return result.log_abs() + std::log(pi(eqn)) + log_norm_prod(eqn);
	}

	double logK_FT_conv(Eqn const* eqn) const
	{
		assert(eqn->links().size() > 0);
		std::vector<Bounded_peak_log_func> n_y;
		for (Link const* link : eqn->links())
		{
			double A = msg_n_A_y(link), B = msg_n_B_y(link);
			assert(A < B);

			std::function<Log_Num(double)> f = [this, link](double y) {
				return Log_Num(log_n(link, y / link->s) - msg_n(link).log_norm_trunc(), PositiveSign); };
			Peak p = Peak::find_log_peak(f, A, B, std::log(1e-5));
			n_y.push_back({f, p, A, B});
		}

		Log_Num conv = ftconv.conv(n_y, 0.0, 1e-7, 1e-10);
		return conv.log_abs() + log_norm_prod(eqn) - log_sto_prod(eqn) + std::log(pi(eqn));
	}

	double logK_logconv(Eqn const* eqn) const
	{
		assert(eqn->links().size() > 0);

		std::vector<Bounded_log_func> n_y;

		for (Link const* link : eqn->links())
		{
			double A = msg_n(link).left() / link->s;
			double B = msg_n(link).right() / link->s;

			if (link->s < 0)
				std::swap(A, B);
			assert(A < B);
			std::function<double(double)> log_f = [this, link](double y) { return log_n(link, y / link->s); };
			n_y.push_back(Bounded_log_func(log_f, A, B));
		}

		log_conv conv(n_y.begin(), n_y.end());
		return conv.eval(0.0, 1e-5) + std::log(pi(eqn)) - log_sto_prod(eqn);
	}

	double msg_n_A_y(Link const* link) const
	{ return link->s > 0 ? msg_n(link).left() / link->s : msg_n(link).right() / link->s; }

	double msg_n_B_y(Link const* link) const
	{ return link->s > 0 ? msg_n(link).right() / link->s : msg_n(link).left() / link->s; }

	double log_sto_prod(Eqn const* eqn) const
	{
		double result = 0.0;
		for (Link const* link : eqn->links())
			result += log(std::abs(link->s));
		return result;
	}

	double log_norm_prod(Eqn const* eqn) const
	{
		double result = 0.0;
		for (Link const* link : eqn->links())
			result += msg_n(link).log_norm_trunc();
		return result;
	}

	// surface element constant
	double pi(Eqn const* eqn) const
	{
		double sum = 0;
		for (Link* e : eqn->links())
			sum += e->s * e->s;
		return std::sqrt(sum);
	}

	double log_m(Link const* link, double x) const
	{
		MsgM const& m = msg_m(link);
		return m.A() < x && x < m.B() ? (m.alpha() - 1) * log(x - m.A()) + (m.beta() - 1) * log(m.B() - x) : -Coss::INF;
	}

	double log_n(Link const* link, double x) const
	{
		double result = 0;
		for (Link const* e : link->var->links())
			if (e != link)
				result += log_m(e, x);
		return result;
	}

	double D_log_m(Link const* link, double x) const
	{
		MsgM const& m = msg_m(link);
		return x <= m.A() ? Coss::INF : x >= m.B() ? -Coss::INF : (m.alpha() - 1) / (x - m.A()) - (m.beta() - 1) * (m.B() - x);
	}

	double D_log_n(Link const* link, double x) const
	{
		double result = 0.0;
		for (Link const* e : link->var->links())
			if (e != link)
				result += D_log_m(link, x);
		return result;
	}

	double log_b_unnorm(Var const* var, double x) const
	{
		double result = 0;
		for (Link const* e : var->links())
			result += log_m(e, x);
		return result;
	}

	std::complex<double> msg_n_FT(Link const* link, double w) const
	{ return beta_characteristic_function(msg_n(link).alpha_trunc(), msg_n(link).beta_trunc(), msg_n(link).left(), msg_n(link).right(), w); }

	static std::complex<double> beta_characteristic_function(double alpha, double beta, double A, double B, double w)
	{ return std::polar(1.0, A * w) * Coss::py::hyp1f1(alpha, alpha + beta, Coss::I * (B - A) * w); }
};

constexpr double BP_B::min_ab;
constexpr double BP_B::max_lambda;

}}	// namespace Coss::bpB