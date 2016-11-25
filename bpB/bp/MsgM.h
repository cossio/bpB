#pragma once
#include <cassert>
#include <Coss/num/num.h>
#include <Coss/Beta.h>

namespace Coss { namespace bpB
{

class MsgM
{
public:

	MsgM(double lb, double ub, double A, double B, double alpha, double beta)
			: _lb(lb), _ub(ub), _A(A), _B(B), _alpha(alpha), _beta(beta)
	{
		assert(Coss::is_finite(alpha) && Coss::is_finite(beta) && alpha > 0 && beta > 0);
		assert(alpha > 0 && beta > 0);
		assert(alpha >= 1 || beta >= 1);
		assert(Coss::is_finite(A) && Coss::is_finite(B));
		assert(A < B);
		assert(Coss::is_finite(lb) && Coss::is_finite(ub));
		assert(lb < ub);
		assert(left() < right());
		assert(left_01_scaled() < right_01_scaled());
		assert(0 <= left_01_scaled() && right_01_scaled() <= 1);
	}

	MsgM(double lb, double ub) : MsgM(lb, ub, lb, ub, 2, 2)
	{ }

	bool updates_next = true;

	double A() const
	{ return _A; }

	double B() const
	{ return _B; }

	double alpha() const
	{ return _alpha; }

	double beta() const
	{ return _beta; }

	double lb() const
	{ return _lb; }

	double ub() const
	{ return _ub; }

	double L() const
	{ return B() - A(); }

	double lambda() const
	{ return alpha() + beta(); }

	double mean() const
	{ return (alpha() * B() + beta() * A()) / lambda(); }

	double var() const
	{ return L() * L() * alpha() * beta() / (lambda() * lambda() * (1 + lambda())); }

	double mode_01_scaled() const
	{ return (alpha() - 1) / (lambda() - 2); }

	double mode() const
	{ return (B() - A()) * mode_01_scaled() + A(); }

	double mode_trunc() const
	{ return lb() <= mode() ? std::min(ub(), mode()) : lb(); }

	double mode_trunc_01_scaled() const
	{ return left_01_scaled() <= mode_01_scaled() ? std::min(right_01_scaled(), mode_01_scaled()) : left_01_scaled(); }

	double left() const
	{ return std::max(A(), lb()); }

	double right() const
	{ return std::min(B(), ub()); }

	double left_01_scaled() const
	{ return (left() - A()) / L(); }

	double right_01_scaled() const
	{ return (right() - A()) / L(); }

	double w1() const
	{
		return exp(alpha() * std::log(left_01_scaled()) + beta() * std::log(1 - left_01_scaled()) -
				   Beta::log_beta(alpha(), beta(), left_01_scaled(), right_01_scaled()));
	}

	double w2() const
	{
		return exp(alpha() * std::log(right_01_scaled()) + beta() * std::log(1 - right_01_scaled()) -
				   Beta::log_beta(alpha(), beta(), left_01_scaled(), right_01_scaled()));
	}

	double first_moment() const
	{ return L() * (alpha() - w2() + w1()) / lambda(); };

	double second_moment() const
	{
		return L() * L() / (lambda() + 1) *
			   ((alpha() + 1) / lambda() * (alpha() - w2() + w1()) - right_01_scaled() * w2() + left_01_scaled() * w1());
	};

	double log_norm() const
	{ return (lambda() - 1) * std::log(L()) + Beta::log_beta(alpha(), beta()); }

	double log_norm_trunc() const
	{ return (lambda() - 1) * std::log(L()) + Beta::log_beta(alpha(), beta(), left_01_scaled(), right_01_scaled()); }

	double mean_trunc() const
	{ return first_moment() + A(); }

	double var_trunc() const
	{ return second_moment() - first_moment() * first_moment(); }

	double lambda_trunc() const
	{ return (mean_trunc() - left()) * (right() - mean_trunc()) / var_trunc() - 1; }

	void regularize()
	{
		// forces A >= lb, B <= ub

		double mean_tr = mean_trunc();
		double var_tr = var_trunc();
		double lambda_tr = (mean_tr - left()) * (right() - mean_tr) / var_tr - 1;
		lambda_tr = std::max(1.0, lambda_tr);

		_A = left();
		_B = right();

		_alpha = lambda_tr * (mean_tr - A()) / L();
		_beta = lambda_tr * (B() - mean_tr) / L();

		assert(_alpha > 0 && _beta > 0);
		assert(_alpha >= 1 || _beta >= 1);

//		_alpha = std::max(_alpha, 1.0);
//		_beta = std::max(_beta, 1.0);
	}

private:
	double _lb, _ub, _A, _B, _alpha, _beta;

};

}}	// namespace Coss::bpB