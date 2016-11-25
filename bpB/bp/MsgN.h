#pragma once
#include <cassert>
#include <complex>
#include <Coss/num/num.h>
#include <Coss/Peak.h>
#include <Coss/sf/hyp1F1.h>

namespace Coss { namespace bpB {

class MsgN
{
private:
	double _left, _right, _mean_trunc, _var_trunc, _log_norm_trunc, _alpha_trunc, _beta_trunc;
	Coss::Peak _log_peak;

public:

	MsgN(double left, double right, double mean_trunc, double var_trunc, double log_norm_trunc, Coss::Peak log_peak)
			: _left(left), _right(right),
			  _mean_trunc(mean_trunc), _var_trunc(var_trunc),
			  _log_norm_trunc(log_norm_trunc), _log_peak(log_peak),
			  _alpha_trunc(lambda_trunc() * (mean_trunc - left) / L()),
			  _beta_trunc(lambda_trunc() * (right - mean_trunc) / L())
	{
		assert(left < mean_trunc && mean_trunc < right);
		assert(var_trunc > 0);
		assert(var_trunc < (right - mean_trunc) * (mean_trunc - left));
		assert(Coss::is_finite(log_norm_trunc));
	}

	MsgN(double left, double right)
			: MsgN(left, right,
				   (left + right) / 2.0, (right - left) * (right - left) / 12.0, std::log(right - left),
				   Coss::Peak((left + right) / 2.0, left, right, 1.0 / (right - left), 1.0 / (right - left)))
	{ }

	bool needs_update = true;

	double left() const
	{ return _left; }

	double right() const
	{ return _right; }

	double L() const
	{ return right() - left(); }

	double log_norm_trunc() const
	{ return _log_norm_trunc; }

	double mean_trunc() const
	{ return _mean_trunc; }

	double var_trunc() const
	{ return _var_trunc; }

	Coss::Peak log_peak() const
	{ return _log_peak; }

	double lambda_trunc() const
	{ return (mean_trunc() - left()) * (right() - mean_trunc()) / var_trunc() - 1; }

	double alpha_trunc() const
	{ return _alpha_trunc; }

	double beta_trunc() const
	{ return _beta_trunc; }
};

}}	// namespace Coss::bpB