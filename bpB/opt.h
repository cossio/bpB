#pragma once
#include <functional>
#include <vector>
#include <iostream>
#include <Coss/Vec.h>
#include <Coss/func.h>
#include <Coss/num/num.h>
#include <Coss/golden_search.h>
#include <Coss/bissect.h>

namespace Coss
{

struct Opt
{
	// Finds a point x = P + lambda * dir that minimizes F, for a <= lambad <= b.
	// Returns F evaluated at this point. Replaces P with P + lambda * dir.
	static double linmin(VecD& P, VecD const& dir, Vec_func const& F, double a, double b)
	{
		assert(a <= b);
		assert(Coss::is_finite(a) && Coss::is_finite(b));
		size_t N = P.size();
		assert(P.size() == N && dir.size() == N);

		VecD P0 = P;
		Re_func linF = [N, &F, &P0, &P, &dir] (double lambda) {
			P = P0 + lambda * dir;
			return -F(P);
		};

		double lambda = a < b ? Coss::golden_search(linF, a, b, 1e-7, 0) : a;
		assert(a <= lambda && lambda <= b);

		return -linF(lambda);
	}

	/*
	Minimization of a function func of n variables. Input consists of an initial starting point
	p[1..n] ; an initial matrix xi[1..n][1..n] , whose columns contain the initial set of di-
	rections (usually the n unit vectors); and ftol , the fractional tolerance in the function value
	such that failure to decrease by more than this amount on one iteration signals doneness. On
	output, p is set to the best point found, xi is the then-current direction set, fret is the returned
	function value at p, and iter is the number of iterations taken. The routine linmin is used.
	 */
	/* The above paragraph is from Numerical Recipes, but we did some modifications.
	 * A and B define the rectangular domain of the function, A <= x <= B. The iterations
	 * stay within these bounds.
	 * In the NR implementation, the columns of xi are the directions. In our implementation
	 * the rows are the directions. This allows us to use any number of directions, instead of
	 * being restricted to n directions. This is useful when the objective function is restrained
	 * by linear constrains in the coordinates. In this case, you want to choose only the directions
	 * that satisfy these linear constrains.
	 */
	static void powell(VecD& p, Vec<VecD>& xi, size_t& iter, double& fret,
					   VecD const& A, VecD const& B, Vec_func const& func)
	{
		const double FTOL = 1e-13;
		const size_t ITMAX = 200;

		size_t n = p.size();

		/* in linmin:
		 * p: point
		 * xi: direction
		 * fret: value of function
		 */
//    void linmin(double p[], double xi[], int n, double *fret,
//        double (*func)(double []));

		size_t i, ibig;
		double del, fp, fptt, lambda_min, lambda_max;

		VecD pt = p, ptt(n, 1.0), xit(n, 1.0);

		fret = func(p);

		// save the initial point
		pt = p;

		for (iter = 0; ; ++iter) // this is the loop that goes over iterations
		{
			fp = fret;
			ibig = 0;   // will be index of direction with biggest function decrease
			del = 0.0;  // will be the biggest function decrease

			for (i = 0; i < xi.size(); i++) // loop over all directions
			{
				xit = xi[i];    // copy the direction
				assert(xit.norm2() > 0);

				fptt = fret;

				// line minimization along current direction
				// minimum and maximum lambdas such that A <= P + lambda*xit <= B
				lambda_minmax(A, B, p, xit, lambda_min, lambda_max);
				fret = linmin(p, xit, func, lambda_min, lambda_max);

				// to decide if we update directions. Record direction with biggest change
				if (fabs(fptt - fret) > del)
				{
					del = fabs(fptt - fret);
					ibig = i;
				}
			}

			// stopping criterion, check function change
			// Note that if the matrix of directions is empty, we will
			// exit through this stopping criterion in the first iteration (because fp == fret).
			if (2.0 * fabs(fp - fret) <= FTOL * (fabs(fp) + fabs(fret)))
				return;

			// error if maximum iterations reached
			if (iter == ITMAX)
				nrerror("powell exceeding maximum iterations.");

			// calculate average direction moved
			xit = p - pt;

			// save starting point
			pt = p;

			// extrapolate along average direction
			lambda_minmax(A, B, p, xit, lambda_min, lambda_max);
			ptt = p + std::min(1.0, 0.9 * lambda_max) * xit;
			fptt = func(ptt);

			if (fptt < fp)
			{
				if (2.0 * (fp - 2.0 * fret + fptt) * Coss::square(fp - fret - del) < del * Coss::square(fp - fptt))
				{
					// line minimization
					fret = linmin(p, xit, func, lambda_min, lambda_max);

					// update directions
					xi[ibig] = xi.back();  // replace direction of biggest change (it is spent out!) with last direction
					xi.back() = xit;    // place average change direction at the end of direction vector
				}
			}
		}
	}

	// Compute the min and max values of lambda such that A <= P + lambda * n <= B
	static void lambda_minmax(VecD const& A, VecD const& B, VecD const& P, VecD const& n,
							  double& lambda_min, double& lambda_max)
	{
		size_t N = P.size();
		assert(A.size() == N && B.size() == N && n.size() == N);
		assert(N > 0);

		assert(n.norm2() > 0);

		// initialize lambda_min, lambda_max
		lambda_min = -Coss::INF;
		lambda_max = Coss::INF;

		for (size_t i = 0; i < N; ++i)
		{
			assert(A[i] < B[i]);
			assert(A[i] <= P[i] && P[i] <= B[i]);
			assert(is_finite(n[i]));

			if (n[i] > 0)
			{
				lambda_min = std::max(lambda_min, (A[i] - P[i]) / n[i]);
				lambda_max = std::min(lambda_max, (B[i] - P[i]) / n[i]);
				assert(lambda_min <= lambda_max);
			}
			else if (n[i] < 0)
			{
				lambda_min = std::max(lambda_min, (B[i] - P[i]) / n[i]);
				lambda_max = std::min(lambda_max, (A[i] - P[i]) / n[i]);
				assert(lambda_min <= lambda_max);
			}
		}

		assert(Coss::is_finite(lambda_min) && Coss::is_finite(lambda_max));
	}

	// Numerical Recipes standard error handler
	static void nrerror(std::string const& error_text)
	{
		std::cerr << "Numerical Recipes run-time error..." << std::endl;
		std::cerr << error_text << std::endl;
		std::cerr << "...now exiting to system..." << std::endl;
		exit(1);
	}

	// minimizes F(x). x contains initial point! Returns F evaluated at this point
	static double multimin(Vec_func const& F, VecD const& A, VecD const& B, VecD& x)
	{
		assert(A.size() == x.size() && B.size() == x.size());
		size_t N = x.size();
		assert(N > 0);

		if (N == 1) // problem is 1-dimensional
		{
			Re_func f = [&F] (double x) { return -F({x}); };
			x[0] = Coss::golden_search(f, A[0], B[0], 1e-7, 1e-7);
			return F(x);
		}

		assert(N > 1);
		// problem is multi-dimensional

		// Initialize matrix of directions with Cartesian unit vectors
		Vec<VecD> dirs(N, VecD(N, 0.0));
		for (size_t i = 0; i < N; ++i)
		{
			dirs[i][i] = 1.0;
			assert(A[i] < B[i]);
		}

		// Rotate last two directions, so that the direction associated
		// with the Lagrange multiplier isn't a Cartesian unit vector.
		// (This would bring problems because one of the bounds of the
		// Lagrange multiplier is infinity, and lambda_minmax would
		// not be finite if a Cartesian unit vector were used.
		double const diag = sqrt(2.0) / 2.0;
		dirs[N - 2][N - 2] = dirs[N - 2][N - 1] = dirs[N - 1][N - 1] = diag;
		dirs[N - 1][N - 2] = -diag;

		size_t iter;
		double fret;

		powell(x, dirs, iter, fret, A, B, F);
		return fret;
	}

	static double multimax(Vec_func const& F, VecD const& A, VecD const& B, VecD& x)
	{
		assert(A.size() == x.size() && B.size() == x.size());
		Vec_func mF = [&F] (VecD const v) { return -F(v); };
		return -multimin(mF, A, B, x);
	}

	/*
	 * Finds min/max xj such that:
	 * \sum_i x_i = 0
	 * F(x) >= c,
	 * A <= x <= B,
	 * If k < j, x_k, x_k+1, ..., x_j-1 are fixed.
	 * If j < k, x_k, x_k+1, ..., x_N-1, x_0, x_1, ..., x_j-1 are fixed.
	 * The fixed values are found in X, which is also the initial point. We make
	 * sure that the fixed values stay fixed by appropriate choice of directions in
	 * Powell's method.
	 * We assume that the initial point X satisfies all the constrains.
	 */
	static void xj_minmax(Vec_func const& F, VecD const& A, VecD const& B,
						  size_t k, size_t j, double c, VecD& X, double& xj_min, double& xj_max)
	{
		size_t N = X.size();
		assert(A.size() == N && B.size() == N);
		assert(N > 1);
		assert(j < N && k < N);

		size_t fixed = (N - k + j) % N;  // number of fixed coordinates (x_k, ..., x_j-1)
		size_t free  = (N - j + k - 1) % N;  // number of free coordinates (x_j+1, ..., x_k-1)
		assert(fixed + free + 1 == N);

		if (free == 0)  // in this case, xj is determined by its constrains, and the range is singular
		{
			xj_min = xj_max = X[j];
			return;
		}

		double sum_fixed = 0;
		for (size_t i = k; i != j; i = (i + 1) % N)
			sum_fixed += X[i];

		// The free coordinates are all the coordinates that are not fixed except
		// the j'th coordinate.

		double sum_free = 0, freeA = 0, freeB = 0;
		for (size_t i = (j + 1) % N; i != k; i = (i + 1) % N)
		{
			sum_free += X[i];
			freeA += A[i];
			freeB += B[i];
			assert(freeA < freeB);
		}

		assert(fabs(sum_fixed + sum_free + X[j]) < 1e-7);

		/* Powell's initial directions are chosen such that the constrains
		 * \sum_i x_i = 0 and the fixed values are enforced throughout
		 * Powell's iterations.
		 *
		 * Powell's direction only move the free coordinates. These are restrained to sum:
		 *
		 *    \sum_i\in free x_i = - sum_fixed - X[j]
		 *
		 */

		double constexpr SQR2 = 1.0 / sqrt(2.0);
		auto Powell_dirs = Vec<VecD>(free - 1, VecD(N, 0));
		for (size_t i = 0; i < free - 1; ++i)
		{
			// to enforce \sum_i x_i = 0, we use directions of the type:
			// n_i = (e_i - e_(i+1))/sqrt(2)
			Powell_dirs[i][(j + i + 1) % N] =  SQR2;
			Powell_dirs[i][(j + i + 2) % N] = -SQR2;
			assert((j + i + 1) % N != k);
			assert((j + i + 2) % N != k);
			/* note that these directions are independent,
			 * but they are not orthogonal
			 * hopefully that is not too bad
			 */
		}

		for (VecD const& dir : Powell_dirs)
			assert(dir.norm2() > 0);

		double xjA = std::max(A[j], -sum_fixed - freeB);
		double xjB = std::min(B[j], -sum_fixed - freeA);
		assert(xjA < xjB);

		/* returns the maximum of F(X) for all X satisfying the constrains
		 * and additionally with X[j] = xj. If this maximum is >= c, returns
		 * c + 1 instead. This introduces a discontinuity that is easy to detect
		 * by bissection below.
		 * X[] contains an initial point and it is assumed that
		 * A <= X <= B and \sum_i X[i] = 0
		 *
		 * Assumes that X satisfies all the constrains, that xj is
		 * within its range of allowed values.
		 */
		Re_func fj = [j, k, &F, c, &X, N, &A, &B, &Powell_dirs, free, xjA, xjB, sum_fixed] (double xj)
		{
			assert(free > 0);
			assert(xjA <= xj && xj <= xjB);

			double delta = xj - X[j];
			X[j] = xj;

			/*
			 * Subtract a piece of "delta" from the free coordinates,
			 * to maintain the constrain
			 * \sum X[i] == 0
			 */
			for (size_t i = (j + 1) % N; i != k && fabs(delta) > 1e-7; i = (i + 1) % N)
			{
				assert(A[i] <= X[i] && X[i] <= B[i]);
				double Xi_new = std::max(A[i], std::min(X[i] - delta, B[i]));
				delta += Xi_new - X[i];
				X[i] = Xi_new;
				assert(A[i] <= X[i] && X[i] <= B[i]);
			}
			assert(fabs(delta) < 1e-7);

			// Copy Powell directions (because the call to powell() will change them,
			// and we want to preserve the original directions for later calls to this function.
			auto Powell_dirs_cp = Powell_dirs;
			size_t iter;
			double fret;

			// before Powell runs, sum(X) is zero
			assert(fabs(X.sum()) < 1e-7);

			powell(X, Powell_dirs_cp, iter, fret, A, B, F);

			// sum(X) is still zero after Powell has run
			assert(fabs(X.sum()) < 1e-7);

			return fret >= c ? c + 1 : fret;
		};


		double xj_c = xjA < xjB ? Coss::golden_search(fj, xjA, xjB, 1e-5, 1e-5) : xjA;
		assert(fj(xj_c) >= c + 1);
		assert(fj(xjA) < c);
		assert(fj(xjB) < c);

		xj_min = fj(xjA) < c ? Coss::bissect(fj, xjA, xj_c, c + 0.5) : xjA;
		xj_max = fj(xjB) < c ? Coss::bissect(fj, xj_c, xjB, c + 0.5) : xjB;
	}
};

}   // namespace Coss
