#include <stella/metrics.h>
#include <stella/grid.h>
#include <stella/halo.h>

namespace stella
{

	template<>
	void compute_interior<2>(metrics<2> & met, const grid<2> & sgrid)
	{
		using namespace cedar;
		array<len_t, 1> beg(2);
		array<len_t, 1> end(2);

		for (int axis = 0; axis < 2; axis++) {
			beg(axis) = sgrid.ibase(axis);
			end(axis) = sgrid.ibound(axis);
		}

		auto & x = sgrid.vertex(0);
		auto & y = sgrid.vertex(1);

		auto & coeff = met.coeff;

		real_t x_s, x_t, y_s, y_t, J;
		real_t g_22, g_11, g_12, g11, g22, g12;
		real_t x_xh[3], x_yh[3], y_xh[3], y_yh[3];
		double feps = 1e-13;

		for (len_t j = beg(1); j < end(1) + 1; j++) {
			for (len_t i = beg(0); i < end(0); i++) {
				// compute metrics at (i,j-1/2)
				for (int k = -1; k < 2; k = k + 2) {
					x_yh[k] = 0.5 * (x(i+k,j-1) + x(i+k,j));
					y_yh[k] = 0.5 * (y(i+k,j-1) + y(i+k,j));
				}

				x_s = 0.5 * (x_yh[1] - x_yh[-1]);
				y_s = 0.5 * (y_yh[1] - y_yh[-1]);
				x_t = x(i,j) - x(i,j-1);
				y_t = y(i,j) - y(i,j-1);

				J = x_s*y_t - x_t*y_s;
				g_11 = x_s*x_s + y_s*y_s;

				if (std::abs(J) < feps) g22 = 0.0;
				else g22 = 1.0 / (J*J) * g_11;

				coeff(i,j,met2::s) = J * g22;
			}
		}

		for (len_t j = beg(1); j < end(1); j++) {
			for (len_t i = beg(0); i < end(0) + 1; i++) {
				// compute metrics at (i-1/2,j)
				for (int k = -1; k < 2; k = k + 2) {
					x_xh[k] = 0.5 * (x(i-1,j+k) + x(i,j+k));
					y_xh[k] = 0.5 * (y(i-1,j+k) + y(i,j+k));
				}

				x_s = x(i,j) - x(i-1,j);
				y_s = y(i,j) - y(i-1,j);
				x_t = 0.5 * (x_xh[1] - x_xh[-1]);
				y_t = 0.5 * (y_xh[1] - y_xh[-1]);

				J = x_s * y_t - x_t * y_s;
				g_22 = x_t * x_t + y_t * y_t;

				if (std::abs(J) < feps) g11 = 0.0;
				else g11 = 1.0 / (J*J) * g_22;

				coeff(i,j,met2::w) = J * g11;
			}
		}

		for (len_t j = beg(1); j < end(1); j++) {
			for (len_t i = beg(0); i < end(0); i++) {
				x_s = 0.5 * (x(i+1,j) - x(i-1,j));
				y_s = 0.5 * (y(i+1,j) - y(i-1,j));
				x_t = 0.5 * (x(i,j+1) - x(i,j-1));
				y_t = 0.5 * (y(i,j+1) - y(i,j-1));

				J = x_s*y_t - x_t*y_s;
				met.J(i,j) = J;

				g_12 = x_s * x_t + y_s * y_t;

				if (std::abs(J) < feps) g12 = 0.0;
				else g12 = -1.0 / (J*J) * g_12;
				met.g12(i,j) = g12;
			}
		}

		exchange_halo(sgrid, met.g12);

		for (len_t j = beg(1); j < end(1); j++) {
			for (len_t i = beg(0); i < end(0); i++) {
				coeff(i,j,met2::sw) = 0.25 * J * (met.g12(i-1,j) + met.g12(i,j-1));
				coeff(i,j,met2::se) = -0.25 * J * (met.g12(i,j-1) + met.g12(i+1,j));
			}
		}
	}


	template<>
	void compute_boundary<2>(metrics<2> & met, const grid<2> & sgrid)
	{

	}

}
