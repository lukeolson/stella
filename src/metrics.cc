#include <stella/metrics.h>
#include <stella/grid.h>
#include <stella/halo.h>

#include <float.h>
#define FDIV(num, den) (fabs(den) < DBL_MIN ? 0.0 : num / den)

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
				J = met.J(i,j);
				coeff(i,j,met2::sw) = 0.25 * J * (met.g12(i-1,j) + met.g12(i,j-1));
				coeff(i,j,met2::se) = -0.25 * J * (met.g12(i,j-1) + met.g12(i+1,j));
			}
		}
	}


	template<>
	void compute_boundary<2>(metrics<2> & met, const grid<2> & sgrid)
	{

	}


	template<>
	void compute_interior<3>(metrics<3> & met, const grid<3> & sgrid)
	{
		using namespace cedar;
		len_t beg[3], end[3];

		for (int axis = 0; axis < 3; axis++) {
			beg[axis] = sgrid.ibase(axis);
			end[axis] = sgrid.ibound(axis);
		}

		auto & x = sgrid.vertex(0);
		auto & y = sgrid.vertex(1);
		auto & z = sgrid.vertex(2);

		auto & coors = sgrid.get_vertex();

		auto & coeff = met.coeff;

		real_t x_r, x_s, x_t, y_r, y_s, y_t, z_r, z_s, z_t;
		real_t g, J;
		real_t g_22, g_11, g_33, g_12, g_13, g_23;
		real_t g11, g22, g33, g12, g13, g23;
		real_t pt0[3], pt1[3], pt2[3];

		for (len_t k = beg[2]; k < end[2]; k++) {
			for (len_t j = beg[1]; j < end[1]; j++) {
				for (len_t i = beg[0]; i < end[0] + 1; i++) {
					x_r = x(i,j,k) - x(i-1,j,k);
					y_r = y(i,j,k) - y(i-1,j,k);
					z_r = z(i,j,k) - z(i-1,j,k);

					pt0[0] = 0.5 * (x(i-1,j-1,k) + x(i,j-1,k));
					pt1[0] = 0.5 * (x(i-1,j+1,k) + x(i,j+1,k));
					pt0[1] = 0.5 * (y(i-1,j-1,k) + y(i,j-1,k));
					pt1[1] = 0.5 * (y(i-1,j+1,k) + y(i,j+1,k));
					pt0[2] = 0.5 * (z(i-1,j-1,k) + z(i,j-1,k));
					pt1[2] = 0.5 * (z(i-1,j+1,k) + z(i,j+1,k));
					x_s = 0.5 * (pt1[0] - pt0[0]);
					y_s = 0.5 * (pt1[1] - pt0[1]);
					z_s = 0.5 * (pt1[2] - pt0[2]);

					pt0[0] = 0.5 * (x(i-1,j,k-1) + x(i,j,k-1));
					pt1[0] = 0.5 * (x(i-1,j,k+1) + x(i,j,k+1));
					pt0[1] = 0.5 * (y(i-1,j,k-1) + y(i,j,k-1));
					pt1[1] = 0.5 * (y(i-1,j,k+1) + y(i,j,k+1));
					pt0[2] = 0.5 * (z(i-1,j,k-1) + z(i,j,k-1));
					pt1[2] = 0.5 * (z(i-1,j,k+1) + z(i,j,k+1));
					x_t = 0.5 * (pt1[0] - pt0[0]);
					y_t = 0.5 * (pt1[1] - pt0[1]);
					z_t = 0.5 * (pt1[2] - pt0[2]);

					g_11 = x_r*x_r + y_r*y_r + z_r*z_r;
					g_12 = x_r*x_s + y_r*y_s + z_r*z_s;
					g_13 = x_r*x_t + y_r*y_t + z_r*z_t;
					g_22 = x_s*x_s + y_s*y_s + z_s*z_s;
					g_23 = x_s*x_t + y_s*y_t + z_s*z_t;
					g_33 = x_t*x_t + y_t*y_t + z_t*z_t;

					g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));
					g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);
					coeff(i,j,k,met3::w) = sqrt(g) * g11;
				}
			}
		}

		std::array<real_t*, 3> deriv;
		for (len_t k = beg[2]; k < end[2]; k++) {
			for (len_t j = beg[1]; j < end[1] + 1; j++) {
				for (len_t i = beg[0]; i < end[1]; i++) {
					deriv[0] = &x_r; deriv[1] = &y_r; deriv[2] = &z_r;
					for (int axis = 0; axis < 3; axis++) {
						pt0[axis] = 0.5 * (coors[axis](i-1,j-1,k) + coors[axis](i-1,j,k));
						pt1[axis] = 0.5 * (coors[axis](i+1,j-1,k) + coors[axis](i+1,j,k));
						*deriv[axis] = 0.5 * (pt1[axis] - pt0[axis]);
					}

					deriv[0] = &x_s; deriv[1] = &y_s; deriv[2] = &z_s;
					for (int axis = 0; axis < 3; axis++) {
						*deriv[axis] = coors[axis](i,j,k) - coors[axis](i,j-1,k);
					}

					deriv[0] = &x_t; deriv[1] = &y_t; deriv[2] = &z_t;
					for (int axis = 0; axis < 3; axis++) {
						pt0[axis] = 0.5 * (coors[axis](i,j-1,k-1) + coors[axis](i,j,k-1));
						pt1[axis] = 0.5 * (coors[axis](i,j-1,k+1) + coors[axis](i,j,k+1));
						*deriv[axis] = 0.5 * (pt1[axis] - pt0[axis]);
					}

					g_11 = x_r*x_r + y_r*y_r + z_r*z_r;
					g_12 = x_r*x_s + y_r*y_s + z_r*z_s;
					g_13 = x_r*x_t + y_r*y_t + z_r*z_t;
					g_22 = x_s*x_s + y_s*y_s + z_s*z_s;
					g_23 = x_s*x_t + y_s*y_t + z_s*z_t;
					g_33 = x_t*x_t + y_t*y_t + z_t*z_t;

					g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));
					g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);
					coeff(i,j,k,met3::s) = sqrt(g) * g22;
				}
			}
		}

		for (len_t k = beg[2]; k < end[2] + 1; k++) {
			for (len_t j = beg[1]; j < end[1]; j++) {
				for (len_t i = beg[0]; i < end[0]; i++) {
					deriv[0] = &x_r; deriv[1] = &y_r; deriv[2] = &z_r;
					for (int axis = 0; axis < 3; axis++) {
						pt0[axis] = 0.5 * (coors[axis](i-1,j,k-1) + coors[axis](i-1,j,k));
						pt1[axis] = 0.5 * (coors[axis](i+1,j,k-1) + coors[axis](i+1,j,k));
						*deriv[axis] = 0.5 * (pt1[axis] - pt0[axis]);
					}

					deriv[0] = &x_s; deriv[1] = &y_s; deriv[2] = &z_s;
					for (int axis = 0; axis < 3; axis++) {
						pt0[axis] = 0.5 * (coors[axis](i,j-1,k-1) + coors[axis](i,j-1,k));
						pt1[axis] = 0.5 * (coors[axis](i,j+1,k-1) + coors[axis](i,j+1,k));
						*deriv[axis] = 0.5 * (pt1[axis] - pt0[axis]);
					}

					deriv[0] = &x_t; deriv[1] = &y_t; deriv[2] = &z_t;
					for (int axis = 0; axis < 3; axis++) {
						*deriv[axis] = coors[axis](i,j,k) - coors[axis](i,j,k-1);
					}

					g_11 = x_r*x_r + y_r*y_r + z_r*z_r;
					g_12 = x_r*x_s + y_r*y_s + z_r*z_s;
					g_13 = x_r*x_t + y_r*y_t + z_r*z_t;
					g_22 = x_s*x_s + y_s*y_s + z_s*z_s;
					g_23 = x_s*x_t + y_s*y_t + z_s*z_t;
					g_33 = x_t*x_t + y_t*y_t + z_t*z_t;

					g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));
					g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);
					coeff(i,j,k,met3::b) = sqrt(g) * g33;
				}
			}
		}

		for (len_t k = beg[2]; k < end[2]; k++) {
			for (len_t j = beg[1]; j < end[1]; j++) {
				for (len_t i = beg[0]; i < end[0]; i++) {
					x_r = 0.5 * (x(i+1,j,k) - x(i-1,j,k));
					y_r = 0.5 * (y(i+1,j,k) - y(i-1,j,k));
					z_r = 0.5 * (z(i+1,j,k) - z(i-1,j,k));
					x_s = 0.5 * (x(i,j+1,k) - x(i,j-1,k));
					y_s = 0.5 * (y(i,j+1,k) - y(i,j-1,k));
					z_s = 0.5 * (z(i,j+1,k) - y(i,j-1,k));
					x_t = 0.5 * (x(i,j,k+1) - x(i,j,k-1));
					y_t = 0.5 * (y(i,j,k+1) - y(i,j,k-1));
					z_t = 0.5 * (z(i,j,k+1) - z(i,j,k-1));

					g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));
					met.J(i,j,k) = sqrt(g);

					met.g12(i,j,k) = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);
					met.g13(i,j,k) = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);
					met.g23(i,j,k) = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);
				}
			}
		}

		exchange_halo(sgrid, met.g12);
		exchange_halo(sgrid, met.g13);
		exchange_halo(sgrid, met.g23);

		for (len_t k = beg[2]; k < end[2]; k++) {
			for (len_t j = beg[1]; j < end[1]; j++) {
				for (len_t i = beg[0]; i < end[0]; i++) {
					J = met.J(i,j,k);
					coeff(i,j,k,met3::sw) = 0.25 * J * (met.g12(i-1,j,k) + met.g12(i,j-1,k));
					coeff(i,j,k,met3::se) = -0.25 * J * (met.g12(i+1,j,k) + met.g12(i,j-1,k));
					coeff(i,j,k,met3::wb) = 0.25 * J * (met.g13(i-1,j,k) + met.g13(i,j,k-1));
					coeff(i,j,k,met3::eb) = -0.25 * J * (met.g13(i+1,j,k) + met.g13(i,j,k-1));
					coeff(i,j,k,met3::sb) = 0.25 * J * (met.g23(i,j-1,k) + met.g23(i,j,k-1));
					coeff(i,j,k,met3::nb) = -0.25 * J * (met.g23(i,j+1,k) + met.g23(i,j,k-1));
				}
			}
		}
	}


	template<>
	void compute_boundary<3>(metrics<3> & met, const grid<3> & sgrid)
	{

	}

}
