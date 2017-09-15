#include <stella/op.h>
#include <stella/metrics.h>
#include <stella/grid.h>

namespace stella
{

template<>
void apply_interior<2>(op<2> & sop, const grid<2> & sgrid)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	auto & met = sgrid.metrics().coeff;
	auto & coeff = sop.get_coefficients();
	int beg[2], end[2];
	for (int axis = 0; axis < 2; axis++) {
		beg[axis] = sgrid.ibase(axis);
		end[axis] = sgrid.ibound(axis);
	}

	for (len_t j = beg[1]; j < end[1]; j++) {
		for (len_t i = beg[0]; i < end[0] + 1; i++) {
			coeff(i,j,nine_pt::w) = -1 * met(i,j,met2::w);
		}
	}

	for (len_t j = beg[1]; j < end[1] + 1; j++) {
		for (len_t i = beg[0]; i < end[0]; i++) {
			coeff(i,j,nine_pt::s) = -1 * met(i,j,met2::s);
		}
	}

	for (len_t j = beg[1]; j < end[1] + 1; j++) {
		for (len_t i = beg[0]; i < end[0] + 1; i++) {
			coeff(i,j,nine_pt::sw) = -1 * met(i,j,met2::sw);
		}
	}

	for (len_t j = sgrid.ibase(1); j < sgrid.ibound(1); j++) {
		for (len_t i = sgrid.ibase(0); i < sgrid.ibound(0); i++) {
			coeff(i,j,nine_pt::c) = -1 * (coeff(i  ,j  ,nine_pt::w ) +
			                              coeff(i+1,j  ,nine_pt::w ) +
			                              coeff(i  ,j  ,nine_pt::s ) +
			                              coeff(i  ,j+1,nine_pt::s ) +
			                              coeff(i  ,j+1,nine_pt::nw) +
			                              coeff(i+1,j+1,nine_pt::sw) +
			                              coeff(i  ,j  ,nine_pt::sw) +
			                              coeff(i+1,j  ,nine_pt::nw));
		}
	}
}

}
