#include <functional>

#include <stella/op.h>
#include <stella/metrics.h>
#include <stella/grid.h>
#include <stella/bc/dirichlet.h>

namespace stella
{

template<>
void op<2>::contrib_rhs_metrics(grid_func<2> & rhs)
{
	auto & met = sgrid->metrics();
	auto & classify = this->get_classify();

	for (auto j : rhs.range(1)) {
		for (auto i : rhs.range(0)) {
			if (classify(i,j) == ptype::interior)
				rhs(i,j) *= met.J(i,j);
		}
	}
}


template<>
void op<2>::register_bcs()
{
	bc_op_callbacks.push_back([](op<2> & sop, const grid<2> & sgrid) {
			apply_dirichlet<2>(sop, sgrid);
		});
	bc_rhs_callbacks.push_back([](grid_func<2> & rhs, const op<2> & sop, const grid<2> & sgrid) {
			apply_dirichlet<2>(rhs, sop, sgrid);
		});
}


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
			coeff(i,j,nine_pt::w) = met(i,j,met2::w);
		}
	}

	for (len_t j = beg[1]; j < end[1] + 1; j++) {
		for (len_t i = beg[0]; i < end[0]; i++) {
			coeff(i,j,nine_pt::s) = met(i,j,met2::s);
		}
	}

	for (len_t j = beg[1]; j < end[1] + 1; j++) {
		for (len_t i = beg[0]; i < end[0] + 1; i++) {
			coeff(i,j,nine_pt::sw) = met(i,j,met2::sw);
		}
	}

	for (len_t j = beg[1]; j < end[1] + 1; j++) {
		for (len_t i = beg[0] + 1; i < end[0]; i++) {
			coeff(i+1,j,nine_pt::nw) = met(i,j,met2::se);
		}
	}

	for (len_t j = sgrid.ibase(1); j < sgrid.ibound(1); j++) {
		for (len_t i = sgrid.ibase(0); i < sgrid.ibound(0); i++) {
			coeff(i,j,nine_pt::c) = (coeff(i  ,j  ,nine_pt::w) +
			                         coeff(i+1,j  ,nine_pt::w) +
			                         coeff(i  ,j  ,nine_pt::s) +
			                         coeff(i  ,j+1,nine_pt::s));
		}
	}
}


template<>
void apply_interior<3>(op<3> & sop, const grid<3> & sgrid)
{
	using namespace cedar;
	using namespace cedar::cdr3;

	auto & met = sgrid.metrics().coeff;
	auto & coeff = sop.get_coefficients();
	int beg[3], end[3];
	for (int axis = 0; axis < 3; axis++) {
		beg[axis] = sgrid.ibase(axis);
		end[axis] = sgrid.ibound(axis);
	}

	for (len_t k = beg[2]; k < end[2]; k++) {
		for (len_t j = beg[1]; j < end[1]; j++) {
			for (len_t i = beg[0]; i < end[0] + 1; i++) {
				coeff(i,j,k,xxvii_pt::pw) = met(i,j,k,met3::w);
			}
		}
	}

	for (len_t k = beg[2]; k < end[2]; k++) {
		for (len_t j = beg[1]; j < end[1] + 1; j++) {
			for (len_t i = beg[0]; i < end[0]; i++) {
				coeff(i,j,k,xxvii_pt::ps) = met(i,j,k,met3::s);
			}
		}
	}

	for (len_t k = beg[2]; k < end[2] + 1; k++) {
		for (len_t j = beg[1]; j < end[1]; j++) {
			for (len_t i = beg[0]; i < end[0]; i++) {
				coeff(i,j,k,xxvii_pt::b) = met(i,j,k,met3::b);
			}
		}
	}
}

}
