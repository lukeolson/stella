#include <stella/bc/dirichlet.h>

namespace stella
{

template<>
void apply_dirichlet<2>(op<2> & sop, const grid<2> & sgrid)
{
	using namespace cedar;
	using namespace cedar::cdr2;
	len_t beg[2], end[2];
	for (int axis = 0; axis < 2; axis++) {
		beg[axis] = sgrid.ibase(axis);
		end[axis] = sgrid.ibound(axis);
	}

	auto & coeff = sop.get_coefficients();
	auto & classify = sop.get_classify();

	for (len_t j = beg[1]; j < end[1]; j++) {
		for (len_t i = beg[0]; i < end[0]; i++) {
			if (classify(i,j) == ptype::dirichlet) {
				coeff(i,j,nine_pt::c) = 1.0;
				// Pass 1. delete connections to neighbors
				coeff(i  ,j,  nine_pt::w ) = 0.0;
				coeff(i  ,j  ,nine_pt::s ) = 0.0;
				coeff(i+1,j  ,nine_pt::w ) = 0.0;
				coeff(i  ,j+1,nine_pt::s ) = 0.0;
				coeff(i  ,j  ,nine_pt::sw) = 0.0;
				coeff(i  ,j+1,nine_pt::nw) = 0.0;
				coeff(i+1,j+1,nine_pt::sw) = 0.0;
				coeff(i+1,j  ,nine_pt::nw) = 0.0;
			}
			// pass 2. delete connections to dirichlet points
			if (classify(i-1,j) == ptype::dirichlet)
				coeff(i,j,nine_pt::w) = 0.0;
			if (classify(i+1,j) == ptype::dirichlet)
				coeff(i+1,j,nine_pt::w) = 0.0;
			if (classify(i,j-1) == ptype::dirichlet)
				coeff(i,j,nine_pt::s) = 0.0;
			if (classify(i,j+1) == ptype::dirichlet)
				coeff(i,j+1,nine_pt::s) = 0.0;
			if (classify(i-1,j-1) == ptype::dirichlet)
				coeff(i,j,nine_pt::sw) = 0.0;
			if (classify(i+1,j+1) == ptype::dirichlet)
				coeff(i+1,j+1,nine_pt::sw) = 0.0;
			if (classify(i+1,j-1) == ptype::dirichlet)
				coeff(i+1,j,nine_pt::nw) == 0.0;
			if (classify(i-1,j+1) == ptype::dirichlet)
				coeff(i,j+1,nine_pt::nw) = 0.0;
		}
	}
}


template<>
void apply_dirichlet<2>(grid_func<2> & rhs, const op<2> & sop, const grid<2> & sgrid)
{
	auto & met = sgrid.metrics().coeff;
	auto & classify = sop.get_classify();
	auto & dirichlet = sop.get_bcval();

	using namespace cedar;
	using namespace cedar::cdr2;
	len_t beg[2], end[2];
	for (int axis = 0; axis < 2; axis++) {
		beg[axis] = sgrid.ibase(axis);
		end[axis] = sgrid.ibound(axis);
	}

	for (len_t j = beg[1]; j < end[1]; j++) {
		for (len_t i = beg[0]; i < end[0]; i++) {
			if (classify(i, j) != ptype::dirichlet) {
				if (classify(i-1, j) == ptype::dirichlet)
					rhs(i, j) += met(i,j,met2::w) * dirichlet(i-1, j);
				if (classify(i+1, j) == ptype::dirichlet)
					rhs(i, j) += met(i+1,j,met2::w) * dirichlet(i+1, j);
				if (classify(i, j-1) == ptype::dirichlet)
					rhs(i, j) += met(i,j,met2::s) * dirichlet(i, j-1);
				if (classify(i, j+1) == ptype::dirichlet)
					rhs(i, j) += met(i,j+1,met2::s) * dirichlet(i, j+1);
				if (classify(i+1, j+1) == ptype::dirichlet)
					rhs(i, j) += met(i+1,j+1, met2::sw) * dirichlet(i+1, j+1);
				if (classify(i-1, j-1) == ptype::dirichlet)
					rhs(i, j) += met(i,j,met2::sw) * dirichlet(i-1, j-1);
				if (classify(i+1, j-1) == ptype::dirichlet)
					rhs(i, j) += met(i,j,met2::se) * dirichlet(i+1, j-1);
				if (classify(i-1, j+1) == ptype::dirichlet)
					rhs(i, j) += met(i-1,j+1,met2::se) * dirichlet(i-1, j+1);
			} else {
				rhs(i,j) = dirichlet(i,j);
			}
		}
	}
}

}
