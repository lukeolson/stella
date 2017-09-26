#ifndef STELLA_SOLVER_H
#define STELLA_SOLVER_H

#include <stella/op.h>
#include <cedar/2d/mpi/solver.h>
namespace stella
{
	template<unsigned int nd>
		void solve(op<nd> & sop, grid_func<nd> & x, const grid_func<nd> b)
	{
		using namespace cedar::cdr2;
		auto & cdr_op = sop.get_coefficients();

		mpi::solver<nine_pt> slv(cdr_op);
		slv.solve(b, x);

		// cedar grid d.s. have changed (need to update halo_exchanger)
		auto kreg = slv.kernel_registry();
		auto sgrid = sop.get_grid();
		sgrid->set_halo_exchanger(kreg->get_halo_exchanger());
	}
}

#endif
