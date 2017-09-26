#ifndef STELLA_BC_DIRICHLET_H
#define STELLA_BC_DIRICHLET_H

#include <stella/op.h>
#include <stella/grid.h>

namespace stella
{

	template<unsigned int nd>
		void apply_dirichlet(op<nd> & sop, const grid<nd> & sgrid);

	template<unsigned int nd>
		void apply_dirichlet(grid_func<nd> & rhs, const op<nd> & sop, const grid<nd> & sgrid);

}

#endif
