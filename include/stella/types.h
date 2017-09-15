#ifndef STELLA_TYPES_H
#define STELLA_TYPES_H

#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/mpi/stencil_op.h>

namespace stella {

template<unsigned int nd>
	using grid_func = typename std::conditional<nd == 2,
		cedar::cdr2::mpi::grid_func,
		cedar::cdr3::mpi::grid_func>::type;

}
#endif
