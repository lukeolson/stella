#ifndef STELLA_METRICS_H
#define STELLA_METRICS_H

#include <cedar/halo_exchanger.h>
#include <cedar/mpi/grid_topo.h>
#include <stella/types.h>

namespace stella
{

	template<unsigned int nd, class halo_exchanger=cedar::halo_exchanger<2>>
	class grid;

enum class met2 {c=0, w, s, sw, se, ndirs};
enum class met3 {c=0, w, s, sw, se, b, wb, eb, sb, nb, ndirs};

template<unsigned int nd>
using metric_stencil = typename std::conditional<nd==2, met2, met3>::type;

template<unsigned int nd>
struct metrics
{
	template<class sten>
	using stencil_op = cedar::cdr2::mpi::stencil_op<sten>;

metrics(cedar::topo_ptr topo): coeff(topo), J(topo), g12(topo) {}

	stencil_op<metric_stencil<nd>> coeff;
	grid_func<2> J, g12;
};

template<unsigned int nd>
void compute_interior(metrics<nd> & met, const grid<nd> & sgrid);

template<unsigned int nd>
void compute_boundary(metrics<nd> & met, const grid<nd> & sgrid);

}

#endif
