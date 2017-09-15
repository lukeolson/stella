#ifndef STELLA_OP_H
#define STELLA_OP_H

#include <type_traits>

#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/3d/mpi/stencil_op.h>

#include <stella/metrics.h>

namespace stella
{
	using namespace cedar;
	template<unsigned int nd>
	class op
	{
		using stencil_op = typename std::conditional<nd==2,
			cdr2::mpi::stencil_op<cdr2::nine_pt>,
			cdr3::mpi::stencil_op<cdr3::xxvii_pt>>::type;
		using met_op = typename std::conditional<nd==2,
			cdr2::mpi::stencil_op<met2>,
			cdr3::mpi::stencil_op<met3>>::type;
	public:
	op(std::shared_ptr<grid<nd>> grd) : coeff(grd), mask(grd), sgrid(grd)
		{
			init_mask();
		}

		stencil_op & get_coefficients() { return coeff; }

		void init_mask()
		{
			mask.set(1);
		}

		void assemble() {
			apply_interior(*this, *sgrid);
		}

	protected:
		stencil_op coeff;
		grid_func<nd> mask;
		std::shared_ptr<grid<nd>> sgrid;
	};

	template<unsigned int nd>
		void apply_interior(op<nd> & sop, const grid<nd> & sgrid);
}

#endif
