#ifndef STELLA_OP_H
#define STELLA_OP_H

#include <type_traits>

#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/3d/mpi/stencil_op.h>

#include <stella/metrics.h>
#include <stella/halo.h>

namespace stella
{
	enum class ptype { interior, dirichlet };

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
	op(std::shared_ptr<grid<nd>> grd) : coeff(grd),
			mask(grd->nlocal(0), grd->nlocal(1)), sgrid(grd), bcval(grd)
		{
			init_mask();
			bcval.set(0.0);
			register_bcs();
		}

		stencil_op & get_coefficients() { return coeff; }
		const stencil_op & get_coefficients() const { return coeff; }

		void init_mask()
		{
			for (len_t j = 0; j < mask.len(1); j++) {
				for (len_t i = 0; i < mask.len(0); i++) {
					mask(i,j) = ptype::interior;
				}
			}
		}

		void assemble()
		{
			apply_interior(*this, *sgrid);
			exchange_halo(*sgrid, bcval);
			apply_boundary(*this, *sgrid);
		}

		void call_bcs(const grid<nd> & sgrid)
		{
			for (auto & bc_callback : bc_op_callbacks) {
				bc_callback(*this, sgrid);
			}
		}

		void contrib_rhs(grid_func<nd> & rhs)
		{
			contrib_rhs_metrics(rhs);
			for (auto & bc_callback : bc_rhs_callbacks) {
				bc_callback(rhs, *this, *sgrid);
			}
		}

		void contrib_rhs_metrics(grid_func<nd> & rhs);

		cedar::array<ptype, nd> & get_classify() { return mask; }
		const cedar::array<ptype, nd> & get_classify() const { return mask; }

		grid_func<nd> & get_bcval() { return bcval; }
		const grid_func<nd> & get_bcval() const { return bcval; }

		std::shared_ptr<grid<nd>> get_grid() { return sgrid; }

		//void assign_bc(unsigned short dim, int side, bctype btype);

	protected:
		stencil_op coeff;
		cedar::array<ptype, nd> mask;
		std::shared_ptr<grid<nd>> sgrid;
		grid_func<nd> bcval;
		std::vector<std::function<void(op<nd>&, const grid<nd> &)>> bc_op_callbacks;
		std::vector<std::function<void(grid_func<nd> &, const op<nd> &, const grid<nd> &)>> bc_rhs_callbacks;

		void register_bcs();
	};


	template<unsigned int nd>
		void apply_interior(op<nd> & sop, const grid<nd> & sgrid);
	template<unsigned int nd>
		void apply_boundary(op<nd> & sop, const grid<nd> & sgrid)
	{
		sop.call_bcs(sgrid);
	}
}

#endif
