#ifndef STELLA_GRID_H
#define STELLA_GRID_H

#include <array>

#include <cedar/kernel_params.h>
#include <cedar/halo_exchanger.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/array.h>
#include <cedar/types.h>
#include <cedar/mpi/grid_topo.h>

#include <stella/metrics.h>

namespace stella
{
	template<unsigned int nd,
		class halo_exchanger>
	class grid : public cedar::grid_topo
	{
		using real_t = cedar::real_t;
		using len_t = cedar::len_t;
	public:
	grid(std::array<len_t, nd> len): grid_topo(std::make_shared<std::vector<len_t>>(NBMG_pIGRD),0,1)
		{
			for (len_t i = 0; i < nd; i++) {
				nlocal(i) = len[i] + 2;
				is(i) = 1;
				if (nd == 2)
					verts[i] = cedar::array<real_t, nd>(len[0]+2, len[1]+2);
				else
					verts[i] = cedar::array<real_t, nd>(len[0]+2, len[1]+2, len[2]+2);
			}
		}
		cedar::array<real_t, nd> & vertex(int dim) { return verts[dim]; }
		const cedar::array<real_t, nd> & vertex(int dim) const { return verts[dim]; }
		len_t ie(int dim) const
		{
			return is(dim) + (nlocal(dim) - 2) - 1;
		}

		len_t ibase(int dim) const {
			return 1; // number of ghosts
		}

		len_t ibound(int dim) const {
			return nlocal(dim) - 1;
		}

		len_t isize(int dim) const {
			return nlocal(dim) - 2;
		}

		int mpi_rank() const
		{
			int rank;
			MPI_Comm_rank(this->comm, &rank);
			return rank;
		}

		int mpi_size() const
		{
			int size;
			MPI_Comm_size(this->comm, &size);
			return size;
		}

		bool mpi_root() const
		{
			return (mpi_rank() == 0);
		}

		void init_halo()
		{
			if (halof == nullptr) {
				cedar::config::reader conf("config.json");
				auto params = cedar::build_kernel_params(conf);
				halof = std::make_unique<cedar::halo_exchanger<nd>>(*params, *this);
			}
		}

		halo_exchanger & get_halo_exchanger()
		{
			return *halof;
		}

		template<class T>
			void exchange_halo(T&& q) const
		{
			halof->exchange(std::forward<T>(q));
		}

		void compute_metrics()
		{
			if (met == nullptr) {
				met = std::make_unique<stella::metrics<nd>>(this->get_reference());
				compute_interior(*met, *this);
				compute_boundary(*met, *this);
			}
		}

		stella::metrics<nd> & metrics() {
			return *met;
		}

		const stella::metrics<nd> & metrics() const {
			return *met;
		}

		void store_reference(std::shared_ptr<cedar::grid_topo> ptr)
		{
			ref = ptr;
		}

		std::shared_ptr<cedar::grid_topo> get_reference()
		{
			return ref;
		}


	protected:
		std::array<cedar::array<real_t, nd>, nd> verts;
		std::unique_ptr<stella::metrics<nd>> met;
		std::unique_ptr<halo_exchanger> halof;
		std::shared_ptr<cedar::grid_topo> ref;
	};
}

#endif

