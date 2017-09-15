#ifndef STELLA_HALO_H
#define STELLA_HALO_H

#include <stella/grid.h>

namespace stella
{
	template<unsigned int nd, class T>
		void exchange_halo(const grid<nd> & sgrid, T&& q)
	{
		sgrid.exchange_halo(std::forward<T>(q));
	}
}

#endif
