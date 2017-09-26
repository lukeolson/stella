#include <stella/grid.h>

namespace stella
{

void apply_func(grid_func<2> & gf,
                const grid<2> & sgrid,
                std::function<cedar::real_t(cedar::real_t x, cedar::real_t y)> func)
{
	auto & x = sgrid.vertex(0);
	auto & y = sgrid.vertex(1);

	for (auto j : gf.range(1)) {
		for (auto i : gf.range(0)) {
			gf(i,j) = func(x(i,j), y(i,j));
		}
	}
}


}
