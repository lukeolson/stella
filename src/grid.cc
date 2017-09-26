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


void apply_func(grid_func<3> & gf,
                const grid<3> & sgrid,
                std::function<cedar::real_t(cedar::real_t x, cedar::real_t y, cedar::real_t z)> func)
{
	auto & x = sgrid.vertex(0);
	auto & y = sgrid.vertex(1);
	auto & z = sgrid.vertex(2);

	for (auto k : gf.range(2)) {
		for (auto j : gf.range(1)) {
			for (auto i : gf.range(0)) {
				gf(i,j,k) = func(x(i,j,k), y(i,j,k), z(i,j,k));
			}
		}
	}
}

}
