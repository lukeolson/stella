#include <stella/grid.h>
#include <stella/op.h>

int main(int argc, char *argv[])
{
	using namespace cedar;
	std::cout << "Hello World!" << std::endl;
	auto grd = std::make_shared<stella::grid<2>>(std::array<len_t,2>({4,4}));

	stella::op<2> sop(grd);

	std::cout << sop.get_coefficients() << std::endl;

	return 0;
}
