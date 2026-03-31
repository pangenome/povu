#include "./shared.hpp"

namespace shared
{

void fill_random(meza::pool::split::ov_mat_t &mat)
{
	qt::u32 I = mat.rows();
	qt::u32 J = mat.cols();
	for (qt::u32 i = 0; i < I; i++)
		for (qt::u32 j = 0; j < J; j++)
			mat.set_value(i, j, static_cast<qt::u8>(rand() % 2));
}

void fill_ref_row(const meza::pool::split::ov_mat_t &filter_mat,
		  meza::pool::split::ov_mat_t &ref_mat)
{
	qt::u32 I = ref_mat.rows();
	qt::u32 J = ref_mat.cols();

	for (qt::u32 i = 0; i < I; i++) {
		for (qt::u32 j = 0; j < J; j++) {
			qt::u8 val = filter_mat.get_value(0, j);
			ref_mat.set_value(i, j, val);
		}
	}
}

} // namespace shared
