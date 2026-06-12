#ifndef POVU_STAGE_COST_HPP
#define POVU_STAGE_COST_HPP

#include <cstddef>
#include <cstdint>

namespace povu::stage_cost
{

enum class Stage : std::size_t {
	component_decomposition = 0,
	tip_dummy_augmentation,
	traversal_frame_construction,
	cycle_class_assignment,
	boundary_emission,
	hierarchy_construction,
	count
};

struct Snapshot {
	std::uint64_t calls;
	std::uint64_t input_items;
	std::uint64_t output_items;
	std::uint64_t elapsed_ns;
};

[[nodiscard]] bool enabled() noexcept;
[[nodiscard]] const char *stage_contract_name(Stage stage) noexcept;
[[nodiscard]] const char *stage_trace_name(Stage stage) noexcept;

void record(Stage stage, std::uint64_t input_items,
	    std::uint64_t output_items, std::uint64_t elapsed_ns) noexcept;
[[nodiscard]] Snapshot snapshot(Stage stage) noexcept;
void reset_for_tests() noexcept;

class Scope
{
public:
	Scope(Stage stage, std::uint64_t input_items) noexcept;
	~Scope() noexcept;

	Scope(const Scope &) = delete;
	Scope &operator=(const Scope &) = delete;
	Scope(Scope &&) = delete;
	Scope &operator=(Scope &&) = delete;

	void set_output_items(std::uint64_t output_items) noexcept;

private:
	Stage stage_;
	std::uint64_t input_items_;
	std::uint64_t output_items_;
	std::uint64_t start_ns_;
	bool active_;
};

} // namespace povu::stage_cost

#endif // POVU_STAGE_COST_HPP
