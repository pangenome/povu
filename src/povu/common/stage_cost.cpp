#include "povu/common/stage_cost.hpp"

#include <array>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <mutex>

namespace povu::stage_cost
{
namespace
{

struct Counter {
	std::atomic<std::uint64_t> calls{0};
	std::atomic<std::uint64_t> input_items{0};
	std::atomic<std::uint64_t> output_items{0};
	std::atomic<std::uint64_t> elapsed_ns{0};
};

constexpr std::size_t STAGE_COUNT = static_cast<std::size_t>(Stage::count);

std::array<Counter, STAGE_COUNT> counters{};
std::once_flag report_once;

constexpr std::array<const char *, STAGE_COUNT> CONTRACT_NAMES{
	"componentDecomposition",
	"tipDummyAugmentation",
	"traversalFrameConstruction",
	"cycleClassAssignment",
	"boundaryEmission",
	"hierarchyConstruction",
};

constexpr std::array<const char *, STAGE_COUNT> TRACE_NAMES{
	"component-decomposition",
	"tip-dummy-augmentation",
	"traversal-frame-construction",
	"cycle-class-assignment",
	"flubble-boundary-extraction",
	"flubble-hierarchy-construction",
};

[[nodiscard]] std::uint64_t now_ns() noexcept
{
	using clock = std::chrono::steady_clock;
	return static_cast<std::uint64_t>(
		std::chrono::duration_cast<std::chrono::nanoseconds>(
			clock::now().time_since_epoch())
			.count());
}

[[nodiscard]] std::size_t idx(Stage stage) noexcept
{
	return static_cast<std::size_t>(stage);
}

void write_report() noexcept
{
	try {
		std::cerr << "povu-stage-cost-trace schema=povu.stage-cost.v1"
			  << '\n';
		for (std::size_t i{}; i < STAGE_COUNT; ++i) {
			const auto stage = static_cast<Stage>(i);
			const Snapshot value = snapshot(stage);
			std::cerr << "povu-stage-cost"
				  << " contract=" << stage_contract_name(stage)
				  << " trace=" << stage_trace_name(stage)
				  << " calls=" << value.calls
				  << " input_items=" << value.input_items
				  << " output_items=" << value.output_items
				  << " elapsed_ns=" << value.elapsed_ns << '\n';
		}
	}
	catch (...) {
	}
}

void register_report() noexcept
{
	std::call_once(report_once, []
		      {
			      std::atexit(write_report);
		      });
}

} // namespace

bool enabled() noexcept
{
	const char *value = std::getenv("POVU_STAGE_COST_TRACE");
	return value != nullptr && value[0] != '\0' && value[0] != '0';
}

const char *stage_contract_name(Stage stage) noexcept
{
	const std::size_t i = idx(stage);
	return i < STAGE_COUNT ? CONTRACT_NAMES[i] : "unknown";
}

const char *stage_trace_name(Stage stage) noexcept
{
	const std::size_t i = idx(stage);
	return i < STAGE_COUNT ? TRACE_NAMES[i] : "unknown";
}

void record(Stage stage, std::uint64_t input_items,
	    std::uint64_t output_items, std::uint64_t elapsed_ns) noexcept
{
	const std::size_t i = idx(stage);
	if (i >= STAGE_COUNT)
		return;

	register_report();
	Counter &counter = counters[i];
	counter.calls.fetch_add(1, std::memory_order_relaxed);
	counter.input_items.fetch_add(input_items, std::memory_order_relaxed);
	counter.output_items.fetch_add(output_items, std::memory_order_relaxed);
	counter.elapsed_ns.fetch_add(elapsed_ns, std::memory_order_relaxed);
}

Snapshot snapshot(Stage stage) noexcept
{
	const std::size_t i = idx(stage);
	if (i >= STAGE_COUNT)
		return {};

	Counter &counter = counters[i];
	return {
		counter.calls.load(std::memory_order_relaxed),
		counter.input_items.load(std::memory_order_relaxed),
		counter.output_items.load(std::memory_order_relaxed),
		counter.elapsed_ns.load(std::memory_order_relaxed),
	};
}

void reset_for_tests() noexcept
{
	for (Counter &counter : counters) {
		counter.calls.store(0, std::memory_order_relaxed);
		counter.input_items.store(0, std::memory_order_relaxed);
		counter.output_items.store(0, std::memory_order_relaxed);
		counter.elapsed_ns.store(0, std::memory_order_relaxed);
	}
}

Scope::Scope(Stage stage, std::uint64_t input_items) noexcept
    : stage_(stage), input_items_(input_items), output_items_(0),
      start_ns_(0), active_(enabled())
{
	if (active_)
		start_ns_ = now_ns();
}

Scope::~Scope() noexcept
{
	if (!active_)
		return;

	const std::uint64_t end_ns = now_ns();
	record(stage_, input_items_, output_items_, end_ns - start_ns_);
}

void Scope::set_output_items(std::uint64_t output_items) noexcept
{
	output_items_ = output_items;
}

} // namespace povu::stage_cost
