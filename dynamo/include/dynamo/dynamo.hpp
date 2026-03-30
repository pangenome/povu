#ifndef DYNAMO_IT_HPP
#define DYNAMO_IT_HPP

#include <utility>
#include <vector>

#include "dynamo/ekg_interval_tree.hpp"

namespace dynamo
{
template <class Scalar, class Value>
struct dynamic_interval_tree {
	using T = ekg::interval_tree::IntervalTree<Scalar, Value>;
	using I = ekg::interval_tree::Interval<Scalar, Value>;
	using interval_vector = std::vector<I>;
	using region = std::pair<Scalar, Scalar>;
	using region_vector = std::vector<region>;

	static_assert(std::is_integral_v<Scalar>,
		      "covered_regions() adjacency merge uses +1; Scalar must "
		      "be integral.");

	T tree;
	interval_vector base;
	interval_vector delta;
	bool dirty = false;

	void add(I v)
	{
		delta.emplace_back(std::move(v));
		dirty = true;
	}

	void commit()
	{
		if (!dirty)
			return;
		base.insert(base.end(), std::make_move_iterator(delta.begin()),
			    std::make_move_iterator(delta.end()));
		delta.clear();
		dirty = false;

		T t{interval_vector(base)};
		tree = std::move(t);
	}

	interval_vector findOverlapping(const Scalar &start,
					const Scalar &stop) const
	{
		interval_vector out = tree.findOverlapping(start, stop);
		for (auto const &iv : delta)
			if (iv.stop >= start && iv.start <= stop)
				out.push_back(iv);
		return out;
	}

	// Disjoint union coverage, merging overlaps AND adjacency (touching by
	// 1)
	region_vector covered_regions() const
	{
		interval_vector v;
		v.reserve(base.size() + delta.size());
		v.insert(v.end(), base.begin(), base.end());
		v.insert(v.end(), delta.begin(), delta.end());

		if (v.empty())
			return {};

		std::sort(v.begin(), v.end(),
			  [](const I &a, const I &b)
			  {
				  if (a.start != b.start)
					  return a.start < b.start;
				  return a.stop < b.stop;
			  });

		region_vector out;
		out.reserve(v.size());

		Scalar cur_s = v[0].start;
		Scalar cur_e = v[0].stop;

		for (std::size_t i = 1; i < v.size(); ++i) {
			const I &next = v[i];

			// Merge if overlapping OR adjacent: next.start <= cur_e
			// + 1 Guard overflow on cur_e + 1
			const bool can_inc =
				(cur_e != std::numeric_limits<Scalar>::max());
			const Scalar merge_limit =
				can_inc ? Scalar(cur_e + 1) : cur_e;

			if (next.start <= merge_limit) {
				if (next.stop > cur_e)
					cur_e = next.stop;
			}
			else {
				out.emplace_back(cur_s, cur_e);
				cur_s = next.start;
				cur_e = next.stop;
			}
		}
		out.emplace_back(cur_s, cur_e);
		return out;
	}

	// Return all Values whose intervals contain pos
	std::vector<Value> values_at(const Scalar &pos) const
	{
		std::vector<Value> out;

		// From indexed tree
		tree.visit_overlapping(pos, [&](const I &iv)
				       { out.push_back(iv.value); });

		// From unindexed delta
		for (auto const &iv : delta) {
			if (iv.start <= pos && pos <= iv.stop) {
				out.push_back(iv.value);
			}
		}

		return out;
	}

	// Optional: query covered regions intersecting [start, stop]
	region_vector covered_regions_overlapping(const Scalar &start,
						  const Scalar &stop) const
	{
		auto cov = covered_regions();
		region_vector out;
		for (auto const &r : cov)
			if (r.second >= start && r.first <= stop)
				out.push_back(r);
		return out;
	}

	// Optional: rebuild the tree on the merged coverage (turn it into a
	// "coverage index")
	void rebuild_as_coverage()
	{
		auto cov = covered_regions();
		base.clear();
		base.reserve(cov.size());
		for (auto const &[s, e] : cov)
			base.emplace_back(s, e,
					  Value{}); // only if Value is
						    // default-constructible

		delta.clear();
		dirty = false;

		T t{interval_vector(base)};
		tree = std::move(t);
	}
};

}; // namespace dynamo

#endif // DYAMO_IT_HPP
