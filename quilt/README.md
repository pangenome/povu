# quilt

`quilt` is a small header-only C++17 utility library for projects that want a few modern conveniences without moving beyond C++17.

It provides:

- convenient aliases for common fixed-width and project-level types
- a few lightweight reusable helper types
- small standard library shims for functionality that became more ergonomic in C++20 and later

C++17 is still a practical baseline, but there are a handful of features from newer standards that are nice to have in day-to-day code. `quilt` gathers a small set of those conveniences in one place so code can stay shorter, clearer, and more consistent.

## Features

### `quilt::types`

Convenience aliases and helper types:

- `u8`, `u32` — fixed-width unsigned integer aliases
- `status_t` — small signed status / return type
- `Time` — alias for `std::chrono::high_resolution_clock`
- `id_t`, `idx_t` — identifier and index aliases
- `op_t<T>` — shorthand for `std::pair<T, T>`
- `unordered_pair<T>` / `up_t<T>` — pair-like type that always stores values as `(min, max)`
- `slice_t` / `slice` — compact `(start, len)` range representation

### `quilt::shim`

Small C++17-friendly compatibility helpers:

- `format` — uses `std::format` when available, otherwise falls back to `fmt::format`
- `erase_if` — provides a simple fallback for pre-C++20 code
- `contains` — convenience wrapper for containers with `.find()`
