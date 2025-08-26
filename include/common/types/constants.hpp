#ifndef POVU_TYPES_CONSTANTS_HPP
#define POVU_TYPES_CONSTANTS_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <sys/types.h>

#include "./core.hpp"
#include "../../graph/types.hpp"


namespace povu::constants {

// colors for DOT format
inline constexpr std::string_view GRAY{"gray"};
inline constexpr std::string_view BLACK{"black"};
inline constexpr std::string_view RED{"red"};
inline constexpr std::string_view BLUE{"blue"};

// numeric
inline constexpr std::size_t SIZE_T_MIN = std::numeric_limits<size_t>::min();
inline constexpr std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
inline constexpr int UNDEFINED_INT = std::numeric_limits<int>::min();
inline constexpr std::size_t UNDEFINED_SIZE_T = SIZE_T_MAX;
inline constexpr pt::idx_t MAX_ID = std::numeric_limits<pt::idx_t>::max();
inline constexpr pt::idx_t MAX_IDX = std::numeric_limits<pt::idx_t>::max();

inline constexpr pt::idx_t UNDEFINED_IDX = MAX_IDX; // TODO: replace with INVALID
inline constexpr pt::id_t UNDEFINED_ID = MAX_ID;    // TODO: replace with INVALID
inline constexpr pt::id_t DUMMY_VTX_ID = UNDEFINED_ID;
inline constexpr pt::id_t INVALID_ID = MAX_ID;
inline constexpr pt::idx_t INVALID_IDX = MAX_IDX;
inline constexpr pt::idx_t INVALID_CLS = MAX_IDX; // equivalence class

//
inline constexpr ptg::id_or_t INVALID_ID_OR{INVALID_ID, ptg::or_e::forward};

// strings
inline constexpr std::string_view EMPTY_SET = "\u2205";
inline constexpr std::string_view UNDEFINED_VALUE = EMPTY_SET;
inline constexpr std::string_view WAVY_ARROW = "\u2933";
inline constexpr std::string_view INF = "\u221E"; //infinity

// genomics
inline constexpr std::string_view UNDEFINED_PATH_LABEL{"undefined"};
inline constexpr std::size_t UNDEFINED_PATH_ID{INVALID_ID};
inline constexpr std::size_t UNDEFINED_PATH_POS{INVALID_ID};

// PVST

inline constexpr char PVST_HEADER_SYMBOL = 'H';
inline constexpr char PVST_FLUBBLE_SYMBOL = 'F';
inline constexpr char PVST_CONCEALED_SYMBOL = 'C';
inline constexpr char PVST_SMOTHERED_SYMBOL = 'S';
inline constexpr char PVST_DUMMY_SYMBOL = 'D';
inline constexpr char PVST_TINY_SYMBOL = 'T';
inline constexpr char PVST_OVERLAP_SYMBOL = 'O';
inline constexpr char PVST_MIDI_SYMBOL = 'M';

inline constexpr std::string_view PVST_VERSION = "0.0.3";

// expected number of columns in a PVST file
// size_t instead of casting before comparison with vec size
inline constexpr std::size_t EXPECTED_PVST_COL_NUMS = 5;

// double space for formatting
inline constexpr std::string_view SHORT_TAB = "  ";

// VCF
inline constexpr char COL_SEP = '\t'; // column separator
inline constexpr char NO_VALUE = '.'; // null character
} // namespace povu::constants

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pc = povu::constants;

#endif
