#ifndef POVU_TYPES_CONSTANTS_HPP
#define POVU_TYPES_CONSTANTS_HPP

#include <cstddef>
#include <string>
#include <sys/types.h>

#include "./core.hpp"

namespace povu::constants {

// colors
const std::string GRAY{"gray"};
const std::string BLACK{"black"};
const std::string RED{"red"};
const std::string BLUE{"blue"};

// numeric
const std::size_t SIZE_T_MIN = std::numeric_limits<size_t>::min();
const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
const int UNDEFINED_INT = std::numeric_limits<int>::min();
const std::size_t UNDEFINED_SIZE_T = SIZE_T_MAX;
const pt::idx_t MAX_ID = std::numeric_limits<pt::idx_t>::max();
const pt::idx_t MAX_IDX = std::numeric_limits<pt::idx_t>::max();

const pt::idx_t UNDEFINED_IDX = MAX_IDX; // TODO: replace with INVALID
const pt::id_t UNDEFINED_ID = MAX_ID;    // TODO: replace with INVALID
const pt::id_t DUMMY_VTX_ID = UNDEFINED_ID;
const pt::id_t INVALID_ID = MAX_ID;
const pt::idx_t INVALID_IDX = MAX_IDX;
const pt::idx_t INVALID_CLS = MAX_IDX; // equivalence class

// strings
const std::string EMPTY_SET = "\u2205";
const std::string UNDEFINED_VALUE = EMPTY_SET;
const std::string WAVY_ARROW = "\u2933";
const std::string INF = "\u221E"; //infinity

// genomics constants
const std::string UNDEFINED_PATH_LABEL{"undefined"};
const std::size_t UNDEFINED_PATH_ID{INVALID_ID};
const std::size_t UNDEFINED_PATH_POS{INVALID_ID};

const char PVST_HEADER_SYMBOL = 'H';
const char PVST_FLUBBLE_SYMBOL = 'F';
const char PVST_CONCEALED_SYMBOL = 'C';
const char PVST_SMOTHERED_SYMBOL = 'S';
const char PVST_DUMMY_SYMBOL = 'D';
const char PVST_TINY_SYMBOL = 'T';
const char PVST_OVERLAP_SYMBOL = 'O';
const char PVST_MIDI_SYMBOL = 'M';

const std::string PVST_VERSION = "0.2";

// VCF
const char COL_SEP = '\t'; // column separator
const char NO_VALUE = '.'; // null character
} // namespace povu::constants

#endif
