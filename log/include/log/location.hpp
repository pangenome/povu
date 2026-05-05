#pragma once

#include <string>

namespace logc
{

struct location {
	const char *file;
	u_int32_t line;
	const char *function;
};

inline std::string to_string(const location &loc)
{
	return std::string{loc.file} + ":" + std::to_string(loc.line) + " in " +
	       loc.function;
}

} // namespace logc

#if defined(__GNUC__) || defined(__clang__)
#define LOG_FUNCTION_NAME __PRETTY_FUNCTION__
#elif defined(_MSC_VER)
#define LOG_FUNCTION_NAME __FUNCSIG__
#else
#define LOG_FUNCTION_NAME __func__
#endif

#define LOG_LOCATION                                                           \
	::logc::location                                                       \
	{                                                                      \
		__FILE__, __LINE__, LOG_FUNCTION_NAME                          \
	}

#define LOG_HERE ::logc::to_string(LOG_LOCATION)
