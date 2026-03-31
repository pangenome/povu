/**
* @file log.c
* @brief Log module main header file
*
* SPDX-FileCopyrightText: 2020 rxd
* SPDX-FileCopyrightText: 2023 BDeliers
*
* SPDX-License-Identifier: MIT
*
*/

#ifndef LOG_H
#define LOG_H

#ifdef __cplusplus
extern "C"
{
#endif


#include <stdarg.h>
#include <stdbool.h>
#include <time.h>

/// @remark Compile with LOGC__USER_SETTINGS defined to include log_conf.h and be able to modify the above settings
#ifdef LOGC__USER_SETTINGS
    #include "../src/log_conf.h"
#endif

#define LOGC_VERSION_MAJOR 0U
#define LOGC_VERSION_MINOR 3U
#define LOGC_VERSION_TEST  0U

/// @brief Log event struct
typedef struct
{
    va_list ap;         /// @var Variable arguments sent to the log_log functions
    const char *fmt;    /// @var printf-style format string
    const char *file;   /// @var Full file path where the log event was raised
    const char *fn;   /// @var Function name where the log event was raised
    time_t time;        /// @var Time of the log event
    void *stream;       /// @var Stream to log to
    int line;           /// @var Line in the source file where the log event was raised
    int level;          /// @var Log level
} log_Event;

/// @brief      Log event callback funtion type
/// @param ev   Log event structure pointer
typedef void (*log_LogFn)(log_Event *ev);

/// @brief          Lock function type
/// @param lock     True if the resource has to be locked, false otherwise
/// @param lock_ptr Pointer to the mutex resource
typedef void (*log_LockFn)(bool lock, void *lock_ptr);

/// @brief  Log levels enum type
typedef enum
{
    LOGC_FATAL,
    LOGC_ERROR,
    LOGC_WARN,
    LOGC_INFO,
    LOGC_DEBUG,
    LOGC_TRACE
}
log_LogLevel;

#define log_trace(...) log_log(LOGC_TRACE, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_debug(...) log_log(LOGC_DEBUG, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_info(...) log_log(LOGC_INFO, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_warn(...) log_log(LOGC_WARN, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_error(...) log_log(LOGC_ERROR, __FILE__, __func__,__LINE__, __VA_ARGS__)
#define log_fatal(...) log_log(LOGC_FATAL, __FILE__, __func__,__LINE__, __VA_ARGS__)

/// @brief          Get the log level as string
/// @param level    Log level
/// @return         String representing the corresponding log level
const char *log_level_string(log_LogLevel level);


/// @brief          Set the process lock function
/// @remark         The lock function is passed true to acquire the resource, false to release
/// @param fn       Function to be called to lock/unlock the log module
/// @param lock     Mutex pointer
void log_set_lock(log_LockFn fn, void *lock);

/// @brief          Set the standard output log level
/// @remark         Callbacks and files are not concerned by this setting
///                 Default to LOGC_TRACE
/// @param level    Desired log level
void log_set_level(log_LogLevel level);

/// @brief          Set the standard output to quiet mode
/// @remark         Callbacks and files are not concerned by this setting
///                 Default to false
/// @param enable   True to enable quiet mode, false to disable
void log_set_quiet(bool enable);

/// @brief          Add a callback to the corresponding log events
/// @remark         The callback function will be passed a log_Event structure when called
/// @param fn       Callback function
/// @param stream   Pointer to the output stream
/// @param level    Corresponding log level
/// @return         0 if success, -1 on error
int log_add_callback(log_LogFn fn, void *stream, log_LogLevel level);

/// @brief          Add a file pointer to the corresponding log events
/// @param fp       File pointer
/// @param level    Corresponding log level
/// @return         0 if success, -1 on error
int log_add_fp(void *fp, log_LogLevel level);

/// @brief          Log function
/// @param level    Level to log
/// @param file     File from where the log is raised
/// @param line     Line in the file where the log is raised
/// @param fmt      printf-style format
/// @param va       Variable arguments
void log_log(log_LogLevel level, const char *file, const char *fn, int line, const char *fmt, ...);

#ifdef __cplusplus
}
#endif
#endif // LOG_H
