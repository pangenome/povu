/**
* @file log.c
* @brief Log module main source file
*
* SPDX-FileCopyrightText: 2020 rxd
* SPDX-FileCopyrightText: 2023 BDeliers
*
* SPDX-License-Identifier: MIT
*
*/

#include <stdio.h>
#include <string.h>

#include "../include/log.h"

/// @brief Maximum number of possible callbacks
#ifndef LOGC__MAX_CALLBACKS
    #define LOGC__MAX_CALLBACKS 10
#endif

/// @brief Default standard output
#ifndef LOGC__DEFAULT_STDOUT
    #define LOGC__DEFAULT_STDOUT stderr
#endif

/// @brief Define the used time format. 0= none, 1= timestamp, 2= local time formatted
#ifndef LOGC__TIME_FORMAT
    #define LOGC__TIME_FORMAT 2
#endif

/// @brief Log event callback structure
typedef struct
{
    log_LogFn fn;        /// @var Callback function to be called
    void *stream;        /// @var Stream to write to
    log_LogLevel level;  /// @var Log level
} Callback;

/// @brief Log management structure
struct LogManagemnt
{
    void *lock_ptr;                                /// @var Mutex used in the lock function
    log_LockFn lock_fn;                            /// @var Lock function to be called before logging
    log_LogLevel std_level;                        /// @var Standard output log level
    bool std_quiet;                                /// @var Standard output quiet mode
    Callback callbacks[LOGC__MAX_CALLBACKS];       /// @var Log event callbacks list
} LogManagemnt;
static struct LogManagemnt L = {NULL, NULL, LOGC_TRACE, false, {0}};

static const char *level_strings[] = { "FATAL", "ERROR", "WARN", "INFO", "DEBUG", "TRACE" };

#ifdef LOGC__STDOUT_COLOR
static const char *level_colors[] = {
    "\x1b[35m", "\x1b[31m", "\x1b[33m", "\x1b[32m", "\x1b[36m", "\x1b[94m"};
#endif

static void stdout_callback(log_Event *ev)
{
#if LOGC__TIME_FORMAT == 2
	char buf[16];
	buf[strftime(buf, sizeof(buf), "%H:%M:%S", localtime(&ev->time))] = '\0';
	fprintf(ev->stream, "%s ",  buf);
#elif LOGC__TIME_FORMAT == 1
	fprintf(ev->stream, "%d ",  ev->time);
#endif

#ifdef LOGC__STDOUT_COLOR
	fprintf(
		ev->stream, "%s%s\x1b[0m ",
		level_colors[ev->level], level_strings[ev->level]);

#ifndef LOGC__STDOUT_NO_FILEINFO
	fprintf(
		ev->stream, "\x1b[90m[%s (%s:%d)]\x1b[0m ",
		ev->file, ev->fn, ev->line);
#endif
#else
	fprintf(
		ev->stream, "%s ",
		level_strings[ev->level]);

#ifndef LOGC__STDOUT_NO_FILEINFO
	fprintf(
		ev->stream, "[%s (%s:%d)] ",
		ev->file, ev->fn, ev->line);
#endif
#endif
	vfprintf(ev->stream, ev->fmt, ev->ap);
	fprintf(ev->stream, "\n");
	fflush(ev->stream);
}

static void file_callback(log_Event *ev)
{
#if LOGC__TIME_FORMAT == 2
    char buf[64];
    buf[strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", localtime(&ev->time))] = '\0';
    fprintf(ev->stream, "%s ",  buf);
#elif LOGC__TIME_FORMAT == 1
    fprintf(ev->stream, "%d ",  ev->time);
#endif
    fprintf(
	ev->stream, "%-5s %s:%d: ",
	level_strings[ev->level], ev->file, ev->line);
    vfprintf(ev->stream, ev->fmt, ev->ap);
    fprintf(ev->stream, "\n");
    fflush(ev->stream);
}

static void lock(void)
{
	if (L.lock_fn)
	{
		L.lock_fn(true, L.lock_ptr);
	}
}

static void unlock(void)
{
	if (L.lock_fn)
	{
		L.lock_fn(false, L.lock_ptr);
	}
}

const char *log_level_string(log_LogLevel level)
{
    return level_strings[level];
}

void log_set_lock(log_LockFn fn, void *lock)
{
    L.lock_fn = fn;
    L.lock_ptr = lock;
}

void log_set_level(log_LogLevel level)
{
    L.std_level = level;
}

void log_set_quiet(bool enable)
{
    L.std_quiet = enable;
}

int log_add_callback(log_LogFn fn, void *stream, log_LogLevel level)
{
    for (int i = 0; i < LOGC__MAX_CALLBACKS; i++)
    {
	if (!L.callbacks[i].fn)
	{
	    L.callbacks[i] = (Callback){fn, stream, level};
	    return 0;
	}
    }
    return -1;
}

int log_add_fp(void *fp, log_LogLevel level)
{
    return log_add_callback(file_callback, fp, level);
}

static void init_event(log_Event *ev, void *stream)
{
    // Store event time
    ev->time = time(NULL);
    ev->stream = stream;
}

void log_log(log_LogLevel level, const char *file, const char *fn, int line, const char *fmt, ...)
{
#ifdef LOGC__FULL_FILE_NAME
	const char *f = file;
#else
	const char *f = strrchr(file, '/');
	if (f)
		f++;
	else
		f = file;
#endif

  // Make an event struct from the specified data
	log_Event ev = {
	    .fmt = fmt,
	    .file = f,
	    .fn = fn,
	    .line = line,
	    .level = level,
	};

	// Acquire mutex to avoid concurrent logs
	lock();

	// Log to standard output
	if (!L.std_quiet && level <= L.std_level) {
	  init_event(&ev, LOGC__DEFAULT_STDOUT);
	  va_start(ev.ap, fmt);
	  stdout_callback(&ev);
	  va_end(ev.ap);
	}

    // Call callbacks applicable to the log event
    for (int i = 0; i < LOGC__MAX_CALLBACKS && L.callbacks[i].fn; i++)
    {
	Callback *cb = &L.callbacks[i];
	if (level <= cb->level)
	{
	    init_event(&ev, cb->stream);
	    va_start(ev.ap, fmt);
	    cb->fn(&ev);
	    va_end(ev.ap);
	}
    }

    // Release mutex
    unlock();
}
