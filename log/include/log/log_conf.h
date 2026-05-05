/**
* @file log.c
* @brief Log module settings file
*
* SPDX-FileCopyrightText: 2023 BDeliers
*
* SPDX-License-Identifier: MIT
*
*/

#ifndef LOG_CONF_H
#define LOG_CONF_H

#ifdef LOG_USE_COLOR
/// @brief Enable colors in standard output
#define LOGC__STDOUT_COLOR
#endif

// fallback if not provided
#ifndef LOG_DEFAULT_STREAM
#define LOG_DEFAULT_STREAM stderr
#endif

/// @brief Standard output to write to (stdout/stderr)
#define LOGC__DEFAULT_STDOUT LOG_DEFAULT_STREAM

#ifndef LOG_FILE_INFO
/// @brief Hide file path and line number for standard output
#define LOGC__STDOUT_NO_FILEINFO
#endif

#ifdef LOG_FULL_FILE_NAME
/// @brief Use full file path in logs
///        undefine to use only the file name
#define LOGC__FULL_FILE_NAME
#endif

/// @brief Time format used in the logs.
///        2 for local time string
///        1 for epoch
///        0 for none
#define LOGC__TIME_FORMAT 0

/// @brief  Maximum number of callback functions
/// @remark Will impact memory footprint
#define LOGC__MAX_CALLBACKS 5

#endif // LOG_CONF_H
