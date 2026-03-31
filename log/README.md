# log.c
A simple logging library implemented in C99.

A hard fork of https://github.com/Balth-D/log.c which is a fork of https://github.com/rxi/log.c

![screenshot](https://cloud.githubusercontent.com/assets/3920290/23831970/a2415e96-0723-11e7-9886-f8f5d2de60fe.png)



## Usage
**[log.c](src/log.c?raw=1)** and **[log.h](src/log.h?raw=1)** should be dropped
into an existing project and compiled along with it. The library provides 6
function-like macros for logging:

```c
log_trace(const char *fmt, ...);
log_debug(const char *fmt, ...);
log_info(const char *fmt, ...);
log_warn(const char *fmt, ...);
log_error(const char *fmt, ...);
log_fatal(const char *fmt, ...);
```

Each function takes a printf format string followed by additional arguments:

```c
log_trace("Hello %s", "world")
```

Resulting in a line with the given format printed to standard output:

```
20:18:26 TRACE src/main.c:11: Hello world
```

And to the output files:

```
2047-03-11 20:18:26 TRACE src/main.c:11: Hello world
```

**Detailed documentation is provided in the 'log.h' file**

### Configure library
To alter the default settings, the library has to be compiled with the 'LOGC__USER_SETTINGS' flag.
This way, a 'log_conf.h' file can be specified to alter the default settings.

#### Using CMake
You can use the library using CMake and toggle some of the configs at compile time
The current defaults are:
```
-DLOG_USE_COLOR=ON -DLOG_FILE_INFO=ON -DLOG_DEFAULT_STREAM=stderr
````

## Differences with the [original project from rxi](https://github.com/rxi/log.c)
The most interesting pull-requests (to me at least) have been integrated.
Thanks to the following for their indirect contribution:
- [**flemingoo**, Log levels re-ordering](https://github.com/rxi/log.c/pull/9)
- [**dianjixz**, C++ support](https://github.com/rxi/log.c/pull/36)
- [**chiefnoah**, Log levels renaming to avoid conflicts with Syslog](https://github.com/rxi/log.c/pull/17)
- [**Allenhe123**, Example with threads](https://github.com/rxi/log.c/pull/27)

My changes:
- Code indentation changed to 4 spaces
- Added doxygen-style documentation for all functions and structures
- Log levels enum converted to a typedef
- Ability to configure library with the 'log_conf.h' file
- Renamed 'udata' variables in code to more meaningful
- Added sample source files
    - **basic** Is a basic example, to log to stdout and to a file
    - **callbacks** Is a minimal example using callbacks with different log levels
    - **threading** Shows a threading example


## License
This library is free software; you can redistribute it and/or modify it under
the terms of the MIT license. See [LICENSE](LICENSE) for details.
