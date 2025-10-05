# Development Guide

While not mandatory, running your changes through the following tools is highly encouraged. If not, contributions are still welcome—a maintainer will review them in due time.

## Toolchain and Compiler Preferences

Although `povu` does not enforce the use of a specific compiler, the development toolchain primarily relies on **Clang/LLVM** for code style and linting. GCC would be preferavble for the sake of freedom but the there seems to be little development tooling around it.

`povu` specifically uses the following tools:
- [clang-format](https://clang.llvm.org/docs/ClangFormat.html)
- [clang-tidy](https://clang.llvm.org/extra/clang-tidy/)
- [Include What You Use (IWYU)](https://github.com/include-what-you-use/include-what-you-use)



### Code Style

`povu` follows a **C++ adaptation** of the [Linux kernel coding style](https://www.kernel.org/doc/html/v4.12/process/coding-style.html). The style configuration is defined in the `.clang-format` file located in the project root.

For more details on configuring `clang-format` style options, see [ClangFormatStyleOptions](https://clang.llvm.org/docs/ClangFormatStyleOptions.html).

---

## Linting: `clang-tidy`

`clang-tidy` is used for linting your code, catching potential bugs, and enforcing consistent style. Be sure to integrate it during development for high-quality contributions.

---

## Formatting: `clang-format`

To ensure consistent formatting, run `clang-format` on your code using the following command:

```bash
find tests include src app -name '*.cpp' -o -name '*.hpp' -o -name '*.cc' | xargs clang-format -i
```

This will apply the defined formatting style to all source files.

---

## Include What You Use (IWYU)

**IWYU** helps manage header file includes by ensuring that each file directly includes what it uses. Follow these steps to use IWYU:

### Running with Clang and IWYU
The following commands enable IWYU during the compilation process:
```bash
CC="clang" CXX="clang++" cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE=include-what-you-use -Bbuild -S.
cmake --build build
```

---

### Running IWYU with GCC or without Clang Integration
If you're using GCC or another compiler, you can still use IWYU manually. First, make sure to enable the export of compile commands:

```bash
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -Bbuild -S.
```

Then, run IWYU by pointing it to the build directory as follows:
```bash
iwyu-tool -p build -- -Xiwyu --no_fwd_decls > iwyu.log
iwyu-fix-includes --ignore '(^|.*/)_deps/' -b --comments --nosafe_headers < iwyu.log
```

This generates an IWYU report (`iwyu.log`) and optionally modifies the includes for you.
