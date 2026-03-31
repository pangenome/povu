# liteseq

A variation graph implementation that takes [GFA](https://gfa-spec.github.io/GFA-spec/) as input.

## Features
- [x] **Configurable Parsing**: Easily toggle the inclusion of vertex labels, references, and more.
- [x] **Minimal Graph**: Optimize performance by retrieving only essential data.
- **GFA Version Support**:
  - [x] [GFA v1.0](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) (`VN:Z:1.0`)
  - [ ] Planned: [GFA v1.1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#gfa-11) (`VN:Z:1.1`)

## Limitations
- The GFA format represents a superset of a variation graph. Therefore, liteseq does not support every feature in the GFA specification.
- Sequence IDs in the GFA must be numerical. Non-numeric IDs are not supported and may result in parsing errors.

## Configuration

To ensure that we parse only the essential subset of GFA data, liteseq uses a
`gfa_config` struct to determine which graph properties to read from the GFA file.
The fields of the `gfa_config` struct are described in the table below.


| Field            | Type     | Description                                                           |
|:-----------------|:---------|:----------------------------------------------------------------------|
| `gfa_file_path`  | `char *` | The path to the GFA file to be parsed. Ensure the file path is valid. |
| `inc_vtx_labels` | `bool`   | If set to `true`, vertex labels will be read.                         |
| `inc_refs`       | `bool`   | If set to `true`, reference fields in the GFA file will be parsed**.    |


**Note**
To verify successful parsing, check if `g->status == 0` in the `gfa_props` returned by `gfa_new`.

## Building liteseq

Prerequisites:
  - CMake (3.0+ recommended)
  - C compiler (e.g., GCC or Clang)


### Fetch the source code
```
git clone https://github.com/urbanslug/liteseq.git

cd liteseq
```


### Release mode
```
cmake -DCMAKE_BUILD_TYPE=Release  -H. -Bbuild && cmake --build build -- -j 3
```

### Debug mode
```
cmake -DCMAKE_BUILD_TYPE=Debug  -H. -Bbuild && cmake --build build -- -j 3
```


### Examples

To compile the examples binary set `LITESEQ_BUILD_EXAMPLE` `ON` when configuring the build

```
cmake -DLITESEQ_BUILD_EXAMPLE=ON -DCMAKE_BUILD_TYPE=Debug  -H. -Bbuild && cmake --build build -- -j 3
```
Run examples with `./bin/liteseq-example <path/to/gfa>`.

## Usage and Examples

1. Include Headers:

```
#include <liteseq/gfa.h>
```

2. Set Up Configuration:

```
gfa_config config = {
    .inc_vtx_labels = false,
    .inc_refs = false,
    // ...
};
```
3. Initialize, Parse, and Use:

```
gfa_props *gfa = gfa_new(config);
// ...
// parse, process, etc.
```

4. Cleanup:

```
gfa_free(gfa);
```

See [`examples/example.c`](./examples/example.c) for a more detailed example.

## Contributing
Contributions, bug reports, and feature requests are welcome.
Please open an issue or a pull request here on GitHub.

## License
MIT

Happy parsing!
