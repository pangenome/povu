# povu

Explore regions of variation

## Building povu


1. Fetch the source code
```
git clone --recursive https://github.com/urbanslug/povu.git
```

2. Compile
```
cmake -H. -Bbuild && cmake --build build -- -j 3
```

3. The binary should be in `./bin/povu`


### Building specific target

Building only the povu library

```
cmake -H. -DCMAKE_BUILD_TYPE=Debug -Bbuild && cmake --build build --target povulib -- -j 8
```

Building only the povu binary

```
cmake -H. -DCMAKE_BUILD_TYPE=Debug -Bbuild && cmake --build build --target povu -- -j 8
```

## Usage and Examples

For general help text run `./bin/povu -h` or just `./bin/povu`

 - **deconstruct**: finds flubbles, reports hairpin inversion boundaries
 - **call**: call variants 


## Development

To compile povu with debug symbols and with address sanitizer

```
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER=address  -H. -Bbuild && cmake --build build -- -j 3
```

## Name

The etymology of the name is rooted in profound philosophy 🤔. "Povu," is [Kiswahili](https://en.wikipedia.org/wiki/Swahili_language) for "foam." Foam, by nature, comprises countless flubbles.
