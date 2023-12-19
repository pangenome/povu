# povu
Variant caller based on cycle equivalence

# Input
Input GFA

Expect a sorted graph in GFA
Expect first node to have id 1

Compile
```
cmake -H. -Bbuild && cmake --build build -- -j 3
```

Run
```
./bin/povu 
```

Example
```
./bin/povu -v 2 call -i test_data/LPA.max120.gfa  --  HG02572__LPA__tig00000001
```

## Development

Compile with debug symbols
```
cmake -DCMAKE_BUILD_TYPE=Debug -H. -Bbuild && cmake --build build -- -j 3
```

A release version (default)
```
cmake -DCMAKE_BUILD_TYPE=Release -H. -Bbuild && cmake --build build -- -j 3
```
