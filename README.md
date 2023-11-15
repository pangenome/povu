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

## Development

Compile with debug symbols
```
cmake -DCMAKE_BUILD_TYPE=Debug -H. -Bbuild && cmake --build build -- -j 3
```

A release version (default)
```
cmake -DCMAKE_BUILD_TYPE=Release -H. -Bbuild && cmake --build build -- -j 3
```
