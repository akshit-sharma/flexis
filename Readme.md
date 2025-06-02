# FLEXIS

## FLEXIS: FLEXible Frequent Subgraph Mining using Maximal Independent Sets

## Prerequisties

- Install [CMake](https://github.com/Kitware/CMake)
- Install [Ninja](https://github.com/ninja-build/ninja)
- Install [vcpkg](https://learn.microsoft.com/en-us/vcpkg/get_started/get-started?pivots=shell-bash) and set VCPKG_ROOT environment

## Generating

### Debug build

```cmake --preset debug```

### Release build

```cmake --preset release```

## Compiling the executable

### Debug build

```cmake --build --preset debug```

### Release build

```cmake --build --preset release```


## Download and process Datasets

```bash ./scripts/processDatasets.sh --grami```


## Running the flexis system

### Flexis command skeleton
```./build/release/bin/flexis --input <vf3.bin file> --support <int> --slider [0.0...1.0]```


### Flexis executable --help

```$ ./build/release/bin/flexis --help
flexis


./build/release/bin/flexis [OPTIONS]


OPTIONS:
  -h,     --help              Print this help message and exit
          --input TEXT:FILE REQUIRED
                              Input bin file
          --output TEXT       Output json file
  -s,     --support INT [10]  support
  -k,     --max-size INT:INT in [2 - 10] [6]
                              max size of candidate pattern
  -a,     --slider FLOAT:FLOAT in [0 - 1] [0.98]
                              slider
          --verbose [0]       verbose mode
          --print-output [0]  print output
          --output-json [0]   output json
          --print-patterns [0]
                              print patterns
```

#### slider argument
Tells the slider used for Maximal Independent Set metric

#### max-size argument
Used to tell the FLEXIS the max-size vertex pattern to explore>


### Command example

#### Example command 1
```./build/release/bin/flexis --input ./data/processed/p2p-Gnutella08.vf3.bin --support 55 --slider 1.0```

#### Example command 2
```./build/release/bin/flexis --input ./data/processed/mico.vf3.bin --support 8000 --slider 0.9```

#### Example with get the time and patterns count for frequent and merging
```./build/release/bin/flexis --input ./data/processed/p2p-Gnutella08.vf3.bin --support 55 --slider 1.0 --print-output```

```
Initial generation took 1 ms, generated 125 patterns
Step 2: starting with 125 patterns
Matching took 3 ms, reduced 125 patterns to 125 frequent patterns
Generation took 1230 ms, extended 125 patterns to 27225 new patterns
Step 3: starting with 27225 patterns
Matching took 2298 ms, reduced 27225 patterns to 2 frequent patterns
Generation took 0 ms, extended 2 patterns to 6 new patterns
Step 4: starting with 6 patterns
Matching took 0 ms, reduced 6 patterns to 0 frequent patterns
Generation took 0 ms, extended 0 patterns to 0 new patterns
Finished mining 4 steps out of 6 steps, 0 patterns left
```

#### Example with verbose flag
Prints the threashold and filename
```./build/release/bin/flexis --input ./data/processed/p2p-Gnutella08.vf3.bin --support 55 --slider 1.0 --verbose```

```
Filename ./data/processed/p2p-Gnutella08.vf3.bin
threshold is set to 55
Initial generation took 1 ms, generated 125 patterns
Step 2: starting with 125 patterns
Matching took 5 ms, reduced 125 patterns to 125 frequent patterns
Generation took 1272 ms, extended 125 patterns to 27225 new patterns
Step 3: starting with 27225 patterns
Matching took 2427 ms, reduced 27225 patterns to 2 frequent patterns
Generation took 0 ms, extended 2 patterns to 6 new patterns
Step 4: starting with 6 patterns
Matching took 0 ms, reduced 6 patterns to 0 frequent patterns
Generation took 0 ms, extended 0 patterns to 0 new patterns
Finished mining 4 steps out of 6 steps, 0 patterns left
```

#### Example with printing patterns
```./build/release/bin/flexis --input ./data/processed/p2p-Gnutella08.vf3.bin --support 55 --slider 1.0 --print-output --print-patterns```

It prints the patterns with the following format:
```
p <num_vertices> <num_edges>
v <vid> <vlabel>
...
e <src> <dst> <elabel>
...
```

The candidate and frequent patterns are printed.

- Candidate patterns are printed as:
Merged patterns: <num_of_patterns>
followed by each candidate pattern

- Frequent patterns are printed as:
Frequent patterns: <num_of_patterns>
followed by each frequent pattern

```
Initial generation took 1 ms, generated 125 patterns
Merged patterns: 125
p 2 1
v 0 2
v 1 1
e 0 1 1
p 2 1
v 0 2
v 1 0
e 0 1 0
p 2 1
v 0 4
v 1 2
e 0 1 2
p 2 1
v 0 1
v 1 3
e 0 1 0
p 2 1
v 0 0
v 1 3
e 0 1 1
...

p 2 1
v 0 3
v 1 1
e 0 1 1
p 2 1
v 0 1
v 1 2
e 0 1 2
Step 2: starting with 125 patterns
Matching took 3 ms, reduced 125 patterns to 125 frequent patterns
Frequent patterns: 125
p 2 1
v 0 2
v 1 1
e 0 1 1
...

e 0 1 0
p 2 1
v 0 3
v 1 1
e 0 1 1
p 2 1
v 0 1
v 1 2
e 0 1 2
Generation took 1275 ms, extended 125 patterns to 27225 new patterns
Merged patterns: 27225
p 3 3
v 0 0
v 1 3
v 2 3
e 2 1 2
e 0 1 4
e 0 2 3
...

p 3 3
v 0 4
v 1 3
v 2 0
e 2 0 3
e 1 2 1
e 0 1 4
Step 3: starting with 27225 patterns
Matching took 2389 ms, reduced 27225 patterns to 2 frequent patterns
Frequent patterns: 2
p 3 2
v 0 1
v 1 3
v 2 0
e 0 1 2
e 0 2 4
p 3 2
v 0 0
v 1 3
v 2 2
e 0 1 4
e 0 2 0
Generation took 0 ms, extended 2 patterns to 6 new patterns
Merged patterns: 6
p 4 3
v 0 1
v 1 3
v 2 0
v 3 0
e 0 1 2
e 0 2 4
e 0 3 4
```

