# FLEXIS

## Prerequisties

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


### Command example
```./build/release/bin/flexis --input ./data/processed/p2p-Gnutella08.vf3.bin --support 55 --slider 1.0```
```./build/release/bin/flexis --input ./data/processed/mico.vf3.bin --support 8000 --slider 0.9```

