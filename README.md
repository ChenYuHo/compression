# compression

## Build
```bash
rm -rf build && mkdir build && cd build
# -DNORANDOM=1 will do deterministic rounding for intml
cmake .. [-DDEBUG=1] [-DNORANDOM=1]
make
```

## Test
```bash
# after cmake
make test
./test
```

## Run
```
Usage:
  ./compress [OPTION...]

  -s, --size arg                 number of float32 elements (default: 26214400 if DEBUG flag is off, else 100)
  -i, --input arg                input file name, randomly generate elements if not provided
  -o, --data-output arg          output file name to save original_data (no-op if given input file)
  -r, --compressed-output arg    output file name to save compressed compressed_data
      --decompressed-output arg  output file name to save decompressed compressed_data
  -m, --method arg               compress method (default: intml)
      --repeat arg               repeat compression and/or decompression (default: 10)
  -p, --print                    print original, compressed, and decompressed data
  -c, --compress                 do and measure compression
  -d, --decompress               do and measure decompression
  -h, --help                     print usage
```

## CUDA
```
nvcc -O3 cnat.cu -o cnat_cuda -lcurand
./cnat_cuda
```
