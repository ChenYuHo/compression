# compression

## Build
```bash
rm -rf build && mkdir build && cd build
# -DNORANDOM=1 will do deterministic rounding for intml
cmake .. [-DDEBUG=1] [-DNORANDOM=1]
make
```

## Run
```
Usage:
  ./compress [OPTION...]

  -s, --size arg           number of float32 elements (default: 26214400 if DEBUG flag is off, else 100)
  -i, --input arg          input file name, randomly generate elements if not provided
  -o, --data-output arg    output file name to save data (no-op if given input file)
  -r, --result-output arg  output file name to save compressed result
  -m, --method arg         compress method (default: intml)
  --repeat arg             repeat compression and/or decompression (default: 10)
  -p, --print              print data and result
  -c, --compress           do and measure compression
  -d, --decompress         do and measure decompression
  -h, --help               Print help
```
