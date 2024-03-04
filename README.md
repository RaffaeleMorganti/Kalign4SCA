# Power traces alignment using Kalign

Sources and scripts to use [Kalign](https://github.com/TimoLassmann/kalign) for power trace alignment.

## Usage

1. Launch `setup.sh` to install customized version of kalign.
2. Use `tracesToFile(traces, "misaligned.afa")` from python scripts to convert power traces to sequences.
3. Run `kalign -i "misaligned.afa" -o "aligned.afa"` to align sequences.
4. Use `tracesFromFile(traces, "aligned.afa")` from python scripts to get aligned power traces.
