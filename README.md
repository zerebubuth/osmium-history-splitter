# (Lib)Osmium History Splitter

This is a sort-of port or rewrite of [MaZderMind's history splitter](https://github.com/MaZderMind/osm-history-splitter) to use [libosmium](https://github.com/osmcode/libosmium) rather than the older, [apparently deprecated](https://github.com/joto/osmium/blob/master/README#L5-L11), library.

It's kinda working at the moment, and it has a few basic unit tests. But still needs more work before you should use it.

## Building

It uses [CMake](https://cmake.org/), and has no extra dependencies beyond [libosmium](https://github.com/osmcode/libosmium) and what it needs. You can build using something like this:

```
mkdir build
cd build
cmake ..
```

There are tests, which you can run using `ctest`, or by running the individual tests from the `build/test/` directory. If you do this, you should set `OSMIUM_HISTORY_SPLITTER_TEST_DATA_DIR` to be the source test subdirectory, e.g:

```
OSMIUM_HISTORY_SPLITTER_TEST_DATA_DIR=../test test/unit_tests
```

## Bugs

Please report any bugs to the issues page. Reports are very welcome, and those with test cases and patches even more so!

