#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}
find ../src/. -iname *.h -o -iname *.c -o -iname *.cpp -o -iname *.hpp \
    | xargs clang-format -style=file -i -fallback-style=none
