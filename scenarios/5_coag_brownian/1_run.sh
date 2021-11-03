#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

mkdir -p out

../../build/partmc example.spec
