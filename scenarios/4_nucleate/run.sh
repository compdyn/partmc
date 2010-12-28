#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

mkdir -p out

../../build/partmc nucleate_with_coag.spec
../../build/partmc nucleate_no_coag.spec
