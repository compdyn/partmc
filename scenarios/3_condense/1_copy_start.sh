#!/bin/sh

echo "This script copies the output of 2_urban_plume2 to serve as the starting state"
echo "You must have run scenarios/2_urban_plume2 first"

# exit on error
set -e
# turn on command echoing
set -v

cp -r ../2_urban_plume2/out start
