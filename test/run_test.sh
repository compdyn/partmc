#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v

((counter = 1))
while [ true ]
do
  echo Attempt $counter

  if ! $($1) &> /dev/null; then
	  echo Failure "$counter"
	  if [ "$counter" -gt $(($2)) ]
	  then
		  echo FAIL
		  exit 1
	  fi
	  echo retrying...
  else
	  echo PASS
	  exit 0
  fi
  ((counter++))
done
