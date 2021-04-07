#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

MAX_FAILURES=10

if [ "$#" -ne 1 ]; then
  echo "USAGE: $0 test_subdir" >& 2
  echo "Example: $0 bidisperse" >& 2
  echo "This will run all scripts like bidisperse/test_bidisperse_XX.sh" >& 2
  echo "It will retry failures up to a maxium of $MAX_FAILURES times" >& 2
  exit 1
fi

TEST_DIRECTORY=$1

cd $TEST_DIRECTORY

FAILURE_COUNT=0
while true ; do
    ANY_FAILURE=0
    for f in test_${TEST_DIRECTORY}_*.sh ; do
        echo "###############################################################################"
        echo "Running test $f..."
        if ./$f ; then
            echo "Test $f succeeded"
        else
            echo "Test $f failed"
            ANY_FAILURE=1
        fi
    done
    echo "###############################################################################"
    if (( ANY_FAILURE )) ; then
        echo "One or more tests failed"
        FAILURE_COUNT=$((FAILURE_COUNT + 1))
        echo "Total number of failures: $FAILURE_COUNT"
        if (( FAILURE_COUNT >= MAX_FAILURES )) ; then
            echo "Failed $MAX_FAILURES times in a row, giving up..."
            exit 1
        fi
        echo "Trying all tests again..."
    fi
    if (( ! ANY_FAILURE )) ; then break ; fi
done

echo "All tests succeeded"
