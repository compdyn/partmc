#!/bin/bash

for d in r?? ; do
  cd $d;
  cp ../start_state0800_large.d .
  ../run_sedi_fix_hybrid
  for f in state_????.d ; do
    ../process_state $f
    rm $f
  done
  cd ..
done

