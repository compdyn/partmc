#!/bin/bash

# exit on error
set -e
# turn on command echoing
#set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

((counter = 1))
while [ true ]
do
  echo Attempt $counter

#todo change .txt to .csv and save with the .c code instead of bash

#Arguments are netcdf read parameters: i, j, k, t, i_count, j_count, k_count, t_count
i=(136 136)
j=(134 134)
k=(1 1)
t=(14 14)
i_n=(1 5)
j_n=(1 5)
k_n=(1 5)
t_n=(1 1)

#echo "Input arguments: \
#i=(136 136)
#j=(134 134)
#k=(1 1)
#t=(14 14)
#i_n=(1 5)
#j_n=(1 5)
#k_n=(1 10)
#t_n=(1 1)
#" > ../../../../../conf_test.csv

#todo pass profile_stats string and num_time_steps as a parameter

#Save test config for stats
echo "i,j,k,t,i_n,j_n,k_n,t_n" > ../../../../../profile_stats.csv

#Test with different parameters
for id in ${!i_n[@]};
do

#todo each iteration create a new line on the profile_stats.csv file

    echo "${i[$id]},${j[$id]},${k[$id]},\
${t[$id]},${i_n[$id]},${j_n[$id]},${k_n[$id]},${t_n[$id]}" >> ../../../../../profile_stats.csv

    exec_str="../../../test_chemistry_cb05cl_ae5_big ${i[$id]} ${j[$id]} ${k[$id]} \
${t[$id]} ${i_n[$id]} ${j_n[$id]} ${k_n[$id]} ${t_n[$id]}"

    if ! $exec_str; then
        echo Failure "$counter"
         if [ "$counter" -gt 1 ]
        then
          echo FAIL
          exit 1
        fi
        ((counter++))
    fi

    #if ! ../../../test_chemistry_cb05cl_ae5_big; then
    #if ! $exec_str; then
    #      echo Failure "$counter"
    #      if [ "$counter" -gt 1 ]
    #      then
    #          echo FAIL
    #         exit 1
    #      fi
    #      echo retrying...
    #  else
    #      echo PASS
    #      exit 0
    #  fi
    #  ((counter++))
    #done





    done
echo PASS
exit 0
done
