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

#Config input cases
READ_NETCDF=false

#todo change .txt to .csv and save with the .c code instead of bash
#maybe most clear option is here only repeat X times (maybe option for ALL columns), and then read
#the file from csv inside code?
#well I prefer put all in bash because I can set options like READ_NETCDF which means one config. or another
#if all in excel it would be difficult to identify this options of testing, and multiple csv are not the way
#el problema es que bueno, esto tmb lo deberia entender el codigo no? xDDD FAIL
#well atleast read all kind of parameters I think is not bad.. so it could still be in the code but not used

if $READ_NETCDF ; then

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
  echo "i,j,k,t,i_n,j_n,k_n,t_n,test_tsteps" > ../../../../../profile_stats.csv

#check total repeats and which must be repeated
#vale si lo suyo seria un for por cada parametro

#Test with different parameters
#for id in ${!i_n[@]}; #good
  for id in "${!i_n[@]-1}"; #test only 1 execution
  do

  #todo each iteration create a new line on the profile_stats.csv file
      echo "${i[$id]},${j[$id]},${k[$id]},\
  ${t[$id]},${i_n[$id]},${j_n[$id]},${k_n[$id]},${t_n[$id]},1" >> ../../../../../profile_stats.csv


  # Original
      exec_str="../../../test_chemistry_cb05cl_ae5_big ${i[$id]} ${j[$id]} ${k[$id]} \
  ${t[$id]} ${i_n[$id]} ${j_n[$id]} ${k_n[$id]} ${t_n[$id]}"

  #need cd to directory, and the info isn't completed (missing cvode and cuda)
  #    exec_str="gprof ../../../test_chemistry_cb05cl_ae5_big ${i[$id]} ${j[$id]} ${k[$id]} \
  #${t[$id]} ${i_n[$id]} ${j_n[$id]} ${k_n[$id]} ${t_n[$id]} > gprof.txt"
  #-I $pwd

  #--print-gpu-summary
  #--analysis-metrics -f -o ../../../../mock_monarch_1000.nvprof
      #exec_str="nvprof --print-gpu-summary \
      #../../../test_chemistry_cb05cl_ae5_big ${i[$id]} ${j[$id]} ${k[$id]} \
      #${t[$id]} ${i_n[$id]} ${j_n[$id]} ${k_n[$id]} ${t_n[$id]}"

  #module load ddt #not working
  #    exec_str="ddt --connect ../../../test_chemistry_cb05cl_ae5_big ${i[$id]} ${j[$id]} ${k[$id]} \
  #${t[$id]} ${i_n[$id]} ${j_n[$id]} ${k_n[$id]} ${t_n[$id]}"

  # MPI (todo improve personalization)
      #exec_str="mpirun -v -np 2 ../../../test_chemistry_cb05cl_ae5_big ${i[$id]} ${j[$id]} ${k[$id]} \
  #${t[$id]} ${i_n[$id]} ${j_n[$id]} ${k_n[$id]} ${t_n[$id]}"

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

else

  #The only that I dont like is that if I want to delete/ºadd one cell I need a lot of work for the others variables...
  #tendre que poner que el string repita si no encuentra xD
  #supongo que habra algun size o asi

  #Arguments are netcdf read parameters: i, j, k, t, i_count, j_count, k_count, t_count
  n_cells=(200) #(100 1125 3375 5625 7875 10800)
  offset_conc=(0) #0.1
  offset_temp=(0) #0.0006
  pmc_multicells=(1) #do multicells? 0=false, 1=true

  echo "Test configuration:" > ../../../../../profile_stats.csv

  #Test with different parameters
  for id_pmc_multicells in "${!pmc_multicells[@]}";
  do
    echo "multi_cells ${pmc_multicells[$id_pmc_multicells]}"
    for id_offset_conc in "${!offset_conc[@]}";
    do
      echo "offset_conc ${offset_conc[$id_offset_conc]}"
      for id_offset_temp in "${!offset_temp[@]}";
      do
        echo "offset_temp ${offset_temp[$id_offset_temp]}"
        for id_n_cells in "${!n_cells[@]}";
        do

          echo "n_cells ${n_cells[id_n_cells]}"

          #todo in the future save all the interesting profile stats in a file apart from the general output
          #{
            #echo "n_cells:${n_cells[id_n_cells]}" >> ../../../../../profile_stats.csv
            #echo "offset_conc:${offset_conc[$id_offset_conc]}"
          #} >> ../../../../../profile_stats.csv

          # Original
          exec_str="../../../test_chemistry_cb05cl_ae5_big ${n_cells[id_n_cells]} \
          ${offset_conc[$id_offset_conc]} ${offset_temp[$id_offset_temp]} ${pmc_multicells[$id_pmc_multicells]}"

          #--print-gpu-summary
          #--analysis-metrics -f -o ../../../../mock_monarch_1000.nvprof
          #exec_str="nvprof --analysis-metrics -f -o ../../../../../test_cb05_10800.nvprof \
          #../../../test_chemistry_cb05cl_ae5_big ${n_cells[id_n_cells]} \
          #${offset_conc[$id_offset_conc]} ${offset_temp[$id_offset_temp]} ${pmc_multicells[$id_pmc_multicells]}"

          if ! $exec_str; then
              echo Failure "$counter"
               if [ "$counter" -gt 1 ]
              then
                echo FAIL
                exit 1
              fi
              ((counter++))
          fi
        done
      done
    done
  done
fi


echo PASS
exit 0
done
