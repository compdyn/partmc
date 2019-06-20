#!/bin/bash

increment=5
executions=2

for (( i=0; i<$executions; i++ ))
do

    cat cb05cl_ae5_mechanism_header.json > tmp.json

    repeat=$(($increment ** $i))
    for (( j=0; j<$repeat; j++ ))
    do
        cat cb05cl_ae5_mechanism_template.json >> tmp.json #body
        #cat cb05cl_ae5_mechanism_template.json >> tmp.json #body
    done

    sed -i '$ s/.$//' tmp.json #delete last character
    cat cb05cl_ae5_mechanism_footer.json >> tmp.json

    cp tmp.json ../../../mechanisms/cb05cl_ae5/cb05cl_ae5_mechanism_big.json

    #../../../build/test_run/chemistry/cb05cl_ae5/test_chemistry_cb05cl_ae5_big.sh
    cd ../../../../
    ./compile.partmc.marenostrum4.sh
    cd partmc/test/chemistry/cb05cl_ae5_scalability

done

#TTODO: For some reason, some rxn sizes (repetitions of mechanism data) makes the program crush
#Example: copying mechanism data 10 times don'work, but 28 times and 3 times works.
#Tested also copying the data manually, both gives segmentation fault

rm tmp.json

