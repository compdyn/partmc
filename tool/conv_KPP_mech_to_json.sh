#!/bin/bash
MECH_FILE_ORIG=$1
MECH_FILE=$2
NL='
'
echo "Processing KPP mechanism..."
echo "this is still a work in progress..."

# copy to new file
cp $MECH_FILE_ORIG $MECH_FILE 

# reformat comment lines
sed -i '' -E "s/^#[[:space:]]*(.*)$/\/\/ \1;/g" $MECH_FILE
sed -i '' -E "s/;[[:space:]]*([^[:space:]]+.*)$/;\\$NL\/\/ \1;/g" $MECH_FILE

# remove blank lines
sed -i '' '/^[[:space:]]*$/ d' $MECH_FILE

# remove excess whitespace and linebreaks;
sed -i '' -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' $MECH_FILE
sed -i '' -E "s/;[[:space:]]*/;\\$NL/g" $MECH_FILE
sed -i '' -E 's/[[:space:]]+/ /g' $MECH_FILE
sed -i '' -E 's/^[[:space:]]+//g' $MECH_FILE

# id rate params and end of reactions
sed -i '' -E "s/:([^;]*);/\\$NL<RATE_PARAM>\1<\/RATE_PARAM><\/RXN>/g" $MECH_FILE

# id beginning of reactions, reactants, and products
sed -i '' -E "s/^{[[:space:]]*([a-zA-Z]+[0-9]+)[[:space:]]*}[[:space:]]*([^\=]+)\=[[:space:]]*(.*)$/<RXN\1><REACT>\2<\/REACT><PROD>\3<\/PROD>/g" $MECH_FILE

# remove brackets surrounding constant reactants
for run in {1..20}
do
  sed -i '' -E "s/(<REACT>[^{}]*)[{}]+(.*<\/REACT>)/\1\2/g" $MECH_FILE
  sed -i '' -E "s/(<PROD>[^{}]*)[{}]+(.*<\/PROD>)/\1\2/g" $MECH_FILE
done

# split reactants
for run in {1..20}
do
  sed -i '' -E "s/(<REACT>[^\+]*)\+[[:space:]]*(.+<\/REACT>)/\1<REACT_SPLIT>\2/g" $MECH_FILE
done

# split products
for run in {1..20}
do
  sed -i '' -E "s/(<PROD>[^\+]*)\+[[:space:]]*(.+<\/PROD>)/\1<PROD_SPLIT>\2/g" $MECH_FILE
  sed -i '' -E "s/(<PROD>[^\-]*)\-(.+<\/PROD>)/\1<PROD_SPLIT><MINUS>\2/g" $MECH_FILE
done
sed -i '' -E "s/<MINUS>[[:space:]]*([[:alpha:]])/\-1.0\1/g" $MECH_FILE
sed -i '' -E "s/<MINUS>[[:space:]]*/\-/g" $MECH_FILE

# get yields
for run in {1..20}
do
  sed -i '' -E "s/(<PROD[^\*]*>)[[:space:]]*([0-9\.\-]+)\*[[:space:]]*([[:alnum:]]+)(.*<\/PROD>)/\1\3<YIELD\2>\4/g" $MECH_FILE
done

###################
# create the json #
###################

# header
sed -i '' -E "1s/^/{ \"pmc-data\" : [\\$NL  {\\$NL    \"name\" : \"\\$MECH_FILE\",\\$NL    \"type\" : \"MECHANISM\",\\$NL    \"reactions\" : [\\$NL/" $MECH_FILE

# footer
sed -i '' -E "\$a\\$NL    ]\\$NL  }\\$NL]}" $MECH_FILE

# reaction wrapper
sed -i '' -E "s/<RXN([a-zA-Z0-9]+)>/      {\\$NL        \"rxn id\" : \"\1\",\\$NL/g" $MECH_FILE
sed -i '' -E "s/<\/RXN>/      },/g" $MECH_FILE

# reactants
sed -i '' -E "s/<REACT>[[:space:]]*([^<[:space:]]+)/        \"reactants\" : {\\$NL          \"\1\" : {}/g" $MECH_FILE
sed -i '' -E "s/<REACT_SPLIT>[[:space:]]*([^<[:space:]]+)/,\\$NL          \"\1\" : {}/g" $MECH_FILE
sed -i '' -E "s/<\/REACT>/\\$NL        },\\$NL/g" $MECH_FILE

# products
sed -i '' -E "s/<PROD>[[:space:]]*([^<[:space:]]+)/        \"products\" : {\\$NL          \"\1\" : {}/g" $MECH_FILE
sed -i '' -E "s/<PROD>/        \"products\" : {/g" $MECH_FILE
sed -i '' -E "s/<PROD_SPLIT>[[:space:]]*([^<[:space:]]+)/,\\$NL          \"\1\" : {}/g" $MECH_FILE
sed -i '' -E "s/{}<YIELD([^>]+)>/{ \"yield\" : \1 }/g" $MECH_FILE
sed -i '' -E "s/<\/PROD>/\\$NL        },/g" $MECH_FILE

# rate params
sed -i '' -E "s/<RATE_PARAM>[[:space:]]*/        \"orig params\" : \"/g" $MECH_FILE
sed -i '' -E "s/<\/RATE_PARAM>/\"\\$NL/g" $MECH_FILE

# remove final comma at end of reaction list
sed -i '' -e ':a' -e 'N' -e '$!ba' -e "s/},[[:space:]]*\n[[:space:]]\]/}\\$NL    ]/g" $MECH_FILE


###############################
# process reaction parameters #
###############################

# prep parameters
sed -i '' -E "s/(\"orig params\" : \")([^\"]+)\"$/\1\2\",\\$NL<RC>\2/g" $MECH_FILE

# remove multiplication by contant species
for run in {1..20}
do
  sed -i '' -E "s/^<RC>[[:space:]]*[a-zA-Z][[:alnum:]]*[[:space:]]*[\*]+[0-9]*[[:space:]]*[\*]*/<RC>/g" $MECH_FILE
done

# photolysis
sed -i '' -E "s/^<RC>.*TUV_J.*$/        \"type\" : \"PHOTOLYSIS\"/g" $MECH_FILE

ARG_MATCH='[[:space:]]*([^,\)[:space:]]+)[[:space:]]*'

# CMAQ_1to4 (Arrhenius)
# sed -i '' -E "s/^<RC>[[:space:]]*CMAQ_1to4\([[:space:]]*([^,[:space:]]+)[[:space:]]*,[[:space:]]*([^,[:space:]]+)[[:space:]]*,[[:space:]]*([^\)[:space:]]+)[[:space:]]*\).*$/        \"type\" : \"ARRHENIUS\",\\$NL        \"A\" : \1,\\$NL        \"B\" : \2,\\$NL        \"C\" : -\3/g" $MECH_FILE
sed -i '' -E "s/^<RC>[[:space:]]*CMAQ_1to4\($ARG_MATCH,$ARG_MATCH,$ARG_MATCH\).*$/        \"type\" : \"ARRHENIUS\",\\$NL        \"A\" : \1,\\$NL        \"B\" : \2,\\$NL        \"C\" : -\3/g" $MECH_FILE

# CMAQ_8 (CMAQ_OH_HNO3)
sed -i '' -E "s/^<RC>[[:space:]]*CMAQ_8\($ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH\).*$/        \"type\" : \"CMAQ_OH_HNO3\",\\$NL        \"k0_A\" : \1,\\$NL        \"k0_C\" : \2,\\$NL        \"k2_A\" : -\3,\\$NL        \"k2_C\" : \4,\\$NL        \"k3_A\" : \5,\\$NL        \"k3_C\" : \6/g" $MECH_FILE

# CMAQ_9 (CMAQ_H2O2)
sed -i '' -E "s/^<RC>[[:space:]]*CMAQ_9\($ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH\).*$/        \"type\" : \"CMAQ_H2O2\",\\$NL        \"k1_A\" : \1,\\$NL        \"k1_C\" : \2,\\$NL        \"k2_A\" : -\3,\\$NL        \"k2_C\" : \4/g" $MECH_FILE

# CMAQ_10 (Troe)
sed -i '' -E "s/^<RC>[[:space:]]*CMAQ_10\($ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH,$ARG_MATCH\).*$/        \"type\" : \"TROE\",\\$NL        \"k0_A\" : \1,\\$NL        \"k0_B\" : \2,\\$NL        \"k0_C\" : -\3,\\$NL        \"kinf_A\" : \4,\\$NL        \"kinf_B\" : \5,\\$NL        \"kinf_C\" : -\6,\\$NL        \"Fc\" : \7,\\$NL        \"N\" : \8/g" $MECH_FILE

# Zeroed reactions
sed -i '' -E "s/^<RC>[[:space:]]*0[[:space:]]*$/        \"type\" : \"ARRHENIUS\",\\$NL        \"preproc flag\" : \"REMOVE\"/g" $MECH_FILE

############
# clean up #
############

# remove double negatives
sed -i '' -E "s/\-\-//g" $MECH_FILE

# change all exponent flags to E
sed -i '' -E "s/([0-9]+\.[0-9]*)[D]([\+\-][0-9]+)/\1E\2/g" $MECH_FILE

# add trailing zeros
sed -i '' -E "s/([0-9]+\.)([^0-9]+)/\10\2/g" $MECH_FILE
sed -i '' -E "s/([0-9]+\.)$/\10/g" $MECH_FILE

# remove comments
sed -i '' '/\/\// d' $MECH_FILE
