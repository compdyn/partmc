#!/bin/bash
SPEC_FILE_ORIG=$1
SPEC_FILE=$2
NL='
'
echo "Processing KPP species..."
echo "this is still a work in progress..."

# copy to new file
cp $SPEC_FILE_ORIG $SPEC_FILE 

# reformat comment lines
sed -i '' -E "s/^#[[:space:]]*(.*)$/\/\/ \1;/g" $SPEC_FILE
sed -i '' -E "s/;[[:space:]]*([^[:space:]]+.*)$/;\\$NL\/\/ \1;/g" $SPEC_FILE

# remove blank lines
sed -i '' '/^[[:space:]]*$/ d' $SPEC_FILE

# remove excess whitespace and linebreaks;
sed -i '' -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' $SPEC_FILE
sed -i '' -E "s/;[[:space:]]*/;\\$NL/g" $SPEC_FILE
sed -i '' -E 's/[[:space:]]+/ /g' $SPEC_FILE
sed -i '' -E 's/^[[:space:]]+//g' $SPEC_FILE

# find species and convert to json
sed -i '' -E "s/^[[:space:]]*([[:alpha:]][[:alnum:]]*)[[:space:]]*\=.*$/  {\\$NL    \"name\" : \"\1\",\\$NL    \"type\" : \"CHEM_SPEC\"\\$NL  },/g" $SPEC_FILE

# header
sed -i '' -E "1s/^/{ \"pmc-data\" : [\\$NL/g" $SPEC_FILE

# footer
sed -i '' -E "\$a\\$NL]}" $SPEC_FILE

# remove final comma at end of species list
sed -i '' -e ':a' -e 'N' -e '$!ba' -e "s/},[[:space:]]*\n[[:space:]]]/}\\$NL]/g" $SPEC_FILE

# remove comments
sed -i '' '/\/\// d' $SPEC_FILE
