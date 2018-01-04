#!/bin/bash
MECH_FILE_ORIG=$1
MECH_FILE=$2
NL='
'
echo "Processing CMAQ mechanism..."
echo "this is still a work in progress..."

# copy to new file
cp $MECH_FILE_ORIG $MECH_FILE 

# remove comment lines
sed -i '' '/^!/ d' $MECH_FILE

# remove blank lines
sed -i '' '/^[[:space:]]*$/ d' $MECH_FILE

# remove special lines
sed -i '' '/^REACTIONS/ d' $MECH_FILE

# remove excess whitespace and linebreaks;
sed -i '' -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' $MECH_FILE
sed -i '' -E "s/;[[:space:]]*/;\\$NL/g" $MECH_FILE
sed -i '' -E 's/[[:space:]]+/ /g' $MECH_FILE
sed -i '' -E 's/^[[:space:]]+//g' $MECH_FILE

# id rate params and end of reactions
sed -i '' -E "s/#([^;]*);/\\$NL<RATE_PARAM>\1<\/RATE_PARAM><\/RXN>/g" $MECH_FILE

# id beginning of reactions, reactants, and products
sed -i '' -E "s/<[[:space:]]*([A-Z]+[0-9]+)[[:space:]]*>[[:space:]]*([^\=]+)\=[[:space:]]*(.*)$/<RXN\1><REACT>\2<\/REACT><PROD>\3<\/PROD>/g" $MECH_FILE

# split reactants
for run in {1..20}
do
  sed -i '' -E "s/(<REACT>[^\+]*)\+[[:space:]]*(.+<\/REACT>)/\1<REACT_SPLIT>\2/g" $MECH_FILE
done

# split products
for run in {1..20}
do
  sed -i '' -E "s/(<PROD>[^\+]*)\+[[:space:]]*(.+<\/PROD>)/\1<PROD_SPLIT>\2/g" $MECH_FILE
  sed -i '' -E "s/(<PROD>[^\-]*)(\-.+<\/PROD>)/\1<PROD_SPLIT>\2/g" $MECH_FILE
done

# get yields
for run in {1..20}
do
  sed -i '' -E "s/(<PROD[^\*]*>)[[:space:]]*([0-9\.\-]+)\*([[:alnum:]]+)(.*<\/PROD>)/\1\3<YIELD\2>\4/g" $MECH_FILE
done

###################
# create the json #
###################

# header
sed -i '' -E "1s/^[[:space:]]*([^[:space:]]+)/{ \"pmc-data\" : [\\$NL  {\\$NL    \"name\" : \"\1\",\\$NL    \"type\" : \"MECHANISM\",\\$NL    \"reactions\" : [\\$NL/" $MECH_FILE

# footer
sed -i '' -E "\$a\\$NL  }\\$NL]}" $MECH_FILE

# reaction wrapper
sed -i '' -E "s/<RXN([A-Z0-9]+)>/      {\\$NL        \"rxn id\" : \"\1\",\\$NL/g" $MECH_FILE
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

