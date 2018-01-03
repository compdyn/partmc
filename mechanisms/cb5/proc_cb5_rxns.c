//
//
//  proc_cb5_rxns.c
//
//  Created by Matt Dawson 1/3/18
//
//

#include "proc_fortran_chem.h"

#define MAX_STRING 500

int main (int argc, const char *argv[])
{

	FILE *fProdLoss;	// cb5 production and loss file
	char strInput[MAX_STRING];	// line being read in
	char strTemp[MAX_STRING];
	char *strPtr;
	int prodOrLoss = 0;	// 1 for prod 2 for loss
        char strSpecies[20] = "";	// Species name
	int num_rxn = 187;
        float fltTemp;
	int intTemp;

	model_t model;		// model data

	if ((fProdLoss=fopen("ext-hrprodloss.F90", "r")) == NULL) {
		printf("Could not open file");
		return 1;
	}

	model.rxns = (rxns_t*) malloc(sizeof(rxns_t)*num_rxn);

	while (!feof(fProdLoss)) {

		fgets(strInput, MAX_STRING, fProdLoss);
		printf("%s", strInput);

		// Make sure something is there
		if (strlen(strInput) > 5) {
			
			// don't process comments
			if (isCode(strInput)==1) {

				// if this is not a continuation line, reset the species and prodOrLoss
				if (strInput[5]==' ') {
					
					// get word in first position
					strTemp[0] = '\0';
					sscanf(strInput, "      %s", strTemp);

					// trim keyword of non-alphanumeric characters
 					keywordTrim(strTemp);

					if (strcicmp(strTemp, "PROD", -1)==0) {
						prodOrLoss = 1;
						strPtr = strstr(strInput, "(");
						strPtr++;
						sscanf(strPtr, "%s", strSpecies);
					} else if (strcicmp(strTemp, "LOSS", -1)==0) {
						prodOrLoss = 2;
						strPtr = strstr(strInput, "(");
						strPtr++;
						sscanf(strPtr, "%s", strSpecies);
					} else if (prodOrLoss == 0) {
						continue;
					} else {
						break;
					}
					printf("\nReading species '%s' type %d\n", strSpecies, prodOrLoss);
				}

				// check for a yield or qty
				if (strInput[26] != ' ') {
					strPtr = &strInput[26];
					sscanf(strPtr, "%f", &fltTemp);
				} else {
					fltTemp = 1.0;
				}

				// get the reaction number
				strPtr = &strInput[40];
				if (sscanf(strPtr, "%d", &intTemp)==0) continue;

 				printf("\nrxn %d\n", intTemp);

				if (prodOrLoss == 1) {
					rate_info_t rate;
					rate.yield = (char*) malloc(sizeof(char*)*20);
					sprintf(rate.yield, "%f", fltTemp);
					addRxnProduct(&(model.rxns[intTemp-1]), &rate, strSpecies);
				} else {
					int i_react = 0;
					for(;model.rxns[intTemp-1].reactants[i_react]!=NULL;i_react++);
					printf("\nAdding spec '%s' %d to rxn %d\n", strSpecies, i_react, intTemp);
					model.rxns[intTemp-1].reactants[i_react] = (char*) malloc(sizeof(char)*20);
					strcpy(model.rxns[intTemp-1].reactants[i_react++], strSpecies);
					if (fltTemp == 2.0) {
						model.rxns[intTemp-1].reactants[i_react] = (char*) malloc(sizeof(char)*20);
						strcpy(model.rxns[intTemp-1].reactants[i_react], strSpecies);
					}
				}
			}
		}
	}

	printRxns(model.rxns, num_rxn);

}
