//
//  utils.c
//  
//
//  Created by Matthew Dawson on 10/1/14.
//
//

#include <stdio.h>
#include "proc_fortran_chem.h"


// Compare string case-insensitively
//    (if size = -1 will compare whole strings
//     otherwise, will only compare first 'size' characters)
int strcicmp(char const *a, char const *b, int size)
{
    
    int output_diag = 0;
    
    if (a==NULL || b==NULL) {
        return -1;
    }
    
    if (output_diag) {
        printf("strcicmp - comparing strings: %s AND %s\n", b, b);
    }
    
    for (;; a++, b++, size--) {
        if (size==0)
            return 0;
        int d = tolower(*a) - tolower(*b);
        if (d != 0 || !*a)
            return d;
    }
}

// Determine if a string contains real code to process
int isCode(char const *a)
{
    int i;
    
    if (a==NULL) {
        return 0;
    }
    
    if (strlen(a)<6) {
        return 0;
    }
    
    for (i=0; i<5; i++) {
        if (a[i] == '\t') {
            return 1;
        }
        if (a[i] != ' ') {
            return 0;
        }
    }
    return 1;
}


// trim keyword of non-alphanumeric characters
int keywordTrim(char *a)
{
    do {
        if (a[0]=='\0') {
            return 0;
        }
        if (a[0]>='a' && a[0]<='z') {
            a++;
            continue;
        }
        if (a[0]>='A' && a[0]<='Z') {
            a++;
            continue;
        }
        if (a[0]>='0' && a[0]<='9') {
            a++;
            continue;
        }
        if (a[0]=='_') {
            a++;
            continue;
        }
        a[0] = '\0';
        return 0;
    } while(1);
}



// Reset a rate variable
int resetRate(rate_info_t *currRate      // Rate being reset
              )
{

    (*currRate).intRateConstIndex = -1;
    (*currRate).reactants[0] = NULL;
    (*currRate).reactants[1] = NULL;
    (*currRate).reactants[2] = NULL;
    (*currRate).yield = NULL;
    (*currRate).intNomDenom = EXP_UNKNOWN;

    return 0;
}

// Reset a reaction variable
int resetRxn(rxns_t *currRxn                // Reaction being reset
             )
{
    int i;
    
    (*currRxn).intRateConstIndex = -1;
    (*currRxn).dblConst = 1.0;

    for (i=0; i<3; i++) {
        (*currRxn).reactants[i] = NULL;
        (*currRxn).confReact[i] = 0;
    }
    
    for (i=0; i<20; i++) {
        (*currRxn).products[i] = NULL;
        (*currRxn).yield[i] = NULL;
    }
    
    return 0;
    
}

// Add a reactant to a rate variable in the next available spot
int addReactant(rate_info_t *currRate,      // Rate being modified
                char *strVarName            // Variable name of reactant being added
                )
{
    int i;

    for (i=0; i<3; i++) {
        if ((*currRate).reactants[i]==NULL) {
            (*currRate).reactants[i] = strdup(strVarName);
            return 0;
        }
    }
    
    return 1;
    
}

// Check for a specific reactant in a rate
//  and return it's index if found, otherwise return -1
int isReactant(rate_info_t *currRate,      // Rate being modified
               char *strVarName            // Variable name of reactant being added
               )
{
    int i;
    
    for (i=0; i<3; i++) {
        if ((*currRate).reactants[i]==NULL) {
            continue;
        }
        if (strcicmp((*currRate).reactants[i], strVarName, -1)==0) {
            return i;
        }
    }
    
    return -1;
    
}


// Check if a rate element is empty
int isRateEmpty(rate_info_t *currRate       // Rate being evaluated
)
{
    int retVal = 1;

    if ((*currRate).intRateConstIndex != -1 ||
        (*currRate).reactants[0] != NULL ||
        (*currRate).reactants[1] != NULL ||
        (*currRate).reactants[2] != NULL ||
        (*currRate).yield != NULL ) {
        retVal = 0;
    }
    
    return retVal;
    
}

// Copy a rate variable
int copyRates(rate_info_t *destRate,        // Destination rate
              rate_info_t *srcRate,         // Source rate
              int intReplace                // (0=append; 1=replace; 2=update)
              )
{
    int i;
    
    // Copy rate constant
    if ((*destRate).intRateConstIndex==-1 || intReplace!=0) {
        (*destRate).intRateConstIndex = (*srcRate).intRateConstIndex;
    }

    // Copy reactants
    for (i=0; i<3; i++) {
        if (intReplace==0 && (*srcRate).reactants[i]!=NULL) {
            if (addReactant(destRate, (*srcRate).reactants[i])!=0) {
                printf("copyRates - Error copying reactant.\n");
                return 1;
            }
        }
        if (intReplace==1) {
            (*destRate).reactants[i] = (*srcRate).reactants[i];
        }
        if (intReplace==2) {
            if (isReactant(destRate, (*srcRate).reactants[i])==-1) {
                if (addReactant(destRate, (*srcRate).reactants[i])!=0) {
                    printf("copyRates - Error copying reactant.\n");
                    return 1;
                }
            }
        }
    }

    // Copy yield
    if (((*destRate).yield==NULL || intReplace!=0) && (*srcRate).yield!=NULL) {
        (*destRate).yield = strdup((*srcRate).yield);
    }
    
    // Copy intNomDenom
    if ((*destRate).intNomDenom==-1 || intReplace!=0) {
        (*destRate).intNomDenom = (*srcRate).intNomDenom;
    }

    return 0;
    
}

// Print out a set of reactions
int printRxns(rxns_t *rxns,                 // Set of reactions to print
              int intNumRxns                // Number of reactions to print
              )
{
    int i;
    int nReact, nProd;
    
    printf("\t");
    
    for (i=0; i<intNumRxns; i++) {
        printf("K(%3d)", rxns[i].intRateConstIndex);
    
        // Reactants
        nReact=0;
        printf("\t%s(%d)", rxns[i].reactants[nReact], rxns[i].confReact[nReact]);
        nReact++;
        while (rxns[i].reactants[nReact] != NULL && nReact<3) {
            printf(" + %s(%d)", rxns[i].reactants[nReact], rxns[i].confReact[nReact]);
            nReact++;
        }
    
        // Products
        nProd=0;
        printf (" -->> ");
        if (rxns[i].yield[nProd]!=NULL)
            printf("%s*", rxns[i].yield[nProd]);
        printf("%s", rxns[i].products[nProd++]);
        while (rxns[i].products[nProd] != NULL && nProd<20) {
            printf (" + ");
            if (rxns[i].yield[nProd]!=NULL)
                printf("%s*", rxns[i].yield[nProd]);
            printf("%s", rxns[i].products[nProd++]);
        }
        
        printf("\n");
    }
    
    return 0;
}


// Print out a clump of rates
int printClump(rate_info_t *rateClump,      // Clump of rates
               int intNumRates              // Number of rates to print
               )
{

    int i;
    
    printf("     Rate C.\tReact 1\tReact 2\tReact 3\tYield\tNomDenom\n");
    for (i=0; i<intNumRates; i++) {
        printf("     K(%d)\t%s\t%s\t%s\t%s\t%d\n", rateClump[i].intRateConstIndex,
               rateClump[i].reactants[0], rateClump[i].reactants[1], rateClump[i].reactants[2],
               rateClump[i].yield, rateClump[i].intNomDenom);
    }
    
    return 0;
    
}


// Update/Verify a reaction with rate info
int updateRxnRate(rxns_t *rxn,          // Pointer to reaction
                  rate_info_t *rate,    // Pointer to rate info
                  double dblConst,      // Constant rate coefficient
                  char *strVerify       // Name of reactant species to verify (NULL if none)
                  )
{

    int intIndex;
    
    int output_diag = OUTPUT_DIAG;
    
    if (rxn==NULL || rate==NULL) {
        printf("updateRxnRate - Error: received null pointer\n");
        return 1;
    }

    if (output_diag) {
        printf("\nupdateRxnRate - Updating K(%d)\n", (*rate).intRateConstIndex);
        printf("   Existing Rxn:\n");
        printRxns(rxn, 1);
        printf("   Updating Rate Info:\n");
        printClump(rate, 1);
    }
    
    if ((*rate).intRateConstIndex != (*rxn).intRateConstIndex) {
        printf("updateRxnRate - Error: rate constant mismatch\n");
        return 1;
    }

    if (isRxnEmpty(rxn)) {
        if ((*rate).reactants[0] !=NULL) {
            (*rxn).reactants[0] = strdup((*rate).reactants[0]);
            (*rxn).confReact[0] = 1;
        }
        if ((*rate).reactants[1] !=NULL) {
            (*rxn).reactants[1] = strdup((*rate).reactants[1]);
            (*rxn).confReact[1] = 1;
        }
        if ((*rate).reactants[2] !=NULL) {
            (*rxn).reactants[2] = strdup((*rate).reactants[2]);
            (*rxn).confReact[2] = 1;
        }
        if (dblConst<=0.0000001 && dblConst>=-0.0000001) {
            (*rxn).dblConst = 1.0;
        } else {
            (*rxn).dblConst = dblConst;
        }
    } else {
        if ((*rate).reactants[0] != NULL) {
            if (isRxnReactant(rxn, (*rate).reactants[0])==-1) {
                printf ("updateRxnRate - reactant mismatch: %s\n", (*rate).reactants[0]);
                return 1;
            }
        }
        if ((*rate).reactants[1] != NULL) {
            if (isRxnReactant(rxn, (*rate).reactants[1])==-1) {
                printf ("updateRxnRate - reactant mismatch: %s\n", (*rate).reactants[1]);
                return 1;
            }
        }
        if ((*rate).reactants[2] != NULL) {
            if (isRxnReactant(rxn, (*rate).reactants[2])==-1) {
                printf ("updateRxnRate - reactant mismatch: %s\n", (*rate).reactants[2]);
                return 1;
            }
        }
    }
    
    if (strVerify != NULL) {
        intIndex = isRxnReactant(rxn, strVerify);
        if (intIndex==-1) {
            printf("updateRxnRate - unexpected error verifying reactant %s\n", strVerify);
            return 1;
        } else {
            (*rxn).confReact[intIndex] = 2;
        }
    }

    if (dblConst > 0.0000001 || dblConst < -0.0000001) {
        if ((*rxn).dblConst != 0.0 && (*rxn).dblConst != dblConst) {
            printf("Warning - conflicting rate coefficients: %lf, %lf\n", dblConst, (*rxn).dblConst);
            return 1;
        } else {
            (*rxn).dblConst = dblConst;
        }
    }
    
    return 0;
}

// Determine if a reaction element is empty
int isRxnEmpty(rxns_t *rxn          // Pointer to reaction
               )
{
    int intRetVal = 1;
    
    int output_diag = OUTPUT_DIAG;
    
    if ((*rxn).reactants[0]!=NULL)
        intRetVal = 0;
    if ((*rxn).products[0]!=NULL)
        intRetVal = 0;
    if ((*rxn).yield[0]!=NULL)
        intRetVal = 0;
    
    if (output_diag) {
        if (intRetVal) {
            printf("\nisRxnEmpty - Yes!\n");
        } else {
            printf("\nisRxnEmpty - No!\n");
        }
    }
    
    return intRetVal;
    
}

// Determine is a reaction element contains a specific species
// Returns -1 if not found, otherwise returns index of reactant
int isRxnReactant(rxns_t *rxn,      // Pointer to reaction
                  char *specName    // Variable name for species
                  )
{
    int i = 0;
    
    // Don't match a NULL string
    if (specName==NULL) {
        return -1;
    }

    while ((*rxn).reactants[i]!=NULL && i<3) {
        if (strcicmp((*rxn).reactants[i++], specName, -1)==0) {
            return i-1;
        }
    }
    
    return -1;
}

// Determine if a reaction element contains a specific product
int isRxnProduct(rxns_t *rxn,       // Pointer to reaction
                 char *specName     // Variable name for species
                 )
{
    int i=0;
    
    // Don't match a NULL string
    if (specName==NULL) {
        return -1;
    }

    while ((*rxn).products[i]!=NULL && i<20) {
        if (strcicmp((*rxn).products[i++], specName, -1)==0) {
            return i-1;
        }
    }
    
    return -1;
    
    
}

// Add a product to a reaction
int addRxnProduct(rxns_t *rxn,          // Pointer to reaction
                  rate_info_t *rate,    // Pointer to rate
                  char *specName        // Variable name for species
                  )
{
    int i=0;
    
    while ((*rxn).products[i]!=NULL && i<20) {
        i++;
    }
    
    if (i==20) {
        printf("addRxnProduct - Error: products full: rxn %d\n", (*rxn).intRateConstIndex);
        return 1;
    }

    (*rxn).products[i] = strdup(specName);
    if ((*rate).yield != NULL)
        (*rxn).yield[i] = strdup((*rate).yield);
    
    return 0;
}


// Verify a reaction product has a certain yield
int checkRxnProductYield(rxns_t *rxn,           // Pointer to reaction
                         rate_info_t *rate,     // Pointer to rate
                         char *specName         // Variable name for species
                         )
{
    int i=0;
    int intCompVal;
    
    while (i<20) {
        if (strcicmp((*rxn).products[i], specName, -1)==0) {
            // Found product in rxn!
            
            // Compare the existing and new yields
            if (strcicmp((*rxn).yield[i], (*rate).yield, -1)==0) {
                // They're the same, so continue
                return 0;
            } else {
                // They're different...
                
                // Check if the existing yield is NULL
                if ((*rxn).yield[i] == NULL) {
                    
                    // Check if the new yield is NULL
                    if ((*rate).yield == NULL) {
                        
                        // They must both be 1
                        (*rxn).yield[i] = strdup("2");
                        
                        return 0;
                        
                    } else {
                        
                        // Assume existing yield is 1 and add the new yield
                        char *strTemp;
                        if ((strTemp=(char*)malloc(sizeof(char)*
                                                   (strlen((*rate).yield) +
                                                    2)))==NULL) {
                                                       printf("checkRxnProductYield - Error allocating space for new yield.\n");
                                                       return 1;
                                                   }
                        strcpy(strTemp, "1");
                        if ((*rate).yield[0] != '-')
                            strcat(strTemp, "+");
                        strcat(strTemp, (*rate).yield);
                        (*rxn).yield[i] = strTemp;
                        return 0;
                    }

                } else {
                
                    // Check if new yield is NULL
                    if ((*rate).yield == NULL) {

                        // Assume the new yield is 1 and add it to the existing
                        char *strTemp;
                        if ((strTemp=(char*)malloc(sizeof(char)*
                                                   (strlen((*rxn).yield[i]) +
                                                    2)))==NULL) {
                                                       printf("checkRxnProductYield - Error allocating space for new yield.\n");
                                                       return 1;
                                                   }
                        strcpy(strTemp, (*rxn).yield[i]);
                        strcat(strTemp, "+1");
                        free((*rxn).yield[i]);
                        (*rxn).yield[i] = strTemp;
                        return 0;

                    } else {
                        
                        // They both have values, so add them together
                        char *strTemp;
                        if ((strTemp=(char*)malloc(sizeof(char)*
                                                   (strlen((*rxn).yield[i]) +
                                                    strlen((*rate).yield) +
                                                    1)))==NULL) {
                                                       printf("checkRxnProductYield - Error allocating space for new yield.\n");
                                                       return 1;
                                                   }
                        strcpy(strTemp, (*rxn).yield[i]);
                        if ((*rate).yield[0] != '-')
                            strcat(strTemp, "+");
                        strcat(strTemp, (*rate).yield);
                        free((*rxn).yield[i]);
                        (*rxn).yield[i] = strTemp;
                        return 0;
                    }
                }
            }
        }
        i++;
    }
    
    printf("checkRxnProductYield - Error: cannot find %s in:\n", specName);
    printRxns(rxn,1);
    return 1;
    
}



// Evaluate a comment to see if we are moving to a new section
int evalComment(char *strInput,                 // String containing comment to evaluate
                int *intRoutine,                // Pointer to a variable holding the current routine
                int *intLoc                     // Pointer to a variable holding the current location in the routine
                )
{
    
    // Entering species descriptions
    if (*intRoutine==SR_DIFFR) {
        if (strcicmp(strInput, "C...PROGRAM VARIABLES...", 24)==0) {
            printf("\n*** Entering Program Variable Description section ***\n\n");
            *intLoc = LOC_SPEC_DESC;
        }
    }
    // Entering coefficient desctriptions
    if (*intRoutine==SR_DIFFR) {
        if (strcicmp(strInput, "C...VARIABLE STOICHIOMETRIC COEFFICIENTS...", 43)==0) {
            printf("\n*** Entering Coefficient Description section ***\n\n");
            *intLoc = LOC_COEFF;
        }
    }
    // Leaving steady state species descriptions
    if (*intRoutine==SR_DIFFR && *intLoc==LOC_COEFF) {
        if (strcicmp(strInput, "C**********************************************************************", 70)==0) {
            printf("\n*** Entering PSSA Species Variable section ***\n\n");
            *intLoc = LOC_SPEC_VARS;
        }
    }
    // Entering set concentration area
    if (*intRoutine==SR_DIFFR) {
        if (strcicmp(strInput, "C     SET THE CONSTANT SPECIES", 30)==0) {
            printf("\n*** Entering the Initial Concentration Setting section ***\n\n");
            *intLoc = LOC_SET_CONC;
        }
    }
    // Entering psuedo steady state area
    if (*intRoutine==SR_DIFFR) {
        if (strcicmp(strInput, "C     CALCULATE THE EXPLICIT STEADY STATE EXPRESSIONS", 53)==0) {
            printf("\n*** Entering PSSA Calculation section ***\n\n");
            *intLoc = LOC_PSSA;
        }
    }
    // Entering rate setting area
    if (*intRoutine==SR_DIFFR) {
        if (strcicmp(strInput, "C     CALCULATE FORWARD RATES FOR EACH REACTION", 47)==0) {
            printf("\n*** Entering Rate section ***\n\n");
            *intLoc = LOC_SET_FOR_RATE;
        }
    }
    // Entering rate setting area
    if (*intRoutine==SR_DIFFR) {
        if (strcicmp(strInput, "C     CALCULATE THE DIFFERENTIAL RATES", 38)==0) {
            printf("\n*** Entering Differential Rates section ***\n\n");
            *intLoc = LOC_DIFF_RATE;
        }
    }

    return 0;
    
}


// Transform the model according to this custom transformation code
int transformModel(params_t *params,      // Pointer to model parameters
                   model_t *model         // Pointer to model data
                   )
{
    char **strPSSAOldName;      // Old name of PSSA species
    char **strPSSANewName;      // New name of PSSA species
    int *intPSSANewIdx;         // New index for PSSA species
    int numPSSAChanges = 52;    // Number of PSSA species name/index changes
    int numSpecChanges = 76;    // Number of regular species name changes
    
    int i, j, k;
    
    // *** Allocate space for arrays ***
    
    if((strPSSAOldName = (char**)malloc(sizeof(char*)*(numPSSAChanges+numSpecChanges)))==NULL) {
        printf("\nError allocating space for old PSSA names.\n");
        return 1;
    } else {
        for (i=0; i<(numPSSAChanges+numSpecChanges); i++) {
            if ((strPSSAOldName[i] = (char*)malloc(sizeof(char)*10))==NULL) {
                printf("\nError allocating space for old PSSA names.\n");
                return 1;
            }
        }
    }

    if((strPSSANewName = (char**)malloc(sizeof(char*)*(numPSSAChanges+numSpecChanges)))==NULL) {
        printf("\nError allocating space for new PSSA names.\n");
        return 1;
    } else {
        for (i=0; i<(numPSSAChanges+numSpecChanges); i++) {
            if ((strPSSANewName[i] = (char*)malloc(sizeof(char)*10))==NULL) {
                printf("\nError allocating space for new PSSA names.\n");
                return 1;
            }
        }
    }

    if((intPSSANewIdx = (int*)malloc(sizeof(int)*numPSSAChanges))==NULL) {
        printf("\nError allocating space for new PSSA indices.\n");
        return 1;
    }
    
    // *** Set up changes to PSSA species ***
    
    // RO275 --> RO259
    strcpy(strPSSAOldName[0], "RO259");
    strcpy(strPSSANewName[0], "XO259");
    intPSSANewIdx[0] = 0;
    strcpy(strPSSAOldName[1], "RO275");
    strcpy(strPSSANewName[1], "RO259");
    intPSSANewIdx[1] = 67;

    // RO276 --> RO260
    strcpy(strPSSAOldName[2], "RO260");
    strcpy(strPSSANewName[2], "XO260");
    intPSSANewIdx[2] = 0;
    strcpy(strPSSAOldName[3], "RO276");
    strcpy(strPSSANewName[3], "RO260");
    intPSSANewIdx[3] = 68;

    // RO277 --> RO261
    strcpy(strPSSAOldName[4], "RO261");
    strcpy(strPSSANewName[4], "XO261");
    intPSSANewIdx[4] = 0;
    strcpy(strPSSAOldName[5], "RO277");
    strcpy(strPSSANewName[5], "RO261");
    intPSSANewIdx[5] = 69;

    // RO278 --> RO262
    strcpy(strPSSAOldName[6], "RO262");
    strcpy(strPSSANewName[6], "XO262");
    intPSSANewIdx[6] = 0;
    strcpy(strPSSAOldName[7], "RO278");
    strcpy(strPSSANewName[7], "RO262");
    intPSSANewIdx[7] = 70;

    // RO279 --> RO263
    strcpy(strPSSAOldName[8], "RO263");
    strcpy(strPSSANewName[8], "XO263");
    intPSSANewIdx[8] = 0;
    strcpy(strPSSAOldName[9], "RO279");
    strcpy(strPSSANewName[9], "RO263");
    intPSSANewIdx[9] = 71;

    // RO280 --> RO264
    strcpy(strPSSAOldName[10], "RO264");
    strcpy(strPSSANewName[10], "XO264");
    intPSSANewIdx[10] = 0;
    strcpy(strPSSAOldName[11], "RO280");
    strcpy(strPSSANewName[11], "RO264");
    intPSSANewIdx[11] = 72;

    // RO281 --> RO265
    strcpy(strPSSAOldName[12], "RO265");
    strcpy(strPSSANewName[12], "XO265");
    intPSSANewIdx[12] = 0;
    strcpy(strPSSAOldName[13], "RO281");
    strcpy(strPSSANewName[13], "RO265");
    intPSSANewIdx[13] = 73;

    // RO282 --> RO266
    strcpy(strPSSAOldName[14], "RO266");
    strcpy(strPSSANewName[14], "XO266");
    intPSSANewIdx[14] = 0;
    strcpy(strPSSAOldName[15], "RO282");
    strcpy(strPSSANewName[15], "RO266");
    intPSSANewIdx[15] = 74;

    // RO283 --> RO267
    strcpy(strPSSAOldName[16], "RO267");
    strcpy(strPSSANewName[16], "XO267");
    intPSSANewIdx[16] = 0;
    strcpy(strPSSAOldName[17], "RO283");
    strcpy(strPSSANewName[17], "RO267");
    intPSSANewIdx[17] = 75;

    // RO291 --> RO268
    strcpy(strPSSAOldName[18], "RO268");
    strcpy(strPSSANewName[18], "XO268");
    intPSSANewIdx[18] = 0;
    strcpy(strPSSAOldName[19], "RO291");
    strcpy(strPSSANewName[19], "RO268");
    intPSSANewIdx[19] = 76;

    // RO292 --> RO269
    strcpy(strPSSAOldName[20], "RO269");
    strcpy(strPSSANewName[20], "XO269");
    intPSSANewIdx[20] = 0;
    strcpy(strPSSAOldName[21], "RO292");
    strcpy(strPSSANewName[21], "RO269");
    intPSSANewIdx[21] = 77;
    
    // RO293 --> RO270
    strcpy(strPSSAOldName[22], "RO270");
    strcpy(strPSSANewName[22], "XO270");
    intPSSANewIdx[22] = 0;
    strcpy(strPSSAOldName[23], "RO293");
    strcpy(strPSSANewName[23], "RO270");
    intPSSANewIdx[23] = 78;
    
    // RO294 --> RO271
    strcpy(strPSSAOldName[24], "RO271");
    strcpy(strPSSANewName[24], "XO271");
    intPSSANewIdx[24] = 0;
    strcpy(strPSSAOldName[25], "RO294");
    strcpy(strPSSANewName[25], "RO271");
    intPSSANewIdx[25] = 79;
    
    // RO295 --> RO272
    strcpy(strPSSAOldName[26], "RO272");
    strcpy(strPSSANewName[26], "XO272");
    intPSSANewIdx[26] = 0;
    strcpy(strPSSAOldName[27], "RO295");
    strcpy(strPSSANewName[27], "RO272");
    intPSSANewIdx[27] = 80;
    
    // RO296 --> RO273
    strcpy(strPSSAOldName[28], "RO273");
    strcpy(strPSSANewName[28], "XO273");
    intPSSANewIdx[28] = 0;
    strcpy(strPSSAOldName[29], "RO296");
    strcpy(strPSSANewName[29], "RO273");
    intPSSANewIdx[29] = 81;
    
    // RO297 --> RO274
    strcpy(strPSSAOldName[30], "RO274");
    strcpy(strPSSANewName[30], "XO274");
    intPSSANewIdx[30] = 0;
    strcpy(strPSSAOldName[31], "RO297");
    strcpy(strPSSANewName[31], "RO274");
    intPSSANewIdx[31] = 82;
    
    // RO298 --> RO275
    strcpy(strPSSAOldName[32], "RO275");
    strcpy(strPSSANewName[32], "XO275");
    intPSSANewIdx[32] = 0;
    strcpy(strPSSAOldName[33], "RO298");
    strcpy(strPSSANewName[33], "RO275");
    intPSSANewIdx[33] = 83;
    
    // RO299 --> RO276
    strcpy(strPSSAOldName[34], "RO276");
    strcpy(strPSSANewName[34], "XO276");
    intPSSANewIdx[34] = 0;
    strcpy(strPSSAOldName[35], "RO299");
    strcpy(strPSSANewName[35], "RO276");
    intPSSANewIdx[35] = 84;
    
    i=36;
    // RAD31 --> RAD18
    strcpy(strPSSAOldName[i], "RAD18");
    strcpy(strPSSANewName[i], "XRD18");
    intPSSANewIdx[i++] = 0;
    strcpy(strPSSAOldName[i], "RAD31");
    strcpy(strPSSANewName[i], "RAD18");
    intPSSANewIdx[i++] = 0;

    // RAD32 --> RAD19
    strcpy(strPSSAOldName[i], "RAD19");
    strcpy(strPSSANewName[i], "XRD19");
    intPSSANewIdx[i++] = 0;
    strcpy(strPSSAOldName[i], "RAD32");
    strcpy(strPSSANewName[i], "RAD19");
    intPSSANewIdx[i++] = 0;
    
    // RAD33 --> RAD20
    strcpy(strPSSAOldName[i], "RAD20");
    strcpy(strPSSANewName[i], "XRD20");
    intPSSANewIdx[i++] = 0;
    strcpy(strPSSAOldName[i], "RAD33");
    strcpy(strPSSANewName[i], "RAD20");
    intPSSANewIdx[i++] = 0;
    
    // RAD34 --> RAD21
    strcpy(strPSSAOldName[i], "RAD21");
    strcpy(strPSSANewName[i], "XRD21");
    intPSSANewIdx[i++] = 0;
    strcpy(strPSSAOldName[i], "RAD34");
    strcpy(strPSSANewName[i], "RAD21");
    intPSSANewIdx[i++] = 0;
    
    // RAD35 --> RAD22
    strcpy(strPSSAOldName[i], "RAD22");
    strcpy(strPSSANewName[i], "XRD22");
    intPSSANewIdx[i++] = 0;
    strcpy(strPSSAOldName[i], "RAD35");
    strcpy(strPSSANewName[i], "RAD22");
    intPSSANewIdx[i++] = 0;
    
    // RAD36 --> RAD23
    strcpy(strPSSAOldName[i], "RAD23");
    strcpy(strPSSANewName[i], "XRD23");
    intPSSANewIdx[i++] = 0;
    strcpy(strPSSAOldName[i], "RAD36");
    strcpy(strPSSANewName[i], "RAD23");
    intPSSANewIdx[i++] = 0;
    
    // RAD37 --> RAD24
    strcpy(strPSSAOldName[i], "RAD24");
    strcpy(strPSSANewName[i], "XRD24");
    intPSSANewIdx[i++] = 0;
    strcpy(strPSSAOldName[i], "RAD37");
    strcpy(strPSSANewName[i], "RAD24");
    intPSSANewIdx[i++] = 0;
    
    // RAD38 --> RAD25
    strcpy(strPSSAOldName[i], "RAD25");
    strcpy(strPSSANewName[i], "XRD25");
    intPSSANewIdx[i++] = 0;
    strcpy(strPSSAOldName[i], "RAD38");
    strcpy(strPSSANewName[i], "RAD25");
    intPSSANewIdx[i++] = 0;
    
    // *** Swap out PSSA names and indices ***

    // Change PSSA species
    for (i=0; i<MAX_UNKNOWN_SPEC; i++) {
     
        // Skip empty records
        if ((*model).unknown_spec_var_name[i]==NULL) {
            continue;
        }
        
        // Check each PSSA species for possible update
        for (j=0; j<numPSSAChanges; j++) {
            if (strcicmp((*model).unknown_spec_var_name[i], strPSSAOldName[j], -1)==0) {
                strcpy((*model).unknown_spec_var_name[i], strPSSANewName[j]);
                if ((*model).unknown_spec_name[i] == NULL) {
                    (*model).unknown_spec_name[i] = strdup(strPSSANewName[j]);
                } else {
                    strcpy((*model).unknown_spec_name[i], strPSSANewName[j]);
                }
//                printf("oldname = %s %s new = %s\n", (*model).unknown_spec_var_name[i], (*model).unknown_spec_name[i], strPSSANewName[j]);
                if (intPSSANewIdx[j] != 0) {
                    (*model).unknown_PSSA_idx[i] = intPSSANewIdx[j];
                }
                break;
            }
        }
        
    }

    // Update reactions involving PSSA species
    for (i=0; i<(*params).nrmax; i++) {
        
        // Check reactants first
        for (j=0; j<3; j++) {
            // Skip empty records
            if ((*model).rxns[i].reactants[j]==NULL) {
                continue;
            }
            // Swap out new PSSA names where appropriate
            for (k=0; k<numPSSAChanges; k++) {
                if (strcicmp((*model).rxns[i].reactants[j], strPSSAOldName[k], -1)==0) {
                    strcpy((*model).rxns[i].reactants[j], strPSSANewName[k]);
                    break;
                }
            }
        }

        // ... then check products
        for (j=0; j<20; j++) {
            // Skip empty records
            if ((*model).rxns[i].products[j]==NULL) {
                continue;
            }
            // Swap out new PSSA names where appropriate
            for (k=0; k<numPSSAChanges; k++) {
                if (strcicmp((*model).rxns[i].products[j], strPSSAOldName[k], -1)==0) {
                    strcpy((*model).rxns[i].products[j], strPSSANewName[k]);
                    break;
                }
            }
        }

    
    }
    
    // *** Update Regular Gas-phase Names ***
    
    i=0;
    // RP24 --> RP20
    strcpy(strPSSAOldName[i], "RP20");
    strcpy(strPSSANewName[i++], "XRP20");
    strcpy(strPSSAOldName[i], "RP24");
    strcpy(strPSSANewName[i++], "RP20");

    // RP25 --> RP21
    strcpy(strPSSAOldName[i], "RP21");
    strcpy(strPSSANewName[i++], "XRP21");
    strcpy(strPSSAOldName[i], "RP25");
    strcpy(strPSSANewName[i++], "RP21");

    // RP26 --> RP22
    strcpy(strPSSAOldName[i], "RP22");
    strcpy(strPSSANewName[i++], "XRP22");
    strcpy(strPSSAOldName[i], "RP26");
    strcpy(strPSSANewName[i++], "RP22");

    // RP27 --> RP23
    strcpy(strPSSAOldName[i], "RP23");
    strcpy(strPSSANewName[i++], "XRP23");
    strcpy(strPSSAOldName[i], "RP27");
    strcpy(strPSSANewName[i++], "RP23");

    // RP28 --> RP24
    strcpy(strPSSAOldName[i], "RP28");
    strcpy(strPSSANewName[i++], "RP24");

    // RP29 --> RP25
    strcpy(strPSSAOldName[i], "RP29");
    strcpy(strPSSANewName[i++], "RP25");
    
    // RP30 --> RP26
    strcpy(strPSSAOldName[i], "RP30");
    strcpy(strPSSANewName[i++], "RP26");
    
    // RP31 --> RP27
    strcpy(strPSSAOldName[i], "RP31");
    strcpy(strPSSANewName[i++], "RP27");
    
    // RP32 --> RP28
    strcpy(strPSSAOldName[i], "RP32");
    strcpy(strPSSANewName[i++], "RP28");
    
    // RP33 --> RP29
    strcpy(strPSSAOldName[i], "RP33");
    strcpy(strPSSANewName[i++], "RP29");
    
    // RP34 --> RP30
    strcpy(strPSSAOldName[i], "RP34");
    strcpy(strPSSANewName[i++], "RP30");
    
    // RP35 --> RP31
    strcpy(strPSSAOldName[i], "RP35");
    strcpy(strPSSANewName[i++], "RP31");
    
    // RP36 --> RP32
    strcpy(strPSSAOldName[i], "RP36");
    strcpy(strPSSANewName[i++], "RP32");
    
    // RP37 --> RP33
    strcpy(strPSSAOldName[i], "RP37");
    strcpy(strPSSANewName[i++], "RP33");
    
    // RP38 --> RP34
    strcpy(strPSSAOldName[i], "RP38");
    strcpy(strPSSANewName[i++], "RP34");
    
    // RP39 --> RP35
    strcpy(strPSSAOldName[i], "RP39");
    strcpy(strPSSANewName[i++], "RP35");

    // RP40 --> RP36
    strcpy(strPSSAOldName[i], "RP40");
    strcpy(strPSSANewName[i++], "RP36");

    // RP41 --> RP37
    strcpy(strPSSAOldName[i], "RP41");
    strcpy(strPSSANewName[i++], "RP37");

    // RP91 --> RP38
    strcpy(strPSSAOldName[i], "RP91");
    strcpy(strPSSANewName[i++], "RP38");

    // RP92 --> RP39
    strcpy(strPSSAOldName[i], "RP92");
    strcpy(strPSSANewName[i++], "RP39");
    
    // RP93 --> RP40
    strcpy(strPSSAOldName[i], "RP93");
    strcpy(strPSSANewName[i++], "RP40");
    
    // RP94 --> RP41
    strcpy(strPSSAOldName[i], "RP94");
    strcpy(strPSSANewName[i++], "RP41");
    
    // RP95 --> RP42
    strcpy(strPSSAOldName[i], "RP95");
    strcpy(strPSSANewName[i++], "RP42");
    
    // RP96 --> RP43
    strcpy(strPSSAOldName[i], "RP96");
    strcpy(strPSSANewName[i++], "RP43");
    
    // RP97 --> RP44
    strcpy(strPSSAOldName[i], "RP97");
    strcpy(strPSSANewName[i++], "RP44");
    
    // RP98 --> RP45
    strcpy(strPSSAOldName[i], "RP98");
    strcpy(strPSSANewName[i++], "RP45");
    
    // RP99 --> RP46
    strcpy(strPSSAOldName[i], "RP99");
    strcpy(strPSSANewName[i++], "RP46");
    
    // R100 --> RP47
    strcpy(strPSSAOldName[i], "R100");
    strcpy(strPSSANewName[i++], "RP47");
    
    // R101 --> RP48
    strcpy(strPSSAOldName[i], "R101");
    strcpy(strPSSANewName[i++], "RP48");
    
    // R102 --> RP49
    strcpy(strPSSAOldName[i], "R102");
    strcpy(strPSSANewName[i++], "RP49");
    
    // AP25 --> AP13
    strcpy(strPSSAOldName[i], "AP25");
    strcpy(strPSSANewName[i++], "AP13");
    
    // AP41 --> AP14
    strcpy(strPSSAOldName[i], "AP41");
    strcpy(strPSSANewName[i++], "AP14");
    
    
    
    // UR55 --> UR36
    strcpy(strPSSAOldName[i], "UR36");
    strcpy(strPSSANewName[i++], "XUR36");
    strcpy(strPSSAOldName[i], "UR55");
    strcpy(strPSSANewName[i++], "UR36");

    // UR56 --> UR37
    strcpy(strPSSAOldName[i], "UR37");
    strcpy(strPSSANewName[i++], "XUR37");
    strcpy(strPSSAOldName[i], "UR56");
    strcpy(strPSSANewName[i++], "UR37");

    // UR57 --> UR38
    strcpy(strPSSAOldName[i], "UR38");
    strcpy(strPSSANewName[i++], "XUR38");
    strcpy(strPSSAOldName[i], "UR57");
    strcpy(strPSSANewName[i++], "UR38");
    
    // UR58 --> UR39
    strcpy(strPSSAOldName[i], "UR39");
    strcpy(strPSSANewName[i++], "XUR39");
    strcpy(strPSSAOldName[i], "UR58");
    strcpy(strPSSANewName[i++], "UR39");
    
    // UR59 --> UR40
    strcpy(strPSSAOldName[i], "UR40");
    strcpy(strPSSANewName[i++], "XUR40");
    strcpy(strPSSAOldName[i], "UR59");
    strcpy(strPSSANewName[i++], "UR40");
    
    // UR60 --> UR41
    strcpy(strPSSAOldName[i], "UR41");
    strcpy(strPSSANewName[i++], "XUR41");
    strcpy(strPSSAOldName[i], "UR60");
    strcpy(strPSSANewName[i++], "UR41");
    
    // UR61 --> UR42
    strcpy(strPSSAOldName[i], "UR42");
    strcpy(strPSSANewName[i++], "XUR42");
    strcpy(strPSSAOldName[i], "UR61");
    strcpy(strPSSANewName[i++], "UR42");
    
    // UR62 --> UR43
    strcpy(strPSSAOldName[i], "UR43");
    strcpy(strPSSANewName[i++], "XUR43");
    strcpy(strPSSAOldName[i], "UR62");
    strcpy(strPSSANewName[i++], "UR43");
    
    // UR63 --> UR44
    strcpy(strPSSAOldName[i], "UR44");
    strcpy(strPSSANewName[i++], "XUR44");
    strcpy(strPSSAOldName[i], "UR63");
    strcpy(strPSSANewName[i++], "UR44");
    
    // UR64 --> UR45
    strcpy(strPSSAOldName[i], "UR45");
    strcpy(strPSSANewName[i++], "XUR45");
    strcpy(strPSSAOldName[i], "UR64");
    strcpy(strPSSANewName[i++], "UR45");
    
    // UR65 --> UR46
    strcpy(strPSSAOldName[i], "UR46");
    strcpy(strPSSANewName[i++], "XUR46");
    strcpy(strPSSAOldName[i], "UR65");
    strcpy(strPSSANewName[i++], "UR46");
    
    // UR66 --> UR47
    strcpy(strPSSAOldName[i], "UR47");
    strcpy(strPSSANewName[i++], "XUR47");
    strcpy(strPSSAOldName[i], "UR66");
    strcpy(strPSSANewName[i++], "UR47");
    
    // UR69 --> UR48
    strcpy(strPSSAOldName[i], "UR48");
    strcpy(strPSSANewName[i++], "XUR48");
    strcpy(strPSSAOldName[i], "UR69");
    strcpy(strPSSANewName[i++], "UR48");
    
    // UR71 --> UR49
    strcpy(strPSSAOldName[i], "UR49");
    strcpy(strPSSANewName[i++], "XUR49");
    strcpy(strPSSAOldName[i], "UR71");
    strcpy(strPSSANewName[i++], "UR49");
    
    // UR72 --> UR50
    strcpy(strPSSAOldName[i], "UR50");
    strcpy(strPSSANewName[i++], "XUR50");
    strcpy(strPSSAOldName[i], "UR72");
    strcpy(strPSSANewName[i++], "UR50");
    
    // UR73 --> UR51
    strcpy(strPSSAOldName[i], "UR51");
    strcpy(strPSSANewName[i++], "XUR51");
    strcpy(strPSSAOldName[i], "UR73");
    strcpy(strPSSANewName[i++], "UR51");
    
    // UR74 --> UR52
    strcpy(strPSSAOldName[i], "UR52");
    strcpy(strPSSANewName[i++], "XUR52");
    strcpy(strPSSAOldName[i], "UR74");
    strcpy(strPSSANewName[i++], "UR52");

    // UR75 --> UR53
    strcpy(strPSSAOldName[i], "UR53");
    strcpy(strPSSANewName[i++], "XUR53");
    strcpy(strPSSAOldName[i], "UR75");
    strcpy(strPSSANewName[i++], "UR53");
    
    // UR76 --> UR54
    strcpy(strPSSAOldName[i], "UR54");
    strcpy(strPSSANewName[i++], "XUR54");
    strcpy(strPSSAOldName[i], "UR76");
    strcpy(strPSSANewName[i++], "UR54");
    
    // UR77 --> UR55
    strcpy(strPSSAOldName[i], "UR55");
    strcpy(strPSSANewName[i++], "XUR55");
    strcpy(strPSSAOldName[i], "UR77");
    strcpy(strPSSANewName[i++], "UR55");

    
    // Update reactions involving these species
    for (i=0; i<(*params).nrmax; i++) {
        
        // Check reactants first
        for (j=0; j<3; j++) {
            // Skip empty records
            if ((*model).rxns[i].reactants[j]==NULL) {
                continue;
            }
            // Swap out new species names where appropriate
            for (k=0; k<numSpecChanges; k++) {
                if (strcicmp((*model).rxns[i].reactants[j], strPSSAOldName[k], -1)==0) {
                    strcpy((*model).rxns[i].reactants[j], strPSSANewName[k]);
                    break;
                }
            }
        }
        
        // ... then check products
        for (j=0; j<20; j++) {
            // Skip empty records
            if ((*model).rxns[i].products[j]==NULL) {
                continue;
            }
            // Swap out new species names where appropriate
            for (k=0; k<numSpecChanges; k++) {
                if (strcicmp((*model).rxns[i].products[j], strPSSAOldName[k], -1)==0) {
                    strcpy((*model).rxns[i].products[j], strPSSANewName[k]);
                    break;
                }
            }
        }
        
        
    }
    
    // Move updated species to correct indices for CON() statement
    // (this will mess stuff up if the species that were in these spots
    //   are used in any of the reactions shown)

    // Remove aromatic species from current locations
    (*model).spec_var_name[124] = NULL;
    (*model).spec_name[124] = NULL;
    (*model).spec_var_name[125] = NULL;
    (*model).spec_name[125] = NULL;
    (*model).spec_var_name[126] = NULL;
    (*model).spec_name[126] = NULL;
    (*model).spec_var_name[127] = NULL;
    (*model).spec_name[127] = NULL;
    (*model).spec_var_name[128] = NULL;
    (*model).spec_name[128] = NULL;
    (*model).spec_var_name[129] = NULL;
    (*model).spec_name[129] = NULL;
    (*model).spec_var_name[130] = NULL;
    (*model).spec_name[130] = NULL;
    (*model).spec_var_name[132] = NULL;
    (*model).spec_name[132] = NULL;
    (*model).spec_var_name[133] = NULL;
    (*model).spec_name[133] = NULL;
    (*model).spec_var_name[134] = NULL;
    (*model).spec_name[134] = NULL;

    // Add in species to new locations
    i = 136;
    (*model).spec_var_name[i] = strdup("TOLU");
    (*model).spec_name[i++]   = strdup("TOLU");

    (*model).spec_var_name[i] = strdup("OCRE");
    (*model).spec_name[i++]   = strdup("OCRE");
    
    (*model).spec_var_name[i] = strdup("MCRE");
    (*model).spec_name[i++]   = strdup("MCRE");
    
    (*model).spec_var_name[i] = strdup("PCRE");
    (*model).spec_name[i++]   = strdup("PCRE");
    
    (*model).spec_var_name[i] = strdup("SD1");
    (*model).spec_name[i++]   = strdup("SD1");
    
    (*model).spec_var_name[i] = strdup("MXYL");
    (*model).spec_name[i++]   = strdup("MXYL");
    
    (*model).spec_var_name[i] = strdup("DMP1");
    (*model).spec_name[i++]   = strdup("DMP1");
    
    (*model).spec_var_name[i] = strdup("DMP2");
    (*model).spec_name[i++]   = strdup("DMP2");
    
    (*model).spec_var_name[i] = strdup("PN11");
    (*model).spec_name[i++]   = strdup("PN11");
    
    (*model).spec_var_name[i] = strdup("PN12");
    (*model).spec_name[i++]   = strdup("PN12");
    
    (*model).spec_var_name[i] = strdup("PN13");
    (*model).spec_name[i++]   = strdup("PN13");
    
    (*model).spec_var_name[i] = strdup("PN14");
    (*model).spec_name[i++]   = strdup("PN14");
    
    (*model).spec_var_name[i] = strdup("RP20");
    (*model).spec_name[i++]   = strdup("RP20");
    
    (*model).spec_var_name[i] = strdup("RP21");
    (*model).spec_name[i++]   = strdup("RP21");
    
    (*model).spec_var_name[i] = strdup("RP22");
    (*model).spec_name[i++]   = strdup("RP22");
    
    (*model).spec_var_name[i] = strdup("RP23");
    (*model).spec_name[i++]   = strdup("RP23");
    
    (*model).spec_var_name[i] = strdup("RP24");
    (*model).spec_name[i++]   = strdup("RP24");
    
    (*model).spec_var_name[i] = strdup("RP25");
    (*model).spec_name[i++]   = strdup("RP25");
    
    (*model).spec_var_name[i] = strdup("RP26");
    (*model).spec_name[i++]   = strdup("RP26");

    (*model).spec_var_name[i] = strdup("RP27");
    (*model).spec_name[i++]   = strdup("RP27");
    
    (*model).spec_var_name[i] = strdup("RP28");
    (*model).spec_name[i++]   = strdup("RP28");
    
    (*model).spec_var_name[i] = strdup("RP29");
    (*model).spec_name[i++]   = strdup("RP29");
    
    (*model).spec_var_name[i] = strdup("RP30");
    (*model).spec_name[i++]   = strdup("RP30");
    
    (*model).spec_var_name[i] = strdup("RP31");
    (*model).spec_name[i++]   = strdup("RP31");
    
    (*model).spec_var_name[i] = strdup("RP32");
    (*model).spec_name[i++]   = strdup("RP32");
    
    (*model).spec_var_name[i] = strdup("RP33");
    (*model).spec_name[i++]   = strdup("RP33");
    
    (*model).spec_var_name[i] = strdup("RP34");
    (*model).spec_name[i++]   = strdup("RP34");
    
    (*model).spec_var_name[i] = strdup("RP35");
    (*model).spec_name[i++]   = strdup("RP35");
    
    (*model).spec_var_name[i] = strdup("RP36");
    (*model).spec_name[i++]   = strdup("RP36");
    
    (*model).spec_var_name[i] = strdup("RP37");
    (*model).spec_name[i++]   = strdup("RP37");
    
    (*model).spec_var_name[i] = strdup("RP38");
    (*model).spec_name[i++]   = strdup("RP38");
    
    (*model).spec_var_name[i] = strdup("RP39");
    (*model).spec_name[i++]   = strdup("RP39");
    
    (*model).spec_var_name[i] = strdup("RP40");
    (*model).spec_name[i++]   = strdup("RP40");
    
    (*model).spec_var_name[i] = strdup("RP41");
    (*model).spec_name[i++]   = strdup("RP41");
    
    (*model).spec_var_name[i] = strdup("RP42");
    (*model).spec_name[i++]   = strdup("RP42");
    
    (*model).spec_var_name[i] = strdup("RP43");
    (*model).spec_name[i++]   = strdup("RP43");
    
    (*model).spec_var_name[i] = strdup("RP44");
    (*model).spec_name[i++]   = strdup("RP44");
    
    (*model).spec_var_name[i] = strdup("RP45");
    (*model).spec_name[i++]   = strdup("RP45");

    (*model).spec_var_name[i] = strdup("RP46");
    (*model).spec_name[i++]   = strdup("RP46");
    
    (*model).spec_var_name[i] = strdup("RP47");
    (*model).spec_name[i++]   = strdup("RP47");
    
    (*model).spec_var_name[i] = strdup("RP48");
    (*model).spec_name[i++]   = strdup("RP48");
    
    (*model).spec_var_name[i] = strdup("RP49");
    (*model).spec_name[i++]   = strdup("RP49");
    
    (*model).spec_var_name[i] = strdup("AP13");
    (*model).spec_name[i++]   = strdup("AP13");
    
    (*model).spec_var_name[i] = strdup("AP14");
    (*model).spec_name[i++]   = strdup("AP14");
    
    (*model).spec_var_name[i] = strdup("UR36");
    (*model).spec_name[i++]   = strdup("UR36");
    
    (*model).spec_var_name[i] = strdup("UR37");
    (*model).spec_name[i++]   = strdup("UR37");
    
    (*model).spec_var_name[i] = strdup("UR38");
    (*model).spec_name[i++]   = strdup("UR38");
    
    (*model).spec_var_name[i] = strdup("UR39");
    (*model).spec_name[i++]   = strdup("UR39");
    
    (*model).spec_var_name[i] = strdup("UR40");
    (*model).spec_name[i++]   = strdup("UR40");
    
    (*model).spec_var_name[i] = strdup("UR41");
    (*model).spec_name[i++]   = strdup("UR41");
    
    (*model).spec_var_name[i] = strdup("UR42");
    (*model).spec_name[i++]   = strdup("UR42");
    
    (*model).spec_var_name[i] = strdup("UR43");
    (*model).spec_name[i++]   = strdup("UR43");
    
    (*model).spec_var_name[i] = strdup("UR44");
    (*model).spec_name[i++]   = strdup("UR44");
    
    (*model).spec_var_name[i] = strdup("UR45");
    (*model).spec_name[i++]   = strdup("UR45");
    
    (*model).spec_var_name[i] = strdup("UR46");
    (*model).spec_name[i++]   = strdup("UR46");
    
    (*model).spec_var_name[i] = strdup("UR47");
    (*model).spec_name[i++]   = strdup("UR47");
    
    (*model).spec_var_name[i] = strdup("UR48");
    (*model).spec_name[i++]   = strdup("UR48");
    
    (*model).spec_var_name[i] = strdup("UR49");
    (*model).spec_name[i++]   = strdup("UR49");
    
    (*model).spec_var_name[i] = strdup("UR50");
    (*model).spec_name[i++]   = strdup("UR50");
    
    (*model).spec_var_name[i] = strdup("UR51");
    (*model).spec_name[i++]   = strdup("UR51");
    
    (*model).spec_var_name[i] = strdup("UR52");
    (*model).spec_name[i++]   = strdup("UR52");
    
    (*model).spec_var_name[i] = strdup("UR53");
    (*model).spec_name[i++]   = strdup("UR53");
    
    (*model).spec_var_name[i] = strdup("UR54");
    (*model).spec_name[i++]   = strdup("UR54");
    
    (*model).spec_var_name[i] = strdup("UR55");
    (*model).spec_name[i++]   = strdup("UR55");
    
    for (; i<225; i++) {
        (*model).spec_var_name[i] = strdup("reloc");
        (*model).spec_name[i] = strdup("reloc");
    }
    
    // *** Update Reaction Indices ***

    int addToIdx = 76;
    for (i=400; i<536; i++) {
        switch (i+1) {
            case 454:
            case 472:
            case 473:
            case 481:
            case 488:
            case 514:
                (*model).rxns[i].intRateConstIndex = 0;
                addToIdx-=1;
                break;
                
            default:
                (*model).rxns[i].intRateConstIndex += addToIdx;
                break;
        }
    }
    
    return 0;
}








