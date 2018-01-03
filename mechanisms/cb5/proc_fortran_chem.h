//
//  proc_fortran_chem.h
//  
//
//  Created by Matthew Dawson on 9/25/14.
//
//

#ifndef ____proc_fortran_chem__
#define ____proc_fortran_chem__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// ****************************
// *** Constant Definitions ***
// ****************************

#define OUTPUT_MECH         1     // Print out mechanism?
#define OUTPUT_DIAG         0     // Print out diagnostic info?

#define MAX_UNKNOWN_SPEC 5000     // maximum number of unknown species
#define MAX_AR_BR         500     // maximum number of AR() or BR() elements
#define MAX_LINE        10000     // maximum length of a line of code or math expression

// Constants for setting current subroutine
#define SR_NONE    0
#define SR_DOCUMS  1
#define SR_DIFFR   2
#define SR_RATE    3
#define SR_RURALSP 4
#define SR_PHORAT  5
#define SR_CHMFLX  6
#define SR_SETPAN  7

// Constants for describing location within subroutine
#define LOC_NONE            0
#define LOC_SPEC_DESC       1
#define LOC_COEFF           2
#define LOC_SPEC_VARS       3
#define LOC_SET_CONC        4
#define LOC_PSSA            5
#define LOC_SET_FOR_RATE    6
#define LOC_DIFF_RATE       7

// Variables declared with DATA or REAL statements
#define DATA_SNAMS     1
#define DATA_AERONAME  2
#define DATA_SMWS      3
#define DATA_CNAMS     4
#define REAL_SPEC_VARS 5

// Constants describing chemical species type
#define SPEC_GAS     1
#define SPEC_AEROSOL 2

// Constants describing math expressions
#define EXP_UNKNOWN     0
#define EXP_VARIABLE    1
#define EXP_ARRAY       2
#define EXP_DIVIDE      3
#define EXP_MULTIPLY    4
#define EXP_ADD         5
#define EXP_SUBTRACT    6
#define EXP_PAREN       7
#define EXP_EQUAL       8
#define EXP_NUMBER      9
#define EXP_NOM         10
#define EXP_DENOM       11
#define EXP_END         12


// *****************************
// *** Structure Definitions ***
// *****************************

// model parameters and output options
typedef struct  {
    int nasect;     // Number of aerosol size bins
    int ngspec;     // Number of gas phase species
    int naspec;     // Number of aerosol species
    int nrmax;      // Number of gas-phase reactions
    int ncmax;      // Number of special coefficients
    int minRxn;     // Minimum reaction to include in A(), B(), PSSA output or -1 for all rxns
    int maxRxn;     // Maximum reaction to include in A(), B(), PSSA output
    int offsetRxn;      // Offset reaction indices starting with this index
    int offsetRxnBy;    // Amount to offset reactions by
    int offsetSpec;     // Offset species indices starting with this index
    int offsetSpecBy;   // Amount to offset species by
    int printCode;      // (1 = print generated code; 0 = only print mechanism)
    int printDetail;    // (2 = print simplified mechanism; 1 = print detailed species and array info; 0 = don't)
    int noMW;           // (1 = exclude MW info from all printouts; 0 = don't)
    int runTransform;   // (1 = run the custom transformation; 0 = don't)
} params_t;

// reaction data
typedef struct {
    int intRateConstIndex;      // Index for reaction rate constant in K()
    double dblConst;            // Constant rate coefficient
    char *reactants[3];         // Name of each reactant
    int confReact[3];           // Confirmation for each reactant (0=no reactant; 1=unconfirmed; 2=confirmed)
    char *products[20];         // Name of each product
    char *yield[20];            // Yield of each product
    char *rxnEqn;               // String containing chemical reaction
} rxns_t;

// reaction rate info
typedef struct {
    int intRateConstIndex;      // Index of reaction rate constant in K()
    char *reactants[3];         // Name of each reactant
    char *yield;                // Product yield from this reaction
    int intNomDenom;            // (1=nominator; 2=denominator)
} rate_info_t;

// model data
typedef struct {
    char **spec_name;           // Name of each chemical species (gas & aerosol)
    int *spec_type;             // Type of species (gas or aerosol phase)
    double *spec_MW;            // Molecular weight of each species
    char **spec_desc;           // Text description of variable
    char **spec_var_name;       // Variable name of each chemical species (gas & aerosol)
    char **spec_aer_name;       // Name of aerosol species
    char **spec_init_conc;      // Initial concentration set
    char **coeff_name;          // Name of the special coefficients
    char **unknown_spec_name;   // Names of unknown species
    char **unknown_spec_desc;   // Descriptions of unknown species
    char **unknown_spec_var_name;       // Names of unknown species variables
    char **unknown_spec_init_conc;      // Initial concentration of an unknown species set
    int unknown_PSSA_idx[MAX_UNKNOWN_SPEC];     // Index of PSSA species in AR() and BR()
    char *unknown_PSSA_str[MAX_UNKNOWN_SPEC];   // PSSA expression for unknown species
    double unknown_const_conc[MAX_UNKNOWN_SPEC];// Constant concentration set for unknown species
    int unknown_AR[MAX_UNKNOWN_SPEC];   // Corresponding index for this species in AR
    int unknown_BR[MAX_UNKNOWN_SPEC];   // Corresponding index for this species in BR
    int nUnknownSpec;           // Number of unknown species
    char *AR[MAX_AR_BR];        // AR array
    char *BR[MAX_AR_BR];        // BR array
    int nARBR;                  // Highest index used in AR() or BR() arrays
    rxns_t *rxns;               // Set of reactions
    char **strForRates;         // Forward Rates
    char **strRateConst;        // Rate constant expressions
} model_t;


// *****************************
// *** Function Declarations ***
// *****************************

// *** utils.c ***

// Compare string case-insensitively
//    (if size = -1 will compare whole strings
//     otherwise, will only compare first 'size' characters)
int strcicmp(char const *a, char const *b, int size);

// Determine if a string contains real code to process
int isCode(char const *a);

// trim keyword of non-alphanumeric characters
int keywordTrim(char *a);

// Reset a rate variable
int resetRate(rate_info_t *currRate      // Rate being reset
);

// Reset a reaction variable
int resetRxn(rxns_t *currRxn                // Reaction being reset
);

// Add a reactant to a rate variable in the next available spot
int addReactant(rate_info_t *currRate,      // Rate being modified
                char *strVarName            // Variable name of reactant being added
);

// Check for a specific reactant in a rate
int isReactant(rate_info_t *currRate,      // Rate being modified
               char *strVarName            // Variable name of reactant being added
);

// Check if a rate element is empty
int isRateEmpty(rate_info_t *currRate       // Rate being evaluated
);

// Copy a rate variable
int copyRates(rate_info_t *destRate,        // Destination rate
              rate_info_t *srcRate,         // Source rate
              int intReplace                // (0=append; 1=replace; 2=update)
);

// Print out a set of reactions
int printRxns(rxns_t *rxns,                 // Set of reactions to print
              int intNumRxns                // Number of reactions to print
);

// Print out a clump of rates
int printClump(rate_info_t *rateClump,      // Clump of rates
               int intNumRates              // Number of rates to print
);

// Update/Verify a reaction with rate info
int updateRxnRate(rxns_t *rxn,          // Pointer to reaction
                  rate_info_t *rate,    // Pointer to rate info
                  double dblConst,      // Constant rate coefficient
                  char *strVerify       // Name of reactant species to verify (NULL if none)
);

// Determine if a reaction element is empty
int isRxnEmpty(rxns_t *rxn          // Pointer to reaction
);

// Determine is a reaction element contains a specific species
int isRxnReactant(rxns_t *rxn,      // Pointer to reaction
                  char *specName    // Variable name for species
);

// Determine if a reaction element contains a specific product
int isRxnProduct(rxns_t *rxn,       // Pointer to reaction
                 char *specName     // Variable name for species
);

// Add a product to a reaction
int addRxnProduct(rxns_t *rxn,          // Pointer to reaction
                  rate_info_t *rate,    // Pointer to rate
                  char *specName        // Variable name for species
);

// Verify a reaction product has a certain yield
int checkRxnProductYield(rxns_t *rxn,           // Pointer to reaction
                         rate_info_t *rate,     // Pointer to rate
                         char *specName         // Variable name for species
);

// Evaluate a comment to see if we are moving to a new section
int evalComment(char *strInput,                 // String containing comment to evaluate
                int *intRoutine,                // Pointer to a variable holding the current routine
                int *intLoc                     // Pointer to a variable holding the current location in the routine
);

// Transform the model according to this custom transformation code
int transformModel(params_t *params,      // Pointer to model parameters
                   model_t *model         // Pointer to model data
);

// DIAGNOSTIC
int checkDummy(rxns_t *rxns,
               int intNumRxns
               );



#endif /* defined(____proc_fortran_chem__) */
