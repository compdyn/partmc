! Copyright (C) 2019 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \page camp_tutorial Boot CAMP: The CAMP Tutorial
!!
!! In the Boot CAMP tutorial, we build a simple box model that
!! incorporates all the CAMP functionality to demonstrate how CAMP can be
!! included in your favorite model.
!!
!! - \ref camp_tutorial_part_1
!! - \ref camp_tutorial_part_2
!! - \ref camp_tutorial_part_3
!! - \ref camp_tutorial_part_4
!!
!! Model code described in Boot CAMP can be found in
!! \c doc/camp_tutorial .

!> \page camp_tutorial_part_1 Boot CAMP: Part 1 - Box Model
!!
!! Prior to beginning this tutorial, the PartMC library should be
!! installed on your system with PartMC-CAMP enabled. Installation
!! instructions can be found \ref camp_chem "here".
!!
!! The purpose of this tutorial is to demonstrate how to integrate CAMP
!! into a new host model, and use it to build and solve a unique
!! chemical system. Over the the course of this tutorial,
!! we will build a simple box model, add CAMP for the chemistry solving,
!!  and develop
!! a simple chemical mechanism that makes use of all the CAMP
!! functionality. Key features of CAMP will be demonstrated, including
!! it's use with various aerosol representations (e.g., bins, modes,
!! single particles), it's combined gas- and aerosol-phase solving, and
!! it's run-time chemical mechanism configuration.
!!
!! The box-model code in this tutorial is designed to make
!! clear how to interact with CAMP, but it is not meant to represent
!! a well-designed model.
!! The code used in the tutorial can be found in `doc/camp_tutorial`.
!!
!! ## CAMP box model ##
!!
!! First, let's build our simple box model. In a file named
!! `box_model.F90`, add the following code to load the CAMP modules
!! we'll need to get started:
!!
!! \snippet camp_tutorial/part_1_code/box_model.F90 Modules to use
!!
!! These modules provide two derived types that are needed in any CAMP
!! implementation, \ref pmc_camp_core::camp_core_t "camp_core_t" and
!! \ref pmc_camp_state::camp_state_t "camp_state_t". CAMP employs an
!! object-oriented design in which all of its functionality is exposed
!! through instances of derived types. Most CAMP modules include one
!! public derived type, with some exceptions we'll see later on.
!!
!! Next, declare the needed variables and create and initialize a
!! CAMP core:
!!
!! \snippet camp_tutorial/part_1_code/box_model.F90 Core and state
!! \snippet camp_tutorial/part_1_code/box_model.F90 Initialize core
!!
!! The \ref pmc_camp_core::camp_core_t "camp_core_t()"  constructor
!! takes one argument: the path to a
!! configuration file for the chemical mechanism you would like to run.
!! (We'll create this a little later in the tutorial.) The
!! constructor reads the input data and creates internal objects to
!! describe the system (reactions, species, etc.). The
!! \ref pmc_camp_core::initialize "initialize()"
!! function instructs these model elements to validate their input data
!! and assemble the information they will need during solving.
!!
!! In the first implementation of our box model, we will assume a fixed
!! mechanism, so we will hard-code some species names.
!! Later on we will use some
!! \ref pmc_camp_core::camp_core_t "camp_core_t" functions to
!! retrieve the species present in the model at run time to avoid this.
!!
!! \snippet camp_tutorial/part_1_code/box_model.F90 Chem-spec data module
!! \snippet camp_tutorial/part_1_code/box_model.F90 Species ids
!! \snippet camp_tutorial/part_1_code/box_model.F90 Get species ids
!!
!! The \ref pmc_chem_spec_data::chem_spec_data_t "camp_spec_data_t"
!! object provides access to information about
!! the chemical species present in the system. If there is a problem
!! setting the \ref pmc_chem_spec_data::chem_spec_data_t "chem_spec_data_t"
!! pointer, for example if this function
!! was called before initializing the
!! \ref pmc_camp_core::camp_core_t "camp_core_t" object, this function
!! returns \c false, otherwise it returns \c true. The
!! \ref pmc_chem_spec_data::gas_state_id "gas_state_id()"
!! function returns the index on the state array for a gas-phase
!! species. We will use these indices later on to set the initial
!! conditions and to retrieve the modeled species concentrations.
!! The \ref pmc_chem_spec_data::gas_state_id "gas_state_id()"
!! function returns 0 if a species is not found,
!! so it is important to check the values after calling this function.
!!
!! Now, let's set up some conditions for our box model:
!!
!! \snippet camp_tutorial/part_1_code/box_model.F90 Set initial conditions
!!
!! The \ref pmc_camp_core::solver_initialize "solver_initialize()" function
!! gets the external solver (CVODE) ready to solve the chemical system.
!! The \c camp_state now can be used to describe the state of the
!! chemical system. It contains species concentrations (for which we will
!! provide initial conditions in the next part of the tutorial) and
!! environmental parameters. The temperature and pressure can be
!! updated at any time before or between calls to the CAMP
!! \ref pmc_camp_core::solve "solve()"
!! function, but \ref pmc_camp_state::update_env_state "update_env_state()"
!! must be called after changing
!! either of these properties. Gas-phase species on the state array are
!! in ppm.
!!
!! With that, all that's left is to solve the chemistry and output the
!! results. We'll include a loop over time so we have something
!! interesting to plot.
!!
!! \snippet camp_tutorial/part_1_code/box_model.F90 Time id
!!
!! \snippet camp_tutorial/part_1_code/box_model.F90 Solve and output
!!
!! The \ref pmc_camp_core::solve "solve()" function advances the model
!! state stored in \c camp_state
!! by solving the chemistry over the time step specified in the
!! second argument (in this case, 0.01 s). Deallocating the \c camp_core
!! and \c camp_state pointers releases all the memory associated with
!! CAMP.
!!
!! That's it for the initial box model code. Before we run the model,
!! we'll need to describe the chemical system to solve. We'll do that in
!! \ref camp_tutorial_part_2 "part 2 of Boot CAMP"!
!!
!! The full box model code can be found in
!! `\doc\camp_tutorial\part_1_code`.


! ***********************************************************************
! ***********************************************************************
! ***********************************************************************


!> \page camp_tutorial_part_2 Boot CAMP: Part 2 - Mechanism
!!
!! In the \ref camp_tutorial_part_1 "last installment of Boot CAMP" we
!! wrote the code for a simple box model. This time we'll build a simple
!! chemical mechanism to run in the box model.
!!
!! All CAMP input files are in <a href="http://json.org">json</a>
!! format, a widely used standard for structured data. Many free tools
!! are available online to compose and validate \c json files, and many
!! editors provide \c json syntax highlighting. If you are new to
!! working with \c json files, we recommend writing the code in an online
!! editor/validator like <a href="https://jsonlint.com">JSONLint</a>.
!!
!! There are two types of input files used by CAMP. We'll start with the
!! simplest one. This is the file whose path we passed to the
!! \ref pmc_camp_core::camp_core_t "camp_core_t()" constructor in
!! \ref camp_tutorial_part_1 "part 1" of the tutorial. We named this
!! file \c my_config_file.json and we'll make its contents as follows:
!! \code{.json}
!! {
!!   "pmc-files" : [
!!     "my_simple_mechanism.json"
!!   ]
!! }
!! \endcode
!! CAMP configuration \c json files begin and end with curly brackets
!! ("{}") that
!! enclose an object named \b pmc-files, which is
!! a comma-separated array of file names that make up the chemical
!! mechanism to load into CAMP. The mechanism data can be organized
!! however you like, into as many files as you'd like. Thus, any number of
!! \b pmc-files may be specified and the arrangement of mechanism
!! elements (species, reactions, etc.) within those files is up to you.
!! Also note that \c json ignores most white
!! space, so the code above is equivalent to:
!! \code{.json}
!! { "pmc-files" :     [ "my_simple_mechanism.json" ]   }
!! \endcode
!!
!! One more note about the CAMP \c json files before we move on. CAMP
!! ignores information in the input files that it is not interested in,
!! as long as the file is in valid \c json format, and this information is
!! not included as an element in an array CAMP uses. Thus, our
!! configuration file could be:
!! \code{.json}
!! {
!!   "note" : "Remember to rename 'my simple mechanism' to something more meaningful",
!!   "pmc-files" : [
!!     "my_simple_mechanism.json"
!!   ],
!!   "change log" : [
!!     "030919 md - created file",
!!     "031019 md - revised file"
!!   ]
!! }
!! \endcode
!! As far as CAMP is concerned, these files are equivalent. This is also
!! a way to include comments in your \c json files, as comment
!! flags are not part of the \c json standard. Note however that adding
!! extra information as an element of the \b pmc-files array (an array
!! that CAMP uses) won't work,
!! as CAMP expects these to be valid input file names.
!!
!! The remaining CAMP input files describe the chemical mechanism and
!! use the following format:
!! \code{.json}
!! {
!!   "pmc-data" : [
!!
!!      ...
!!
!!   ]
!! }
!! \endcode
!! Here, \b pmc-data is a comma-separated array of model element
!! objects. There can be any number of these input files, but they must
!! all enclose their model elements with this text.
!!
!! We'll start off wth a single file that describes our mechanism,
!! \c my_simple_mechanism.json. The order of model elements in
!! the \b pmc-data array is arbitrary. We'll start with chemical
!! species. In our first mechanism, we'll just have four: \f$\ce{O3}\f$,
!! \f$\ce{NO}\f$, \f$\ce{NO2}\f$ and \f$\ce{O2}\f$. The input data for
!! these gas-phase species in the \b pmc-data array is:
!! \code{.json}
!!     {
!!       "name" : "O3",
!!       "type" : "CHEM_SPEC",
!!     },
!!     {
!!       "name" : "NO",
!!       "type" : "CHEM_SPEC",
!!     },
!!     {
!!       "name" : "NO2",
!!       "type" : "CHEM_SPEC",
!!     },
!!     {
!!       "name" : "O2",
!!       "type" : "CHEM_SPEC",
!!     }
!! \endcode
!! All CAMP model elements must have a unique name that is chosen by the
!! user and a type that must be one of a set of CAMP data types. For
!! chemical species, this type is \c CHEM_SPEC. Chemical species default
!! to being gas-phase species, but can be specified as being
!! condensed-phase, as we'll see later on.
!!
!! Now, let's build our mechanism. We'll start with just two
!! Arrhenius-type reactions:
!! \code{.json}
!!     {
!!       "name" : "my simple mechanism",
!!       "type" : "MECHANISM",
!!       "reactions" : [
!!         {
!!           "type" : "ARRHENIUS",
!!           "reactants" : {
!!             "NO" : { },
!!             "O3" : { }
!!           },
!!           "products" : {
!!             "NO2" : { },
!!             "O2" : { }
!!           },
!!           "A" : 26.59
!!         },
!!         {
!!           "type" : "PHOTOLYSIS",
!!           "reactants" : {
!!             "NO2" : { }
!!           },
!!          "products" : {
!!            "NO" : { },
!!            "O" : { }
!!          },
!!          "my photo label" : "NO2 photolysis"
!!         },
!!         {
!!           "type" : "ARRHENIUS",
!!           "reactants" : {
!!             "O" : { },
!!             "O2" : { }
!!           },
!!           "products" : {
!!             "O3" : { }
!!           },
!!           "A" : 2.183E-5
!!         }
!!       ]
!!     }
!! \endcode
!! CAMP \b MECHANISM objects are collections of reactions, which are
!! specified in the \b reactions array.
!! For each reaction several elements must be specified. For
!! Arrhenius-like reactions, these include the \b reactants and
!! \b products, as well as the pre-exponential factor \b A. They
!! also typically have some optional parameters, which assume default
!! values unless they are specified in the input files. A description of
!! the format used for each reaction's input data is described
!! \ref camp_rxn "here". The empty curly brackets after the products and
!! reactants allow for the inclusion of information specific to these
!! species, such as reactant quantities (for self reactions) and product
!! yields. For Arrhenius-like reactions, these are described in more
!! detail \ref camp_rxn_arrhenius "here".
!!
!! The only key-value pair not required by CAMP, but that is present in this
!! input data is <b>my photo label</b> in the \f$\ce{NO2}\f$ photolysis
!! reaction. We'll use this label in
!! \ref camp_tutorial_part_3 "part 3 of Boot CAMP" to set the photolysis
!! rate from our box model, and start generating results!
!!
!! The full configuration and mechanism \c json files described in the
!! part of the tutorial can be found in \c /doc/camp_tutorial/part_2_code.


! ***********************************************************************
! ***********************************************************************
! ***********************************************************************


!> \page camp_tutorial_part_3 Boot CAMP: Part 3 - Updating CAMP Parameters
!!
!! So far, we've \ref camp_tutorial_part_1 "built a simple box model"
!! and \ref camp_tutorial_part_2 "set up a simple chemical mechanism"
!! to run in that box model. Now, we're going to make some final
!! adjustments to make this a fully operational box model. In many
!! cases, the host model needs to provide some information to CAMP that
!! changes over the course of a model run. We've seen how to update
!! environmental conditions like temperature and pressure, but other
!! information requires some knowledge of the system being modeled. In
!! the case of our simple mechanism, this is \f$\ce{NO2}\f$ photolysis. CAMP
!! photolysis reactions need to know the photolysis rate at a given
!! moment in the model run. Other reactions that need some external help
!! include \ref pmc_rxn_emission::rxn_emission_t "emissions",
!! \ref pmc_rxn_first_order_loss::rxn_first_order_loss_t
!! "first order loss", and
!! \ref pmc_rxn_wet_deposition::rxn_wet_deposition_t "wet deposition".
!!
!! Let's update our box model to set the \f$\ce{NO2}\f$ photolysis rate.
!! Before the call to \ref pmc_camp_core::solver_initialize
!! "solver_initialize()", we'll add the following code:
!!
!! \snippet camp_tutorial/part_3_code/box_model.F90 NO2 photolysis modules
!! \snippet camp_tutorial/part_3_code/box_model.F90 NO2 photolysis variables
!! \snippet camp_tutorial/part_3_code/box_model.F90 Find NO2 photolysis
!!
!! We first find the chemical mechanism, which we named <b>my simple
!! mechanism</b> in the last part of the tutorial. Next, we cycle through
!! the reactions in that mechanism looking for
!! \ref pmc_rxn_photolysis::rxn_photolysis_t "rxn_photolysis_t"
!! reactions. Then, we check the properties of the photolysis reactions
!! looking for a key-value pair name <b>my photo label</b> and make sure
!! its value is <b>NO2 photolysis</b>, as we specified in the last
!! section. The \ref pmc_rxn_data::rxn_data_t::property_set "property_set"
!! of reactions gives you direct access to the input data for each reaction.
!! This includes the information CAMP uses, like \b reactants and \b A, as
!! well as those it doesn't use, like our <b>my photo label</b>. There are
!! functions to get string, integers, real numbers, and subsets of
!! properties from
!! \ref pmc_rxn_data::rxn_data_t::property_set "property_set". To see
!! the available functions, see \ref pmc_property. The last step is to
!! fix our \ref pmc_rxn_photolysis::rxn_update_data_photolysis_t
!! "rxn_update_data_photolysis_t" object to the \f$\ce{NO2}\f$
!! photolysis reaction we located. This and similar objects for
!! other reactions that accept external updates allow you to change
!! reaction parameters during a model run.
!!
!! Finally, we'll add the following code before we start to loop over
!! time and solve the chemistry:
!!
!! \snippet camp_tutorial/part_3_code/box_model.F90 Set NO2 photolysis
!!
!! Here, we simply set the reaction property of interest in our update
!! data object, in this case the photolysis rate, and then pass this data
!! to CAMP. This can be done before any call to
!! \ref pmc_camp_core::solve "solve()". This rate will remain the same
!! until we change it again.
!!
!! Now, our box model code and our input files are complete. To compile
!! the code, try something like:
!! \code{.sh}
!! 
!! \endcode
!!
!! Then, to run to mechanism:
!! \code{.sh}
!! 
!! \endcode
!! If you have <a href="http://www.gnuplot.info/">gnuplot</a> installed,
!! you can use /doc/camp_tutorial/part_3_code/plot.conf to check out the
!! results.
!!
!! In the \ref camp_tutorial_part_4 "next installment" of Boot CAMP,
!! we'll start passing messages!
!!
!! The files described in this part of the tutorial and needed to run the
!! box model and plot the results can be found in
!! \c /doc/camp_tutorial/part_2_code.


! ***********************************************************************
! ***********************************************************************
! ***********************************************************************


!> \page camp_tutorial_part_4 Boot CAMP: Part 4 - Message Passing
!!
!!



