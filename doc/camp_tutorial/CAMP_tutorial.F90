! Copyright (C) 2019 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \page camp_tutorial Boot CAMP: The CAMP Tutorial
!!
!! In the Boot CAMP tutorial, we build a simple box model that
!! incorporates all the CAMP functionality to demonstrate how CAMP can be
!! included in your favorite model.
!!
!! - \ref camp_tutorial_part_0
!! - \ref camp_tutorial_part_1
!! - \ref camp_tutorial_part_2
!! - \ref camp_tutorial_part_3
!! - \ref camp_tutorial_part_4
!! - \ref camp_tutorial_part_5
!! - \ref camp_tutorial_part_6
!! - \ref camp_tutorial_part_7
!!
!! Model code described in Boot CAMP can be found in
!! \c doc/camp_tutorial .
!!
!! If you have Docker installed and want to quickly run the code
!! described in the tutorial, start a container with PartMC:
!! \code{.sh}
!!   docker run -it compdyn/partmc:develop-85-tutorial bash
!! \endcode
!! Then, follow the instructions at the bottom of the
!! sections of the tutorial that include executable code.
!! To remove the containers once you're done:
!! \code{.sh}
!!   docker system prune
!! \endcode

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_tutorial_part_0 Boot CAMP: Part 0 - Why CAMP?
!!
!! We get it, there are CAMPing people and there are people who will
!! never go CAMPing. If you're on the fence, this part of the tutorial
!! describes the novel things CAMP can do, and why you should consider
!! CAMPing. If you're already onboard, you can skip this section, pack
!! your bags and move on to \ref camp_tutorial_part_1.
!!
!! There are three features of CAMP that set it apart from other
!! approaches to solving chemical systems in atmospheric models. These
!! are its \ref bc0_ss_combined_solving "combined solving" of gas- and
!! aerosol-phase chemistry and partitioning, its
!! \ref bc0_ss_runtime_config "run-time configuration" and its
!! \ref bc0_ss_portability "portability" across models with different ways
!! of representing aerosol systems and different aerosol mircophysical
!! schemes. These features are described in more detail below and are
!! followed by an example of how a chemist doing laboratory experiments
!! on atmospheric systems can use CAMP to \ref bc0_ss_example "rapidly deploy"
!! their new chemistry to a suite of atmospheric models.
!!
!! ## \anchor bc0_ss_combined_solving Combined solving ##
!!
!! Typically, chemical systems are spread across several distinct modules
!! within atmospheric models (\ref bc0_fig_typical_model "Fig. 1")
!! These modules
!! have often been developed independently and may require significant
!! modification to incorporate new chemical processes (particularly when
!! these span the gas and condensed phases) or to port the module to a
!! new host model with its own set of existing modules. In addition,
!! when interrelated processes with similar rates are solved in separate
!! modules (like gas-phase reactions whose products partition to the
!! condensed phase and undergo further reaction), artifacts from the
!! separated solving can be introduced.
!!
!! \image html schematic_typical_model.png
!! \anchor bc0_fig_typical_model Fig. 1. A typical configuration for
!! chemistry and chemistry-adjacent process in atmospheric models.
!!
!!
!! CAMP takes a difference approach to solving chemical systems. Instead
!! of solving chemistry across a collection of modules, CAMP accepts
!! rates for processes from modules that would typically directly update
!! the model state (like emissions and deposition), and it uses an
!! object-oriented approach to build an integrated multi-phase chemical
!! system that includes emissions, deposition, gas-phase reactions
!! including photolysis, partitioning between the gas and aerosol phase,
!! and condensed phase reactions in any number of unique aerosol phases
!! (organic, aqueous aerosol, cloud droplets, etc.;
!! \ref bc0_fig_campground "Fig. 2"). The collection of
!! objects representing the chemistry and related processes are then
!! solved as a single kinetic system, avoiding artifacts from operator
!! splitting and greatly easing the processes of incorporating new
!! multi-phase chemical processes.
!!
!! \image html schematic_CAMP_structure.png
!! \anchor bc0_fig_campground Fig. 2. Schematic showing how CAMP
!! interacts with a host model and how integrated chemical systems are
!! described within CAMP.
!!
!! We'll see in later parts of the tutorial how these objects are
!! generated and configured, and how CAMP works with the way your model
!! describes aerosol systems (bins, modes, etc.).
!!
!!
!!
!! ## \anchor bc0_ss_runtime_config Run-time configuration ##
!!
!! Another key feature of CAMP is its run-time configurability. As a
!! stand-alone library, you don't need to modify any of the CAMP source
!! code to incorporate CAMP into a new model, but you do need to tell CAMP
!! about the chemical system you want to solve and how your model
!! describes aerosol systems. This is primarily done at run-time using a
!! collection of `JSON` input files (\ref bc0_fig_json "Fig. 3").
!! (There are also ways to update
!! certain CAMP parameters during a model run, like emissions and
!! photolysis rates.) This opens up many exciting possibilities. To
!! start, you can change the chemical mechanism without any changes to the
!! model source code. Data assimilation and model sensitivity analyses can
!! take advantage of the ability to adjust a wide variety of model
!! parameters (everything from activation energies for Arrhenius reactions
!! to ion-pair interaction parameters in activity calculations) without
!! any changes to source code or recompilation of the model.
!!
!! \image html schematic_json_objects.png
!! \anchor bc0_fig_json Fig. 3. Examples of `JSON` input data used with
!! CAMP.
!!
!!
!!
!! ## \anchor bc0_ss_portability Portability ##
!!
!! The last thing we'll mention that makes CAMP unique is its
!! portability across models with different aerosol representations. We'll
!! describe how to make CAMP work with your model's aerosol representation
!! in \ref camp_tutorial_part_5. Here, we'll just provide an overview of
!! how CAMP makes this work. Similar to other CAMP model elements, your
!! input data can tell CAMP to create any number of
!! \ref camp_aero_phase "aerosol phases". These are simply collections
!! of species, like an "aqueous" phase that includes "water", "sulfate",
!! "nitrate", etc. Condensed-phase and partitioning reactions must specify
!! the aerosol phase they apply to. Thus, you could have a
!! \ref camp_rxn_HL_phase_transfer "Henry's law phase transfer" reaction
!! for nitric acid that applies to the "aqueous" phase.
!!
!! The key to separating the chemistry from the host model's aerosol
!! representation is that the chemistry is the same for every instance of
!! a particular aerosol phase, but the number of instances of these phases,
!! and the physical properties of the aerosols of which they are a part,
!! are determined by the host model's aerosol representation. A
!! description of this representation is provided by another `JSON` input
!! object. This input data specifies not only the type and dimensions of
!! the aerosol representation (e.g., the number and size of bins, or the
!! number and shape of modes) but also which aerosol phases are associated
!! with with which aerosol groups. Some examples of how different aerosol
!! representations implement instances of the same set of aerosol phases in
!! the CAMP state array are shown in \ref bc0_fig_aero_reps "Fig. 4".
!!
!! \image html modal_aero_rep.png
!! \image html binned_aero_rep.png
!! \image html single_particle_aero_rep.png
!! \anchor bc0_fig_aero_reps Fig. 4. How modal (top) binned (middle) and
!! single-particle (bottom) aerosol representations might implement the
!! same three aerosol phases.
!!
!!
!!
!! ## \anchor bc0_ss_example How you might CAMP ##
!!
!! To demonstrate how CAMP functionality can be used in a real-world
!! scenario, imagine you are a reasearch chemist who just discovered a new
!! gas-phase reaction:
!! \f[\ce{
!!   A + OH -> B
!! }\f]
!! Species B can partition to an aqueous phase where it reacts with
!! nitrate:
!! \f[\ce{
!!   B + NO3- -> C
!! }\f]
!! Species C can then repartition back to the gas-phase.
!!
!! You first try this out in a box model running CAMP with a standard
!! gas and aerosol phase mechansim by adding a
!! single `JSON` file to the CAMP input file list. Your new file has three
!! \ref input_format_species "chemical species", one
!! \ref camp_rxn_arrhenius "gas-phase Arrhenius" reaction, one
!! \ref camp_rxn_HL_phase_transfer "Henry's law phase transfer" reaction,
!! and one \ref camp_rxn_condensed_phase_arrhenius "condensed-phase Arrhenius"
!! reaction. The only modification to the existing mechanism input files
!! you make is to add species C to the "aequeous" aerosol phase.
!! You evaluate the model results and tweak the reaction parameters
!! until you're able to fit your experimental results.
!!
!! Next, you take your mechanism input files and run them in a PartMC
!! particle-resolved urban plume scenario to evaluate the effects of
!! mixing state on your newly discovered system. Then you try out your
!! mechanism in the MONARCH chemical weather prediction system to see the
!! regional and global scale impacts of your new chemistry. What is most
!! important about this process, and what makes CAMP so unique, is that
!! you use the \a the \a same \a input \a files in all three models and
!! you never
!! modify source code or recompile a model. Your chemistry is running in
!! the exact same version of CAMP in all three models. If that doesn't
!! make you want to try CAMPing, probably nothing will.
!!
!! So lace-up your hiking boots and stay hydrated. It's time to start
!! \ref camp_tutorial_part_1 "the first stage of Boot CAMP"!
!!
!! <hr>
!! \image{inline} html icon_trail.png
!! \ref camp_tutorial "Index"
!! \image{inline} html icon_trail.png <b> Next: </b>
!! \ref camp_tutorial_part_1 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_tutorial_part_1 Boot CAMP: Part 1 - Box Model
!!
!! Prior to beginning this tutorial, the PartMC library should be
!! installed on your system with PartMC-CAMP enabled. Installation
!! instructions can be found \ref camp_chem "here". Alternatively,
!! you can run PartMC in Docker following the instructions
!! \ref camp_tutorial "here".
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
!! The \ref pmc_chem_spec_data::chem_spec_data_t "chem_spec_data_t"
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
!! gets the external solver
!! (<a href="https://computing.llnl.gov/projects/sundials/cvode">CVODE</a>)
!! ready to solve the chemical system.
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
!! second argument. We're keeping the time step small (\f$10^{-15}\f$ s)
!! for now so we have something interesting to plot from these very fast
!! reactions. Deallocating the \c camp_core
!! and \c camp_state pointers releases all the memory associated with
!! CAMP.
!!
!! That's it for the initial box model code. Before we run the model,
!! we'll need to describe the chemical system to solve. We'll do that in
!! \ref camp_tutorial_part_2 "part 2 of Boot CAMP"!
!!
!! The full box model code can be found in
!! `\doc\camp_tutorial\part_1_code`.
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_tutorial_part_0
!! \image{inline} html icon_trail.png
!! \ref camp_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_tutorial_part_2 <b> > </b>

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
!! species. In our first mechanism, we'll just have five: \f$\ce{O3}\f$,
!! \f$\ce{NO}\f$, \f$\ce{NO2}\f$, \f$\ce{O2}\f$ and \f$\ce{O}\f$.
!! The input data for
!! these gas-phase species in the \b pmc-data array is:
!! \code{.json}
!!     {
!!       "name" : "O3",
!!       "type" : "CHEM_SPEC"
!!     },
!!     {
!!       "name" : "NO",
!!       "type" : "CHEM_SPEC"
!!     },
!!     {
!!       "name" : "NO2",
!!       "type" : "CHEM_SPEC"
!!     },
!!     {
!!       "name" : "O2",
!!       "type" : "CHEM_SPEC"
!!     },
!!     {
!!       "name" : "O",
!!       "type" : "CHEM_SPEC"
!!     },
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
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_tutorial_part_1
!! \image{inline} html icon_trail.png
!! \ref camp_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_tutorial_part_3 <b> > </b>

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
!! photolysis reaction we located, then at the end just make sure we
!! found the reaction we were looking for. This and similar objects for
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
!!   gfortran -o run_box_model box_model.F90 -lpartmc -I/usr/local/include/partmc
!! \endcode
!! Where the include path points to where the PartMC library \c .mod
!! files are installed. If you have trouble compiling or running because
!! of missing libraries, make sure your `LD_LIBRARY_PATH` and `PATH`
!! include the directories where the PartMC, json-fortran, SUNDIALS,
!! netCDF, and SuiteSparse libraries are installed.
!!
!! Then, to run to mechanism:
!! \code{.sh}
!!   ./run_box_model > output.txt
!! \endcode
!! If you have <a href="http://www.gnuplot.info/">gnuplot</a> installed,
!! you can copy /doc/camp_tutorial/part_3_code/plot.conf to the
!! directory with your results and check out them out:
!! \code{.sh}
!!   gnuplot plot.conf
!! \endcode
!! Our results look like this:
!!
!! \image html BootCAMP_part_3_results.png
!!
!! In the \ref camp_tutorial_part_4 "next installment" of Boot CAMP,
!! we'll start passing messages!
!!
!! The files described in this part of the tutorial and needed to run the
!! box model and plot the results can be found in
!! \c /doc/camp_tutorial/part_3_code.
!!
!! <hr>
!! ### Docker Instructions ###
!! Inside the container:
!! \code{.sh}
!!   dnf install -y gnuplot
!!   mkdir boot-camp
!!   cd boot-camp
!!   cp ../partmc/doc/camp_tutorial/part_3_code/* .
!!   gfortran -o run_box_model box_model.F90 -lpartmc -I/usr/local/include/partmc
!!   ./run_box_model > output.txt
!!   gnuplot plot.conf
!!   exit
!! \endcode
!! Back outside the container:
!! \code{.sh}
!!   docker cp pmc:/boot-camp/results.png .
!!   docker container rm pmc
!!   open results.png
!! \endcode
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_tutorial_part_2
!! \image{inline} html icon_trail.png
!! \ref camp_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_tutorial_part_4 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************


!> \page camp_tutorial_part_4 Boot CAMP: Part 4 - Message Passing
!!
!! This part of \ref camp_tutorial "Boot CAMP" shows how to use CAMP's
!! message passing functions. If you're only interested in using CAMP on
!! a single processor, you can skip this part and move on to
!! \ref camp_tutorial_part_5.
!!
!! We'll wrap our MPI code with a compiler flag named `USE_MPI` to make
!! sure our box model can be built with or without MPI. The order of
!! operations is important for MPI runs and is summarized in the following
!! table.
!!
!! | Process   | Operation                                                      |
!! |-----------|----------------------------------------------------------------|
!! | primary   | `camp_core => camp_core_t( input_files )`                      |
!! | primary   | `call camp_core%%initialize( )`                                |
!! | primary   | access `camp_core_t` properties/set up `update_data_t` objects |
!! | primary   | pack all objects on a buffer                                   |
!! | all       | pass the buffer                                                |
!! | secondary | `camp_core => camp_core_t( )`                                  |
!! | secondary | unpack the `camp_core_t` and other objects from the buffer     |
!! | all       | `call camp_core%%solver_initialize( )`                         |
!! | all       | use `update_data_t` objects to update rates, etc.              |
!! | all       | `call camp_core%%solve( camp_state, time_step )`               |
!! | all       | deallocate all objects                                         |
!!
!! We'll go through this step-by-step, update our box model and discuss
!! why each process in done when and where it is.
!!
!! Note that the PartMC MPI functions use `MPI_WORLD_COMM` by default,
!! but they accept an optional `comm` argument if you would like to use a
!! different communicator. See the specific function documentation for
!! details.
!!
!! First, let's add the modules we need for MPI. We'll use the standard
!! mpi module and the PartMC mpi module, with some custom functions.
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 MPI modules
!!
!! Now we'll declare a buffer, a position index, and a pack size
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 MPI variables
!!
!! Next, let's initialize MPI and wrap some of our existing code
!! in a conditional statement
!! that ensures we load the input data and initialize CAMP on the
!! primary process only (we're including the existing call to the
!! \ref pmc_camp_core::camp_core_t "camp_core_t" constructor and
!! `camp_core_t::initialize()` to show the
!! placement of the start of our new conditional block):
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 wrap initialization
!!
!! The `camp_core_t::initialize()` subroutine instructs the internal
!! model elements to take their input data and condense it down into a
!! small data block containing only the information they need to solve the
!! chemical system during calls to `camp_core_t::solve()`. The
!! \ref pmc_camp_core::camp_core_t "camp_core_t" MPI functions pass only
!! this condensed data to other processes. So, after the core is passed,
!! you will not have access
!! to the raw input data or model \ref pmc_property::property_t "property_t"
!! objects that we used to set up the
!! \ref pmc_rxn_data::rxn_update_data_t "rxn_update_data_t" objects in
!! \ref camp_tutorial_part_3 "part 3".
!! Thus, all the setup of
!! \ref pmc_rxn_data::rxn_update_data_t "rxn_update_data_t"
!! objects must be done on the
!! primary process, before passing the core and update objects to the
!! other processes.
!!
!! So, let's end our first MPI conditional block after we setup
!! the \f$\ce{NO2}\f$ photolysis
!! \ref pmc_rxn_data::rxn_update_data "rxn_update_data_t" object and
!! before the call to `camp_core_t::solver_initialize()`.
!! The first step is to get the size of the buffer to be used to pass
!! the objects
!! (the existing check that the \f$\ce{NO2}\f$
!! photolysis update data object is attached is included to show the
!! placement of the following code block):
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 get pack size
!!
!! After we allocate the buffer on the primary process, we'll pack it
!! with the object data:
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 pack objects
!!
!! Next, we'll pass the species indexes we looked up. (Remember, we
!! won't be able to do this on the secondary processes.)
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 pass indices
!!
!! After we pack the objects and exit the primary process block, we'll
!! pass the buffer to the other processes:
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 pass the buffer
!!
!! Next, we'll unpack the objects on the secondary processes:
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 unpack the objects
!!
!! Note that we call the \ref pmc_camp_core::camp_core_t "camp_core_t"
!! constructor without passing the input file list. This creates an
!! empty core on the secondary processes that we can fill with the packed
!! data from the buffer.
!! After unpacking the objects and deallocating the buffer, our message
!! passing is complete, and the rest of the code remains the same, beginning
!! with the call to `solver_initialize()`.
!!
!! This is not a very useful parallelization of our box model, as we're
!! just solving the same system on every process, but it demonstrates how
!! to initialize and pass the `camp_core_t` and `update_data_t` objects.
!! The `camp_state_t::state_var(:)` array can be accessed directly and
!! passed however your model passes double-precision
!! floating-point arrays, or you can use the
!! `pmc_mpi_pack_size_real_array()`, `pmc_mpi_pack_real_array()`,
!! and `pmc_mpi_unpack_real_array()` functions.
!!
!! To finish up, let's add a conditional block around the output to
!! print the results from the first secondary process, just
!! to make sure our message passing is working, and finalize MPI.
!!
!! \snippet camp_tutorial/part_4_code/box_model.F90 output
!!
!! To compile the model code with mpi, be sure to include the `USE_MPI`
!! flag definition:
!! \code{.sh}
!!   mpif90 -o run_box_model box_model.F90 -DUSE_MPI -lpartmc -I/usr/local/include/partmc
!!   mpirun -v -np 2 run_box_model > output.txt
!! \endcode
!!
!! In later installments of \ref camp_tutorial "Boot CAMP" we'll include
!! a section towards the end that describes any MPI-related code needed to
!! run the updates described.
!!
!! Now that our messages are passed, it's aerosol time. That's the
!! topic of the \ref camp_tutorial_part_5 "next installment of Boot CAMP"!
!!
!! <hr>
!! ### Docker Instructions ###
!! To run a Docker container with MPI support, we'll need to build the
!! image locally. So, we'll clone the PartMC repo, build the container
!! with MPI and then run it:
!! \code{.sh}
!!   git clone https://github.com/compdyn/partmc.git
!!   cd partmc
!!   git checkout develop-85-tutorial
!!   docker build -f Dockerfile.mpi -t partmc-test-mpi .
!!   docker run ---name pmc -it partmc-test-mpi bash
!! \endcode
!! Inside the container:
!! \code{.sh}
!!   sudo dnf install -y gnuplot
!!   mkdir boot-camp
!!   cd boot-camp
!!   cp ../partmc/doc/camp_tutorial/part_4_code/* .
!!   mpif90 -o run_box_model box_model.F90 -DUSE_MPI -lpartmc -I/usr/local/include/partmc
!!   mpirun -v -np 2 run_box_model > output.txt
!!   gnuplot plot.conf
!!   exit
!! \endcode
!! Back outside the container:
!! \code{.sh}
!!   docker cp pmc:/boot-camp/results.png .
!!   docker container rm pmc
!!   open results.png
!! \endcode
!! You should get the same results as described in
!! \ref camp_tutorial_part_3
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_tutorial_part_3
!! \image{inline} html icon_trail.png
!! \ref camp_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_tutorial_part_5 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_tutorial_part_5 Boot CAMP: Part 5 - Aerosol Representations
!!
!! \todo finish
!!
!!
!!
!!
!!
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_tutorial_part_4
!! \image{inline} html icon_trail.png
!! \ref camp_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_tutorial_part_6 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_tutorial_part_6 Boot CAMP: Part 6 - Aerosol Representation Input Data
!!
!! \todo finish
!!
!!
!!
!!
!!
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_tutorial_part_5
!! \image{inline} html icon_trail.png
!! \ref camp_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_tutorial_part_7 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_tutorial_part_7 Boot CAMP: Part 7 - Sub Models
!!
!! \todo finish
!!
!!
!!
!!
!!
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_tutorial_part_6
!! \image{inline} html icon_trail.png
!! \ref camp_tutorial "Index"
!! \image{inline} html icon_trail.png
