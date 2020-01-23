! Copyright (C) 2019 Matt Dawson and Christian Guzman
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \page camp_dev_tutorial Dev CAMP: The CAMP developer guide
!!
!! Welcome developer!
!!
!! This is the Dev CAMP tutorial.
!!
!! The tutorial aims for understanding the logic behind CAMP, his
!! main components and the interactions between them. We hope this
!! guide will help you to understand the insides of CAMP and allows you
!! to contribute to the CAMP development.
!!
!! Before starting, it's recommended have finalized the \ref camp_tutorial
!! Boot CAMP tutorial.
!!
!! - \ref camp_dev_tutorial_part_0
!! - \ref camp_dev_tutorial_part_1
!! - \ref camp_dev_tutorial_part_2
!! - \ref camp_dev_tutorial_part_3
!! - \ref camp_dev_tutorial_part_4
!! - \ref camp_dev_tutorial_part_5
!! - \ref camp_dev_tutorial_part_6
!! - \ref camp_dev_tutorial_part_7
!!
!! Model code described in Dev CAMP can be found in
!! \c doc/camp_tutorial/dev_camp.
!!
!! This tutorial focuses in main concepts. Some of the concepts explained can differ
!! in the future if better implementations are found, but the tutorial should be
!! enough to replicate a consistent CAMP skeleton.

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_0 Dev CAMP: Part 0 - Introduction
!!
!! CAMP can be divided into the following main components:
!!
!! - Input user data: All input data received by CAMP.
!! It includes the user-defined JSON files and the variables coming from
!! the host model such as chemical concentrations and environmental variables
!! like temperature or pressure. The \ref camp_tutorial "Boot CAMP tutorial"
!! teaches how to enter these variables in CAMP.
!!
!! - CAMP_core interface: User interface.
!! Responsible for adapting the input data to efficent data structures.
!! Includes all the files that prepare the data for the solving system.
!!
!! - Solving system: Corresponds to all the files related to the resolution of the
!! chemical mechanism. Includes the solving of the chemical reactions equations and the
!! ODE solver. At present, the external package CVODE is used as the CAMP ODE solver.
!!
!! \image html very_simple_CAMP_components.jpg
!!
!! <hr>
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png <b> Next: </b>
!! \ref camp_dev_tutorial_part_1 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_1 Dev CAMP: Part 1 - Mathematical considerations
!!
!! In general, CAMP tries to predict the concentrations of a chemical species
!! at some point in the future by solving a set of ordinary differential
!! equations (ODE) that represent the reactions that compose a chemical
!! \ref camp_mechanism "mechanism".
!!
!! \ref camp_rxn "Reactions" take the general form:
!!
!! \f[
!! c_1y_1+c_2y_2 + \dots + c_my_m \leftrightarrow c_{m+1}y_{m+1}+c_{m+2}y_{m+2} + \dots + c_ny_n,
!! \f]
!!
!! where species \f$y_i\f$ is a participant in the reaction with stoichiometric
!! coefficient \f$c_i\f$. The rate of change for each participating species \f$y_i\f$ with
!! respect to reaction \f$j\f$ is given by
!!
!! \f[
!!   \left(\frac{dy_i}{dt}\right)_j =
!!  \begin{cases}
!!    \quad -c_ir_j(\mathbf{y},T,P,\dots) & \quad \text{for } i \le m; \\
!!    \quad c_ir_j(\mathbf{y},T,P,\dots) & \quad \text{for } m < i \le n \\
!!  \end{cases},
!! \f]
!!
!! where the rate \f$r_j\f$ of reaction \f$j\f$ is an often complex function of the entire
!! model state (including species concentrations \f$y\f$, environmental conditions,
!! such as temperature, \f$T\f$, and pressure, \f$P\f$, physical aeorosl properties,
!! such as surface area density and number concentration, etc.). The
!! overall rate of change for each species \f$y_i\f$ at any given time is thus,
!!
!! \f[
!! f_i \equiv \frac{dy_i}{dt} = \sum_j \left(\frac{dy_i}{dt}\right)_j,
!! \f]
!!
!! where \f$f\f$ is referred to as the derivative of the system.
!!
!! Each reaction type calculates its rate differently.
!! The simplest example of a typical atmospheric reaction is an
!! \ref camp_rxn_arrhenius "Arrhenius-like reaction". For more information,
!! check the other \ref camp_rxn "reaction" types pages.
!!
!! ## ODE properties ##
!!
!! Solving an ODE implies selecting the solving method that better fills
!! its properties.
!!
!! For the chemical mechanism scope, the system to solve can be either
!! a non-stiff problem as a stiff problem.
!! Since CAMP aims to solve all the types of mechanisms, the ODE method
!! must be implicit to secure stability and a
!! reasonable number of solving time-steps for stiff mechanism cases.
!!
!! All solver methods require the system derivative in each solver iteration.
!! But implicit methods tipically needs also the system Jacobian.
!!
!! Solving the derivative and the Jacobian also includes a wide range
!! of solver options that depends on the system properties. Common mechanisms
!! consist of a great number of reactions (from tens to hundreds) and a
!! high number of total different parameters. These properties translate
!! into a small Jacobian matrix with few non-zero values, suitable for
!! direct sparse solving methods.
!!
!! ## CVODE configuration ##
!!
!! To fullfill the properties explained, CAMP uses the
!! <a href=https://computing.llnl.gov/projects/sundials/cvode>CVODE</a>
!! solving package
!! with the following configuration:
!!
!! - Stiffness: Backward Differentiation Formula (BDF) method.
!! - Direct sparse linear solver: KLU sparse solver.
!! - Level of parallelization: No parallelization, serial execution.
!!
!! Overral structure diagram of the CVODE package with the options
!! configured inside red circles:
!!
!! \image html cvode_configuration.png
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_0
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_dev_tutorial_part_2 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_2 Dev CAMP: Part 2 - CAMP_core interface
!!
!! This section explains the camp_core_interface. We recommend read the
!! \ref camp_tutorial "Boot CAMP tutorial" to understand the definitions of
!! JSON files, input data and the user interaction with the CAMP_core
!! interface.
!!
!! The CAMP_core interface follows an structure based on classes, where the
!! more relevants are:
!!
!!  - \ref pmc_camp_core::camp_core_t "camp_core":
!!  Main class, includes and manage all the program classes.
!!  - \ref pmc_camp_state::camp_state_t "camp_state":
!!  Includes concentration array and the env_state class.
!!  - \ref pmc_env_state::env_state_t "env_state": Includes all the
!!  environmental program variables. The most remarkables are temperature
!!  and pressure.
!!  - \ref pmc_camp_solver_data "camp_solver_data":
!!  Allows the CAMP_core interface (programmed in fortran)
!!  use the functions from the solving system (programmed in C as a
!!  requisite of the CVODE package).
!!  - \ref pmc_chem_spec_data "chem_spec_data" : Reads JSON data.
!!  The user can found here the species data defined in JSON files
!!  by searching the reactants names like "O3" or "NO2".
!!  - \ref pmc_rxn_data::rxn_data_t "rxn_data" : Also read JSON data by calling
!!  \ref pmc_property::property_t "property_t".
!!  Includes declarations of the arrays where JSON data is stored: one array for
!!  integer type data and other for floating data (for example, qty from \ref
!!  camp_rxn_arrhenius "arrhenius reaction" is an integer, but \f$A\f$ is floating data).
!!  todo: why more than 1 class read JSON data?
!!  - \ref pmc_rxn_arrhenius::rxn_arrhenius_t "rxn_REACTIONTYPE" : Set of classes,
!!  one for each reaction type defined in CAMP,
!!  where REACTIONTYPE is the reaction name (arrhenius, photolysis, etc.)
!!  Extends from rxn_data, meaning that also has the same variables
!!  and functions as rxn_data. Receive the data read by chem_spec_data
!!  and stores them in the integer and float arrays defined in rxn_data. It also
!!  adds to the integer array some data necesary for future calculations and
!!  iterate the arrays, like the number of reactants and products.
!!  todo: rxn_factory? differences between rxn data, aero_rep_data, aero_phase_data and submodel_data?
!!
!! The workflow corresponds to: (todo: diagram)
!!
!! - <b> From the user point of view </b>: initialize camp_core, create camp_state,
!!  call chem_spec_data through camp_core, set camp_state chemical data and solve.
!! - <b> For the code part </b>: camp_core object is created, camp_core calls
!!  chem_spec_data to read JSON files, chem_spec_data calls rxn_data and
!!  rxn_REACTIONTYPE to save the data read in an integer and a float array,
!!  camp_core receives a camp_state object with the environmental and
!!  concentrations data, and then camp_core calls the \ref
!!  camp_dev_tutorial_part_3 "solving system" through
!!  the camp_solver_data interface, updating the chemical concentrations.
!!
!! Notice that all the reaction data read in the JSON files will be stored
!! in the integer and the float aray. We will name these arrays as int_data and
!! float_data from now on. Only the environmental data array (env_data) and the
!! concentrations array (state_var) are outside off int_data and float_data.
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_1
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_dev_tutorial_part_3 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_3 Dev CAMP: Part 3 - Solving system
!!
!!
!! This section describes the CAMP solving process. This process receive as
!! an input data the arrays explained in the previous section
!! (int_data, float_data, state_var and env_data), solve the derivative and
!! Jacobian functions defined in section 1 and send the results to the CVODE
!! package. The CVODE package will handle the rest to solve the system ODE.
!!
!! The CAMP_core interface interacts only (or should it be) with the
!! \ref camp_solver.c "camp_solver" file, the main file of the solving system.
!!
!! ## Pre-solver operations ##
!!
!! Before calling the CVODE solving process, some operations must be done:
!!
!! 1- Mount the sparse Jacobian structure. Save the index to the
!! non-zero values in the Jacobian. Each index should correspond to a reaction
!! (Jacobian row) and the species (Jacobian species), given that the Jacobian
!! value for a position is zero if the reactant do not exists for this reaction,
!! and non-zero otherwise. Remember that the species present
!! in a reaction are defined in the JSON files, and the Jacobian must be
!! set following the sparse structure.
!!
!! 2- Configure CVODE with the optimal solver options defined in section 1
!! (BDF, KLU sparse, etc.).
!!
!! 3- Calculate the possible rate constants or reaction parameters for each
!! reaction. There are reactions with rate constant/s that only depends
!! in the initial values, like the env_data or initial concentrations.
!! Computing them before the solving will reduce the amount of
!! solver calculations.
!!
!! ## Derivative and Jacobian##
!!
!! CAMP needs to define the function f and Jac, which computes the derivative
!! and the Jacobian respectively. The declaration of these functions is
!! defined in the CVODE package (including the input parameters). We only
!! need to fill the content.
!!
!! These functions will be called in the middle of the CVODE solving process,
!! depending of the system will be called less or more times. CVODE will
!! pass as a function parameters the current concentration array (state_var),
!! the pointer to the future concentration array or Jacobian (our function
!! result must be stored here), the current model time (s) and the pointer
!! to the main CAMP data structure (solver_data). solver_data contains all
!! the CAMP data used during the solving. For example, this struct
!! contains the int_data and float_data arrays, necessaries for solving
!! both derivative and Jacobian.

Now, both derivative and Jacobian functions must calculate the
the rate of change \f$r_j\f$ (mentioned in section 1) for each participating species
\f$y_i\f$ in each reaction, and sum the reaction contributions for the
same species. Example: specie 1 is present in reaction 2 and 3, so the rate
calculated in these reactions corresponding to specie 1 must be added. Remember
also that the sign changes if the species act as a reactant (-) or as a product (+).
Some reaction parameters like the yield can affect also the rate.









!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_2
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_dev_tutorial_part_4 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_4 Dev CAMP: Part 4 - GPU Interface
!!
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_3
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_dev_tutorial_part_5 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_5 Dev CAMP: Part5 - Testing
!!
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_4
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png
!! <b> Next: </b>
!! \ref camp_dev_tutorial_part_6 <b> > </b>

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_6 Dev CAMP: Part 6 -
!!
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_5
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png
!!
!!