! Copyright (C) 2019 Matt Dawson and Christian Guzman
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!! todo: Guide the tutorial with images with clearly identified shapes and colours
!! for each module. Basically regarding the user where this part comes from the initial
!! basic diagram, but with a more specific diagram for the part explained.

!! todo: algorithmic part about operations (multiplying this variable or
!! getting this index), workflow, and this kind of things must be clear.
!! add more diagrams or data structures image with okay i have this jacobian
!! with this indices, something like that. or this struct data image
!! with all of this necessary data (apart from writing it).

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
!! Before starting, it's recommended to have finalized the \ref camp_tutorial
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
!!
!! After complete the tutorial, we recommend take a look at the complete
!! sequence diagram of CAMP version 1.0. It contains most of the variable and functions
!! explained during the tutorial and some more.
!!
!!  <a href=https://www.draw.io/#G1XbNwQZ58d3RWbOK_6_aVvFr9_vW0hwBc>CAMP
!!  sequence diagram</a>
!!
!! Q: develop using the link?
!!
!! \image html CAMP_complete_sequence.jpg
!!

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_0 Dev CAMP: Part 0 - Introduction
!!
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
!! \image html very_simple_CAMP_components_2.jpg
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
!! \ref camp_rxn_arrhenius "Arrhenius-like reaction".
!!
!! \f[
!! r_j = k_j \prod_{i=1}^my_i^{c_i}
!! \f]
!!
!! where \f$k_j\f$ is referred to as the rate constant.
!! For more information about the rate constant, check the
!! \ref camp_rxn "reaction" types pages.
!!
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
!!  Allows the CAMP_core interface (programmed in Fortran)
!!  use the functions from the solving system (programmed in C as a
!!  requisite of the CVODE package).
!!  - \ref pmc_chem_spec_data "chem_spec_data" : Reads JSON data.
!!  The user can found here the species data defined in JSON files
!!  by searching the reactants names like "O3" or "NO2".
!!  - \ref pmc_rxn_data::rxn_data_t "rxn_data" : Also read JSON data by calling
!!  \ref pmc_property::property_t "property_t".
!!  Includes declarations of the arrays where JSON data is stored: one array for
!!  integer type data and other for floating data (for example, qty from \ref
!!  camp_rxn_arrhenius "arrheniufs reaction" is an integer, but \f$A\f$ is floating data).
!!  Q: why more than 1 class read JSON data?
!!  - \ref pmc_rxn_arrhenius::rxn_arrhenius_t "rxn_REACTIONTYPE" : Set of classes,
!!  one for each reaction type defined in CAMP,
!!  where REACTIONTYPE is the reaction name (arrhenius, photolysis, etc.)
!!  Extends from rxn_data, meaning that also has the same variables
!!  and functions as rxn_data. Receive the data read by chem_spec_data
!!  and stores them in the integer and float arrays defined in rxn_data. It also
!!  adds to the integer array some data necesary for future calculations and
!!  iterate the arrays, like the number of reactants and products.
!!  - \ref pmc_mpi "pmc_mpi" : Interface to MPI library. Apart from including a MPI
!!  call on each function, these functions checks if the MPI call success. The
!!  purpose is to compress complex MPI implementations in one call (for example,
!!  broadcasting to all nodes a 2d matrix).
!!  - \ref pmc_rxn_factory::rxn_factory_t "rxn_factory": Interface to create and
!!  add to the system a new reaction in runtime. Apart from JSON files, this can be
!!  another way to enter the reaction data. Check \ref camp_rxn_add for more information.
!!  - \ref pmc_util "pmc_util": Treat an array of characters as a string_t object. Used to store string
!!   informations like species names in chem_spec_data or camp_core.
!!
!!  todo: identify what partmc files are used and which ones not:
!! Used (explain before): mpi, util,
!! Not used: stats, run_coagulation, run_part, run_sect,  aero_binned, aero_state,
!! spec_file (integrated in env_state but only used by partmc)
!!   ...
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
!! concentrations array (concs_data) are outside off int_data and float_data.
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
!! This section describes the CAMP solving process. This process receive as
!! an input data the arrays explained in the previous section
!! (int_data, float_data, concs_data and env_data), solve the Derivative and
!! Jacobian functions defined in \ref camp_dev_tutorial_part_1
!! and send the results to the ODE solver
!! package.
!!
!! The \ref camp_solver.c "CAMP_solver" file have to be the only file that interacts
!! with the external solving parts (CAMP_core interface and ODE solver). The following
!! image illustrates the overall solving system structure:
!!
!! \image html camp_solving_system.jpg
!!
!! ## Choosing the ODE solver: CVODE ##
!!
!! Solving an ODE implies selecting the solving method that better fills
!! its properties.
!!
!! For the chemical mechanism scope, the system to solve can be either
!! a non-stiff problem as a stiff problem.
!! Since CAMP aims to solve all the types of mechanisms, the ODE method
!! have to be implicit to secure a
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
!! To fullfill the properties explained, CAMP uses the
!! <a href=https://computing.llnl.gov/projects/sundials/cvode>CVODE</a>
!! solving package
!! with the following configuration:
!!
!! - Stiffness: Backward Differentiation Formula (BDF) method.
!! - Direct sparse linear solver: KLU sparse solver.
!! - Level of parallelization: No parallelization, serial execution.
!!
!! Note that <b> CAMP is not locked to CVODE </b> and can be adapted to any ODE
!! solver availabe with the proper implementation. The only requirements
!! to achieve an optimal behaviour is select a solver designed for stiff systems.
!!
!! We recommend check the sparse matrix structure before continuing the tutorial.
!! Good documentation about the KLU sparse solver and his structure can be found at
!! <a href=https://fossies.org/linux/SuiteSparse/KLU/Doc/palamadai_e.pdf>
!! "KLU-A HIGH PERFORMANCE SPARSE LINEAR SOLVER FOR CIRCUIT SIMULATION
!! PROBLEMS", by Ekanathan Palamadai Natarajan</a>
!!
!! ## Pre-solver operations ##
!!
!! Before calling the CVODE solving process, some operations must be done:
!!
!! 1- Mount the sparse Jacobian structure. A sparse matrix has the data
!! values array and two arrays of indices that indicates the matrix rows
!! and columns. For our case, rows represents the system reactions and
!! columns represents the species. Each row has a value for all the species
!! present, but if the reaction does not affect the specie the value will be
!! zero everytime and won't be stored in the sparse data array. In addition
!! to this, an extra index structure is created to link reactions with Jacobian
!! in CAMP (\ref JacMap "JacMap").
!!
!! 2- Configure CVODE with the optimal solver options (BDF, KLU sparse, etc.).
!!
!! 3- Calculate the possible rate constants or reaction parameters for each
!! reaction. There are reactions with rate constant/s that only depends
!! in the initial values, like the env_data or initial concentrations.
!! Computing them before the solving will reduce the amount of
!! solver calculations.
!!
!! ## Derivative and Jacobian##
!!
!! CAMP needs to define the function f and Jac, which computes the Derivative
!! and the Jacobian respectively. The declaration of these functions is
!! defined in the CVODE package (including the input parameters). We only
!! need to fill the content.
!!
!! These functions will be called in the middle of the CVODE solving process,
!! depending of the system will be called less or more times. CVODE will
!! pass as a function parameters the current concentration array (concs_data),
!! the pointer to the future concentration array or Jacobian (deriv or J
!! variables), the current model time (s) and the pointer
!! to the main CAMP data structure (solver_data struct). solver_data contains all
!! the CAMP data used during the solving. For example, this struct
!! contains the int_data and float_data arrays, necessaries for solving
!! both Derivative and Jacobian.
!!
!! Once Derivative and Jacobian finish, CVODE will try to solve the ODE. The implicit
!! solving method configured (BDF) will aproximate the result. If the
!! aproximation is near enough to the correct result (convergence), the ODE is solved and
!! CAMP returns to the CAMP solver interface part. Otherwise, CVODE will
!! repeat the iteration using the last result obtained.
!!
!! \image html workflow_deriv_jac.jpg
!!
!! These functions must solve the operations defined at
!! \ref camp_dev_tutorial_part_1. Remember that they can be more
!! or less complex, depending if the reaction rate constant was calculated
!! before the solving or it needs to be calculated in these function.
!! Some reaction parameters like the \ref camp_rxn_arrhenius "arrhenius yield"
!! can also affect the rate.
!!
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

!> \page camp_dev_tutorial_part_4 Dev CAMP: Part 4 - MPI
!!
!!
!! All the calls are done throught the MPI class, which derives from the PartMC original
!! module. This class adapt more of the MPI functions with a safe call, checking
!! the flag MPI_SUCCESS and exiting the program if the call do not success.
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

!> \page camp_dev_tutorial_part_5 Dev CAMP: Part 5 - GPU interface
!!
!!
!!
!! In order to accelerate the CAMP module, a GPU module was developed.
!!
!! This module is developed in the CUDA language. The idea behind this
!! module is translate the most computational expensive part (Solving
!! system) to the GPU. The user can choose to compute all the module in
!! the CPU (set as default) or enable the GPU execution.
!!
!! At the present, only the Derivative and Jacobian functions has a GPU
!! version. However, our team is currently investigating alternatives for
!! a complete GPU ODE solver. ODE solvers, (including CVODE) are releasing
!! increasingly more options to compute expensive operations in GPU,
!! related mostly with matrix operations and linear solving (search cuBLAS and
!! cuSolver for more information). In the future, more GPU solving parts
!! will be implemented to reduce the high computational time of big experiments.
!!
!! ## Parallelization strategy ##
!!
!! The parallelization strategy can be translated briefly into the
!! following question: How we divide the work?
!!
!! Reminding the concepts learned in this tutorial, both Derivative and Jacobian
!! solves a reaction by computing the rate of change \f$r_j\f$ for each
!! participating species \f$y_i\f$. The resulting rate of change is added to the
!! next iteration specie concentration. In this way, each reaction represents
!! a contribution to the related species, and can increase or decrease its value.
!!
!! This concept means that there are no relation with each reaction till we
!! sum the resulting rates of change. In other words, the calculations done
!! in each reaction before the sum can be computed in parallel.
!!
!! So, in conclusion and in the code meanings, instead of a loop iterating
!! all the reactions we will have GPU threads doing this work in parallell.
!!
!! ## GPU pre-solving ##
!!
!! At the start of the solving system, before calling the CVODE solving process,
!! some preparations must to be done:
!!
!!  - Allocate in the GPU all the CAMP variables necessaries for the solving. The
!!  solver_data struct, int_data and float_data arrays, etc.
!!  - Check GPU properties. Check the maximum capacities of the GPU available
!!  to adjust the parallelism accordingly (enough threads, blocks, etc). This
!!  can be done in runtime thanks to the CUDA library.
!!  - Update GPU variables computed in the CPU, like the rate constants and the
!!  Jacobian indices.
!!
!! During the updating of GPU variables, an extra preparation can be done
!! to improve the GPU execution. Taking into account that all the reaction
!! information (int_data and float_data) is stored consecutively in memory and
!! the GPU threads acces simultaneously each reaction first value, we have to
!! set all the first reaction values consecutively in memory. In other words,
!! instead of having in the firsts values of the int_data array the values
!! for the first reaction, we will have as first value the first value of
!! the first reaction, as second the first value of the second reaction, and
!! so on.
!!
!! \image html reverse_rxn_data.jpg
!!
!! Note that this optimization will change the way that we read the data.
!! We name that implementation as "reverse optimization" to refer to it later.
!!
!! ## Derivative and Jacobian ##
!!
!! To launch a kernel in CUDA is necessary define the number of threads and
!! blocks required by the program. For our case, the number of threads is
!! directly the number of reactions, meanwhile the number of blocks can be
!! the minimum (assign the maximum number of threads per block).
!!
!! \code{.c}
!! int n_threads = n_reactions
!! int n_blocks = (n_threads + max_threads_per_block - 1) / max_threads_per_block);
!! \endcode
!!
!! After this, we can launch a kernel passing as input parameter all the pointers
!! to the data previously stored during the GPU pre-solving. Now the kernel function
!! must realize two tasks: divide the work for each thread and call the appropiate
!! Derivative or Jacobian function. Note that starting from this step, all the
!! code will be computed in parallel for all the threads enable.
!!
!! Divide the work translates to assign to each thread a diffent reaction pointer
!! for each input data variable. In the same manner that we update that pointer
!! to point the first value of each reaction in the base version reaction loop,
!! now we use the threads unique id as a loop iterator.
!! Note that with the reverse optimization, the first value of each reaction
!! is directly the next value in the reaction parameter arrays. After that,
!! only remains call the appropiate reaction function.
!!
!! \code{.c}
!! int index = blockIdx.x * blockDim.x + threadIdx.x;
!!
!! double *rxn_float_data = &(float_data)[index]);
!! int *rxn_int_data = &(int_data)[index]);
!! double *rxn_rate_constants = &(rate_constants[index*n_reaction_rate_constants]
!!
!! int rxn_type = int_data[0];
!! rxn_int_data = &(int_data[1*n_rxn]);
!! switch (rxn_type){
!!   case ARRHENIUS:
!!      rxn_arrhenius_calc_deriv_contrib(solver_data, deriv_data, rxn_int_data,
!!      rxn_float_data, rxn_rate_constants)
!!  }
!!
!! \endcode

!! division in threads, defines, get n_rxn, atomic, err idk (but this should be named before)


!! ## Multi-cells approach ##

!! Basically:

!! - Amplify state deriv and rate constants array with each indepedent value.
!! - Multiply the number of threads also for the number of cells
!! (apart from the number of reactions).
!! - Assign to the block of reactions of each cell, the corresponding
!! state_array, deriv array and rate constants for that cell.

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

!> \page camp_dev_tutorial_part_6 Dev CAMP: Part 6 - Testing
!!
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_5
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png
!!
!!

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_6 Dev CAMP: Part 7 -
!!
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_5
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png
!!
!!

!! todo: notes from Mario: Ensure to write the part of optimization
!! changes, like this was the speedup of base version and now this optimize
!! this part, etc. NOTE that normally for the thesis people ask if the CPU
!! version was optimized already, thing that I can say NO, because I optimize
!! this, but the sparse atleast is optimized blablabla. Maybe something like
!! "patch notes" related to optimization, but with a more extensive section.
!! Don't know if let this in the tutorial, or MENTION it in the end or somewhere
!! and link to the documentation. With this done, the phd thesis should be
!! more or less copy paste from the documentation.


!! todo: create some section to write current development, like a historial. Example:
!! I defined an issue explaining a problem. OKAY, this problem should be
!! yes or yes in documentation as recent changes or something like that
!! Once finish, if is an optimization work, check the last optimization work
!! time result and compare the improvement to said okay this optimization
!! changes the time measured last time with the last optimization
!! from 1s to 0.8, so 20% improvement so nice.

!! todo: Create sections to fill as much kind of documentation as possible
!! (like the one written in papers, manuals, or technical reports)
!! (maybe state of the art can suits too?)
!!
!!
!!
