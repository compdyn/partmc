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
!! Before starting, it's recommended to have finalized \ref camp_tutorial
!! Boot CAMP tutorial.
!!
!! - \ref camp_dev_tutorial_part_0
!! - \ref camp_dev_tutorial_part_1
!! - \ref camp_dev_tutorial_part_2
!! - \ref camp_dev_tutorial_part_3
!! - \ref camp_dev_tutorial_part_4
!! - \ref camp_dev_tutorial_part_5
!!
!! Model code described in Dev CAMP can be found in
!! \c doc/camp_tutorial/dev_camp.
!!
!! This tutorial focuses in main concepts. Some of the concepts explained can differ
!! in the future if better implementations are found, but the tutorial should be
!! enough to replicate a consistent CAMP skeleton.
!!
!! ## Testing through CMake##
!! CAMP includes multiple tests with different scenarios (testing each reaction type,
!! chemical mechanisms like Carbond Bond 05, etc.). CMake leads this test execution
!! through the CMakeList file. The developer should ensure his implementation passes all these
!! tests available to execute with the CMake command "make test". The Github CAMP project
!! also check these tests on each commit.
!!
!! ## Documenting your contribution ##
!!
!! CAMP generates all his documentation through Doxygen. Ensure to document properly your
!! apportation to CAMP with commentaries inside the code and using Doxygen if necessary.
!! Don't hesitate to contact the main developers for doubts or help in your contribution.
!! Remember also that the project on Github includes multiple issues with different
!! ideas in development, so ensure to take a look before developing something that is
!! already in progress.

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page camp_dev_tutorial_part_0 Dev CAMP: Part 0 - Introduction
!!
!!
!! CAMP can be divided into the following main components:
!!
!! - <b> Input user data </b>: All input data received by CAMP.
!! It includes the user-defined JSON files and the variables coming from
!! the host model such as chemical concentrations and environmental variables
!! like temperature or pressure. The \ref camp_tutorial "Boot CAMP tutorial"
!! teaches how to enter these variables in CAMP.
!!
!! - <b> CAMP_core interface </b>: User interface.
!! Responsible for adapting the input data to efficent data structures.
!! Includes all the files that prepare the data for the solving system.
!!
!! - <b> Solving system </b>: Corresponds to all the files related to the resolution of the
!! chemical mechanism (and ODE system). Includes the solving of the chemical reactions equations and the
!! ODE solver to advance in time to the next time-step. In the default version of CAMP 1.0, the ODE solving is
!! solved iteratively using the Backward Differentiation Formula method from the CVODE package, alongside
!! the SuiteSparse KLU method as the internal linear solver. These algorithms can differ for other versions
!! or enabling some options in CAMP, for example enabling the GPU flag can change the linear solving for
!! performance reasons.
!!
!! The following images describes the general workflow of CAMP.
!!
!! The first image summarize the main CAMP components and interactions. The chemistry API
!! refers to the JSON files, chemistry core to CAMP_core interface and the integrated
!! chemical mechanism and solver part represents the solving system.
!!
!! \image html camp_main_components.jpg
!!
!! In the second image, the preparations, user and host
!! model parts corresponds to the category "Input user data". The CAMP_core interface is
!! inside CAMP, and the solving system is composed by the ODE solver and some calculations
!! from CAMP like solving the system equations \f$f(y)\f$.
!!
!! \image html preparations_camp.jpg
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
!! ## Derivative workflow ##
!!
!! The next image summarize the workflow on computing the Derivative.
!!
!! \image html Deriv_one_cell.jpg
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
!!  - \ref pmc_util "pmc_util": Treat an array of characters as a string_t object.
!!   Used to store string information like species names in chem_spec_data or camp_core.
!!
!! The workflow corresponds to:
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
!! concentrations array (concs_data) are outside of int_data and float_data.
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
!! with the external solving parts (CAMP_core interface and ODE solver). It's
!! recommended to check the documentation generated for this file, including
!! the main functions described during this section and specially their
!! call graphs.
!!
!! The next image is a simple workflow representation to take an overview of the components
!! interactions. The chemistry core represents the camp core interface and the integrated
!! chemical mechanism represents the solving system.
!!
!! \image html camp_and_solving.jpg
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
!! All solver methods require the solution to the system equations (\f$f(y)\f$)
!! in each ODE solver iteration.
!! But implicit methods tipically needs also the system Jacobian (\f$J = df/dy\f$).
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
!! Note that CAMP is not locked to CVODE and can be adapted to any ODE
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
!! These functions will be called in the middle of the CVODE solving process. To give an idea,
!! in the default CVODE version used the Derivative is called a minimum 
!! of 4 times in different parts of the BDF algorithm, and the Jacobian 1 time in the 
!! setup of the KLU Sparse solver. The number of calls differs depending on the input system.
!!  
!! CVODE will pass as a function parameters the current concentration array (concs_data),
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
!! These functions must solve the operations defined at
!! \ref camp_dev_tutorial_part_1. Remember that they can be more
!! or less complex, depending on the dependencies of the rate constants.
!! Maybe all these rate constants can be calculated at the start of these functions
!! or they have a dependance with the current concentration of some species and need
!! to be calculated during the solving.
!!
!! ## Workflow ##
!!
!! The next image summarize the workflow of the solving system. Note this workflow
!! only represent the main interactions and variables and the last CAMP version can differ.
!!
!! \image html camp_solving_system.jpg
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
!! Most of the classes of the CAMP_core interface has some MPI functions to enable
!! send the class object to another MPI thread. The common approach to communicate
!! the class data is pack all the components (e.g. integers and floats) in the same
!! MPI communicator and unpack the data in the receiving node/s.
!!
!! Also, the whole module CAMP can be executed in independent MPI threads. This approach
!! is used when integrating a CAMP into a host model. The host model will replicate
!! all his modules (including CAMP), divide the work to each node and eventually communicate
!! the results obtained in each module between the rest of the nodes.
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
!! A GPU module was developed in the CUDA language to speedup the solving
!! for large chemical systems. The user can choose to compute all the module in
!! the CPU (set as default) or enable the GPU execution.
!!
!! We will explain the idea of this GPU version with an example: The 
!! parallelization of the Derivative function. This was the first function
!! translated to the GPU during the GPU development to test the efficiency of the
!! GPU, and will be serve also as our introduction point to understand the general
!! concept of how to parallellize a chemical module like CAMP in the GPU.
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
!! This concept means that there are no relation with each reaction until we
!! sum the resulting rates of change. In other words, the calculations done
!! in each reaction before the sum can be computed in parallel.
!!
!! So, instead of a loop iterating
!! all the reactions we will have GPU threads doing this work in parallell.
!!
!! \image html Deriv_GPU_one_cell.jpg
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
!!  - Update GPU variables computed in the CPU, like the rate constants.
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
!! In the next image we consider that the reaction data is organized in
!! the form of a matrix, where each row represent a different reaction and the
!! columns the parameters for the reaction. With this implementation, the threads
!! access the first parameter of each reaction simultaneously. In the original CPU
!! data structure, these access are not consecutive in memory. So, we apply
!! this optimization which is similar to invert the original "matrix": Now,
!! the first parameters accessed are consecutive in memory, improving the memory
!! access. We call this optimization "reverse matrix" to reference it later.
!!
!! \image html reverse_matrix.jpg
!!
!! ## Derivative##
!!
!! To launch a kernel in CUDA is necessary define the number of threads and
!! blocks required by the program. For our case, the number of threads is
!! directly the number of reactions, meanwhile the number of blocks can be
!! the minimum (assign the maximum number of threads per block).
!!
!! \code{.c}
!! int n_threads = n_reactions;
!! int n_blocks = (n_threads + max_threads_per_block - 1) / max_threads_per_block);
!! \endcode
!!
!! After this, we can launch a kernel passing as input parameter all the pointers
!! to the data previously stored during the GPU pre-solving. Now the kernel function
!! must realize two tasks: divide the work for each thread and call the appropiate
!! Derivative function. Note that starting from this step, all the
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
!!
!! Warning: Take into account the reverse optimization also when accessing the code
!! in the reaction type Derivative functions. E.g. in the CPU version accesing the
!! next reaction parameter is sum 1 in the float_array, but now in the reverse
!! optimization we need to sum the lenght of the "matrix rows" (n_rxn a.k.a number of
!! reactions.
!!
!! Also notice that in this reaction type Derivative functions we should add the
!! contribution of each reaction type to the Derivative throught the CUDA instruction
!! "atomicAdd", since the reactions can affect the same specie and the threads will
!! try to access simultaneously the same data.
!!
!! ## Multi-cells in GPU ##
!!
!! In CAMP, the user can compute simultaneously multiple chemical mechanisms (including
!! in the GPU parts). The user only needs to send the number of cells (a.k.a chemical
!! mechanisms) present in the concentration array. Then, this array will have a size of
!! n_reactions * n_cells. Meanwhile in the "original one-cell" implementation the user/host
!! model sends the mechanism (e.g Carbond Bond 05 mechanism) one per one cell, the user
!! group the data in the same data structure. This also translates that the ODE solver
!! will solve the mechanism at once, instead of reinitializing the solving per each cell,
!! improving the general performance.
!!
!! \image html camp_multicells_def.jpg
!!
!! In the GPU, we need to take into account these n_cells by:
!!
!! - Multiply the number of threads also for the number of cells (apart from the number
!! of reactions). Notice that n_cells is set to 1 per default, so 
!! if the user didn't specify n_cells (which mean he don't want
!! to compute multi-cells), the multiplication per n_cells won't affect the computation.
!!
!! \code{.c}
!! int n_threads = n_reactions * n_cells;
!! int n_blocks = (n_threads + max_threads_per_block - 1) / max_threads_per_block);
!! \endcode 
!!
!! - Recalculate sizes on the GPU. Now the concentration array is multiplied by n_cells.
!! The rate constants that depends on environmental data will increase also by
!! the number of cells, since each cell can have different temperature for example (and
!! consider also that each reaction type can need more or less rate constants).
!! The reaction data will remain with the same size, since we are applying the same
!! chemical mechanism with the same coefficients to all the cells.
!!
!! \code{.c}
!! int index = blockIdx.x * blockDim.x + threadIdx.x;
!! int i_cell=index/n_rxn;
!! int i_rxn=index%n_rxn;
!!
!! //Get indices of each reaction
!! double *rxn_float_data = &(float_data)[i_rxn]);
!! int *rxn_int_data = &(int_data)[i_rxn]);
!! //Get indices for concentrations
!! double *rxn_concs = &( concs[n_concs_cell*i_cell]);
!! //Get indices for rates
!! double *rxn_rate_constants = &(rate_constants[n_rate_constants_cell*i_cell+
!! n_rate_constants_reaction[i_rxn]]
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
!!
!! <hr>
!! <b> < Previous: </b> \ref camp_dev_tutorial_part_4
!! \image{inline} html icon_trail.png
!! \ref camp_dev_tutorial "Index"
!! \image{inline} html icon_trail.png