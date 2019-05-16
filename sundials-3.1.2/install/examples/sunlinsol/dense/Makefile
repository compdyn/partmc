# -*- mode: makefile -*-
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
# -----------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# -----------------------------------------------------------------
# Makefile for  serial examples
#
# This file is generated from a template using variables
# set at configuration time. It can be used as a template for
# other user Makefiles.
# -----------------------------------------------------------------

SHELL = sh

prefix       = /home/Earth/cguzman/Downloads/sundials-3.1.2/install
exec_prefix  = /home/Earth/cguzman/Downloads/sundials-3.1.2/install
includedir   = /home/Earth/cguzman/Downloads/sundials-3.1.2/install/include
libdir       = /home/Earth/cguzman/Downloads/sundials-3.1.2/install/lib

CPP      = /usr/bin/cc
CPPFLAGS = -g -pg
CC       = /usr/bin/cc
CFLAGS   = -g -pg
LDFLAGS  = -g -pg
LIBS     =  -lm /usr/lib64/librt.so -lsundials_sunmatrixdense

LINKFLAGS = -Wl,-rpath,/home/Earth/cguzman/Downloads/sundials-3.1.2/install/lib

# -----------------------------------------------------------------------------------------

LIBRARIES_LAPACK = -lsundials_sunlinsollapackdense -lsundials_sunlinsollapackband  
LINKFLAGS_LAPACK = ${LINKFLAGS}::

INCLUDES_KLU  = /include
LIBRARIES_KLU = -lsundials_sunlinsolklu 
LINKFLAGS_KLU = ${LINKFLAGS}:/lib

INCLUDES_SLUMT  = 
LIBRARIES_SLUMT = -lsundials_sunlinsolsuperlumt  
LINKFLAGS_SLUMT = ${LINKFLAGS}::

TMP_INCS  = ${includedir} ${INCLUDES_SLUMT} ${INCLUDES_KLU}
INCLUDES  = $(addprefix -I, ${TMP_INCS})
LIBRARIES = -lsundials_sunlinsoldense -lsundials_nvecserial ${LIBS}

# -----------------------------------------------------------------------------------------

EXAMPLES =  test_sunlinsol_dense
EXAMPLES_DEPENDENCIES =  test_sunlinsol sundials_nvector sundials_matrix sundials_linearsolver

OBJECTS = ${EXAMPLES:=.o}
OBJECTS_DEPENDENCIES = ${EXAMPLES_DEPENDENCIES:=.o}

# LAPACK Examples
EXAMPLES_BL = 
EXAMPLES_DEPENDENCIES_BL = 

OBJECTS_BL = ${EXAMPLES_BL:=.o}
OBJECTS_DEPENDENCIES_BL = ${EXAMPLES_DEPENDENCIES_BL:=.o}

# KLU Examples
EXAMPLES_KLU = 
EXAMPLES_DEPENDENCIES_KLU = 

OBJECTS_KLU = ${EXAMPLES_KLU:=.o}
OBJECTS_DEPENDENCIES_KLU = ${EXAMPLES_DEPENDENCIES_KLU:=.o}

# SuperLU_MT Examples
EXAMPLES_SLUMT = 
EXAMPLES_DEPENDENCIES_SLUMT = 

OBJECTS_SLUMT = ${EXAMPLES_SLUMT:=.o}
OBJECTS_DEPENDENCIES_SLUMT = ${EXAMPLES_DEPENDENCIES_SLUMT:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CPPFLAGS} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------

all: examples examples_bl examples_klu examples_slumt

examples: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

examples_bl: ${OBJECTS_BL}
	@for i in ${EXAMPLES_BL} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_BL} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_LAPACK} ${LINKFLAGS_LAPACK}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_BL} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_LAPACK} ${LINKFLAGS_LAPACK} ; \
	done

examples_klu: ${OBJECTS_KLU}
	@for i in ${EXAMPLES_KLU} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_KLU} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_KLU} ${LINKFLAGS_KLU}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_KLU} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_KLU} ${LINKFLAGS_KLU} ; \
	done

examples_slumt: ${OBJECTS_SLUMT}
	@for i in ${EXAMPLES_SLUMT} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_SLUMT} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_SLUMT} ${LINKFLAGS_SLUMT}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_SLUMT} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_SLUMT} ${LINKFLAGS_SLUMT} ; \
	done

${OBJECTS}: ${OBJECTS_DEPENDENCIES}

${OBJECTS_BL}: ${OBJECTS_DEPENDENCIES_BL}

${OBJECTS_KLU}: ${OBJECTS_DEPENDENCIES_KLU}

${OBJECTS_SLUMT}: ${OBJECTS_DEPENDENCIES_SLUMT}

clean:
	rm -f ${OBJECTS_DEPENDENCIES} ${OBJECTS_DEPENDENCIES_BL} ${OBJECTS_DEPENDENCIES_KLU} ${OBJECTS_DEPENDENCIES_SLUMT}
	rm -f ${OBJECTS} ${OBJECTS_BL} ${OBJECTS_KLU} ${OBJECTS_SLUMT}
	rm -f ${EXAMPLES} ${EXAMPLES_BL} ${EXAMPLES_KLU} ${EXAMPLES_SLUMT}

# -----------------------------------------------------------------------------------------
