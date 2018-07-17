/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * CUDA utilities
 *
*/
/** \file
 * \brief CUDA utilities
*/
#ifndef CUDA_UTIL_H_
#define CUDA_UTIL_H_
#include <cuda.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
  if (err != cudaSuccess) {
    printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
            file, line );
    exit( EXIT_FAILURE );
  }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

#endif
