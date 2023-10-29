/*
  Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton,
    Roberto Sabatini
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef MAGIC_DIV_TAU_H
#define MAGIC_DIV_TAU_H

#include <fclaw_base.h>   /* Needed for FCLAW_F77_FUNC */

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>


#ifdef __cplusplus
extern "C"
{
#endif

void magic_div_tau_src3(fclaw2d_global_t *glob,
                        fclaw2d_patch_t *patch,
                        int blockno,
                        int patchno,
                        double t,
                        double dt);


#define MAGIC_DIV_TAU_MAPPED FCLAW_F77_FUNC(magic_div_tau_mapped, \
    MAGIC_DIV_TAU_MAPPED)

void MAGIC_DIV_TAU_MAPPED(int *blockno, int* mx, int* my, int* mz, 
                          int* mbc, int* meqn, int* maux, 
                          double* dx, double* dy, double* dz, double* dt, 
                          double xp[], double yp[], double zp[], 
                          double xd[], double yd[], double zd[], 
                          double volumes[], double xrot[], 
                          double yrot[], double zrot[], 
                          double faceareas[], double aux[], double q[]);


#ifdef __cplusplus
}
#endif

#endif /* !MAGIC_DIV_TAU_H */

