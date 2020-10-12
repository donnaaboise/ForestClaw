/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef MGTEST_DIAGNOSTICS_H
#define MGTEST_DIAGNOSTICS_H

#include <fclaw2d_include_all.h>

#include <fc2d_multigrid.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


/* --------------------------- Problem dependent functions -----------------------------*/

void mgtest_diagnostics_initialize(fclaw2d_global_t *glob,
                                   void **acc_patch);


void mgtest_diagnostics_reset(fclaw2d_global_t *glob,
                              void* patch_acc);

void mgtest_diagnostics_gather(fclaw2d_global_t *glob,
                               void* patch_acc,
                               int init_flag);

void mgtest_diagnostics_finalize(fclaw2d_global_t *glob,
                                 void** patch_acc);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
