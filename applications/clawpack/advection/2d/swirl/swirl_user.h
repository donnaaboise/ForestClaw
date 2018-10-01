/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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

#ifndef SWIRL_USER_H
#define SWIRL_USER_H

#include <fc2d_cudaclaw.h> /* Needed for cuda_rpn typedef */

#include <fclaw2d_include_all.h>


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    double period;
    int claw_version;
    int cuda;
    int is_registered;

} user_options_t;

void swirl_link_solvers(fclaw2d_global_t *glob);

void swirl_problem_setup(fclaw2d_global_t* glob);

/* ------------------------------------- Options ---------------------------------------*/
user_options_t* swirl_options_register (fclaw_app_t * app,
                                        const char *configfile);

void swirl_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* swirl_get_options(fclaw2d_global_t* glob);


/* --------------------------------------- Cuda ----------------------------------------*/

void swirl_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2);

/* ------------------------------------ Fortran ----------------------------------------*/
#define SWIRL_SETPROB FCLAW_F77_FUNC(swirl_setprob, SWIRL_SETPROB)
void SWIRL_SETPROB(double* tperiod);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
