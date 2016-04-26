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

#ifndef FC2D_GEOCLAW_H
#define FC2D_GEOCLAW_H

#include <fclaw2d_forestclaw.h>
#include <fclaw_package.h>

#include "fc2d_geoclaw_options.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


typedef void (*fc2d_geoclaw_setprob_t)();

typedef void (*fc2d_geoclaw_bc2_t)(const int* meqn, const int* mbc,
                                     const int* mx, const int* my,
                                     const double* xlower, const double* ylower,
                                     const double* dx, const double* dy,
                                     const double q[], const int* maux,
                                     const double aux[], const double* t,
                                     const double* dt, const int mthbc[]);

typedef  void (*fc2d_geoclaw_qinit_t)(const int* meqn,const int* mbc,
                                        const int* mx, const int* my,
                                        const double* xlower, const double* ylower,
                                        const double* dx, const double* dy,
                                        double q[], const int* maux, double aux[]);

typedef void (*fc2d_geoclaw_setaux_t)(const int* mbc,
                                        const int* mx, const int* my,
                                        const double* xlower, const double* ylower,
                                        const double* dx, const double* dy,
                                        const int* maux, double aux[]);

typedef void (*fc2d_geoclaw_b4step2_t)(const int* mbc,
                                         const int* mx, const int* my, const int* meqn,
                                         double q[], const double* xlower,
                                         const double* ylower,
                                         const double* dx, const double* dy,
                                         const double* t, const double* dt,
                                         const int* maux, double aux[]);

typedef void (*fc2d_geoclaw_src2_t)(const int* meqn,
                                      const int* mbc, const int* mx,const int* my,
                                      const double* xlower, const double* ylower,
                                      const double* dx, const double* dy, double q[],
                                      const int* maux, double aux[], const double* t,
                                      const double* dt);


typedef void (*fc2d_geoclaw_rpn2_t)(const int* ixy,const int* maxm, const int* meqn,
                                      const int* mwaves, const int* maux,
                                      const int* mbc,const int* mx,
                                      double ql[], double qr[], double auxl[], double auxr[],
                                      double fwave[], double s[],double amdq[], double apdq[]);

typedef void (*fc2d_geoclaw_rpt2_t)(const int* ixy, const int* imp, const int* maxm, const int* meqn,
                                       const int* mwaves, const int* maux, const int* mbc,const int* mx,
                                       double ql[], double qr[], double aux1[], double aux2[],
                                       double aux3[],  double asdq[],
                                       double bmasdq[], double bpasdq[]);


typedef void (*fc2d_geoclaw_flux2_t)(const int* ixy,const int* maxm, const int* meqn,
                                        const int* maux,const int* mbc,const int* mx,
                                        double q1d[], double dtdx1d[],
                                        double aux1[], double aux2[], double aux3[],
                                        double faddm[],double faddp[], double gaddm[],
                                        double gaddp[],double cfl1d[], double wave[],
                                        double s[], double amdq[],double apdq[],double cqxx[],
                                        double bmasdq[], double bpasdq[],
                                        fc2d_geoclaw_rpn2_t rpn2,
                                        fc2d_geoclaw_rpt2_t rpt2);

typedef void (*fc2d_geoclaw_fluxfun_t)(const int* meqn, double q[], double aux[],
                                          double fq[]);


typedef struct fc2d_geoclaw_vtable
{
    fc2d_geoclaw_setprob_t setprob;
    fc2d_geoclaw_bc2_t bc2;
    fc2d_geoclaw_qinit_t qinit;
    fc2d_geoclaw_setaux_t setaux;
    fc2d_geoclaw_b4step2_t b4step2;
    fc2d_geoclaw_src2_t src2;
    fc2d_geoclaw_rpn2_t rpn2;
    fc2d_geoclaw_rpt2_t rpt2;
    fc2d_geoclaw_fluxfun_t fluxfun;
} fc2d_geoclaw_vtable_t;


void fc2d_geoclaw_init_vtables(fclaw2d_vtable_t* vt,
                               fc2d_geoclaw_vtable_t* geoclaw_vt);

void fc2d_geoclaw_set_vtables(fclaw2d_domain_t *doamin,
                              fclaw2d_vtable_t *vt,
                              fc2d_geoclaw_vtable_t* geoclaw_vt);

#define GEOCLAW_BC2 FCLAW_F77_FUNC(geoclaw_bc2,GEOCLAW_BC2)
void GEOCLAW_BC2(const int* meqn, const int* mbc,
                   const int* mx, const int* my,
                   const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   const double q[], const int* maux,
                   const double aux[], const double* t,
                   const double* dt, const int mthbc[]);

/* --------------------------------------------------------------------
   Classic routines
   - These are provided only for convenience;  these files are not
   compiled into the library, but will be provided by the user.
   -------------------------------------------------------------------- */

/* Macros for C/Fortran portability */
#define SETPROB FCLAW_F77_FUNC(setprob,SETPROB)
#define GEOCLAW_QINIT   FCLAW_F77_FUNC(geoclaw_qinit,GEOCLAW_QINIT)
#define GEOCLAW_SETAUX  FCLAW_F77_FUNC(geoclaw_setaux,GEOCLAW_SETAUX)
#define GEOCLAW_B4STEP2 FCLAW_F77_FUNC(geoclaw_b4step2,GEOCLAW_B4STEP2)
#define GEOCLAW_SRC2    FCLAW_F77_FUNC(geoclaw_src2,GEOCLAW_SRC2)
#define BC2     FCLAW_F77_FUNC(bc2,BC2)
#define GEOCLAW_RPN2    FCLAW_F77_FUNC(geoclaw_rpn2,GEOCLAW_RPN2)
#define GEOCLAW_RPT2    FCLAW_F77_FUNC(geoclaw_rpt2,GEOCLAW_RPT2)

/* These will be converted to MACROS slowly ... */

/* Specific to geoclaw */
#define GEOCLAW_SET_MODULES   FCLAW_F77_FUNC(geoclaw_set_modules, \
                                             GEOCLAW_SET_MODULES)
void GEOCLAW_SET_MODULES(const int* mwaves_in, const int* mcapa_in,
                         const int mthlim_in[], const int method_in[],
                         const double *ax, const double *bx, const double *ay,
                         const double *by);

void SETPROB();

void GEOCLAW_QINIT(const int* meqn,const int* mbc,
                   const int* mx, const int* my,
                   const double* xlower, const double* ylower,
                   const double* dx, const double* dy,
                   double q[], const int* maux, double aux[]);

void GEOCLAW_SETAUX(const int* mbc,
                    const int* mx, const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const int* maux, double aux[]);

void BC2(const int* meqn, const int* mbc,
         const int* mx, const int* my,
         const double* xlower, const double* ylower,
         const double* dx, const double* dy,
         const double q[], const int* maux,
         const double aux[], const double* t,
         const double* dt, const int mthbc[]);

void GEOCLAW_B4STEP2(const int* mbc,
                     const int* mx, const int* my, const int* meqn,
                     double q[], const double* xlower,
                     const double* ylower,
                     const double* dx, const double* dy,
                     const double* t, const double* dt,
                     const int* maux, double aux[]);

void GEOCLAW_SRC2(const int* meqn,
                  const int* mbc, const int* mx,const int* my,
                  const double* xlower, const double* ylower,
                  const double* dx, const double* dy, double q[],
                  const int* maux, double aux[], const double* t,
                  const double* dt);

/* Riemann solvers */
void GEOCLAW_RPN2(const int* ixy,const int* maxm, const int* meqn,
                  const int* mwaves, const int* maux,
                  const int* mbc,const int* mx,
                  double ql[], double qr[], double auxl[], double auxr[],
                  double wave[], double s[],double amdq[], double apdq[]);

void GEOCLAW_RPT2(const int* ixy, const int* imp, const int* maxm, const int* meqn,
                  const int* mwaves, const int* maux, const int* mbc,const int* mx,
                  double ql[], double qr[], double aux1[], double aux2[],
                  double aux3[],  double asdq[],
                  double bmasdq[], double bpasdq[]);

/* --------------------------------------------------------------------
   Time stepping
   -------------------------------------------------------------------- */

#define GEOCLAW_STEP2_WRAP FCLAW_F77_FUNC(geoclaw_step2_wrap,GEOCLAW_STEP2_WRAP)
void GEOCLAW_STEP2_WRAP(const int* maxm, const int* meqn, const int* maux,
                            const int* mbc, const int method[], const int mthlim[],
                            const int* mcapa, const int* mwaves, const int* mx,
                            const int* my, double qold[], double auxold[],
                            const double* dx, const double* dy, const double* dt,
                            const double* cfl, double work[], const int* mwork,
                            const double* xlower, const double* ylower, const int* level,
                            const double* t, double fp[], double fm[], double gp[],
                            double gm[],
                            fc2d_geoclaw_rpn2_t rpn2,
                            fc2d_geoclaw_rpt2_t rpt2,
                            int block_corner_count[]);

#define GEOCLAW_STEP2 FCLAW_F77_FUNC(geoclaw_step2,GEOCLAW_STEP2)
void GEOCLAW_STEP2(const int* maxm, const int* meqn, const int* maux,
                            const int* mbc, const int* mx,
                            const int* my, double qold[], double aux[],
                            const double* dx, const double* dy, const double* dt,
                            const double* cflgrid, double fm[], double fp[], double gm[],
                            double gp[],
                            fc2d_geoclaw_rpn2_t rpn2,
                            fc2d_geoclaw_rpt2_t rpt2);

#define GEOCLAW_FLUX2 FCLAW_F77_FUNC(geoclaw_flux2,GEOCLAW_FLUX2)
void GEOCLAW_FLUX2(const int* ixy,const int* maxm, const int* meqn,
                      const int* maux,const int* mbc,const int* mx,
                      double q1d[], double dtdx1d[],
                      double aux1[], double aux2[], double aux3[],
                      double faddm[],double faddp[], double gaddm[],
                      double gaddp[],double cfl1d[], double wave[],
                      double s[], double amdq[],double apdq[],double cqxx[],
                      double bmasdq[], double bpasdq[],
                      fc2d_geoclaw_rpn2_t rpn2,fc2d_geoclaw_rpt2_t rpt2);
/*
#define geoclaw_FLUX2FW FCLAW_F77_FUNC(geoclaw_flux2fw,geoclaw_FLUX2FW)
void GEOCLAW_FLUX2FW(const int* ixy,const int* maxm, const int* meqn, //
                        const int* maux,const int* mbc,const int* mx,
                        double q1d[], double dtdx1d[],
                        double aux1[], double aux2[], double aux3[],
                        double faddm[],double faddp[], double gaddm[],
                        double gaddp[],double cfl1d[], double fwave[],
                        double s[], double amdq[],double apdq[],double cqxx[],
                        double bmasdq[], double bpasdq[],
                        fc2d_geoclaw_rpn2_t rpn2,fc2d_geoclaw_rpt2_t rpt2,
                        const int* mwaves, const int* mcapa,
                        int method[], int mthlim[]);*/

#define GEOCLAW_SET_CAPACITY FCLAW_F77_FUNC(geoclaw_set_capacity,GEOCLAW_SET_CAPACITY)
void GEOCLAW_SET_CAPACITY(const int* mx, const int *my, const int *mbc,
                             const double *dx, const double* dy, double area[],
                             const int *mcapa, const int* maux, double aux[]);


#define GEOCLAW_SET_BLOCK FCLAW_F77_FUNC(geoclaw_set_block,GEOCLAW_SET_BLOCK)
void GEOCLAW_SET_BLOCK(int* blockno);

#define FC2D_GEOCLAW_GET_BLOCK FCLAW_F77_FUNC(fc2d_geoclaw_get_block, \
                                                 FC2D_GEOCLAW_GET_BLOCK)
int FC2D_GEOCLAW_GET_BLOCK();


#define GEOCLAW_UNSET_BLOCK FCLAW_F77_FUNC(geoclaw_unset_block, \
                                              GEOCLAW_UNSET_BLOCK)
void GEOCLAW_UNSET_BLOCK();

/************************ Regridding ******************************/
#define GEOCLAW_TAG4REFINEMENT FCLAW_F77_FUNC(geoclaw_tag4refinement,GEOCLAW_TAG4REFINEMENT)
void GEOCLAW_TAG4REFINEMENT(int* mx,int* my, int* mbc, int *meqn,int*maux,
                            double * xlower,double * ylower,
                            double* dx,double* dy, double*t, int* blockno,
                            double q[], double aux[], int* level,int* maxlevel,
                            int* init_flag,int* tag_patch);



/***************************** MINIMAL API ******************************/

void fc2d_geoclaw_register_vtable (fclaw_package_container_t *
                                      pkg_container,
                                      fc2d_geoclaw_options_t *
                                      clawopt);

/* -------------------------------------------------------------------------
   New routines
   ------------------------------------------------------------------------- */
void fc2d_geoclaw_define_auxarray(fclaw2d_domain_t* domain,
                                    fclaw2d_patch_t* this_patch);

void fc2d_geoclaw_aux_data(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t *this_patch,
                              double **aux, int* maux);

int fc2d_geoclaw_get_maux(fclaw2d_domain_t* domain);
void fc2d_geoclaw_maux(fclaw2d_domain_t* domain, int* maux);

void fc2d_geoclaw_register (fclaw_app_t* app, const char *configfile);

void fc2d_geoclaw_package_register(fclaw_app_t* app,
                                      fc2d_geoclaw_options_t* clawopt);

int fc2d_geoclaw_get_package_id (void);

fc2d_geoclaw_options_t* fc2d_geoclaw_get_options(fclaw2d_domain_t *domain);

/* -------------------------------------------------------------------------
   Routines that won't change
   ------------------------------------------------------------------------- */
void fc2d_geoclaw_setup(fclaw2d_domain_t *domain);

void
    fc2d_geoclaw_setprob(fclaw2d_domain_t* domain);

void fc2d_geoclaw_patch_setup(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx);
void
    fc2d_geoclaw_setaux(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx);

void
    fc2d_geoclaw_set_capacity(fclaw2d_domain_t *domain,
                                 fclaw2d_patch_t *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx);

void
    fc2d_geoclaw_qinit(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx);

void
    fc2d_geoclaw_b4step2(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            double t,
                            double dt);

void
    fc2d_geoclaw_bc2(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt,
                        fclaw_bool intersects_bc[],
                        fclaw_bool time_interp);

void
    fc2d_geoclaw_src2(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt);


/* A single step method that advances the solution a single step on a single grid
   using a time step dt determined by the subcycle manager */
double
    fc2d_geoclaw_step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

/* Use this ro return only the right hand side of the clawpack algorithm */
double
    fc2d_geoclaw_step2_rhs(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs);

double
fc2d_geoclaw_update(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt);

int fc2d_geoclaw_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *fine_patches,
                                      int blockno, int patchno);

int fc2d_geoclaw_patch_tag4refinement(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno, int this_patch_idx,
                                      int initflag);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC2D_CLAWPACH5_H */
