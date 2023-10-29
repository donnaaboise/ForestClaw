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

#include "fc3d_clawpack46.h"
#include "fc3d_clawpack46_options.h"
#include "fc3d_clawpack46_fort.h"

#include "magic_div_tau.h"

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>

#include <fclaw_pointer_map.h>

#include <fclaw3dx_clawpatch.h>

#include <fclaw3dx_clawpatch_options.h>
#include <fclaw3dx_clawpatch_fort.h>

#include <fclaw3d_metric.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_options.h>

void magic_div_tau_src3(fclaw2d_global_t *glob,
                        fclaw2d_patch_t *patch,
                        int blockno,
                        int patchno,
                        double t,
                        double dt)
{
    //fc3d_clawpack46_vtable_t*  claw46_vt = fc3d_clawpack46_vt(glob);

    int mx,my,mz, mbc;
    double xlower,ylower,zlower,dx,dy,dz;
    fclaw3dx_clawpatch_grid_data(glob,patch, &mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);

    double *q;
    int meqn;
    fclaw3dx_clawpatch_soln_data(glob,patch,&q,&meqn);

    double *aux;
    int maux;
    fclaw3dx_clawpatch_aux_data(glob,patch,&aux,&maux);

    double *xp,*yp,*zp;
    double *xd,*yd,*zd;
    double *volume, *facearea;
    fclaw3d_metric_patch_mesh_data(glob,patch,
                                   &xp,&yp,&zp,&xd,&yd,&zd,
                                   &volume,&facearea);

    double *xrot, *yrot, *zrot;
    fclaw3d_metric_patch_basis(glob,patch,&xrot,&yrot,&zrot);

    MAGIC_DIV_TAU_MAPPED(&blockno, &mx, &my, &mz, &mbc, &meqn, 
    &maux, &dx, &dy, &dz, &dt, xp,yp,zp, xd,yd,zd, 
    volume, xrot,yrot,zrot, facearea, aux, q);
}
