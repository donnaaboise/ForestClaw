/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#define PATCH_DIM 3
#define REFINE_DIM 2
#define fclaw2d_clawpatch_set_refinement_criteria fclaw3dx_clawpatch_set_refinement_criteria
#define fclaw2d_clawpatch_get_refinement_criteria fclaw3dx_clawpatch_get_refinement_criteria
//options.h
#define fclaw2d_clawpatch_options_store fclaw3dx_clawpatch_options_store
#define fclaw2d_clawpatch_get_options fclaw3dx_clawpatch_get_options
#define fclaw2d_clawpatch_options_t fclaw3dx_clawpatch_options_t
#define fclaw2d_clawpatch_options_register fclaw3dx_clawpatch_options_register
//clwaptch.h
#define fclaw2d_clawpatch_vtable_initialize fclaw3dx_clawpatch_vtable_initialize
#define fclaw2d_clawpatch_time_sync_new fclaw3dx_clawpatch_time_sync_new
#define fclaw2d_clawpatch_time_sync_setup fclaw3dx_clawpatch_time_sync_setup
#define fclaw2d_clawpatch_time_sync_delete fclaw3dx_clawpatch_time_sync_delete
#define fclaw2d_clawpatch_timesync_data fclaw3dx_clawpatch_timesync_data
#define fclaw2d_clawpatch_grid_data fclaw3dx_clawpatch_grid_data
#define fclaw2d_clawpatch_soln_data fclaw3dx_clawpatch_soln_data
#define fclaw2d_clawpatch_get_q fclaw3dx_clawpatch_get_q
#define fclaw2d_clawpatch_vt fclaw3dx_clawpatch_vt
#define fclaw2d_clawpatch_diagnostics_vtable_initialize fclaw3dx_clawpatch_diagnostics_vtable_initialize
#define fclaw2d_clawpatch_pillow_vtable_initialize fclaw3dx_clawpatch_pillow_vtable_initialize
#define fclaw2d_clawpatch_diagnostics_error_default fclaw3dx_clawpatch_diagnostics_error_default
#define fclaw2d_clawpatch_diagnostics_cons_default fclaw3dx_clawpatch_diagnostics_cons_default
#define fclaw2d_clawpatch_time_sync_pack_registers  fclaw3dx_clawpatch_time_sync_pack_registers
#define fclaw2d_clawpatch_face_transformation_intra fclaw3dx_clawpatch_face_transformation_intra
#define cb_clawpatch_output_ascii fclaw3dx_clawpatch_output_ascii_cb
#define fclaw2d_clawpatch_time_header_ascii fclaw3dx_clawpatch_time_header_ascii
#define fclaw2d_clawpatch_face_transformation fclaw3dx_clawpatch_face_transformation
#define fclaw2d_clawpatch_transform_init_data fclaw3dx_clawpatch_transform_init_data
#define fclaw2d_clawpatch_time_sync_reset fclaw3dx_clawpatch_time_sync_reset
#define fclaw2d_clawpatch_time_sync_samesize fclaw3dx_clawpatch_time_sync_samesize
#define fclaw2d_clawpatch_time_sync_f2c fclaw3dx_clawpatch_time_sync_f2c
#define fclaw2d_clawpatch_diagnostics_vtable_initialize fclaw3dx_clawpatch_diagnostics_vtable_initialize

#define fclaw2d_clawpatch_t fclaw3dx_clawpatch_t
#define fclaw2d_clawpatch_vtable_t fclaw3dx_clawpatch_vtable_t
#define fclaw2d_clawpatch_registers_t fclaw3dx_clawpatch_registers_t
#define fclaw2d_vtk_patch_data_t  fclaw3dx_vtk_patch_data_t
#define fclaw2d_clawpatch_packmode_t fclaw3dx_clawpatch_packmode_t