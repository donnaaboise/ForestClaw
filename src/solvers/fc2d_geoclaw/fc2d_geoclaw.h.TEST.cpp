/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_clawpatch_options.h>
#include <fc2d_geoclaw.h>
#include <fc2d_geoclaw_options.h>
#include <fclaw_forestclaw.h>
#include <fclaw_convenience.h>
#include <test.hpp>

TEST_CASE("fc2d_geoclaw_solver_initialize stores two seperate vtables in two seperate globs")
{
	fclaw_global_t* glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
	fclaw_global_t* glob2 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;

	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 1);
	fclaw_global_store_domain(glob1, domain);
	fclaw_global_store_domain(glob2, domain);

	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_options_new(2);
	fclaw_clawpatch_options_store(glob1, clawpatch_opt);
	fclaw_clawpatch_options_store(glob2, clawpatch_opt);

	/* create some empty options structures */
	fclaw_options_t *fclaw_opt1 = FCLAW_ALLOC_ZERO(fclaw_options_t,1);
	fc2d_geoclaw_options_t *geoclaw_opt1 = FCLAW_ALLOC_ZERO(fc2d_geoclaw_options_t,1);
	fclaw_options_store(glob1, fclaw_opt1);
	fc2d_geoclaw_options_store(glob1, geoclaw_opt1);

	fclaw_options_t *fclaw_opt2 = FCLAW_ALLOC_ZERO(fclaw_options_t,1);
	fc2d_geoclaw_options_t *geoclaw_opt2 = FCLAW_ALLOC_ZERO(fc2d_geoclaw_options_t,1);
	fclaw_options_store(glob2, fclaw_opt2);
	fc2d_geoclaw_options_store(glob2, geoclaw_opt2);

	fclaw_vtables_initialize(glob1);
	fc2d_geoclaw_solver_initialize(glob1);

	fclaw_vtables_initialize(glob2);
	fc2d_geoclaw_solver_initialize(glob2);

	CHECK_NE(fc2d_geoclaw_vt(glob1), fc2d_geoclaw_vt(glob2));

	fclaw_domain_destroy(domain);
	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
	FCLAW_FREE(fclaw_opt1);
	FCLAW_FREE(geoclaw_opt1);
	FCLAW_FREE(fclaw_opt2);
	FCLAW_FREE(geoclaw_opt2);
	fclaw_clawpatch_options_destroy(clawpatch_opt);
}

TEST_CASE("fc2d_geoclaw_solver_initialize sets is_set flag")
{
	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 1);
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
	fclaw_global_store_domain(glob, domain);

	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_options_new(2);
	fclaw_clawpatch_options_store(glob, clawpatch_opt);

	/* create some empty options structures */
	fclaw_options_t *fclaw_opt = FCLAW_ALLOC_ZERO(fclaw_options_t,1);
	fc2d_geoclaw_options_t *geoclaw_opt = FCLAW_ALLOC_ZERO(fc2d_geoclaw_options_t,1);
	fclaw_options_store(glob, fclaw_opt);
	fc2d_geoclaw_options_store(glob, geoclaw_opt);

	fclaw_vtables_initialize(glob);
	fc2d_geoclaw_solver_initialize(glob);


	CHECK_UNARY(fc2d_geoclaw_vt(glob)->is_set);

	fclaw_domain_destroy(domain);
	fclaw_global_destroy(glob);
	FCLAW_FREE(fclaw_opt);
	FCLAW_FREE(geoclaw_opt);
	fclaw_clawpatch_options_destroy(clawpatch_opt);
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fc2d_geoclaw_vt fails if not intialized")
{
	fclaw_global_t* glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
	fclaw_global_t* glob2 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;

	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 1);
	fclaw_global_store_domain(glob1, domain);
	fclaw_global_store_domain(glob2, domain);

	CHECK_SC_ABORTED(fc2d_geoclaw_vt(glob1));

	CHECK_SC_ABORTED(fc2d_geoclaw_vt(glob2));

	fclaw_domain_destroy(domain);
	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fc2d_geoclaw_vtable_initialize called twice on a glob")
{
	fclaw_global_t* glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
	fclaw_global_t* glob2 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;

	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 1);
	fclaw_global_store_domain(glob1, domain);
	fclaw_global_store_domain(glob2, domain);

	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_options_new(2);
	fclaw_clawpatch_options_store(glob1, clawpatch_opt);
	fclaw_clawpatch_options_store(glob2, clawpatch_opt);

	/* create some empty options structures */
	fclaw_options_t *fclaw_opt = FCLAW_ALLOC_ZERO(fclaw_options_t,1);
	fc2d_geoclaw_options_t *geoclaw_opt = FCLAW_ALLOC_ZERO(fc2d_geoclaw_options_t,1);
	fclaw_options_store(glob1, fclaw_opt);
	fc2d_geoclaw_options_store(glob1, geoclaw_opt);

	fclaw_options_store(glob2, fclaw_opt);
	fc2d_geoclaw_options_store(glob2, geoclaw_opt);

	fclaw_vtables_initialize(glob1);
	fc2d_geoclaw_solver_initialize(glob1);
	fc2d_geoclaw_solver_initialize(glob1);

	fclaw_vtables_initialize(glob2);
	fc2d_geoclaw_solver_initialize(glob2);
	fc2d_geoclaw_solver_initialize(glob2);

	fclaw_domain_destroy(domain);
	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
	fclaw_clawpatch_options_destroy(clawpatch_opt);
	FCLAW_FREE(fclaw_opt);
	FCLAW_FREE(geoclaw_opt);
}

#endif