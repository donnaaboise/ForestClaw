/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "hemisphere_user.h"

static void *
hemisphere_register (user_options_t* user_opt, sc_options_t * opt)
{
    sc_options_add_int (opt, 0, "example", &user_opt->example, 1,
                        "[user] Example number");

    sc_options_add_int (opt, 0, "mapping", &user_opt->mapping, 0,
                        "[user] 0 for five-patch hemisphere, "    \
                        "1 for pillowsphere [0]");

    sc_options_add_double (opt, 0, "alpha", &user_opt->alpha, 0.4,
                           "Ratio of outer square to inner square [0.4]");

    sc_options_add_double (opt, 0, "revs_per_sec", &user_opt->revs_per_second, 0.5,
                           "Revolutions per second [0.5]");

    sc_options_add_int (opt, 0, "claw-version", &user_opt->claw_version, 4,
                           "Clawpack_version (4 or 5) [4]");

    user_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
hemisphere_check (user_options_t *user_opt)
{
    if (user_opt->example < 0 || user_opt->example > 1) {
        fclaw_global_essentialf ("Option --user:example must be 0 or 1\n");
        return FCLAW_EXIT_ERROR;
    }

    return FCLAW_NOEXIT;
}

static void
hemisphere_destroy (user_options_t *user)
{
    /* Nothing to destroy */
}


/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (user_options_t*) package;

    return hemisphere_register(user,opt);
}

static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;

    return hemisphere_check(user);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (user_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    hemisphere_destroy (user);

    FCLAW_FREE (user);
}

static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    NULL,
    options_check,
    options_destroy
};

/* ------------- User options access functions --------------------- */

user_options_t* hemisphere_options_register (fclaw_app_t * app,
                                             const char *configfile)
{
    user_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (user_options_t, 1);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"user",user);
    return user;
}

void hemisphere_options_store (fclaw_global_t* glob, user_options_t* user)
{
    fclaw_global_options_store(glob, "user", user);
}

const user_options_t* hemisphere_get_options(fclaw_global_t* glob)
{
    return (user_options_t*) fclaw_global_get_options(glob, "user");
}

