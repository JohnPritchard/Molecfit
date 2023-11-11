/*
 * This file is part of the ESO Telluric Correction Library
 * Copyright (C) 2001-2018 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef MF_FIT_H
#define MF_FIT_H

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include <cpl.h>

#include "mf_constants.h"
#include "mf_parameters.h"
#include "mf_io.h"

CPL_BEGIN_DECLS

/*----------------------------------------------------------------------------*/
/**
 *                 Typedefs: Enumeration types
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Defines
 */
/*----------------------------------------------------------------------------*/

/* Maximum number of iterations per MF_FIT run */
#define MF_FIT_ITER_MAX               100

/*----------------------------------------------------------------------------*/
/**
 *                 Global variables
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Macros
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Typedefs: Structured types
 */
/*----------------------------------------------------------------------------*/

/* Results structure for fit */
typedef struct {
  int    mpfit_calls;   /* Number of MPFIT fitting function evaluation calls                                       */
  int    lblrtm_calls;  /* Number of LBLRTM calls (wavenumber-restricted sub-spectra are not counted individually) */
  double orignorm;      /* Starting value of chi^2                                                                 */
  int    niter;         /* Number of iterations                                                                    */
  int    status;        /* Fitting status code                                                                     */
  int    npar;          /* Total number of parameters                                                              */
  int    nfunc;         /* Number of residuals (= num. of data points)                                             */
  double *xerror;       /* Final parameter uncertainties (1-sigma) npar-vector, or 0 if not desired                */
} mf_fit_results;

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Execute fit */
MF_INTERNAL cpl_error_code mf_fit(
    const cpl_size           run_execution,
    mf_io_lblrtm_config      *lblrtm_config ,
    mf_io_lnfl_config        *lnfl_config,
    cpl_table                *spec,
    cpl_matrix               *kernel,
    cpl_table                *atm_profile,
    cpl_table                **atm_profile_out,
    mf_parameters            *params,
    mf_fit_results           *results,
    cpl_table                **spec_out,
    cpl_error_code           *range_status);

/* Run fit of telluriccorr in a temporary folder */
MF_INTERNAL mf_fit_results * mf_fit_results_create(void);

/* Delete fit results structure */
MF_INTERNAL void mf_fit_results_delete(
    mf_fit_results           *results);

CPL_END_DECLS


#endif /* MF_FIT_H */
