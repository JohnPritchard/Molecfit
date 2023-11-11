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

#ifndef MF_CALCTRANS_H
#define MF_CALCTRANS_H

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include <cpl.h>

#include "mf_constants.h"
#include "mf_parameters.h"
#include "mf_configuration.h"

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

#define MF_CALCTRANS_MINIMUN_OBS_FLUX       1e-2    /* Minimum transmission criterion for good quality of telluric absorption correction */

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

/* Result structure of mf_calctrans_lblrtm */
typedef struct {
  char                       *tmp_folder;                   /* Temporary directory                                                 */
  cpl_size                   n_range;                       /* Number of ranges                                                    */
  cpl_table                  **spec_out;                    /* Array with the cpl_table calculated in each range                   */
  cpl_error_code             *range_status;                 /* Array with the error_code status    in each range                   */
} mf_calctrans_lblrtm_results;

/* Result structure of mf_calctrans_convolution */
typedef struct {
  cpl_table                  *spec_telluriccorr_format;     /* Spectrum table telluric correction                                  */
  cpl_matrix                 *kernel_resampled_normalized;  /* If user kernel provide -> internal kernel [resample and normalized] */
} mf_calctrans_convolution_results;

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Run calctrans-lblrtm execution in a temporary folder and return the structure results */
MF_EXPORT mf_calctrans_lblrtm_results * mf_calctrans_lblrtm(
    mf_configuration            *config,               /* Config structure with the parameters of telluriccorr                */
    const cpl_table             *spec_telluriccorr,    /* Input Spectrum table                                                */
    const cpl_table             *molecules,            /* Table defining molecules to fit                                     */
    double                      wl_start,              /* Lower wavelength override for LNFL                                  */
    double                      wl_end,                /* Upper wavelength override for LNFL                                  */
    cpl_table                   *atm_parameters,       /* Input:  Atmospheric profile           of mf_model(...)              */
    cpl_table                   *best_fit_parameters); /* Input:  Best fitting model parameters of mf_model(...)              */

/* Run calctrans-convolution in the same temporary folder as calctrans-convolution
   and return the structure results with the input spectrum corrected inside */
MF_EXPORT mf_calctrans_convolution_results * mf_calctrans_convolution(
    mf_parameters_config        *config,               /* Config structure with the parameters of telluriccorr                */
    mf_io_lnfl_config           *lnfl_config,
    mf_calctrans_lblrtm_results *lblrtm_results,       /* Results structure obtained in the calctrans_lblrtm execution        */
    const cpl_propertylist      *header_spec,          /* Header cpl_propertylist of the         input spectrum               */
    const cpl_table             *spec_telluriccorr,    /* Input Spectrum table                                                */
    const cpl_propertylist      *header_kernel,        /* Header cpl_propertylist of the kernel spectrum                      */
    const cpl_matrix            *kernel,               /* Input Kernel matrix, one row per total input spectra pixel          */
    double                      wl_start,              /* Lower wavelength override for LNFL                                  */
    double                      wl_end,                /* Upper wavelength override for LNFL                                  */
    cpl_table                   *best_fit_parameters); /* Input:  Best fitting model parameters of mf_model(...)              */

/* Cleanup the mf_calctrans_lblrtm_results structure */
MF_EXPORT void mf_calctrans_lblrtm_results_delete(
    mf_calctrans_lblrtm_results *results);             /* Calctrans-lblrtm result structure                                   */

/* Cleanup the mf_calctrans_convolution_results structure */
MF_EXPORT void mf_calctrans_convolution_results_delete(
    mf_calctrans_convolution_results *results);        /* Calctrans-convolution result structure                              */

/* Corrects the observed flux (and error if present) for telluric absorption */
MF_EXPORT cpl_error_code mf_calctrans_correct_spectrum(
    cpl_table                *spec,                    /* in/out: CPL table with observed flux and modeled transmission curve */
    const cpl_boolean        trans);                   /* Transmission / Emission                                             */


CPL_END_DECLS


#endif /* MF_CALCTRANS_H */
