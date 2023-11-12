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

#ifndef MF_MODEL_H
#define MF_MODEL_H

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include <cpl.h>

#include "mf_constants.h"
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

/* Result structure of mf_model */
typedef struct {
  char                       *tmp_folder;            /* Temporary directory                                                    */
  cpl_table                  *gdas_before;           /* GDAS file immediately before MDJ-OBS                                   */
  cpl_table                  *gdas_after;            /* GDAS file immediately after  MDJ-OBS                                   */
  cpl_table                  *gdas_interpolate;      /* GDAS interpolate used                                                  */
  cpl_table                  *atm_profile_standard;  /* Atmospheric profile standard. Initial atm_profile                      */
  cpl_table                  *atm_profile_combined;  /* Atmospheric profile merged with GDAS_interpolate and user as reference */
  cpl_table                  *atm_profile_fitted;    /* Atmospheric profile calculated by telluriccorr (After fit)             */
  cpl_table                  *res;                   /* Results table                                                          */
  cpl_table                  *spec;                  /* Spectrum table                                                         */
  cpl_matrix                 *kernel;                /* NULL or telluriccorr user kernel : resampled and normalized            */
  cpl_size                   n_range;                /* Number of ranges                                                       */
  cpl_table                  **spec_out;             /* Array with the cpl_table calculated in each range                      */
  cpl_error_code             *range_status;          /* Array with the error code status    in each range                      */
  const char                 *gdas_src;              /* String to store the GDAS SOURCE                                        */
} mf_model_results;

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Run telluriccorr in temporary folder and return the structure results */
MF_EXPORT mf_model_results * mf_model(
    mf_configuration         *config,                /* Config structure with the parameters of telluriccorr                     */
    const cpl_table          *molecules,             /* Table defining molecules to fit                                          */
    const cpl_propertylist   *header_spec,           /* Header cpl_propertylist of the         input spectrum                    */
    const cpl_table          *spec_telluriccorr,     /* Input Spectrum table                                                     */
    const cpl_table          *inc_wranges,           /* Wavelength inclusion table                                               */
    const cpl_table          *exc_wranges,           /* Wavelength exclusion table                                               */
    const cpl_table          *exc_pranges,           /* Pixel exclusion table                                                    */
    const cpl_propertylist   *header_kernel,         /* Header cpl_propertylist of the kernel spectrum                           */
    const cpl_matrix         *kernel,                /* Input Kernel matrix, one row per total input spectra pixel               */
    const cpl_table          *gdas_user,             /* NULL or user GDAS profile                                                */
    const cpl_table          *atm_profile_standard,  /* NULL or user atm_profile_standard                                        */
    const cpl_table          *atm_profile_combined); /* NULL or atm_profile_combined (atm_profile_standard merged with gdas_user */

/* Deallocate results structure */
MF_EXPORT void mf_model_results_delete(
    mf_model_results         *results);               /* Model result structure                                                   */


CPL_END_DECLS


#endif /* MF_MODEL_H */
