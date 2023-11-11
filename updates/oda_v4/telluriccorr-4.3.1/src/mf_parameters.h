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

#ifndef MF_PARAMETERS_H
#define MF_PARAMETERS_H

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include <cpl.h>

#include "mf_constants.h"
#include "mf_molecules.h"

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

/* Generic string for the KEYWORDS (cpl_propertylist in FITS headers) for put in front of the parameters */
#define MF_PARAMETERS_CONTEX_DEFAULT                  "ESO DRS MF PARAM"

/*** Telluric Correction parameters for keywords ***/
#define MF_PARAMETERS_KEYWORD_NONE                    "NONE"
#define MF_PARAMETERS_NULL                            "NULL"


/*** Telluric Correction parameters : Directories ***/

#define MF_PARAMETERS_TELLURICCORR_PATH               "TELLURICCORR_PATH"
#define MF_PARAMETERS_TELLURICCORR_PATH_DESC          "Installation directory"
#define MF_PARAMETERS_TELLURICCORR_PATH_ENV           "TELLURICCORRDIR"
#define MF_PARAMETERS_TELLURICCORR_PATH_INIT          TELLURICCORR_BASEDIR

#define MF_PARAMETERS_TELLURICCORR_DATA_PATH          "TELLURICCORR_DATA_PATH"
#define MF_PARAMETERS_TELLURICCORR_DATA_PATH_DESC     "Data directory"
#define MF_PARAMETERS_TELLURICCORR_DATA_PATH_ENV      "TELLURICCORRDIR_DATA"
#define MF_PARAMETERS_TELLURICCORR_DATA_PATH_INIT     TELLURICCORR_DATADIR

#define MF_PARAMETERS_TMP_PATH                        "TMP_PATH"
#define MF_PARAMETERS_TMP_PATH_DESC                   "Temporary directory"
#define MF_PARAMETERS_TMP_PATH_MEMORY                 "/dev/shm"
#define MF_PARAMETERS_TMP_PATH_DEFAULT                "/tmp"

/*#define MF_PARAMETERS_OUTPUT_PATH                     "OUTPUT_PATH"
#define MF_PARAMETERS_OUTPUT_PATH_DESC                "Directory for output files"
#define MF_PARAMETERS_OUTPUT_PATH_DEFAULT             "output"
#define MF_PARAMETERS_OUTPUT_PATH_INIT                MF_PARAMETERS_OUTPUT_PATH_DEFAULT

#define MF_PARAMETERS_OUTPUT_NAME                     "OUTPUT_NAME"
#define MF_PARAMETERS_OUTPUT_NAME_DESC                "Root name of output files"
#define MF_PARAMETERS_OUTPUT_NAME_DEFAULT             "out"
#define MF_PARAMETERS_OUTPUT_NAME_INIT                MF_PARAMETERS_OUTPUT_NAME_DEFAULT */


/*** Telluric Correction parameters : Input parameters ***/

// #define MF_PARAMETERS_OMP_NUM_THREADS                 "OPENMP_THREADS"
// #define MF_PARAMETERS_OMP_NUM_THREADS_ENV             "OMP_NUM_THREADS"
/* #define MF_PARAMETERS_OMP_NUM_THREADS_DESC            "Number of OPENMP threads in the executions.\n"\
                                                      " If doesn't exist OMP_NUM_THREAD environment variable\n"\
                                                      " SERIAL execution by default"
*/

#define MF_PARAMETERS_SILENT_EXTERNAL_BINS            "SILENT_EXTERNAL_BINS"
#define MF_PARAMETERS_SILENT_EXTERNAL_BINS_DESC       "Silent the output of the external binaries"
#define MF_PARAMETERS_SILENT_EXTERNAL_BINS_YES        CPL_TRUE
#define MF_PARAMETERS_SILENT_EXTERNAL_BINS_NO         CPL_FALSE
#define MF_PARAMETERS_SILENT_EXTERNAL_BINS_INIT       MF_PARAMETERS_SILENT_EXTERNAL_BINS_YES

#define MF_PARAMETERS_TRANSMISSION                    "TRANSMISSION"
#define MF_PARAMETERS_TRANSMISSION_DESC               "Type of input spectrum : 0 = Emission(radiance); 1 = Transmission"
#define MF_PARAMETERS_TRANSMISSION_FALSE              CPL_FALSE
#define MF_PARAMETERS_TRANSMISSION_TRUE               CPL_TRUE
#define MF_PARAMETERS_TRANSMISSION_INIT               MF_PARAMETERS_TRANSMISSION_TRUE

#define MF_PARAMETERS_COLUMN_LAMBDA                   "COLUMN_LAMBDA"
#define MF_PARAMETERS_COLUMN_LAMBDA_DESC              "Wavelength column ('NULL' can be used if the file is an image and that\n"\
                                                      " the data are in the primary\n"\
                                                      " (data are given by the FITS header keywords [CRVAL1=wave_ini, CD1_1=step])\n"\
                                                      " If CD1_1 is absent, then the DEPRECATED CDELT1 keyword will be used."
#define MF_PARAMETERS_COLUMN_LAMBDA_NULL              MF_PARAMETERS_NULL
#define MF_PARAMETERS_COLUMN_LAMBDA_DEFAULT           MF_COL_IN_LAMBDA
#define MF_PARAMETERS_COLUMN_LAMBDA_INIT              MF_PARAMETERS_COLUMN_LAMBDA_DEFAULT

#define MF_PARAMETERS_COLUMN_FLUX                     "COLUMN_FLUX"
#define MF_PARAMETERS_COLUMN_FLUX_DESC                "Flux column"
#define MF_PARAMETERS_COLUMN_FLUX_NULL                MF_PARAMETERS_NULL
#define MF_PARAMETERS_COLUMN_FLUX_DEFAULT             MF_COL_IN_FLUX
#define MF_PARAMETERS_COLUMN_FLUX_INIT                MF_PARAMETERS_COLUMN_FLUX_DEFAULT

#define MF_PARAMETERS_COLUMN_DFLUX                    "COLUMN_DFLUX"
#define MF_PARAMETERS_COLUMN_DFLUX_DESC               "Flux error column (Avoided by writing 'NULL') : 1-sigma error on the flux"
#define MF_PARAMETERS_COLUMN_DFLUX_NULL               MF_PARAMETERS_NULL
#define MF_PARAMETERS_COLUMN_DFLUX_DEFAULT            MF_COL_IN_DFLUX
#define MF_PARAMETERS_COLUMN_DFLUX_INIT               MF_PARAMETERS_COLUMN_DFLUX_NULL

#define MF_PARAMETERS_COLUMN_MASK                     "COLUMN_MASK"
#define MF_PARAMETERS_COLUMN_MASK_DESC                "Mask column (Avoided by writing 'NULL') : Indicates if a pixel is invalid"
#define MF_PARAMETERS_COLUMN_MASK_NULL                MF_PARAMETERS_NULL
#define MF_PARAMETERS_COLUMN_MASK_DEFAULT             MF_COL_IN_MASK
#define MF_PARAMETERS_COLUMN_MASK_INIT                MF_PARAMETERS_COLUMN_MASK_NULL

#define MF_PARAMETERS_DEFAULT_ERROR                   "DEFAULT_ERROR"
#define MF_PARAMETERS_DEFAULT_ERROR_DESC              "Default error relative to mean for the case that the error column\n"\
                                                      " is not provided"
#define MF_PARAMETERS_DEFAULT_ERROR_INIT              0.01

#define MF_PARAMETERS_WLG_TO_MICRON                   "WLG_TO_MICRON"
#define MF_PARAMETERS_WLG_TO_MICRON_DESC              "Multiplicative factor applied to the wavelength to express is in micron.\n"\
                                                      " E.g.: if wavelength is given in nm, the value should be 0.001."
#define MF_PARAMETERS_WLG_TO_MICRON_INIT              1.

// Remove this comment in the release (with improvements, TODO : NEW from the original improvements DPB proposal, VALUE = 'VAC_RV')
#define MF_PARAMETERS_WAVELENGTH_FRAME               "WAVELENGTH_FRAME"
#define MF_PARAMETERS_WAVELENGTH_FRAME_DESC          "Wavelength in vacuum                                      = 'VAC'.\n"\
                                                      " Wavelength in air    with the observatory reference frame = 'AIR'.\n"\
													  " Wavelength in air    with another         reference frame = 'AIR_RV'.\n"\
                                                      " Wavelength in vacuum with another         reference frame = 'VAC_RV'.\n"\
                                                      "   (typically the  sun or the barycenter of the solar system).\n"\
                                                      " In the latter case, the radial velocity of the observatory relative\n"\
                                                      "   to the external reference frame must be provided in the parameter obs_RV"

#define MF_PARAMETERS_WAVELENGTH_FRAME_VACUUM        "VAC"
#define MF_PARAMETERS_WAVELENGTH_FRAME_VACUUM_RV     "VAC_RV"
#define MF_PARAMETERS_WAVELENGTH_FRAME_AIR           "AIR"
#define MF_PARAMETERS_WAVELENGTH_FRAME_AIR_RV        "AIR_RV"
#define MF_PARAMETERS_WAVELENGTH_FRAME_INIT          MF_PARAMETERS_WAVELENGTH_FRAME_VACUUM

// Remove this comment in the release (with improvements, TODO : NEW from the original improvements DPB proposal, VALUE = 'VAC_RV')
#define MF_PARAMETERS_OBSERVATORY_ERF_RV_DESC         "The radial velocity of the observatory in km/s\n"\
                                                      " relative to the external reference frame;\n"\
                                                      " It is positive if the distance between the science target and the Earth\n"\
                                                      " increases along the line-of-sight to the science target.\n"\
                                                      " It must be provided if MF_PARAMETERS_WAVELENGTH_FRAME = 'VAC_RV'"
#define MF_PARAMETERS_OBSERVATORY_ERF_RV_KEYWORD      "OBS_ERF_RV_KEY"
#define MF_PARAMETERS_OBSERVATORY_ERF_RV_KEYWORD_DEFAULT "OBS_ERF_RV_KEY"
#define MF_PARAMETERS_OBSERVATORY_ERF_RV_KEYWORD_INIT MF_PARAMETERS_KEYWORD_NONE
#define MF_PARAMETERS_OBSERVATORY_ERF_RV_VALUE        "OBS_ERF_RV_VALUE"
#define MF_PARAMETERS_OBSERVATORY_ERF_RV_VALUE_DESC   "If OBS_ERF_RV_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_OBSERVATORY_ERF_RV_VALUE_INIT   0.


#define MF_PARAMETERS_CLEAN_MODEL_FLUX                "CLEAN_MODEL_FLUX"
#define MF_PARAMETERS_CLEAN_MODEL_FLUX_DESC           "Set model flux to 0 for non-fitted pixels"
#define MF_PARAMETERS_CLEAN_MODEL_FLUX_INIT           CPL_FALSE


/* Precision */
#define MF_PARAMETERS_FTOL                            "FTOL"
#define MF_PARAMETERS_FTOL_DESC                       "Relative chi-square convergence criterion"
#define MF_PARAMETERS_FTOL_MIN                        1e-10
#define MF_PARAMETERS_FTOL_INIT                       MF_PARAMETERS_FTOL_MIN

#define MF_PARAMETERS_XTOL                            "XTOL"
#define MF_PARAMETERS_XTOL_DESC                       "Relative parameter convergence criterion"
#define MF_PARAMETERS_XTOL_MIN                        1e-10
#define MF_PARAMETERS_XTOL_INIT                       MF_PARAMETERS_XTOL_MIN


/* Flux units  */
#define MF_PARAMETERS_FLUX_UNIT                       "FLUX_UNIT"
#define MF_PARAMETERS_FLUX_UNIT_DESC                  "Conversion of fluxes from phot/(s*m2*mum*as2) (emission spectrum only)\n"\
                                                      "    to flux unit of observed spectrum:\n"\
                                                      " 0: phot / (s *  m^2 * mum * as^2) [no conversion]\n"\
                                                      " 1:    W / (     m^2 * mum * as^2)\n"\
                                                      " 2:  erg / (s * cm^2 *   A * as^2)\n"\
                                                      " 3:  mJy / (                 as^2)\n"\
                                                      " For other units, the conversion factor has to be considered\n"\
                                                      "  as constant term of the continuum fit."
#define MF_PARAMETERS_FLUX_UNIT_NO_CONVERSION         0
#define MF_PARAMETERS_FLUX_UNIT_1                     1
#define MF_PARAMETERS_FLUX_UNIT_2                     2
#define MF_PARAMETERS_FLUX_UNIT_3                     3
#define MF_PARAMETERS_FLUX_UNIT_INIT                  MF_PARAMETERS_FLUX_UNIT_NO_CONVERSION


/* Telescope background */

#define MF_PARAMETERS_FIT_TELESCOPE_BACK              "FIT_TELESCOPE_BACKGROUND"
#define MF_PARAMETERS_FIT_TELESCOPE_BACK_DESC         "Fit of telescope background -- 1 = yes; 0 = no (emission spectrum only)"
#define MF_PARAMETERS_FIT_TELESCOPE_BACK_INIT         CPL_TRUE

#define MF_PARAMETERS_TELESCOPE_BACK_CONST            "TELESCOPE_BACKGROUND_CONST"
#define MF_PARAMETERS_TELESCOPE_BACK_CONST_DESC       "Initial value for telescope background fit"
#define MF_PARAMETERS_TELESCOPE_BACK_CONST_MIN        0.
#define MF_PARAMETERS_TELESCOPE_BACK_CONST_MAX        1.
#define MF_PARAMETERS_TELESCOPE_BACK_CONST_INIT       0.1


/* Continuum Modelling. (Note Range Specific) */

#define MF_PARAMETERS_FIT_CONTINUUM                   "FIT_CONTINUUM"
#define MF_PARAMETERS_FIT_CONTINUUM_DESC              "Comma deliminated string of flags (1=true, 0=false) for fitting continuum in specific regions.\n"\
                                                      " If set to NULL, check if the TAG[WAVE_INCLUDE] points to a FITS BINTABLE with column CONT_FIT_FLAG provided"
#define MF_PARAMETERS_FIT_CONTINUUM_INIT              "1"
#define MF_PARAMETERS_FIT_CONTINUUM_INIT_SINGLE_VALUE  CPL_TRUE

#define MF_PARAMETERS_CONTINUUM_N                     "CONTINUUM_N"
#define MF_PARAMETERS_CONTINUUM_N_DESC                "Polynomial order for continuum model for each region. Presented as a comma deliminated string.\n"\
                                                      " If set to NULL, check if the TAG[WAVE_INCLUDE] points to a FITS BINTABLE with column CONT_POLY_ORDER provided"
#define MF_PARAMETERS_CONTINUUM_N_INIT                 "0"
#define MF_PARAMETERS_CONTINUUM_N_INIT_SINGLE_VALUE    0

#define MF_PARAMETERS_CONTINUUM_CONST                 "CONTINUUM_CONST"
#define MF_PARAMETERS_CONTINUUM_CONST_DESC            "Initial constant term for continuum fit (valid for all fit ranges)\n"\
                                                      " [emission spectrum: about 1 for correct flux_unit]"
#define MF_PARAMETERS_CONTINUUM_CONST_INIT              1.


/* Mapping of Regions to chip ids (Note Range Specific) */

#define MF_PARAMETERS_MAP_REGIONS_TO_CHIP              "MAP_REGIONS_TO_CHIP"
#define MF_PARAMETERS_MAP_REGIONS_TO_CHIP_DESC         "Comma deliminated string of chip indices that each range is associated with.\n"\
                                                      " If set to NULL, check if the TAG[WAVE_INCLUDE] points to a FITS BINTABLE with column MAPPED_TO_CHIP provided"
#define MF_PARAMETERS_MAP_REGIONS_TO_CHIP_INIT         "1"


/* EXPERT MODE FLAG (To be used only by people who know what they are doing) */
#define MF_PARAMETERS_EXPERT_MODE                       "EXPERT_MODE"
#define MF_PARAMETERS_EXPERT_MODE_DESC                  "If set to true, will check if TAG[INIT_FIT_PARAMETERS] points to a fits file with a bintable of parameter values to use as initial values for the fitting process"
#define MF_PARAMETERS_EXPERT_MODE_INIT                  CPL_FALSE


/* New implementation with 3 new FITS files in ESOTK : CONTINUUM_POLYNOME, CONTINUUM_LIBRARY, CONTINUUM_MAPPING */

// Remove this comment in the release (with improvements, TODO : No chips in the API/refactoring, one call to Molecfict for each detector)
// A new FITS file in the ESOTK sof file (CONTINUUM_POLYNOME)
// CONTINUUM_POLYNOME : FITS file IMAGE1D, one extension for each polynome (Load in a cpl_vector)
// One input in Telluric Correction : cpl_vector where n_pos=degree and pos0-coef0, pos1-coef1, ...
//                                    If provided replace MF_PARAMETERS_CONTINUUM_N and MF_PARAMETERS_CONTINUUM_CONST.
#define MF_PARAMETERS_CONTINUUM_POLYNOME              "CONTINUUM_POLYNOME"
#define MF_PARAMETERS_CONTINUUM_POLYNOME_DESC         "(CONTINUUM_N +1) coefficients of the polynome representing the continuum, separated by a comma.\n"\
                                                      " If provided, replace CONTINUUM_CONST"
#define MF_PARAMETERS_CONTINUUM_POLYNOME_INIT         "0.,1."

// Remove this comment in the release (with improvements, TODO)
// A new FITS file in the ESOTK sof file (CONTINUUM_LIBRARY)
// CONTINUUM_LIBRARY : FITS file BINTABLE, one extension for each detector (Load in a cpl_table - 2 columns WAVELENGHT and FLUX).
// One input in Telluric Correction : cpl_table with 2 columns WAVELENGHT and FLUX
//                                    If provided, replace MF_PARAMETERS_FIT_CONTINUUM and MF_PARAMETERS_CONTINUUM_N.
#define MF_PARAMETERS_CONTINUUM_SPECTRUM              "CONTINUUM_SPECTRUM"
#define MF_PARAMETERS_CONTINUUM_SPECTRUM_DESC         "Path and name of a spectrum to be used as a continuum.\n"\
                                                      " If provided, the spectrum must be a binary FITS table with 2 columns file\n"\
                                                      " giving respectively wavelength and flux;\n"\
                                                      " The flux is then multiplied by the polynome\n"\
                                                      " with parameters given by CONTINUUM_N and CONTINUUM_CONST"
#define MF_PARAMETERS_CONTINUUM_SPECTRUM_INIT         "none"

// CONTINUUM_MAPPING : FITS file BINTABLE, two columns CONTINUUM_POLYNOME_EXT, CONTINUUM_LIBRARY_EXT. Means row1[ext1_RAW] mapping with -> ext[?,CONTINUUM_POLYNOME], ext[?,CONTINUUM_LIBRARY].

// Remove this comment in the release (TODO in the improvements part, now only define with any impact in the code)
/*#define MF_PARAMETERS_OBSERVATORY_BARY_RV             "CONST_BARY_RV"
#define MF_PARAMETERS_OBSERVATORY_BARY_RV_DESC        "Only used if a continuum spectrum is provided, CONST_BARY_RV is needed to take into account\n"\
                                                      " the radial velocity of the observatory relative to the reference frame of the user-provided continuum.\n"\
                                                      " Hence, it represents the radial velocity (in km/s) of the observatory relative\n"\
                                                      " to the barycenter of the solar system along the line-of-sight to the science target;\n"\
                                                      " It is positive if the distance between the science target and the Earth\n"\
                                                      " increases the line-of-sight to the science target."
#define MF_PARAMETERS_OBSERVATORY_BARY_RV_INIT        0.
*/

/* Wavelength solution */

#define MF_PARAMETERS_FIT_WLC                         "FIT_WLC"
#define MF_PARAMETERS_FIT_WLC_DESC                    "Flags for including regions in wavelength corrections.\n"\
                                                      " If set to NULL, check if the TAG[WAVE_INCLUDE] points to\n"\
                                                      " a FITS BINTABLE with column WLC_FIT_FLAG provided."
#define MF_PARAMETERS_FIT_WLC_INIT                     "0"
#define MF_PARAMETERS_FIT_WLC_INIT_SINGLE_VALUE        CPL_FALSE

#define MF_PARAMETERS_WLC_N                           "WLC_N"
#define MF_PARAMETERS_WLC_N_DESC                      "Polynomial degree of the refined wavelength solution"
#define MF_PARAMETERS_WLC_N_INIT                      1

#define MF_PARAMETERS_WLC_CONST                       "WLC_CONST"
#define MF_PARAMETERS_WLC_CONST_DESC                  "Initial constant term for wavelength adjustment\n"\
                                                      " (shift relative to half wavelength range)"
#define MF_PARAMETERS_WLC_CONST_MIN                   -0.05
#define MF_PARAMETERS_WLC_CONST_MAX					  0.05
#define MF_PARAMETERS_WLC_CONST_INIT                  MF_PARAMETERS_WLC_CONST_MIN

/* New implementation with 2 new FITS files in ESOTK : WLC_POLYNOME, WLC_MAPPING */

// Remove this comment in the release (with improvements, TODO : No chips in the API/refactoring, one call to Molecfict for each detector)
// A new FITS file in the ESOTK sof file (WLC_POLYNOME)
// WLC_POLYNOME : FITS file IMAGE1D, one extension for each polynome (Load in a cpl_vector)
// One input in Telluric Correction : cpl_vector where n_pos=degree and pos0-coef0, pos1-coef1, ...
//                                    If provided replace MF_PARAMETERS_WLC_N and MF_PARAMETERS_WLC_CONST.
#define MF_PARAMETERS_WLC_POLYNOME                    "WLC_POLYNOME"
#define MF_PARAMETERS_WLC_POLYNOME_DESC               "(WLC_N +1) coefficients of the polynome representing the wavelength solution, separated by a comma.\n"\
                                                      " If provided, replace CONT_CONST"
#define MF_PARAMETERS_WLC_POLYNOME_INIT               "0"

// WLC_MAPPING : FITS file BINTABLE, one column WLC_POLYNOME_EXT. Means row1[ext1_RAW] mapping with -> ext[?,WLC_POLYNOME]

// Remove this comment in the release (with improvements, TODO : NEW VALUE = 'MODEL'; 'DATA' equal to the previous behavior)
/* #define MF_PARAMETERS_WLC_REF                         "WLC_REF"
#define MF_PARAMETERS_WLC_REF_DESC                    "Indicates that the reference for the wavelength calibration :\n"\
                                                      " - If it is set to 'DATA',  is the input  data.\n"\
                                                      " - If it is set to 'MODEL', is the output model.\n"\
                                                      " In the 'MODEL' case, the wavelength given in the output spectrum\n"\
                                                      " has been corrected with the result of the fit.\n"\
                                                      " If set to 'MODEL', the degree of the polynomial used\n"\
                                                      " for the wavelength correction should be 0 or at most 1,\n"\
                                                      " except if the inclusion regions cover a large\n"\
                                                      " or well-sampled spectral range of the input spectrum"
#define MF_PARAMETERS_WLC_REF_DATA                    "DATA"
#define MF_PARAMETERS_WLC_REF_MODEL                   "MODEL"
#define MF_PARAMETERS_WLC_REF_INIT                    MF_PARAMETERS_WLC_REF_DATA */


/* Resolution */

#define MF_PARAMETERS_FIT_RES_BOX                     "FIT_RES_BOX"
#define MF_PARAMETERS_FIT_RES_BOX_DESC                "Fit the width of a Boxcar LSF"
#define MF_PARAMETERS_FIT_RES_BOX_INIT                CPL_TRUE

#define MF_PARAMETERS_RES_BOX                         "RES_BOX"
#define MF_PARAMETERS_RES_BOX_DESC                    "Initial value for FWHM of Boxcar rel. to slit width\n"\
                                                      " at the centre of the spectrum"
#define MF_PARAMETERS_RES_BOX_MIN                     0.
#define MF_PARAMETERS_RES_BOX_MAX                     2.
#define MF_PARAMETERS_RES_BOX_INIT                    1.


/* Gaussian Kernel */

#define MF_PARAMETERS_FIT_GAUSS                       "FIT_RES_GAUSS"
#define MF_PARAMETERS_FIT_GAUSS_DESC                  "Fit the FWHM of a Gaussian LSF"
#define MF_PARAMETERS_FIT_GAUSS_INIT                  CPL_TRUE

#define MF_PARAMETERS_RES_GAUSS                       "RES_GAUSS"
#define MF_PARAMETERS_RES_GAUSS_DESC                  "Initial value for FWHM of the Gaussian in pixels\n"\
                                                      " at the centre of the spectrum"
#define MF_PARAMETERS_RES_GAUSS_MIN                   0.
#define MF_PARAMETERS_RES_GAUSS_MAX                   100.
#define MF_PARAMETERS_RES_GAUSS_INIT                  1.


/* Lorentzian Kernel */

#define MF_PARAMETERS_FIT_LORENTZ                     "FIT_RES_LORENTZ"
#define MF_PARAMETERS_FIT_LORENTZ_DESC                "Fit the FWHM of a Lorentzian LSF"
#define MF_PARAMETERS_FIT_LORENTZ_INIT                CPL_TRUE

#define MF_PARAMETERS_RES_LORENTZ                     "RES_LORENTZ"
#define MF_PARAMETERS_RES_LORENTZ_DESC                "Initial value for FWHM of the Lorentz in pixels\n"\
                                                      " at the centre of the spectrum"
#define MF_PARAMETERS_RES_LORENTZ_MIN                 0.
#define MF_PARAMETERS_RES_LORENTZ_MAX                 100.
#define MF_PARAMETERS_RES_LORENTZ_INIT                1.


/* Variables for synthetic kernels */

#define MF_PARAMETERS_KERN_MODE                       "KERNMODE"
#define MF_PARAMETERS_KERN_MODE_DESC                  "Voigtian profile approximation instead of independent Gaussian and Lorentzian?"
#define MF_PARAMETERS_KERN_MODE_INIT                  CPL_FALSE

#define MF_PARAMETERS_KERN_FAC                        "KERNFAC"
#define MF_PARAMETERS_KERN_FAC_DESC                   "Size of Voigtian/Gaussian/Lorentzian kernel in FWHM"
#define MF_PARAMETERS_KERN_FAC_MIN                    3.
#define MF_PARAMETERS_KERN_FAC_MAX                    1000.
#define MF_PARAMETERS_KERN_FAC_INIT                   MF_PARAMETERS_KERN_FAC_MIN

#define MF_PARAMETERS_VAR_KERN                        "VARKERN"
#define MF_PARAMETERS_VAR_KERN_DESC                   "Does the kernel size increase linearly with wavelength?"
#define MF_PARAMETERS_VAR_KERN_INIT                   CPL_FALSE


/* If the input data file contains a suitable FITS header, the keyword names of the following parameters will be read,
 * but the corresponding values will not be used.
 * The reading of parameter values from this file can be forced by setting keywords to NONE. */

#define MF_PARAMETERS_OBSERVING_DATE_DESC             "Observing date in years or MJD in days (not string)"
#define MF_PARAMETERS_OBSERVING_DATE_KEYWORD          "OBSERVING_DATE_KEYWORD"
#define MF_PARAMETERS_OBSERVING_DATE_KEYWORD_INIT     "MJD-OBS"
#define MF_PARAMETERS_OBSERVING_DATE_VALUE            "OBSERVING_DATE_VALUE"
#define MF_PARAMETERS_OBSERVING_DATE_VALUE_DESC       "If OBSERVING_DATE_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_OBSERVING_DATE_VALUE_INIT       -1.

#define MF_PARAMETERS_UTC_DESC                        "UTC in s"
#define MF_PARAMETERS_UTC_KEYWORD                     "UTC_KEYWORD"
#define MF_PARAMETERS_UTC_KEYWORD_INIT                "UTC"
#define MF_PARAMETERS_UTC_VALUE                       "UTC_VALUE"
#define MF_PARAMETERS_UTC_VALUE_DESC                  "If UTC_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_UTC_VALUE_INIT                  -1.

/* Altitude is always above 0 (and even 20 for UTs) */
#define MF_PARAMETERS_TELESCOPE_ANGLE_DESC            "Telescope altitude angle in deg"
#define MF_PARAMETERS_TELESCOPE_ANGLE_KEYWORD         "TELESCOPE_ANGLE_KEYWORD"
#define MF_PARAMETERS_TELESCOPE_ANGLE_KEYWORD_INIT    "ESO TEL ALT"
#define MF_PARAMETERS_TELESCOPE_ANGLE_VALUE           "TELESCOPE_ANGLE_VALUE"
#define MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_DESC      "If TELESCOPE_ANGLE_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_MIN       0.                        /* Astronomical Horizon */
#define MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_MAX       90.                       /* Zenith               */
#define MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_INIT      90.

/* TODO : Add range inside checks */
#define MF_PARAMETERS_RELATIVE_HUMIDITY_DESC          "Relative humidity in %"
#define MF_PARAMETERS_RELATIVE_HUMIDITY_KEYWORD       "RELATIVE_HUMIDITY_KEYWORD"
#define MF_PARAMETERS_RELATIVE_HUMIDITY_KEYWORD_INIT  "ESO TEL AMBI RHUM"
#define MF_PARAMETERS_RELATIVE_HUMIDITY_VALUE         "RELATIVE_HUMIDITY_VALUE"
#define MF_PARAMETERS_RELATIVE_HUMIDITY_VALUE_DESC    "If RELATIVE_HUMIDITY_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_RELATIVE_HUMIDITY_MIN           0.
#define MF_PARAMETERS_RELATIVE_HUMIDITY_MAX           100.
#define MF_PARAMETERS_RELATIVE_HUMIDITY_VALUE_INIT    15.

/* TODO : Add range inside checks and test to see if telluriccorr works at 20 hPa, otherwise we should increase the lower limit */
#define MF_PARAMETERS_PRESSURE_DESC                   "Pressure in hPa"
#define MF_PARAMETERS_PRESSURE_KEYWORD                "PRESSURE_KEYWORD"
#define MF_PARAMETERS_PRESSURE_KEYWORD_INIT           "ESO TEL AMBI PRES START"
#define MF_PARAMETERS_PRESSURE_VALUE                  "PRESSURE_VALUE"
#define MF_PARAMETERS_PRESSURE_VALUE_DESC             "If PRESSURE_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_PRESSURE_VALUE_MIN              20
#define MF_PARAMETERS_PRESSURE_VALUE_MAX              1100
#define MF_PARAMETERS_PRESSURE_VALUE_INIT             750.

#define MF_PARAMETERS_TEMPERATURE_DESC                "Ambient temperature in deg C"
#define MF_PARAMETERS_TEMPERATURE_KEYWORD             "TEMPERATURE_KEYWORD"
#define MF_PARAMETERS_TEMPERATURE_KEYWORD_INIT        "ESO TEL AMBI TEMP"
#define MF_PARAMETERS_TEMPERATURE_VALUE               "TEMPERATURE_VALUE"
#define MF_PARAMETERS_TEMPERATURE_VALUE_DESC          "If TEMPERATURE_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_TEMPERATURE_VALUE_INIT          15.

#define MF_PARAMETERS_MIRROR_TEMPERATURE_DESC         "Mirror temperature in deg C"
#define MF_PARAMETERS_MIRROR_TEMPERATURE_KEYWORD      "MIRROR_TEMPERATURE_KEYWORD"
#define MF_PARAMETERS_MIRROR_TEMPERATURE_KEYWORD_INIT "ESO TEL TH M1 TEMP"
#define MF_PARAMETERS_MIRROR_TEMPERATURE_VALUE        "MIRROR_TEMPERATURE_VALUE"
#define MF_PARAMETERS_MIRROR_TEMPERATURE_VALUE_DESC   "If MIRROR_TEMPERATURE_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_MIRROR_TEMPERATURE_VALUE_INIT   15.

#define MF_PARAMETERS_ELEVATION_DESC                  "Elevation above sea level in m (default is Paranal: 2635. m)"
#define MF_PARAMETERS_ELEVATION_KEYWORD               "ELEVATION_KEYWORD"
#define MF_PARAMETERS_ELEVATION_KEYWORD_INIT          "ESO TEL GEOELEV"
#define MF_PARAMETERS_ELEVATION_VALUE                 "ELEVATION_VALUE"
#define MF_PARAMETERS_ELEVATION_VALUE_DESC            "If ELEVATION_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_ELEVATION_VALUE_INIT            2635.

#define MF_PARAMETERS_LONGITUDE_DESC                  "Longitude (default is Paranal: -70.4051 deg)"
#define MF_PARAMETERS_LONGITUDE_KEYWORD               "LONGITUDE_KEYWORD"
#define MF_PARAMETERS_LONGITUDE_KEYWORD_INIT          "ESO TEL GEOLON"
#define MF_PARAMETERS_LONGITUDE_VALUE                 "LONGITUDE_VALUE"
#define MF_PARAMETERS_LONGITUDE_VALUE_DESC            "If LONGITUDE_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_LONGITUDE_VALUE_MIN             -180.
#define MF_PARAMETERS_LONGITUDE_VALUE_MAX              180.
#define MF_PARAMETERS_LONGITUDE_VALUE_INIT            -70.4051

#define MF_PARAMETERS_LATITUDE_DESC                   "Latitude (default is Paranal: -24.6276 deg)"
#define MF_PARAMETERS_LATITUDE_KEYWORD                "LATITUDE_KEYWORD"
#define MF_PARAMETERS_LATITUDE_KEYWORD_INIT           "ESO TEL GEOLAT"
#define MF_PARAMETERS_LATITUDE_VALUE                  "LATITUDE_VALUE"
#define MF_PARAMETERS_LATITUDE_VALUE_DESC             "If LATITUDE_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_LATITUDE_VALUE_MIN              -90.
#define MF_PARAMETERS_LATITUDE_VALUE_MAX               90.
#define MF_PARAMETERS_LATITUDE_VALUE_INIT             -24.6276


/*** Telluric Correction parameters : Instrumental ***/

#define MF_PARAMETERS_SLIT_WIDTH_DESC                 "Slit width in arcsec (taken from FITS header if present)"
#define MF_PARAMETERS_SLIT_WIDTH_KEYWORD              "SLIT_WIDTH_KEYWORD"
#define MF_PARAMETERS_SLIT_WIDTH_KEYWORD_ANY          "ANY"
#define MF_PARAMETERS_SLIT_WIDTH_KEYWORD_DEFAULT      "ESO INS SLIT1 WID"
#define MF_PARAMETERS_SLIT_WIDTH_KEYWORD_INIT         MF_PARAMETERS_SLIT_WIDTH_KEYWORD_DEFAULT
#define MF_PARAMETERS_SLIT_WIDTH_VALUE                "SLIT_WIDTH_VALUE"
#define MF_PARAMETERS_SLIT_WIDTH_VALUE_DESC           "If SLIT_WIDTH_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_SLIT_WIDTH_VALUE_INIT           0.4

#define MF_PARAMETERS_PIXEL_SCALE_DESC                "Pixel scale in arcsec (taken from this file only)"
#define MF_PARAMETERS_PIXEL_SCALE_KEYWORD             "PIX_SCALE_KEYWORD"
#define MF_PARAMETERS_PIXEL_SCALE_KEYWORD_INIT        MF_PARAMETERS_KEYWORD_NONE
#define MF_PARAMETERS_PIXEL_SCALE_VALUE               "PIX_SCALE_VALUE"
#define MF_PARAMETERS_PIXEL_SCALE_VALUE_DESC          "If PIX_SCALE_KEYWORD=='NONE' take this value"
#define MF_PARAMETERS_PIXEL_SCALE_VALUE_INIT          0.086


/*** Telluric Correction parameters : Atmospheric ***/

#define MF_PARAMETERS_REFERENCE_ATMOSPHERIC           "REFERENCE_ATMOSPHERIC"
#define MF_PARAMETERS_REFERENCE_ATMOSPHERIC_DESC      "Reference atmospheric profile. Possible values:\n"\
                                                      " - equ.fits (default; equatorial atmosphere, valid for Paranal);\n"\
                                                      " - tro.fits (tropical atmosphere);\n"\
                                                      " - std.fits (standard atmosphere);\n"\
                                                      " - Other file located in :\n"\
                                                      "   ({TELLURICCORR_DATA_PATH}/profiles/mipas/)"
#define MF_PARAMETERS_REFERENCE_ATMOSPHERIC_EQU_ATM   "equ.fits"  /* Atmospheric profile : Tropical (day) profiles for MIPAS v3 - 03-JAN-2001 - J.J. Remedios - EOS, Space Research Centre, Leicester, U.K.        */
#define MF_PARAMETERS_REFERENCE_ATMOSPHERIC_STD_ATM   "std.fits"  /* Atmospheric profile : Standard atmosphere - Transformed to RFM .atm file format by program USARFM v23-AUG-96 - 50 levels */
#define MF_PARAMETERS_REFERENCE_ATMOSPHERIC_TRO_ATM   "tro.fits"  /* Atmospheric profile : Tropical atmosphere - Transformed to RFM .atm file format by program USARFM v23-AUG-96 - 50 levels */
#define MF_PARAMETERS_REFERENCE_ATMOSPHERIC_INIT      MF_PARAMETERS_REFERENCE_ATMOSPHERIC_EQU_ATM

#define MF_PARAMETERS_GDAS_PROFILE                    "GDAS_PROFILE"
#define MF_PARAMETERS_GDAS_PROFILE_DESC               "Specify which GDAS profile to use. Possible values: \n"\
                                                      " - 'auto',  automatic retrieval of the GDAS profiles\n"\
                                                      " (P[hPa] HGT[km] T[K] RELHUM[%]) close in time to the\n"\
                                                      " observation and in location to the observatory. If \n"\
                                                      " the files are not on disk and there is no internet \n"\
                                                      " connection, the average profile is taken from      \n"\
                                                      " share\/molecfit\/data\/profiles\/lib corresponding to\n"\
                                                      " the month of the observation (GDAS_t0_s1.fits for  \n"\
                                                      " Dec/Jan, GDAS_t0_s2.fits for Feb/Mar, etc)         \n"\
                                                      " See Sec. 8.1.4 of the molecfit manual for more info.\n"\
                                                      " - 'null', use the profile in the SOF with tag GDAS. \n"\
                                                      " If there is no profile in the SOF, the behaviour    \n"\
                                                      " is the same as GDAS_PROF=auto.                      \n"\
                                                      " - 'none', use the average profile taken from       \n"\
                                                      " share/molecfit/data/profiles/lib corresponding to  \n"\
                                                      " the month of observation (see 'auto' description). \n"\
                                                      " - 'directory/file.fits', use the specified path and\n"\
                                                      " filename as the GDAS profile. Either an absolute path\n"\
                                                      " (starting with '/') or a relative path may be used,\n"\
                                                      " however a relative path is preferred, since only  \n"\
                                                      " the first 40 char of the path and filename are    \n"\
                                                      " copied to the FITS header. The file format must be\n"\
                                                      " a FITS binary table with columns having the names: \n"\
                                                      " 'press height temp relhum' and units:              \n"\
                                                      " hPa, km, K and \%, respectively.                   \n" 
#define MF_PARAMETERS_GDAS_PROFILE_NONE               "none"
#define MF_PARAMETERS_GDAS_PROFILE_AUTO               "auto"
#define MF_PARAMETERS_GDAS_PROFILE_INIT               MF_PARAMETERS_GDAS_PROFILE_AUTO

#define MF_PARAMETERS_LAYERS                          "LAYERS"
#define MF_PARAMETERS_LAYERS_DESC                     "Grid of layer heights for merging ref_atm and GDAS profile.\n"\
                                                      " Fixed grid = CPL_TRUE and natural grid = CPL_FALSE"
#define MF_PARAMETERS_LAYERS_INIT                     CPL_TRUE

#define MF_PARAMETERS_EMIX                            "EMIX"
#define MF_PARAMETERS_EMIX_DESC                       "Upper mixing height in km for considering data of a local meteo station.\n"\
                                                      " If emix is below geoelev, rhum, pres, and temp are not used\n"\
                                                      " for modifying the corresponding profiles."
#define MF_PARAMETERS_EMIX_INIT                       5.

#define MF_PARAMETERS_PWV                             "PWV"
#define MF_PARAMETERS_PWV_DESC                        "PWV value in mm for the input water vapor profile.\n"\
                                                      " The merged profile composed of ref_atm, GDAS, and local meteo data\n"\
                                                      " will be scaled to this value if pwv > 0 (default: -1 -> no scaling)"
#define MF_PARAMETERS_PWV_INIT                        -1.

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

#define MF_PARAMETERS_MAXIMUM_NO_OF_RANGES     1001
#define MF_PARAMETERS_RANGE_FLAGS_NOT_DEFINED -1
#define MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS      1001
#define MF_PARAMETERS_CHIP_FLAGS_NOT_DEFINED  -1

/*----------------------------------------------------------------------------*/
/**
 *                 Typedefs: Structured types
 */
/*----------------------------------------------------------------------------*/

/* Structure for a parameter value in different data types and valid for keywords in cpl_paramater and cpl_property    */
typedef struct {
  char                       *key;                     /* Keyword name in the header                                   */
  double                     min;                      /* double minimum                                               */
  double                     max;                      /* double maximum                                               */
  double                     value;                    /* double                                                       */
} mf_parameters_keyword;

/*  */
typedef struct {
  cpl_boolean                fit;                      /* Flag: Polynomial fit                                         */
  cpl_size                   n;                        /* Polynomial degree                                            */
  double                     const_min;                /* Minimum value                                                */
  double                     const_max;                /* Minimum value                                                */
  double                     const_val;                /* Initial value of the constant term of the polynomial fit     */
} mf_parameters_fit;

/******************************************************************************/
/* mf_parameters: CONFIGURATION                                                */
/******************************************************************************/

/*  */
typedef struct {
  char                       *telluric_path;           /*  */
  char                       *telluriccorr_data_path;  /*  */
  char                       *tmp_path;                /*  */
//  char                       *output_path;             /*  */
//  char                       *output_name;             /*  */
} mf_parameters_directories;

/*  */
typedef struct {
  //unsigned int               omp_num_threads;          /*  */
  cpl_boolean                silent_external_bins;     /*  */
  cpl_boolean                transmission;             /*  */
  char                       *column_lambda;           /*  */
  char                       *column_flux;             /*  */
  char                       *column_dflux;            /*  */
  char                       *column_mask;             /*  */
  cpl_boolean                mask_binary;              /*  */
  int                        mask_threshold;           /*  */
  double                     default_error;            /*  */
  double                     wlg_to_micron;            /*  */
  char                       *wavelengths_frame;       /*  */
  mf_parameters_keyword      observing_erv_rv;         /*  */
  cpl_boolean                clean_flux;               /*  */
} mf_parameters_inputs;

/*  */

typedef struct {
  double                     ftol;                     /*  */
  double                     xtol;                     /*  */
  int                        flux_unit;                /*  */
  mf_parameters_fit          fit_telescope_background; /*  */
  mf_parameters_fit          fit_continuum;            /*  */
  /*double                     obs_bary_rv;*/          /*  */
  mf_parameters_fit          fit_wavelenght;           /*  */
//  char                       *wlc_ref;                 /*  */
  mf_parameters_fit          fit_res_box;              /*  */
  mf_parameters_fit          fit_gauss;                /*  */
  mf_parameters_fit          fit_lorentz;              /*  */
  cpl_boolean                kern_mode;                /*  */
  double                     kern_fac;                 /*  */
  cpl_boolean                var_kern;                 /*  */
  cpl_boolean                fit_chips[MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS];
  int                        fit_n_cflags;
  cpl_boolean                fit_ranges[MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];
  int                        fit_n_rflags;
  cpl_boolean                expert_mode;
  int                        cont_poly_order[MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];
  double                     cont_coeffs    [MF_PARAMETERS_MAXIMUM_NO_OF_RANGES][MF_FIT_N_POLYNOME_MAX];
  int                        wlc_poly_order [MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS];
  double                     wlc_coeffs     [MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS ][MF_FIT_N_POLYNOME_MAX];
} mf_parameters_fitting;

/*  */
typedef struct {
  mf_parameters_keyword      observing_date;           /*  */
  mf_parameters_keyword      utc;                      /*  */
  mf_parameters_keyword      telescope_angle;          /*  */
  mf_parameters_keyword      relative_humidity;        /*  */
  mf_parameters_keyword      pressure;                 /*  */
  mf_parameters_keyword      temperature;              /*  */
  mf_parameters_keyword      mirror_temperature;       /*  */
  mf_parameters_keyword      elevation;                /*  */
  mf_parameters_keyword      longitude;                /*  */
  mf_parameters_keyword      latitude;                 /*  */
} mf_parameters_ambient;

/*  */
typedef struct {
  mf_parameters_keyword      slit_width;               /*  */
  mf_parameters_keyword      pixel_scale;              /*  */
} mf_parameters_instrumental;

/*  */
typedef struct {
  char                       *ref_atm;                 /*  */
  char                       *gdas_prof;               /*  */
  cpl_boolean                layers;                   /*  */
  double                     emix;                     /*  */
  double                     pwv;                      /*  */
} mf_parameters_atmospheric;

/*  */
typedef struct {
  char                       *tmp_folder;              /*                                                              */
  cpl_boolean                single_spectrum;          /*                                                              */
  int                        n_range;                  /*                                                              */
  cpl_size                   nchip;                    /* Number of chips                                              */
  cpl_boolean                wlc_lin;                  /* Related with the polynomial degree of the wavelength fit     */
  double                     chi2;                     /*                                                              */
  mf_molecules               molecules;                /* Molecules parameters                                         */
} mf_parameters_internal;

/* Exposed structure : Initially default values for to be modify by the user */
typedef struct {
  mf_parameters_directories  directories;              /* Parameters related to the local directories                  */
  mf_parameters_inputs       inputs;                   /* Input        parameters                                      */
  mf_parameters_fitting      fitting;                  /* Fitting      parameters                                      */
  mf_parameters_ambient      ambient;                  /* Ambient      parameters                                      */
  mf_parameters_instrumental instrumental;             /* Instrumental parameters                                      */
  mf_parameters_atmospheric  atmospheric;              /* Atmospheric  parameters                                      */
  mf_parameters_internal     internal;                 /* Internal     parameters                                      */
} mf_parameters_config;

/* Exposed structure : Initially default values for to be modify by the user */
typedef struct {
  char                       *line_db;               /*  */
  int                        line_db_fmt;            /*  */
  cpl_boolean                use_ODA;                /*  */
} mf_io_lnfl_config;

/******************************************************************************/

/* Calculated variables */
typedef struct {
  double                     time_start;               /* Ini time in the telluriccorr call                            */
  double                     time_fit;                 /* Total FIT        run time                                    */
  double                     time_bin_lnfl;            /* Total LNFL   bin run time                                    */
  double                     time_bin_lblrtm;          /* Total LBLRTM bin run time                                    */
  double                     time_end;                 /* End time in the telluriccorr call                            */
} mf_parameters_timers;

/******************************************************************************/

/* Structure for storing the parameters */
typedef struct {
  mf_parameters_config       *config;                  /*                                                              */
  mf_parameters_timers       timers;                   /* Timer executions                                             */
  cpl_table                  *molectab;                /* CPL table for parameters related to molecules                */
  cpl_table                  *rangetab;                /*                                                              */
  cpl_table                  *chiptab;                 /* CPL table for chip-related parameters                        */
} mf_parameters;

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Create the user parameter structure (with default values) */
MF_EXPORT mf_parameters_config * mf_parameters_config_create(void);

/* Update the mf_parameters_config with the raw cpl_propertylist *primary_header */
MF_EXPORT cpl_error_code mf_parameters_config_update_with_header_keywords(
    mf_parameters_config     *config,
    const cpl_propertylist   *header);

/* Deallocate the user parameter structure */
MF_EXPORT void mf_parameters_config_delete(
    mf_parameters_config     *config);

/* Check the user parameter structure. For check the user modifications */
MF_EXPORT cpl_error_code mf_parameters_config_check(
    mf_parameters_config     *config,
    mf_io_lnfl_config        *lnfl_config);           /* User parameter structure    */

/* Initializes the structure for the parameters with default values */
MF_INTERNAL mf_parameters * mf_parameters_initialize(
    mf_parameters_config     *config,
    const cpl_table          *molecules);

/* Create default mf_parameters configuration */
MF_INTERNAL mf_parameters * mf_parameters_create(
    mf_parameters_config     *config,
    const cpl_table          *molecules,
    const char               *tmp_folder);

/* Delete mf_parameters structure */
MF_INTERNAL void mf_parameters_delete(
    mf_parameters            *params);


CPL_END_DECLS


#endif /* MF_PARAMETERS_H */
