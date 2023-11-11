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

#ifndef MF_LBLRTM_H
#define MF_LBLRTM_H

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include <cpl.h>

#include "mf_constants.h"
#include "mf_io.h"
#include "mf_parameters.h"

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

#define MF_LBLRTM_DELTA_FACTOR            0.14    /* Factor multiplied by delta                         */
#define MF_LBLRTM_DELTA_TOL               1.1     /* Exponent for increase the tolerance of the wl      */
#define MF_LBLRTM_DELTA_MIN             125.      /* Minimum delta                   in cm^{-1}         */
#define MF_LBLRTM_DELTA_MAX            1750.      /* Maximum delta                   in cm^{-1}         */
#define MF_LBLRTM_DELTA_ABS            2020.      /* Maximum delta between v1 and v2 in cm^{-1}         */

/* Continua and Rayleigh extinction
 *
 *  = 0  no continuum calculated
 *  = 1  all continua calculated,             including Rayleigh extinction where applicable
 *  = 2  H2O self not calculated,             all other continua/Rayleigh extinction calculated
 *  = 3  H2O foreign not calculated,          all other continua/Rayleigh extinction calculated
 *  = 4  H2O self and foreign not calculated, all other continua/Rayleigh extinction calculated
 *  = 5  Rayleigh extinction not calculated,  all other continua calculated
 *  = 6  Individual continuum scale factors input (Requires Record 1.2a)
 *
 *  [6 will not be implemented] */
#define MF_LBLRTM_ICNTNM               "LBLRTM_ICNTNM"
#define MF_LBLRTM_ICNTNM_DESC          "Continua and Rayleigh extinction [0,1,2,3,4,5]"
#define MF_LBLRTM_ICNTNM_MIN           0
#define MF_LBLRTM_ICNTNM_MAX           5
#define MF_LBLRTM_ICNTNM_INIT          5

/* Aerosols
 *
 *  = 0  no aerosols used
 *  = 1  internal LOWTRAN aerosol models
 *  = 5  spectral optical depths by layer from file 'in_lblrtm_cld'
 *  = 7  user defined aerosol models
 *  = 9  use precalculated aerosols (TAPE20 from a previous aerosol run)
 *
 *  [5,7,9 will not be implemented]  */
#define MF_LBLRTM_IAERSL               "LBLRTM_IAERSL"
#define MF_LBLRTM_IAERSL_DESC          "Aerosols [0,1]"
#define MF_LBLRTM_IAERSL_MIN           0
#define MF_LBLRTM_IAERSL_MAX           1
#define MF_LBLRTM_IAERSL_INIT          0

/* Optical depth values
 *
 *  Number of optical depth values printed for the beginning and
 *  ending of each panel as a result of convolution for current layer
 *  (for MPTS < O, output printing is suppressed)     */
#define MF_LBLRTM_MPTS                 "LBLRTM_MPTS"
#define MF_LBLRTM_MPTS_DESC            "Number of optical depth values"
#define MF_LBLRTM_MPTS_INIT            5

/* Number of values for each panel
 *
 *  Number of values printed for the beginning and ending of each panel
 *  as result of merge of current layer with previous layers
 *
 *  (optical depth for IEMIT = 0; radiance and transmission for IEMIT = 1)  */
#define MF_LBLRTM_NPTS                 "LBLRTM_NPTS"
#define MF_LBLRTM_NPTS_DESC            "Number of values for each panel"
#define MF_LBLRTM_NPTS_INIT            5

/* Beginning wavenumber value for the calculation [mu m] -> ~2400nm */
#define MF_LBLRTM_V1                   "LBLRTM_V1"
#define MF_LBLRTM_V1_DESC              "Beginning wavenumber value for the calculation"
#define MF_LBLRTM_V1_INIT              1.9

/* Ending wavenumber value for the calculation [mu m]
 *
 * (V2-V1 must be less than 2020 cm-1) -> ~1900nm */
#define MF_LBLRTM_V2                   "LBLRTM_V2"
#define MF_LBLRTM_V2_DESC              "Ending wavenumber value for the calculation"
#define MF_LBLRTM_V2_INIT              2.4

/* Number of sample points per mean halfwidth (between 1 and 4)
 *
 * (default = 4) */
#define MF_LBLRTM_SAMPLE               "LBLRTM_SAMPLE"
#define MF_LBLRTM_SAMPLE_DESC          "Number of sample points per mean halfwidth [between 1 to 4, default=4]"
#define MF_LBLRTM_SAMPLE_MIN           1
#define MF_LBLRTM_SAMPLE_MAX           4
#define MF_LBLRTM_SAMPLE_INIT          4

/* Average collision broadened halfwidth [cm - 1/atm]
 *
 * (default = 0.04) */
#define MF_LBLRTM_ALFAL0               "LBLRTM_ALFAL0"
#define MF_LBLRTM_ALFAL0_DESC          "Average collision broadened halfwidth [cm-1/atm]"
#define MF_LBLRTM_ALFAL0_INIT          0.

/* Average molecular mass [amu] for Doppler halfwidth
 *
 * (default = 36) */
#define MF_LBLRTM_AVMASS               "LBLRTM_AVMASS"
#define MF_LBLRTM_AVMASS_DESC          "Average molecular mass [amu] for Doppler halfwidth"
#define MF_LBLRTM_AVMASS_INIT          0.

/* Minimum molecular optical depth below which lines will be rejected
 *
 * (negative value defaults to DPTMIN = 0.0002) */
#define MF_LBLRTM_DPTMIN               "LBLRTM_DPTMIN"
#define MF_LBLRTM_DPTMIN_DESC          "Minimum molecular optical depth below which lines will be rejected"
#define MF_LBLRTM_DPTMIN_INIT          0.0002

/* Factor multiplying molecular continuum optical depth to determine optical depth below which lines will be rejected
 *
 * (negative value defaults to DPTFAC = 0.001) */
#define MF_LBLRTM_DPTFAC               "LBLRTM_DPTFAC"
#define MF_LBLRTM_DPTFAC_DESC          "Factor multiplying molecular continuum optical depth"
#define MF_LBLRTM_DPTFAC_INIT          0.001

/* Temperature of boundary [K] */
#define MF_LBLRTM_TBOUND               "LBLRTM_TBOUND"
#define MF_LBLRTM_TBOUND_DESC          "Temperature of boundary [K]"
#define MF_LBLRTM_TBOUND_INIT          0.

/* Frequency dependent boundary emissivity coefficients
 *
 * EMISSIVITY = SREMIS1 + SREMIS2*V + SREMIS3*(V**2) */
#define MF_LBLRTM_SREMIS1              "LBLRTM_SREMIS1"
#define MF_LBLRTM_SREMIS2              "LBLRTM_SREMIS2"
#define MF_LBLRTM_SREMIS3              "LBLRTM_SREMIS3"
#define MF_LBLRTM_SREMIS1_DESC         "Emissivity coefficient 1"
#define MF_LBLRTM_SREMIS2_DESC         "Emissivity coefficient 2"
#define MF_LBLRTM_SREMIS3_DESC         "Emissivity coefficient 3"
#define MF_LBLRTM_SREMIS1_INIT         0.
#define MF_LBLRTM_SREMIS2_INIT         0.
#define MF_LBLRTM_SREMIS3_INIT         0.

/* Frequency dependent boundary reflectivity coefficients
 *
 * REFLECTIVITY = SRREFL1 + SRREFL2*V + SRREFL3*(V**2) */
#define MF_LBLRTM_SRREFL1              "LBLRTM_SRREFL1"
#define MF_LBLRTM_SRREFL2              "LBLRTM_SRREFL2"
#define MF_LBLRTM_SRREFL3              "LBLRTM_SRREFL3"
#define MF_LBLRTM_SRREFL1_DESC         "Reflectivity coefficient 1"
#define MF_LBLRTM_SRREFL2_DESC         "Reflectivity coefficient 2"
#define MF_LBLRTM_SRREFL3_DESC         "Reflectivity coefficient 3"
#define MF_LBLRTM_SRREFL1_INIT         0.
#define MF_LBLRTM_SRREFL2_INIT         0.
#define MF_LBLRTM_SRREFL3_INIT         0.

/* Selects atmospheric profile
 *
 *  = 0  user supplied atmospheric profile
 *  = 1  tropical model
 *  = 2  midlatitude summer model
 *  = 3  midlatitude winter model
 *  = 4  subarctic summer model
 *  = 5  subarctic winter model
 *  = 6  U.S. standard 1976     */
#define MF_LBLRTM_MODEL                "LBLRTM_MODEL"
#define MF_LBLRTM_MODEL_DESC           "Atmospheric profile [0,1,2,3,4,5,6]"
#define MF_LBLRTM_MODEL_MIN            0
#define MF_LBLRTM_MODEL_MAX            6
#define MF_LBLRTM_MODEL_INIT           0

/* Selects type of path
 *
 *  = 1  horizontal path (constant pressure, temperature), use RECORD 3.2H
 *  = 2  slant path from H1 to H2, use RECORD 3.2
 *  = 3  slant path from H1 to space (see HSPACE), use RECORD 3.2     */
#define MF_LBLRTM_ITYPE                "LBLRTM_ITYPE"
#define MF_LBLRTM_ITYPE_DESC           "Type of path [1,2,3]"
#define MF_LBLRTM_ITYPE_MIN            1
#define MF_LBLRTM_ITYPE_MAX            3
#define MF_LBLRTM_ITYPE_INIT           3

/* Zeroing of small amounts of absorbers
 *
 *  = 0  zeroes absorber amounts which are less than 0.1 percent of total (default)
 *  = 1  suppresses zeroing of small amounts     */
#define MF_LBLRTM_NOZERO               "LBLRTM_NOZERO"
#define MF_LBLRTM_NOZERO_DESC          "Zeroing of small amounts of absorbers [0,1]"
#define MF_LBLRTM_NOZERO_MIN           0
#define MF_LBLRTM_NOZERO_MAX           1
#define MF_LBLRTM_NOZERO_INIT          0

/* Output
 *
 *  = 0  full printout
 *  = 1  selects short printout     */
#define MF_LBLRTM_NOPRNT               "LBLRTM_NOPRNT"
#define MF_LBLRTM_NOPRNT_DESC          "Do not print output? [0,1]"
#define MF_LBLRTM_NOPRNT_MIN           0
#define MF_LBLRTM_NOPRNT_MAX           1
#define MF_LBLRTM_NOPRNT_INIT          0

/* Write out layer data
 *
 *  = 0  layer data not written (default)
 *  = 1  layer data written to unit ITAPE7)PU (TAPE7)     */
#define MF_LBLRTM_IPUNCH               "LBLRTM_IPUNCH"
#define MF_LBLRTM_IPUNCH_DESC          "Write out layer data to TAPE7 [0,1]"
#define MF_LBLRTM_IPUNCH_MIN           0
#define MF_LBLRTM_IPUNCH_MAX           1
#define MF_LBLRTM_IPUNCH_INIT          0

/* Radius of earth [km]
 *
 *  defaults for RE=0:
 *  a)  MODEL 0,2,3,6    RE = 6371.23 km
 *  b)        1          RE = 6378.39 km
 *  c)        4,5        RE = 6356.91 km     */
#define MF_LBLRTM_RE                   "LBLRTM_RE"
#define MF_LBLRTM_RE_DESC              "Radius of earth [km]"
#define MF_LBLRTM_RE_INIT              0.

/* Altitude definition for space (default = 100 km)
 *
 *  internal models defined to 120 km  */
#define MF_LBLRTM_HSPACE               "LBLRTM_HSPACE"
#define MF_LBLRTM_HSPACE_DESC          "Altitude definition for space [km]"
#define MF_LBLRTM_HSPACE_INIT          120.

/* Latitude of location of calculation [degrees]
 *
 *  defaults for REF_LAT = 0:
 *  a) MODEL 0,2,3,6    REF_LAT = 45.0 degrees
 *  b) MODEL 1          REF_LAT = 15.0
 *  c) MODEL 4,5        REF_LAT = 60.0     */
/*#define MF_LBLRTM_REF_LAT              "LBLRTM_REF_LAT"
#define MF_LBLRTM_REF_LAT_DESC         MF_PARAMETERS_LATITUDE_DESC // "Latitude of location of calculation [degrees] [-90.-90]"
#define MF_LBLRTM_REF_LAT_MIN          MF_PARAMETERS_LATITUDE_VALUE_MIN // -90.
#define MF_LBLRTM_REF_LAT_MAX          MF_PARAMETERS_LATITUDE_VALUE_MAX // 90.
#define MF_LBLRTM_REF_LAT_INIT         MF_PARAMETERS_LATITUDE_VALUE_INIT // -24.63 */

/* Observer altitude [km] */
/*#define MF_LBLRTM_H1                   "LBLRTM_H1"
#define MF_LBLRTM_H1_DESC              "Observer altitude [km]"
#define MF_LBLRTM_H1_INIT              MF_PARAMETERS_ELEVATION_VALUE_INIT/1000.0// 2.64  */

/* Upper height limit
 *
 *  for ITYPE = 2, H2 is the end point altitude [km]
 *      ITYPE = 3, H2 is the tangent height [km] for H2 .GT. 0.
 *                 if H2 = 0. ANGLE determines tangent height */
#define MF_LBLRTM_H2                   "LBLRTM_H2"
#define MF_LBLRTM_H2_DESC              "Upper height limit [km]"
#define MF_LBLRTM_H2_INIT              0.

/* Length of a straight path from H1 to H2 [km] */
#define MF_LBLRTM_RANGE                "LBLRTM_RANGE"
#define MF_LBLRTM_RANGE_DESC           "Length of a straight path from H1 to H2 [km]"
#define MF_LBLRTM_RANGE_INIT           0.

/* Earth centered angle from H1 to H2 [degrees] */
#define MF_LBLRTM_BETA                 "LBLRTM_BETA"
#define MF_LBLRTM_BETA_DESC            "Earth centered angle from H1 to H2 [degrees]"
#define MF_LBLRTM_BETA_INIT            0.

/* Path length
 *
 *  = 0  short path (default)
 *  = 1  long path through a tangent height
 *
 *  LEN is only used for H1 > H2 (ANGLE > 90`)
 *
 *  for ITYPE = 2, only 3 of the first 5 parameters are required to
 *                 specify the path, e.g., H1, H2, ANGLE or H1, H2 and
 *                 RANGE
 *
 *  for ITYPE = 3, H1 = observer altitude must be specified. Either
 *                 H2 = tangent height or ANGLE must be specified.
 *                 Other parameters are ignored.     */
#define MF_LBLRTM_LEN                  "LBLRTM_LEN"
#define MF_LBLRTM_LEN_DESC             "Path length [0,1]"
#define MF_LBLRTM_LEN_MIN              0
#define MF_LBLRTM_LEN_MAX              1
#define MF_LBLRTM_LEN_INIT             0

/* Height of observer
 *
 *  Height of observer, used only for informational purposes in
 *  satellite-type simulations when computing output geometry
 *  above 120 km. */
#define MF_LBLRTM_HOBS                 "LBLRTM_HOBS"
#define MF_LBLRTM_HOBS_DESC            "Height of observer"
#define MF_LBLRTM_HOBS_INIT            0.

/* Maximum Voigt width ratio across a layer (if zero, default = 1.5) */
#define MF_LBLRTM_AVTRAT               "LBLRTM_AVTRAT"
#define MF_LBLRTM_AVTRAT_DESC          "Maximum Voigt width ratio across a layer"
#define MF_LBLRTM_AVTRAT_INIT          2.

/* Maximum layer temperature difference at ALTD1 (if zero, default = 5 K) */
#define MF_LBLRTM_TDIFF1               "LBLRTM_TDIFF1"
#define MF_LBLRTM_TDIFF1_DESC          "Maximum layer temperature difference at ALTD1 [K]"
#define MF_LBLRTM_TDIFF1_INIT          5.

/* Maximum layer temperature difference at ALTD2 (if zero, default = 8 K) */
#define MF_LBLRTM_TDIFF2               "LBLRTM_TDIFF2"
#define MF_LBLRTM_TDIFF2_DESC          "Maximum layer temperature difference at ALTD2 [K]"
#define MF_LBLRTM_TDIFF2_INIT          8.

/* Altitude of TDIFF1 (if zero, default = 0 Km) */
#define MF_LBLRTM_ALTD1                "LBLRTM_ALTD1"
#define MF_LBLRTM_ALTD1_DESC           "Altitude of TDIFF1 [km]"
#define MF_LBLRTM_ALTD1_INIT           0.

/* Altitude of TDIFF2 (if zero, default = 100 Km) */
#define MF_LBLRTM_ALTD2                "LBLRTM_ALTD2"
#define MF_LBLRTM_ALTD2_DESC           "Altitude of TDIFF2 [km]"
#define MF_LBLRTM_ALTD2_INIT           0.

/* Number of wavenumbers [cm-1] per major division. */
#define MF_LBLRTM_DELV                 "LBLRTM_DELV"
#define MF_LBLRTM_DELV_DESC            "Number of wavenumbers [cm-1] per major division"
#define MF_LBLRTM_DELV_INIT            1.

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

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Create the user parameter structure (with default values) */
MF_EXPORT mf_io_lblrtm_config * mf_lblrtm_config_create(void);

/* Deallocate the user parameter structure */
MF_EXPORT void mf_lblrtm_config_delete(
    mf_io_lblrtm_config      *config);

/* Check the user parameter structure. For check the user modifications */
MF_EXPORT cpl_error_code mf_lblrtm_config_check(
    mf_io_lblrtm_config      *config);

/* Execute the LBLRTM binary in all ranges */
MF_INTERNAL cpl_error_code mf_lblrtm(
    const cpl_size           run_execution,
    mf_parameters            *params,
    cpl_boolean              last_call,
    cpl_table                **spec_out,
    cpl_error_code           *range_status,
    cpl_array                *reffitpar,
    int                      mpfit_calls,
    int                      *lblrtm_calls,
    mf_io_lblrtm_config      *lblrtm_config ,
    mf_io_lnfl_config        *lnfl_config,
    const cpl_table          *atm_profile,
    cpl_table                **atm_profile_out,
    cpl_array                *fitpar);


cpl_error_code mf_io_lblrtm_oda(mf_io_lnfl_config   *lnfl_config,
                                mf_io_lblrtm_config *lblrtm_config,
                                mf_parameters       *params,
                                cpl_table           *atm);


CPL_END_DECLS


#endif /* MF_LBLRTM_H */
