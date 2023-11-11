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

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include "mf_io.h"

#include "mf_parameters.h"

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

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Converts MJD into fractional year (return value, double precision) */
static double mf_parameters_convert_obsdate_from_mjd_to_frac_year(
    double                   mjd);

/*  */
static double mf_parameters_slit_width_xshooter(
    const cpl_propertylist   *header);

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_parameters       Telluric Correction parameters and defines : User and internal structures and management
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create the user parameter structure (with default values).
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
mf_parameters_config * mf_parameters_config_create(void)
{
    /*** Initialize ***/
    mf_parameters_config *config = cpl_malloc(sizeof(mf_parameters_config));

    const char *telluriccorr_env = mf_io_getenv(MF_PARAMETERS_TELLURICCORR_PATH_ENV);
    if (telluriccorr_env) {
        cpl_msg_info(cpl_func, "(mf_parameters) Exist environment variable %s --> telluriccorr_path = %s", MF_PARAMETERS_TELLURICCORR_PATH_ENV, telluriccorr_env);
        config->directories.telluric_path              = cpl_sprintf("%s", telluriccorr_env);
    } else {
        cpl_msg_info(cpl_func, "(mf_parameters) Not exist environment variable %s      --> telluriccorr_path = %s", MF_PARAMETERS_TELLURICCORR_PATH_ENV, MF_PARAMETERS_TELLURICCORR_PATH_INIT);
        config->directories.telluric_path              = cpl_sprintf("%s", MF_PARAMETERS_TELLURICCORR_PATH_INIT);
    }

    const char *data_env = mf_io_getenv(MF_PARAMETERS_TELLURICCORR_DATA_PATH_ENV);
    if (data_env) {
        cpl_msg_info(cpl_func, "(mf_parameters) Exist environment variable %s --> telluriccorr_data_path = %s", MF_PARAMETERS_TELLURICCORR_DATA_PATH_ENV, data_env);
        config->directories.telluriccorr_data_path         = cpl_sprintf("%s", data_env);
    } else {
        cpl_msg_info(cpl_func, "(mf_parameters) Not exist environment variable %s --> telluriccorr_data_path = %s", MF_PARAMETERS_TELLURICCORR_DATA_PATH_ENV, MF_PARAMETERS_TELLURICCORR_DATA_PATH_INIT);
        config->directories.telluriccorr_data_path         = cpl_sprintf("%s", MF_PARAMETERS_TELLURICCORR_DATA_PATH_INIT);
    }

    config->directories.tmp_path                       = mf_io_pwd();
/*    config->directories.output_path                    = cpl_sprintf("%s/%s", config->directories.tmp_path, MF_PARAMETERS_OUTPUT_PATH_INIT);
    config->directories.output_name                    = cpl_sprintf("%s", MF_PARAMETERS_OUTPUT_NAME_INIT); */

    /* OMP executions */
    //const char *num_threads = mf_io_getenv(MF_PARAMETERS_OMP_NUM_THREADS_ENV);
    /*if (num_threads) config->inputs.omp_num_threads    = atoi(num_threads); */
    /*else             config->inputs.omp_num_threads    = 1;*/

    config->inputs.silent_external_bins                = MF_PARAMETERS_SILENT_EXTERNAL_BINS_INIT;

    config->inputs.transmission                        = MF_PARAMETERS_TRANSMISSION_INIT;
    config->inputs.column_lambda                       = cpl_sprintf("%s", MF_PARAMETERS_COLUMN_LAMBDA_INIT);
    config->inputs.column_flux                         = cpl_sprintf("%s", MF_PARAMETERS_COLUMN_FLUX_INIT);
    config->inputs.column_dflux                        = cpl_sprintf("%s", MF_PARAMETERS_COLUMN_DFLUX_INIT);
    config->inputs.column_mask                         = cpl_sprintf("%s", MF_PARAMETERS_COLUMN_MASK_INIT);
    config->inputs.default_error                       = MF_PARAMETERS_DEFAULT_ERROR_INIT;
    config->inputs.wlg_to_micron                       = MF_PARAMETERS_WLG_TO_MICRON_INIT;
    config->inputs.wavelengths_frame                   = cpl_sprintf("%s", MF_PARAMETERS_WAVELENGTH_FRAME_INIT);
    config->inputs.observing_erv_rv.key                = cpl_sprintf("%s", MF_PARAMETERS_OBSERVATORY_ERF_RV_KEYWORD_INIT);
    config->inputs.observing_erv_rv.value              = MF_PARAMETERS_OBSERVATORY_ERF_RV_VALUE_INIT;
    config->inputs.clean_flux                          = MF_PARAMETERS_CLEAN_MODEL_FLUX_INIT;

    config->fitting.ftol                               = MF_PARAMETERS_FTOL_INIT;
    config->fitting.xtol                               = MF_PARAMETERS_XTOL_INIT;
    config->fitting.flux_unit                          = MF_PARAMETERS_FLUX_UNIT_INIT;

    config->fitting.fit_telescope_background.fit       = MF_PARAMETERS_FIT_TELESCOPE_BACK_INIT;
    config->fitting.fit_telescope_background.const_min = MF_PARAMETERS_TELESCOPE_BACK_CONST_MIN;
    config->fitting.fit_telescope_background.const_max = MF_PARAMETERS_TELESCOPE_BACK_CONST_MAX;
    config->fitting.fit_telescope_background.const_val = MF_PARAMETERS_TELESCOPE_BACK_CONST_INIT;

    config->fitting.fit_continuum.fit                  = MF_PARAMETERS_CONTINUUM_N_INIT_SINGLE_VALUE;
    config->fitting.fit_continuum.n                    = MF_PARAMETERS_CONTINUUM_N_INIT_SINGLE_VALUE;
    config->fitting.fit_continuum.const_val            = MF_PARAMETERS_CONTINUUM_CONST_INIT;
    /*config->fitting.obs_bary_rv                        = MF_PARAMETERS_OBSERVATORY_BARY_RV_INIT;*/

    config->fitting.fit_wavelenght.fit                 = MF_PARAMETERS_FIT_WLC_INIT_SINGLE_VALUE;
    config->fitting.fit_wavelenght.n                   = MF_PARAMETERS_WLC_N_INIT;
    config->fitting.fit_wavelenght.const_val           = MF_PARAMETERS_WLC_CONST_INIT;
    //config->fitting.wlc_ref                            = cpl_sprintf("%s", MF_PARAMETERS_WLC_REF_INIT);

    config->fitting.fit_res_box.fit                    = MF_PARAMETERS_FIT_RES_BOX_INIT;
    config->fitting.fit_res_box.const_min              = MF_PARAMETERS_RES_BOX_MIN;
    config->fitting.fit_res_box.const_max              = MF_PARAMETERS_RES_BOX_MAX;
    config->fitting.fit_res_box.const_val              = MF_PARAMETERS_RES_BOX_INIT;

    config->fitting.fit_gauss.fit                      = MF_PARAMETERS_FIT_GAUSS_INIT;
    config->fitting.fit_gauss.const_min                = MF_PARAMETERS_RES_GAUSS_MIN;
    config->fitting.fit_gauss.const_max                = MF_PARAMETERS_RES_GAUSS_MAX;
    config->fitting.fit_gauss.const_val                = MF_PARAMETERS_RES_GAUSS_INIT;

    config->fitting.fit_lorentz.fit                    = MF_PARAMETERS_FIT_LORENTZ_INIT;
    config->fitting.fit_lorentz.const_min              = MF_PARAMETERS_RES_LORENTZ_MIN;
    config->fitting.fit_lorentz.const_max              = MF_PARAMETERS_RES_LORENTZ_MAX;
    config->fitting.fit_lorentz.const_val              = MF_PARAMETERS_RES_LORENTZ_INIT;

    config->fitting.kern_mode                          = MF_PARAMETERS_KERN_MODE_INIT;
    config->fitting.kern_fac                           = MF_PARAMETERS_KERN_FAC_INIT;
    config->fitting.var_kern                           = MF_PARAMETERS_VAR_KERN_INIT;
    
    config->fitting.fit_n_cflags                       = MF_PARAMETERS_CHIP_FLAGS_NOT_DEFINED;/*Implies that chip flags
                                                                                                 have not been specified*/
    for (int i=0; i<MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS; i++) {
        config->fitting.fit_chips[i]=CPL_FALSE; /*Not neccessary but makes me feel better (mnb)*/
    }
 
    config->fitting.fit_n_rflags                       = MF_PARAMETERS_RANGE_FLAGS_NOT_DEFINED;/*Implies that range flags
                                                                                                 have not been specified*/
    for (int i=0; i<MF_PARAMETERS_MAXIMUM_NO_OF_RANGES; i++) {
        config->fitting.fit_ranges[i]=CPL_FALSE; /*Not neccessary but makes me feel better (mnb)*/
    }

    config->fitting.expert_mode                        = MF_PARAMETERS_EXPERT_MODE_INIT;

 
    config->ambient.observing_date.key                 = cpl_sprintf("%s", MF_PARAMETERS_OBSERVING_DATE_KEYWORD_INIT);
    config->ambient.observing_date.value               = MF_PARAMETERS_OBSERVING_DATE_VALUE_INIT;

    config->ambient.utc.key                            = cpl_sprintf("%s", MF_PARAMETERS_UTC_KEYWORD_INIT);
    config->ambient.utc.value                          = MF_PARAMETERS_UTC_VALUE_INIT;

    config->ambient.telescope_angle.key                = cpl_sprintf("%s", MF_PARAMETERS_TELESCOPE_ANGLE_KEYWORD_INIT);
    config->ambient.telescope_angle.value              = MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_INIT;
    config->ambient.telescope_angle.min                = MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_MIN;
    config->ambient.telescope_angle.max                = MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_MAX;

    config->ambient.relative_humidity.key              = cpl_sprintf("%s", MF_PARAMETERS_RELATIVE_HUMIDITY_KEYWORD_INIT);
    config->ambient.relative_humidity.value            = MF_PARAMETERS_RELATIVE_HUMIDITY_VALUE_INIT;

    config->ambient.pressure.key                       = cpl_sprintf("%s", MF_PARAMETERS_PRESSURE_KEYWORD_INIT);
    config->ambient.pressure.value                     = MF_PARAMETERS_PRESSURE_VALUE_INIT;

    config->ambient.temperature.key                    = cpl_sprintf("%s", MF_PARAMETERS_TEMPERATURE_KEYWORD_INIT);
    config->ambient.temperature.value                  = MF_PARAMETERS_TEMPERATURE_VALUE_INIT;

    config->ambient.mirror_temperature.key             = cpl_sprintf("%s", MF_PARAMETERS_MIRROR_TEMPERATURE_KEYWORD_INIT);
    config->ambient.mirror_temperature.value           = MF_PARAMETERS_MIRROR_TEMPERATURE_VALUE_INIT;

    config->ambient.elevation.key                      = cpl_sprintf("%s", MF_PARAMETERS_ELEVATION_KEYWORD_INIT);
    config->ambient.elevation.value                    = MF_PARAMETERS_ELEVATION_VALUE_INIT;

    config->ambient.longitude.key                      = cpl_sprintf("%s", MF_PARAMETERS_LONGITUDE_KEYWORD_INIT);
    config->ambient.longitude.value                    = MF_PARAMETERS_LONGITUDE_VALUE_INIT;
    config->ambient.longitude.min                      = MF_PARAMETERS_LONGITUDE_VALUE_MIN;
    config->ambient.longitude.max                      = MF_PARAMETERS_LONGITUDE_VALUE_MAX;

    config->ambient.latitude.key                       = cpl_sprintf("%s", MF_PARAMETERS_LATITUDE_KEYWORD_INIT);
    config->ambient.latitude.value                     = MF_PARAMETERS_LATITUDE_VALUE_INIT;
    config->ambient.latitude.min                       = MF_PARAMETERS_LATITUDE_VALUE_MIN;
    config->ambient.latitude.max                       = MF_PARAMETERS_LATITUDE_VALUE_MAX;

    config->instrumental.slit_width.key                = cpl_sprintf("%s", MF_PARAMETERS_SLIT_WIDTH_KEYWORD_INIT);
    config->instrumental.slit_width.value              = MF_PARAMETERS_SLIT_WIDTH_VALUE_INIT;

    config->instrumental.pixel_scale.key               = cpl_sprintf("%s", MF_PARAMETERS_PIXEL_SCALE_KEYWORD_INIT);
    config->instrumental.pixel_scale.value             = MF_PARAMETERS_PIXEL_SCALE_VALUE_INIT;

    config->atmospheric.ref_atm                        = cpl_sprintf("%s", MF_PARAMETERS_REFERENCE_ATMOSPHERIC_INIT);
    config->atmospheric.gdas_prof                      = cpl_sprintf("%s", MF_PARAMETERS_GDAS_PROFILE_INIT);
    config->atmospheric.layers                         = MF_PARAMETERS_LAYERS_INIT;
    config->atmospheric.emix                           = MF_PARAMETERS_EMIX_INIT;
    config->atmospheric.pwv                            = MF_PARAMETERS_PWV_INIT;

    config->internal.tmp_folder                        = NULL;
    config->internal.single_spectrum                   = CPL_FALSE;
    config->internal.n_range                           = 1;
    config->internal.nchip                             = 1;
    config->internal.chi2                              = 1.e12;

    config->internal.molecules.n_molec                 = MF_MOLECULES_NUMBER_INIT;
    config->internal.molecules.lbl_molecs              = cpl_sprintf("%s", MF_MOLECULES_FIT_ARRAY_FLAGS_INIT);

    return config;
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief Update the mf_parameters_config with the raw cpl_propertylist *primary_header
 *
 * @param config             Initialized mf_parameters_config structure
 * @param header             Input cpl_property *header
 *
 * @return cpl_error_code    CPL_ERROR_NONE if everything is OK or NULL in other case.
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_parameters_config_update_with_header_keywords(
    mf_parameters_config     *config,
    const cpl_propertylist   *header)
{
  /*** Check inputs ***/
  cpl_ensure(config && header, CPL_ERROR_NULL_INPUT, CPL_ERROR_NULL_INPUT);

  /*** Local variables ***/
  cpl_error_code     err   = CPL_ERROR_NONE;
  const cpl_property *prop = NULL;


  /*** CHECK KEYWORDS ***/
  cpl_msg_info(cpl_func, "(mf_parameters) Updating mf_parameters_config with the user provided data header ... !");


  /* MF_PARAMETERS_OBSERVATORY_ERF_RV_KEYWORD */
  if (!err && strcmp(config->inputs.observing_erv_rv.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->inputs.observing_erv_rv.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->inputs.observing_erv_rv.key);

          double value = cpl_property_get_double(prop);

          config->inputs.observing_erv_rv.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_OBSERVATORY_ERF_RV_KEYWORD, config->inputs.observing_erv_rv.key, config->inputs.observing_erv_rv.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->inputs.observing_erv_rv.key);
      }
  }


  /* MF_PARAMETERS_OBSERVING_DATE */
  if (!err) {
      double value = config->ambient.observing_date.value;

      cpl_boolean exist_keyword = CPL_FALSE;
      if (strcmp(config->ambient.observing_date.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {
          exist_keyword = CPL_TRUE;
          err = (cpl_propertylist_has(header, config->ambient.observing_date.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
          if (!err) {
              prop  = cpl_propertylist_get_property_const(header, config->ambient.observing_date.key);
              value = cpl_property_get_double(prop);
          } else {
              cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.observing_date.key);
          }
      }

      if (!err) {

          if (value < 0.)  {
              err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "KEYWORD[%s] : invalid data", config->ambient.observing_date.key);
          } else {

              if(value < 10000.) {
                  /* Years */
                  config->ambient.observing_date.value = value;
              } else {
                  /* MJD */
                  cpl_msg_info(cpl_func, "(mf_parameters) Input MJD-OBS = %lf", value);
                  config->ambient.observing_date.value = mf_parameters_convert_obsdate_from_mjd_to_frac_year(value);
                  cpl_msg_info(cpl_func, "(mf_parameters) Convert MJD into date in years [%lf => %lf]", value, config->ambient.observing_date.value);
              }

              if (exist_keyword) {
                  cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %lf", MF_PARAMETERS_OBSERVING_DATE_KEYWORD, config->ambient.observing_date.key, config->ambient.observing_date.value);
              } else {
                  cpl_msg_info(cpl_func, "(mf_parameters) Set by parameters OBSERVING_DATE = %lf", config->ambient.observing_date.value);
              }
          }
      }
  }

  /* MF_PARAMETERS_UTC */
  if (!err && strcmp(config->ambient.utc.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.utc.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {
          prop  = cpl_propertylist_get_property_const(header, config->ambient.utc.key);

          double value = cpl_property_get_double(prop);
          if (value < 0) {
              err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "KEYWORD[%s] : invalid data", config->ambient.utc.key);
          } else {
              config->ambient.utc.value = value;
              cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_UTC_KEYWORD, config->ambient.utc.key, config->ambient.utc.value);
          }
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.utc.key);
      }
  }

  /* MF_PARAMETERS_TELESCOPE_ANGLE */
  if (!err && strcmp(config->ambient.telescope_angle.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.telescope_angle.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->ambient.telescope_angle.key);

          double value = cpl_property_get_double(prop);

          config->ambient.telescope_angle.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_TELESCOPE_ANGLE_KEYWORD, config->ambient.telescope_angle.key, config->ambient.telescope_angle.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.telescope_angle.key);
      }
  }

  /* MF_PARAMETERS_RELATIVE_HUMIDITY */
  if (!err && strcmp(config->ambient.relative_humidity.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.relative_humidity.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->ambient.relative_humidity.key);

          double value = cpl_property_get_double(prop);

          config->ambient.relative_humidity.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_RELATIVE_HUMIDITY_KEYWORD, config->ambient.relative_humidity.key, config->ambient.relative_humidity.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.relative_humidity.key);
      }
  }

  /* MF_PARAMETERS_PRESSURE */
  if (!err && strcmp(config->ambient.pressure.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.pressure.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->ambient.pressure.key);

          double value = cpl_property_get_double(prop);

          config->ambient.pressure.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_PRESSURE_KEYWORD, config->ambient.pressure.key, config->ambient.pressure.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.pressure.key);
      }
  }

  /* MF_PARAMETERS_TEMPERATURE */
  if (!err && strcmp(config->ambient.temperature.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.temperature.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->ambient.temperature.key);

          double value = cpl_property_get_double(prop);

          config->ambient.temperature.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_TEMPERATURE_KEYWORD, config->ambient.temperature.key, config->ambient.temperature.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.temperature.key);
      }
  }

  /* MF_PARAMETERS_MIRROR_TEMPERATURE */
  if (!err && strcmp(config->ambient.mirror_temperature.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.mirror_temperature.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->ambient.mirror_temperature.key);

          double value = cpl_property_get_double(prop);

          config->ambient.mirror_temperature.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_MIRROR_TEMPERATURE_KEYWORD, config->ambient.mirror_temperature.key, config->ambient.mirror_temperature.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.mirror_temperature.key);
      }
  }

  /* MF_PARAMETERS_ELEVATION */
  if (!err && strcmp(config->ambient.elevation.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.elevation.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->ambient.elevation.key);

          double value = cpl_property_get_double(prop);

          config->ambient.elevation.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_ELEVATION_KEYWORD, config->ambient.elevation.key, config->ambient.elevation.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.elevation.key);
      }
  }

  /* MF_PARAMETERS_LONGITUDE */
  if (!err && strcmp(config->ambient.longitude.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.longitude.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->ambient.longitude.key);

          double value = cpl_property_get_double(prop);

          config->ambient.longitude.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_LONGITUDE_KEYWORD, config->ambient.longitude.key, config->ambient.longitude.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.longitude.key);
      }
  }

  /* MF_PARAMETERS_LATITUDE */
  if (!err && strcmp(config->ambient.latitude.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.latitude.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->ambient.latitude.key);

          double value = cpl_property_get_double(prop);

          config->ambient.latitude.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_LATITUDE_KEYWORD, config->ambient.latitude.key, config->ambient.latitude.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->ambient.latitude.key);
      }
  }


  /* MF_PARAMETERS_SLIT_WIDTH */
  if (!err && strcmp(config->instrumental.slit_width.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      if (cpl_propertylist_has(header, config->instrumental.slit_width.key) == CPL_TRUE) {

          prop  = cpl_propertylist_get_property_const(header, config->instrumental.slit_width.key);

          double value = cpl_property_get_double(prop);

          config->instrumental.slit_width.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_SLIT_WIDTH_KEYWORD, config->instrumental.slit_width.key, config->instrumental.slit_width.value);

      } else {

          /* Check special case, if XSHOOTHER DATA : slitw != -1. */
          double slitw = mf_parameters_slit_width_xshooter(header);
          if (slitw != -1) {
              config->instrumental.slit_width.value = slitw;
          } else {
              //err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND, "KEYWORD[%s] : not found", MF_PARAMETERS_SLITW);
              cpl_msg_warning(cpl_func, "(mf_parameters) KEYWORD[%s] with key = %s : different of 'NONE' but it isn't in the header --> keep 'slitw' previous/default value = %g",
                              MF_PARAMETERS_SLIT_WIDTH_KEYWORD, config->instrumental.slit_width.key, config->instrumental.slit_width.value);
          }
      }
  }

  /* MF_PARAMETERS_PIXEL_SCALE */
  if (!err && strcmp(config->instrumental.pixel_scale.key, MF_PARAMETERS_KEYWORD_NONE) != 0) {

      err = (cpl_propertylist_has(header, config->ambient.observing_date.key) == CPL_TRUE) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_INPUT;
      if (!err) {

          prop  = cpl_propertylist_get_property_const(header, config->instrumental.pixel_scale.key);

          double value = cpl_property_get_double(prop);

          config->instrumental.pixel_scale.value = value;

          cpl_msg_info(cpl_func, "(mf_parameters) Get %10s from header KEYWORD[%s] = %g", MF_PARAMETERS_PIXEL_SCALE_KEYWORD, config->instrumental.pixel_scale.key, config->instrumental.pixel_scale.value);
      } else {
          cpl_error_set_message(cpl_func, err, "KEYWORD[%s] : not found", config->instrumental.pixel_scale.key);
      }
  }

  return err;
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief Deallocate the user parameter structure.
 *
 * @param config             mf_parameters_config structure for delete
 *
 */
/* ---------------------------------------------------------------------------*/
void mf_parameters_config_delete(
    mf_parameters_config     *config)
{
  if (config) {

      if (config->directories.telluric_path     ) cpl_free(config->directories.telluric_path);
      if (config->directories.telluriccorr_data_path) cpl_free(config->directories.telluriccorr_data_path);

      if (config->directories.tmp_path          ) cpl_free(config->directories.tmp_path);
      if (config->internal.tmp_folder           ) cpl_free(config->internal.tmp_folder);

/*      if (config->directories.output_path       ) cpl_free(config->directories.output_path);
      if (config->directories.output_name       ) cpl_free(config->directories.output_name); */

      if (config->inputs.column_lambda          ) cpl_free(config->inputs.column_lambda);
      if (config->inputs.column_flux            ) cpl_free(config->inputs.column_flux);
      if (config->inputs.column_dflux           ) cpl_free(config->inputs.column_dflux);
      if (config->inputs.column_mask            ) cpl_free(config->inputs.column_mask);
      if (config->inputs.wavelengths_frame      ) cpl_free(config->inputs.wavelengths_frame);
      if (config->inputs.observing_erv_rv.key   ) cpl_free(config->inputs.observing_erv_rv.key);

      //if (config->fitting.wlc_ref               ) cpl_free(config->fitting.wlc_ref);

      if (config->ambient.observing_date.key    ) cpl_free(config->ambient.observing_date.key);
      if (config->ambient.utc.key               ) cpl_free(config->ambient.utc.key);
      if (config->ambient.telescope_angle.key   ) cpl_free(config->ambient.telescope_angle.key);
      if (config->ambient.relative_humidity.key ) cpl_free(config->ambient.relative_humidity.key);
      if (config->ambient.pressure.key          ) cpl_free(config->ambient.pressure.key);
      if (config->ambient.temperature.key       ) cpl_free(config->ambient.temperature.key);
      if (config->ambient.mirror_temperature.key) cpl_free(config->ambient.mirror_temperature.key);
      if (config->ambient.elevation.key         ) cpl_free(config->ambient.elevation.key);
      if (config->ambient.longitude.key         ) cpl_free(config->ambient.longitude.key);
      if (config->ambient.latitude.key          ) cpl_free(config->ambient.latitude.key);

      if (config->instrumental.slit_width.key   ) cpl_free(config->instrumental.slit_width.key);
      if (config->instrumental.pixel_scale.key  ) cpl_free(config->instrumental.pixel_scale.key);

      if (config->atmospheric.ref_atm           ) cpl_free(config->atmospheric.ref_atm  );
      if (config->atmospheric.gdas_prof         ) cpl_free(config->atmospheric.gdas_prof);

      if (config->internal.molecules.lbl_molecs ) cpl_free(config->internal.molecules.lbl_molecs);

      cpl_free(config);
  }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Check the user parameter structure. For check the user modifications.
 *
 * @param config             mf_parameters_config structure for check
 *
 * @return cpl_error_code    CPL_ERROR_NONE if everything is OK or NULL in other case.
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_parameters_config_check(
    mf_parameters_config     *config,
    mf_io_lnfl_config        *lnfl_config)
{
  cpl_ensure(config, CPL_ERROR_NULL_INPUT, CPL_ERROR_NULL_INPUT);

  // TODO: Add check for all variables and values
  cpl_error_code err = CPL_ERROR_NONE;

  /* Check LNFL third party binary */
  if (!err) {

      if (!(config->directories.telluric_path)) err = CPL_ERROR_FILE_NOT_FOUND;
      else {

        if (lnfl_config->use_ODA) {
          if (!err) {
            char *bin_lnflODA;
            bin_lnflODA   = cpl_sprintf("%s/%s", config->directories.telluric_path, MF_BIN_LNFLODA);
            err = mf_io_access(bin_lnflODA);
            cpl_free(bin_lnflODA);
            if (err != CPL_ERROR_NONE) {
                cpl_msg_error(cpl_func, "External binary [%s], not found! in the path : %s", MF_BIN_LNFLODA, config->directories.telluric_path);
                err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_BIN_LNFLODA);
            }
          }
        }

        if (!err) {
            char *bin_lnfl = cpl_sprintf("%s/%s", config->directories.telluric_path, MF_BIN_LNFL);
            err = mf_io_access(bin_lnfl);
            cpl_free(bin_lnfl);
            if (err != CPL_ERROR_NONE) {
                cpl_msg_error(cpl_func, "External binary [%s], not found! in the path : %s", MF_BIN_LNFL, config->directories.telluric_path);
                err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_BIN_LNFL);
            }
        }

        if (!err) {
            char *bin_lblrtm = cpl_sprintf("%s/%s", config->directories.telluric_path, MF_BIN_LBLRTM);
            err = mf_io_access(bin_lblrtm);
            cpl_free(bin_lblrtm);
            if (err != CPL_ERROR_NONE) {
                cpl_msg_error(cpl_func, "External binary [%s], not found! in the path : %s", MF_BIN_LBLRTM, config->directories.telluric_path);
                err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_BIN_LBLRTM);
            }
        }
      }
  }

  if (!err && (    strcmp(config->inputs.wavelengths_frame, MF_PARAMETERS_WAVELENGTH_FRAME_VACUUM   )
                && strcmp(config->inputs.wavelengths_frame, MF_PARAMETERS_WAVELENGTH_FRAME_VACUUM_RV)
                && strcmp(config->inputs.wavelengths_frame, MF_PARAMETERS_WAVELENGTH_FRAME_AIR      )
				&& strcmp(config->inputs.wavelengths_frame, MF_PARAMETERS_WAVELENGTH_FRAME_AIR_RV   )) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %s; Invalid value", MF_PARAMETERS_WAVELENGTH_FRAME, config->inputs.wavelengths_frame);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_WAVELENGTH_FRAME);
  }


  if (!err && (   config->fitting.ftol < MF_PARAMETERS_FTOL_MIN) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_FTOL, config->fitting.ftol);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_FTOL);
  }

  if (!err && (   config->fitting.xtol < MF_PARAMETERS_XTOL_MIN) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_XTOL, config->fitting.xtol);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_XTOL);
  }

  if (!err && (   config->fitting.flux_unit != MF_PARAMETERS_FLUX_UNIT_NO_CONVERSION
               && config->fitting.flux_unit != MF_PARAMETERS_FLUX_UNIT_1
               && config->fitting.flux_unit != MF_PARAMETERS_FLUX_UNIT_2
               && config->fitting.flux_unit != MF_PARAMETERS_FLUX_UNIT_3            ) ){
      cpl_msg_info(cpl_func, "Config parameter[%s] = %d; Invalid value", MF_PARAMETERS_FLUX_UNIT, config->fitting.flux_unit);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_FLUX_UNIT);
  }


  if (!err && (   config->fitting.fit_telescope_background.const_val < MF_PARAMETERS_TELESCOPE_BACK_CONST_MIN
               || config->fitting.fit_telescope_background.const_val > MF_PARAMETERS_TELESCOPE_BACK_CONST_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_TELESCOPE_BACK_CONST, config->fitting.fit_telescope_background.const_val);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_TELESCOPE_BACK_CONST);
  }


  if (!err && (   config->fitting.fit_continuum.n < MF_FIT_N_POLYNOME_MIN
               || config->fitting.fit_continuum.n > MF_FIT_N_POLYNOME_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %lld; Invalid value", MF_PARAMETERS_CONTINUUM_N, config->fitting.fit_continuum.n);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_CONTINUUM_N);
  }


  if (!err && (   config->fitting.fit_wavelenght.n < MF_FIT_N_POLYNOME_MIN
               || config->fitting.fit_wavelenght.n > MF_FIT_N_POLYNOME_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %lld; Invalid value", MF_PARAMETERS_WLC_N, config->fitting.fit_wavelenght.n);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_WLC_N);
  }

  if (!err &&     (config->fitting.fit_wavelenght.const_val < MF_PARAMETERS_WLC_CONST_MIN
		  || config->fitting.fit_wavelenght.const_val > MF_PARAMETERS_WLC_CONST_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_WLC_CONST, config->fitting.fit_wavelenght.const_val);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_WLC_CONST);
  }

  /*if (!err && (    strcmp(config->fitting.wlc_ref, MF_PARAMETERS_WLC_REF_DATA )
                && strcmp(config->fitting.wlc_ref, MF_PARAMETERS_WLC_REF_MODEL)) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %s; Invalid value", MF_PARAMETERS_WLC_REF, config->fitting.wlc_ref);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_WLC_REF);
  }*/


  if (!err && (   config->fitting.fit_res_box.const_val < MF_PARAMETERS_RES_BOX_MIN
               || config->fitting.fit_res_box.const_val > MF_PARAMETERS_RES_BOX_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_RES_BOX, config->fitting.fit_res_box.const_val);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_RES_BOX);
  }

  if (!err && (   config->fitting.fit_gauss.const_val < MF_PARAMETERS_RES_GAUSS_MIN
               || config->fitting.fit_gauss.const_val > MF_PARAMETERS_RES_GAUSS_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_RES_GAUSS, config->fitting.fit_gauss.const_val);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_RES_GAUSS);
  }

  if (!err && (   config->fitting.fit_lorentz.const_val < MF_PARAMETERS_RES_LORENTZ_MIN
               || config->fitting.fit_lorentz.const_val > MF_PARAMETERS_RES_LORENTZ_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_RES_LORENTZ, config->fitting.fit_lorentz.const_val);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_RES_LORENTZ);
  }

  if (!err && (   config->fitting.kern_fac < MF_PARAMETERS_KERN_FAC_MIN
               || config->fitting.kern_fac > MF_PARAMETERS_KERN_FAC_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_KERN_FAC, config->fitting.kern_fac);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_KERN_FAC);
  }

  if (!err && (   config->ambient.telescope_angle.value < MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_MIN
               || config->ambient.telescope_angle.value > MF_PARAMETERS_TELESCOPE_ANGLE_VALUE_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_TELESCOPE_ANGLE_VALUE, config->ambient.telescope_angle.value);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_TELESCOPE_ANGLE_VALUE);
  }

  if (!err && (   config->ambient.relative_humidity.value < MF_PARAMETERS_RELATIVE_HUMIDITY_MIN
               || config->ambient.relative_humidity.value > MF_PARAMETERS_RELATIVE_HUMIDITY_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_RELATIVE_HUMIDITY_VALUE, config->ambient.relative_humidity.value);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_RELATIVE_HUMIDITY_VALUE);
  }

  if (!err && (   config->ambient.pressure.value < MF_PARAMETERS_PRESSURE_VALUE_MIN
               || config->ambient.pressure.value > MF_PARAMETERS_PRESSURE_VALUE_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_PRESSURE_VALUE, config->ambient.pressure.value);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_PRESSURE_VALUE);
  }

  if (!err && (   config->ambient.longitude.value < MF_PARAMETERS_LONGITUDE_VALUE_MIN
               || config->ambient.longitude.value > MF_PARAMETERS_LONGITUDE_VALUE_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_LONGITUDE_VALUE, config->ambient.longitude.value);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_LONGITUDE_VALUE);
  }

  if (!err && (   config->ambient.latitude.value < MF_PARAMETERS_LATITUDE_VALUE_MIN
               || config->ambient.latitude.value > MF_PARAMETERS_LATITUDE_VALUE_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %g; Invalid value", MF_PARAMETERS_LATITUDE_VALUE, config->ambient.latitude.value);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_LATITUDE_VALUE);
  }

  if (!err && (   config->internal.molecules.n_molec < MF_MOLECULES_NUMBER_MIN
               || config->internal.molecules.n_molec > MF_MOLECULES_NUMBER_MAX) ){
      cpl_msg_error(cpl_func, "Config parameter[%s] = %d; Invalid value", MF_MOLECULES_NUMBER, config->internal.molecules.n_molec);
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_MOLECULES_NUMBER);
  }


  return err;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Initializes the structure for the parameters with default values.
 *
 * @param config             Initialized mf_parameters_config structure
 * @param molecules          cpl_table defining molecules to fit. If NULL you need to provide the mf_molecules_config parameter != NULL
 *
 * @return mf_parameters     structure if everything is OK or NULL in other case.
 *
 */
/* ---------------------------------------------------------------------------*/
mf_parameters * mf_parameters_initialize(
    mf_parameters_config     *config,
    const cpl_table          *molecules)
{
  /* Initialize structrure */
  mf_parameters *params = cpl_malloc(sizeof(mf_parameters));

  /* Set configuration (default or provided) */
  params->config = config;

  /* Check with the polynomial degree of the wavelength fit */
  if (params->config->fitting.fit_wavelenght.n == 0) {
      params->config->fitting.fit_wavelenght.n = 1;
      params->config->internal.wlc_lin = CPL_FALSE;
  } else {
      params->config->internal.wlc_lin = CPL_TRUE;
  }

  /* Initialize timers */
  params->timers.time_start      = 0.;
  params->timers.time_fit        = 0.;
  params->timers.time_bin_lnfl   = 0.;
  params->timers.time_bin_lblrtm = 0.;
  params->timers.time_end        = 0.;

  /* Initialize the internal parameters based in the input molecules cpl_table */
  if (molecules) {
      params->molectab = cpl_table_duplicate(molecules);
      mf_molecules_lbl_fill(params->molectab, &(params->config->internal.molecules.n_molec), &(params->config->internal.molecules.lbl_molecs));
  } else {
      params->molectab = NULL;  /* Not need params->molecules in the convolution part */
  }

  /* Create table for range-related parameters */
  params->rangetab = cpl_table_new(1);
  cpl_table_new_column(      params->rangetab, MF_COL_CHIP,       CPL_TYPE_INT      );
  cpl_table_new_column(      params->rangetab, MF_COL_FIT_RANGE,  CPL_TYPE_INT      );
  cpl_table_new_column_array(params->rangetab, MF_COL_CONT_COEF,  CPL_TYPE_DOUBLE, 1);
  cpl_table_new_column(      params->rangetab, MF_COL_PIX_RES,    CPL_TYPE_DOUBLE   );
  cpl_table_new_column(      params->rangetab, MF_COL_WN_START,   CPL_TYPE_DOUBLE   );
  cpl_table_new_column(      params->rangetab, MF_COL_WN_END,     CPL_TYPE_DOUBLE   );
  cpl_table_new_column(      params->rangetab, MF_COL_WN_STEP,    CPL_TYPE_DOUBLE   );
  cpl_table_new_column(      params->rangetab, MF_COL_LNFL,       CPL_TYPE_STRING   );


  /* Table for chip-related parameters */
  params->chiptab = cpl_table_new(1);
  cpl_table_new_column(      params->chiptab,  MF_COL_FIT_CHIP,   CPL_TYPE_INT      );
  cpl_table_new_column_array(params->chiptab,  MF_COL_WLC_COEF,   CPL_TYPE_DOUBLE, 2);
  cpl_table_new_column(      params->chiptab,  MF_COL_WL_MIN,     CPL_TYPE_DOUBLE   );
  cpl_table_new_column(      params->chiptab,  MF_COL_WL_MAX,     CPL_TYPE_DOUBLE   );


  /* Check errors */
  if (cpl_error_get_code() != CPL_ERROR_NONE) {
      if (params->molectab) cpl_table_delete(params->molectab);
      if (params->rangetab) cpl_table_delete(params->rangetab);
      if (params->chiptab ) cpl_table_delete(params->chiptab );
      cpl_free(params);
      return NULL;
  }

  return params;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create default mf_parameters configuration.
 *
 * @param config             Initialized mf_parameters_config structure
 * @param molecules          cpl_table defining molecules to fit. If NULL you need to provide the mf_molecules_config parameter != NULL
 * @param tmp_folder         Temporary folder
 *
 * @return mf_parameters     structure if everything is OK or NULL in other case.
 *
 */
/* ---------------------------------------------------------------------------*/
mf_parameters * mf_parameters_create(
    mf_parameters_config     *config,
    const cpl_table          *molecules,
    const char               *tmp_folder)
{
  /* Initialize with the default configuration */
  mf_parameters *params = mf_parameters_initialize(config, molecules);
  if (!params) {
      return NULL;
  }

  /*** CHANGE VARIABLES TO TMP_DIR ***/

  /* Check tmp_folder */
  if (tmp_folder) {

      /* Check if exist previous tmp_folder */
      if (mf_io_access(tmp_folder) != CPL_ERROR_NONE) {

          cpl_msg_error(cpl_func, "(mf_parameters) Previous temporary workspace [%s] doesn't exist!", tmp_folder);
          mf_parameters_delete(params);
          return NULL;

      } else {

          /* Set as reused tmp_dir. The directory won't be deleted when mf_parameters_config will be deleted */
          if (params->config->internal.tmp_folder) cpl_free(params->config->internal.tmp_folder);
          params->config->internal.tmp_folder = cpl_strdup(tmp_folder);
          cpl_msg_info(cpl_func, "(mf_parameters) Reusing temporary workspace [%s]", tmp_folder);
      }

  } else {

      /* Not tmp_folder but check if the config was reused to free the variable */
      if (params->config->internal.tmp_folder) cpl_free(params->config->internal.tmp_folder);

      /* Make a temporary folder and check */
      params->config->internal.tmp_folder = mf_io_mkstemp(params->config->directories.tmp_path, MF_PARAMETERS_TMP_PATH_DEFAULT);
      if (params->config->internal.tmp_folder) {
          cpl_msg_info(cpl_func, "(mf_parameters) Created new temporary workspace in %s", params->config->internal.tmp_folder);
      } else {
          cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                "Temporary directory creation failed!");
          mf_parameters_delete(params);
          return NULL;
      }
  }

  return params;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Delete mf_parameters structure.
 *
 * @param params             mf_parameters structure for delete
 *
 */
/* ---------------------------------------------------------------------------*/
void mf_parameters_delete(
    mf_parameters            *params)
{
  if (params) {

      // NOTE: Not clean params->config is a external pointer and it cannot delete.
      params->config = NULL;

      if (params->molectab) cpl_table_delete(params->molectab);
      if (params->rangetab) cpl_table_delete(params->rangetab);
      if (params->chiptab ) cpl_table_delete(params->chiptab );

      cpl_free(params);
  }
}


/** @cond PRIVATE */

/* ---------------------------------------------------------------------------*/
/**
 * @brief Converts MJD into fractional year (return value, double precision).
 *
 * @param mjd                MJD.
 *
 * @return double            fractional year.
 *
 */
/* ---------------------------------------------------------------------------*/
static double mf_parameters_convert_obsdate_from_mjd_to_frac_year(
    double                   mjd)
{
    double refjd  = mjd - 51544.;
    double ncycle = floor(refjd / 1461.);
    double rest   = refjd - ncycle * 1461.;

    double frac;
    double nyears;
    if (rest < 366.) {
        frac = rest / 366.;
        nyears = 4 * ncycle;
    } else {
        frac = modf((rest - 366.) / 365., &nyears);
        nyears += 4 * ncycle + 1;
    }

    double date = 2000. + nyears + frac;

    return date;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Obtains slit width of XSHOOTER exposures from FITS header.
 *
 * @param header             cpl_property *header of one instrument
 *
 * @return cpl_error_code    Returns 'value' (slit width in arcsec) if XSHOOTER or -1 if the header of another instrument is provided
 *
 */
/* ---------------------------------------------------------------------------*/
static double mf_parameters_slit_width_xshooter(
    const cpl_propertylist   *header)
{
    /* Get properties */
    const cpl_property *prop;

    /* Check instrument */
    prop = cpl_propertylist_get_property_const(header, "INSTRUME");
    if (!prop) return -1;

    if (strstr(cpl_property_get_string(prop), "SHOOT") == NULL) {
        /* Instrument is not XSHOOTER */
        return -1;
    }

    prop = cpl_propertylist_get_property_const(header, "ESO SEQ ARM");
    if (!prop) return -1;

    /* Check XSHOOTER arm (UVB, VIS, or NIR) */
    int opti = 0;
    if (strstr(cpl_property_get_string(prop), "UVB") != NULL) {
        opti = 3;
    } else if (strstr(cpl_property_get_string(prop), "VIS") != NULL) {
        opti = 4;
    } else if (strstr(cpl_property_get_string(prop), "NIR") != NULL) {
        opti = 5;
    } else {
        /* Unexpected XSHOOTER arm */
        return -1;
    }

    /* Get name of fits keyword and write info message */
    char *slitw_key = cpl_sprintf("ESO INS OPTI%d NAME", opti);

    prop = cpl_propertylist_get_property_const(header, slitw_key);
    if (!prop) {
        cpl_free(slitw_key);
        return -1;
    }

    /* Get slit width for identified arm */
    char slitw_str[4];
    strncpy(slitw_str, cpl_property_get_string(prop), 3);
    slitw_str[3] = '\0';

    /* Info */
    cpl_msg_info(cpl_func, "(mf_parameters) Get slitw from key = %s , value = %s", slitw_key, slitw_str);
    cpl_free(slitw_key);

    /* Get value */
    double value = strtod(slitw_str, NULL);

    /* Check and return */
    if (value == 0) return -1;
    else            return value;
}

/** @endcond */


/**@}*/
