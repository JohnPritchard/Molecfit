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

#include "mf_parameters.h"
#include "mf_spectrum.h"
#include "mf_atm_combined_gdas.h"
#include "mf_kernel_user.h"
#include "mf_lnfl.h"
#include "mf_lblrtm.h"
#include "mf_convolution.h"

#include "mf_calctrans.h"

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

/* Initialize parameter mf_cacltrans_lblrtm_results */
static mf_calctrans_lblrtm_results * mf_calctrans_lblrtm_results_create(
    mf_parameters            *params);

/* Initialize parameter mf_cacltrans_convolution_results */
static mf_calctrans_convolution_results * mf_calctrans_convolution_results_create(void);

/* Set array with the best-fit parameters to mf_fit */
static cpl_error_code mf_calctrans_best_fit_parameters(
    cpl_array                *fitpar,
    const mf_parameters      *params,
    cpl_table                *best_fit_parameters);

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_calctrans       Calculate the telluric correction.
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
 * @brief lblrtm calctrans in temporary folder
 *
 * @param config               parameter that contains the configuration for the execution
 * @param spec_telluriccorr    array of input spectra, each entry is considered a chip
 * @param molecules            table defining molecules to fit
 * @param wl_start             lower wavelength override for LNFL -1 to use values from inspec pipeline units, converted via wlgtomicron parameter
 * @param wl_end               upper wavelength override for LNFL -1 to use values from inspec pipeline units, converted via wlgtomicron parameter
 * @param atm_parameters       cpl_table with the calculated atmospheric profile, result of mf_model
 * @param best_fit_parameters  cpl_table with the best_fit parameters, result of mf_model
 *
 * @return mf_calctrans_lblrtm_results   results if OK execution or NULL in error case.
 *
 */
/* ---------------------------------------------------------------------------*/
mf_calctrans_lblrtm_results * mf_calctrans_lblrtm(
    mf_configuration            *config,
    const cpl_table             *spec_telluriccorr,
    const cpl_table             *molecules,
    double                      wl_start,
    double                      wl_end,
    cpl_table                   *atm_parameters,
    cpl_table                   *best_fit_parameters)
{
    /* Check mandatory inputs in the external call */
    cpl_error_ensure(config,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Configuration parameter (config:*mf_configuration) == NULL");
    cpl_error_ensure(spec_telluriccorr,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Input spectrum (spec_telluriccorr:*cpl_table) == NULL");
    cpl_error_ensure(molecules,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Input molecule/s (molectab:*cpl_table) == NULL");
    cpl_error_ensure((wl_start == -1 && wl_end == -1) || (wl_start > 0 && wl_end > 0 && wl_start < wl_end),
                     CPL_ERROR_ILLEGAL_INPUT, return NULL,
                     "Illegal 'wl_start' = [%g] and/or 'wl_end' = [%g] wavelength", wl_start, wl_end);
    cpl_error_ensure(atm_parameters,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Input atmospheric profile (atmprof:*cpl_propertylist) == NULL");
    cpl_error_ensure(best_fit_parameters,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Input Best-fit parameter for the profile (res_table:*cpl_table) == NULL");

    /* Check configuration structures */
    cpl_error_code err = CPL_ERROR_NONE;
    err += mf_lblrtm_config_check(    config->lblrtm    );
    err += mf_lnfl_config_check(      config->lnfl      );
    err += mf_parameters_config_check(config->parameters, config->lnfl);
    cpl_error_ensure(!err,
                     cpl_error_get_code(), return NULL,
                     "Configuration structures failed: %s", cpl_error_get_message());

#ifdef _USE_OPENMP
    /* Set OMP_NUM_THREADS */
    /*cpl_msg_info(cpl_func,    "(mf_calctrans ) Setting telluriccorr OMP_NUM_THREADS = %i", config->parameters->inputs.omp_num_threads);
    omp_set_num_threads(config->parameters->inputs.omp_num_threads); */
#else
    cpl_msg_warning(cpl_func, "(mf_calctrans ) OPENMP unsupported by the compiler!");
#endif

    /* Initialize parameters in a temporary directory */
    mf_parameters *params = mf_parameters_create(config->parameters, molecules, NULL);
    cpl_error_ensure(params,
                     cpl_error_get_code(), return NULL,
                     "mf_parameters structure non-create, configuration problems: %s", cpl_error_get_message());

    /* Update total wavelength range to allow reusing existing lnfl result for multiple input spectra with different kernels [lnfl wave number cm^-1] */
    double wlg_to_micron = params->config->inputs.wlg_to_micron;
    double wn_start      = 1e4 / (wl_end   * wlg_to_micron);
    double wn_end        = 1e4 / (wl_start * wlg_to_micron);

    /* Get initial time */
    //params->timers.time_start = cpl_test_get_walltime();

    /* Read spectrum and header data from data file */
    //cpl_boolean correct_spectrum = CPL_TRUE;
    if (!err) {

        /* Change parameters to calculate transmission curve for full spectrum */
        //if (params->config->inputs.trans == MF_PARAMETERS_TRANS_EMISSION) correct_spectrum = CPL_FALSE;

        /* Update params and discart output cpl_table */
        cpl_table *spec_telluriccorr_format = mf_spectrum_ranges_apply(params, spec_telluriccorr, NULL, NULL, NULL);
        if (!spec_telluriccorr_format) err = CPL_ERROR_ILLEGAL_INPUT;
        else                       cpl_table_delete(spec_telluriccorr_format);

    }

    /* Adapt start and end of rangetab to user provided values used to execute lnfl only once for not fully overlapping spectra */
    if (!err && (wn_start > 0 && wn_end > wn_start)) {

        for (cpl_size i = 0; i < cpl_table_get_nrow(params->rangetab); i++) {

            double o_wnstart = cpl_table_get(params->rangetab, MF_COL_WN_START, i, NULL);
            double o_wnend   = cpl_table_get(params->rangetab, MF_COL_WN_END,   i, NULL);

            if (o_wnstart > 0 && o_wnend > 0) {
              cpl_table_set(params->rangetab, MF_COL_WN_START, i, wn_start);
              cpl_table_set(params->rangetab, MF_COL_WN_END,   i, wn_end  );
            }
        }

        err = cpl_error_get_code();
    }

    /* Read the summary of the fit results and write the fit parameters in a CPL array */
    cpl_array *fitpar = NULL;
    if (!err) {
        fitpar = cpl_array_new(0, CPL_TYPE_DOUBLE);
        err    = mf_calctrans_best_fit_parameters(fitpar, params, best_fit_parameters);
        if (!err && cpl_array_get_size(fitpar) == 0) {
            err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT, "Error read the convert res_table:cpl_table into fitpart:cpl_array");
        }
    }

    /* Create result parameters */
    mf_calctrans_lblrtm_results *results = mf_calctrans_lblrtm_results_create(params);
    err = cpl_error_get_code();
    if (!results && !err) err = CPL_ERROR_ILLEGAL_OUTPUT;

    /* Run LNFL */
    if (!err) err = mf_lnfl(config->lnfl, params);

    /* Run LBLRTM: Calculate model transmission curve : Set counter for mf_lblrtm calls to 1 */
    if (!err) {

        cpl_array *reffitpar   = cpl_array_duplicate(fitpar);
        int       mpfit_calls  = 1;
        int       lblrtm_calls = 0;
        cpl_table *atm_profile_out = NULL;

        /* Call mf_lblrtm */
        cpl_boolean last_call     = CPL_TRUE;
        cpl_size    run_execution = 0;          /* Not fit --> Unique execution */
        err = mf_lblrtm(run_execution,
                        params, last_call, results->spec_out, results->range_status,
                        reffitpar, mpfit_calls, &lblrtm_calls, config->lblrtm, config->lnfl,
                        atm_parameters, &atm_profile_out,
                        fitpar);

        /* Cleanup */
        cpl_array_delete(reffitpar);
        if (atm_profile_out) cpl_table_delete(atm_profile_out);
    }

    /* Cleanup */
    if (fitpar) cpl_array_delete(fitpar);

    if (!err) {

        /* Update state */
        params->timers.time_end = cpl_test_get_walltime();
        double runtime  = params->timers.time_end - params->timers.time_start;

        double time_external = params->timers.time_bin_lnfl + params->timers.time_bin_lblrtm;

        cpl_msg_info(cpl_func, "(mf_calctrans ) time_bin_lnfl      : %.2lf min",  params->timers.time_bin_lnfl   / 60.);
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_bin_lblrtm    : %.2lf min",  params->timers.time_bin_lblrtm / 60.);
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_bin_calls     : %.2lf min",            time_external        / 60.);
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_telluriccorr  : %.2lf min", (runtime - time_external)       / 60.);
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_total         : %.2lf min",  runtime                        / 60.);
    }

    /* Cleanup*/
    mf_parameters_delete(params);

    /* Check and return */
    if (!err) {
        return results;
    } else {
        mf_calctrans_lblrtm_results_delete(results);
        return NULL;
    }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Use the temporary results of lblrtm calctrans execution and apply it to the input spectrum
 *
 * @param config               Parameter that contains the configuration for the execution
 * @param lblrtm_results       Results of the lblrtm calctrans execution
 * @param header_spec          Property list containing values telluriccorr reads (ESO TEL/INS)
 * @param spec_telluriccorr    Array of input spectra, each entry is considered a chip
 * @param header_kernel        Property list containing the keywords of kernel
 * @param kernel               Kernel, one row per total input spectra pixel
 * @param wl_start             Lower wavelength override for LNFL -1 to use values from inspec pipeline units, converted via wlgtomicron parameter
 * @param wl_end               Upper wavelength override for LNFL -1 to use values from inspec pipeline units, converted via wlgtomicron parameter
 * @param best_fit_parameters  cpl_table with the best_fit parameters, result of mf_model
 *
 * @return mf_calctrans_convolution_results   results if OK execution or NULL in error case.
 *
 */
/* ---------------------------------------------------------------------------*/
mf_calctrans_convolution_results * mf_calctrans_convolution(
    mf_parameters_config        *config,
    mf_io_lnfl_config           *lnfl_config,
    mf_calctrans_lblrtm_results *lblrtm_results,
    const cpl_propertylist      *header_spec,
    const cpl_table             *spec_telluriccorr,
    const cpl_propertylist      *header_kernel,
    const cpl_matrix            *kernel,
    double                      wl_start,
    double                      wl_end,
    cpl_table                   *best_fit_parameters)
{
    /* Check mandatory inputs in the external call */
    cpl_error_ensure(config,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Configuration parameters (config:*mf_config_params) == NULL");
    cpl_error_ensure(lblrtm_results,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "result parameters (lblrtm_results:*mf_calctrans_lblrtm_results) of mf_calctrans_lblrtm() == NULL");
    cpl_error_ensure(spec_telluriccorr,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Input spectrum (spec_telluriccorr:*cpl_table) == NULL");
    cpl_error_ensure((wl_start == -1 && wl_end == -1) || (wl_start > 0 && wl_end > 0 && wl_start < wl_end),
                     CPL_ERROR_ILLEGAL_INPUT, return NULL,
                     "Illegal 'wl_start' = [%g] and/or 'wl_end' = [%g] wavelength", wl_start, wl_end);
    cpl_error_ensure(best_fit_parameters,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Input Best-fit parameter for the profile (res_table:*cpl_table) == NULL");

    /* Check configuration structures */
    cpl_error_code err = CPL_ERROR_NONE;
    err += mf_parameters_config_check(config, lnfl_config);
    cpl_error_ensure(!err,
                     cpl_error_get_code(), return NULL,
                     "Configuration structures failed: %s", cpl_error_get_message());

#ifdef _USE_OPENMP
    /* Set OMP_NUM_THREADS */
    /*cpl_msg_info(cpl_func,    "(mf_calctrans ) Setting telluriccorr OMP_NUM_THREADS = %i", config->inputs.omp_num_threads);
    omp_set_num_threads(config->inputs.omp_num_threads); */
#else
    cpl_msg_warning(cpl_func, "(mf_calctrans ) OPENMP unsupported by the compiler!");
#endif

    /* Initialize parameters in a temporary directory */
    mf_parameters *params = mf_parameters_create(config, NULL, lblrtm_results->tmp_folder);
    cpl_error_ensure(params,
                     cpl_error_get_code(), return NULL,
                     "mf_parameters structure non-create, configuration problems: %s", cpl_error_get_message());

    /* Update total wavelength range to allow reusing existing lnfl result for multiple input spectra with different kernels [lnfl wave number cm^-1] */
    double wlg_to_micron = params->config->inputs.wlg_to_micron;
    double wn_start      = 1e4 / (wl_end   * wlg_to_micron);
    double wn_end        = 1e4 / (wl_start * wlg_to_micron);


    /* Get initial time */
    params->timers.time_start = cpl_test_get_walltime();

    /* Create result parameters */
    mf_calctrans_convolution_results *results = mf_calctrans_convolution_results_create();
    if (!results) err = cpl_error_get_code();

    /* Read spectrum */
    cpl_boolean correct_spectrum = CPL_TRUE;
    if (!err && results) {

        /* Change parameters to calculate transmission curve for full spectrum */
        if (params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE) correct_spectrum = CPL_FALSE;

        /* Make output cpl_table telluriccorr spectrum from header and spec and update params */
        results->spec_telluriccorr_format = mf_spectrum_ranges_apply(params, spec_telluriccorr, NULL, NULL, NULL);
        if (!(results->spec_telluriccorr_format)) err = CPL_ERROR_ILLEGAL_INPUT;
    }

    /* Create internal kernel */
    if (!err && kernel) {

        cpl_size spec_n_lambdas = cpl_table_get_nrow(spec_telluriccorr);
        results->kernel_resampled_normalized = mf_kernel_user_create(spec_n_lambdas, header_spec, header_kernel, kernel);
        if (!(results->kernel_resampled_normalized)) err = CPL_ERROR_ILLEGAL_INPUT;
    }

    /* Adapt start and end of rangetab to user provided values used to execute lnfl only once for not fully overlapping spectra */
    if (!err && (wn_start > 0 && wn_end > wn_start)) {

        for (cpl_size i = 0; i < cpl_table_get_nrow(params->rangetab); i++) {

            double o_wnstart = cpl_table_get(params->rangetab, MF_COL_WN_START, i, NULL);
            double o_wnend   = cpl_table_get(params->rangetab, MF_COL_WN_END,   i, NULL);

            if (o_wnstart > 0 && o_wnend > 0) {
              cpl_table_set(params->rangetab, MF_COL_WN_START, i, wn_start);
              cpl_table_set(params->rangetab, MF_COL_WN_END,   i, wn_end);
            }
        }

        err = cpl_error_get_code();
    }

    /* Read the summary of the fit results and write the fit parameters in a CPL array */
    cpl_array *fitpar = NULL;
    if (!err) {
        fitpar = cpl_array_new(0, CPL_TYPE_DOUBLE);
        err    = mf_calctrans_best_fit_parameters(fitpar, params, best_fit_parameters);
        if (!err && cpl_array_get_size(fitpar) == 0) {
            err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT, "Error read the convert res_table:cpl_table into fitpart:cpl_array");
        }
    }

    /* Calculate the Convolution : Set counter for mf_convolution calls to 1 */
    if (!err) {

        /* Kernel spectrum convolution */
        if (results->kernel_resampled_normalized) cpl_msg_info(cpl_func, "(mf_calctrans ) Using user-defined kernel");
        else                 cpl_msg_info(cpl_func, "(mf_calctrans ) Using synthetic kernel"   );

        cpl_boolean last_call   = CPL_TRUE;
        int         mpfit_calls = 1;

        err = mf_convolution(params, correct_spectrum, last_call, lblrtm_results->spec_out, lblrtm_results->range_status,
                             fitpar, mpfit_calls,
                             results->spec_telluriccorr_format,
                             results->kernel_resampled_normalized);
    }

    /* Check kind of spectrum */
    int transmission = params->config->inputs.transmission;
    if (!err) {

        if (transmission == MF_PARAMETERS_TRANSMISSION_TRUE) {
            cpl_msg_info(cpl_func, "(mf_calctrans ) Perform telluric absorption correction");
        } else if (!correct_spectrum) {
            cpl_msg_info(cpl_func, "(mf_calctrans ) Sky emission spectrum: no flux correction");
        } else if (transmission == MF_PARAMETERS_TRANSMISSION_FALSE) {
            err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                        "Invalid object value(s): trans of mf_parameters *params (not transmission as required)");
        }
    }

    /* Check whether there was no continuum correction */
    if (!err) {
        double scale_min = cpl_table_get_column_min(results->spec_telluriccorr_format, MF_COL_MOD_SCALE);
        double scale_max = cpl_table_get_column_max(results->spec_telluriccorr_format, MF_COL_MOD_SCALE);
        if (scale_max - scale_min > MF_TOL) {
            err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                        "Invalid object value(s): column mscal of cpl_table *spec != 1 (unwanted continuum correction)");
        }
    }

    /* Corrects flux for telluric absorption and writes it in new column MF_COL_CFLUX */
    if (!err) {

        /* Rename flux column */
        cpl_table_name_column(results->spec_telluriccorr_format, MF_COL_MOD_FLUX, MF_COL_OUT_TELLURIC_CORR);

        /* Remove unnecessary columns in spectrum table */
        cpl_table_erase_column(results->spec_telluriccorr_format, MF_COL_MOD_RANGE);
        cpl_table_erase_column(results->spec_telluriccorr_format, MF_COL_MOD_SCALE);
        cpl_table_erase_column(results->spec_telluriccorr_format, MF_COL_DEV      );

        /* The quality of the correction is indicated by the column MF_COL_C_MASK (0 = bad, 1 = acceptable) */
        err = mf_calctrans_correct_spectrum(results->spec_telluriccorr_format, transmission);
    }

    if (!err) {

        /* Update state */
        params->timers.time_end = cpl_test_get_walltime();
        double runtime  = params->timers.time_end - params->timers.time_start;

        double time_external = params->timers.time_bin_lnfl + params->timers.time_bin_lblrtm;

        /* Show info */
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_bin_lnfl      : %.2lf min",  params->timers.time_bin_lnfl   / 60.);
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_bin_lblrtm    : %.2lf min",  params->timers.time_bin_lblrtm / 60.);
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_bin_calls     : %.2lf min",            time_external        / 60.);
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_telluriccorr  : %.2lf min", (runtime - time_external)       / 60.);
        cpl_msg_info(cpl_func, "(mf_calctrans ) time_total         : %.2lf min",  runtime                        / 60.);
    }

    /* Cleanup */
    mf_parameters_delete(params);
    if (fitpar) cpl_array_delete(fitpar);

    /* Check and return */
    if (!err) {
        return results;
    } else {
        mf_calctrans_convolution_results_delete(results);
        return NULL;
    }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Cleanup the mf_calctrans_lblrtm_results structure
 *
 * @param results     Parameter structure to clean
 *
 */
/* ---------------------------------------------------------------------------*/
void mf_calctrans_lblrtm_results_delete(
    mf_calctrans_lblrtm_results *results)
{
  if (results) {

      /* Remove telluriccorr temporary files */
      if (results->tmp_folder) {
          mf_io_rm_rf(results->tmp_folder, 10);
          cpl_free(results->tmp_folder  );
      }

      if (results->spec_out) {

          for (cpl_size i = 0; i < results->n_range; i++) {
              if (results->spec_out[i]) {
                  cpl_table_delete(results->spec_out[i]);
              }
          }

          cpl_free(results->spec_out);
      }

      if (results->range_status) cpl_free(results->range_status);

      cpl_free(results);
  }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Cleanup the mf_calctrans_convolution_results structure
 *
 * @param results     Parameter structure to clean
 *
 */
/* ---------------------------------------------------------------------------*/
void mf_calctrans_convolution_results_delete(
    mf_calctrans_convolution_results *results)
{
  if (results) {

      if (results->spec_telluriccorr_format       ) cpl_table_delete( results->spec_telluriccorr_format       );
      if (results->kernel_resampled_normalized) cpl_matrix_delete(results->kernel_resampled_normalized);

      cpl_free(results);
  }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Corrects the observed flux (and error if present) for telluric absorption.
 *
 * @param spec               in/out: CPL table with observed flux and modeled transmission curve.
 * @param trans              Transmission (= 1) or emission (= 0 or 2).
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *
 * @note The corrected fluxes (and errors) are given in column MF_COL_CFLUX (and MF_COL_CDFLUX).
 *              The column MF_COL_C_MASK indicates whether the correction is relatively reliable for a pixel (= 1) or whether it probably failed (= 0).
 *              The minimum transmission mf_TRANS_MIN is used for the classification.
 *              The columns MF_COL_CFLUX and MF_COL_C_MASK are created if they do not exist.
 *              The existence of the input columns MF_COL_FLUX and MF_COL_MOD_TRANS is mandatory.
 *              The column MF_COL_CDFLUX is only created if it does not exist and an input column labelled MF_COL_DFLUX is present.
 *              No flux correction is carried out in the case of a sky emission spectrum (trans != 1).
 *              In this case, the column(s) MF_COL_CFLUX (and MF_COL_CDFLUX) contain(s) the uncorrected input flux (and error respectively).
 *
 * @note The in/out table with columns for corrected flux, (error,) and quality flag
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_calctrans_correct_spectrum(
    cpl_table                *spec,
    const cpl_boolean        trans)
{
    /* Check data points */
    cpl_size nrow = cpl_table_get_nrow(spec);
    if (nrow <= 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "No data: cpl_table *spec");
    }

    /* Check General columns */
    if (   cpl_table_has_column(spec, MF_COL_IN_FLUX          ) != 1
        || cpl_table_has_column(spec, MF_COL_OUT_TELLURIC_CORR) != 1) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                     "Minimun columns not included inside telluriccorr spectrum table (%s,%s)",
                                     MF_COL_IN_FLUX, MF_COL_OUT_TELLURIC_CORR);
    }

    /* Create output columns if they don't exist */
    if (cpl_table_has_column(spec, MF_COL_OUT_FLUX) != 1) {
        cpl_table_duplicate_column(spec, MF_COL_OUT_FLUX, spec, MF_COL_IN_FLUX);
    }
    if (cpl_table_has_column(spec, MF_COL_OUT_MASK) != 1) {
        cpl_table_new_column(spec, MF_COL_OUT_MASK, CPL_TYPE_INT  );
    }


    /* Get FLUX, DFLUX (error flux) and Telluric correction data pointers */
    double *in_flux      = cpl_table_get_data_double(spec, MF_COL_IN_FLUX          );
    double *out_flux     = cpl_table_get_data_double(spec, MF_COL_OUT_FLUX         );
    double *telluric_cor = cpl_table_get_data_double(spec, MF_COL_OUT_TELLURIC_CORR);

    /* Get criteria for minimum observed flux and transmission */
    double flux_min  = MF_CALCTRANS_MINIMUN_OBS_FLUX * cpl_table_get_column_median(spec, MF_COL_IN_FLUX);
    double trans_min = MF_CALCTRANS_MINIMUN_OBS_FLUX;

    for (cpl_size i = 0; i < nrow; i++) {

        if (trans == MF_PARAMETERS_TRANSMISSION_FALSE) {

            out_flux[i]  = in_flux[i];

            if (telluric_cor[i] < trans_min) cpl_table_set(spec, MF_COL_OUT_MASK, i, CPL_BINARY_0);
            else                             cpl_table_set(spec, MF_COL_OUT_MASK, i, CPL_BINARY_1);

        } else {

            /* Only correct spectrum in MF_PARAMETERS_TRANS_TRANSMISSION : Divide model spectrum by transmission curve */
            if (telluric_cor[i] < trans_min) {

                 /* Avoid extreme correction factors */
                 out_flux[i]  = in_flux[i]  / trans_min;

                 cpl_table_set(spec, MF_COL_OUT_MASK, i, 0);

            } else {

                 /* Normal correction */
                 out_flux[i]  = in_flux[i]  / telluric_cor[i];

                 if (in_flux[i] < flux_min) cpl_table_set(spec, MF_COL_OUT_MASK, i, CPL_BINARY_0);
                 else                       cpl_table_set(spec, MF_COL_OUT_MASK, i, CPL_BINARY_1);
            }
        }
    }

    /* DFLUX code correction, if it is needed */
    if (cpl_table_has_column(spec, MF_COL_IN_DFLUX) == 1) {

        /* Create output column if it doesn't exist */
        if (cpl_table_has_column(spec, MF_COL_OUT_DFLUX) != 1) {
            cpl_table_duplicate_column(spec, MF_COL_OUT_DFLUX, spec, MF_COL_IN_DFLUX);
        }

        double *in_dflux  = cpl_table_get_data_double(spec, MF_COL_IN_DFLUX );
        double *out_dflux = cpl_table_get_data_double(spec, MF_COL_OUT_DFLUX);

        for (cpl_size i = 0; i < nrow; i++) {
            if (     trans != MF_PARAMETERS_TRANSMISSION_TRUE) out_dflux[i] = in_dflux[i];
            else if (telluric_cor[i] < trans_min              ) out_dflux[i] = in_dflux[i] / trans_min;
            else                                                out_dflux[i] = in_dflux[i] / telluric_cor[i];
        }
    }

    return cpl_error_get_code();
}


/** @cond PRIVATE */

/* ---------------------------------------------------------------------------*/
/**
 * @brief Initialize parameter mf_cacltrans_lblrtm_results
 *
 * @param params     Parameter structure
 *
 * @return mf_calctrans_lblrtm_results     Initializate structure
 *
 */
/* ---------------------------------------------------------------------------*/
static mf_calctrans_lblrtm_results * mf_calctrans_lblrtm_results_create(
    mf_parameters            *params)
{
  mf_calctrans_lblrtm_results *results = cpl_malloc(sizeof(mf_calctrans_lblrtm_results));

  results->tmp_folder      = cpl_sprintf("%s", params->config->internal.tmp_folder);

  results->n_range         = params->config->internal.n_range;

  results->spec_out        = cpl_calloc(results->n_range, sizeof(cpl_table *));
  for (cpl_size i = 0; i < results->n_range; i++) {
      results->spec_out[i] = NULL;
  }

  results->range_status    = cpl_calloc(results->n_range, sizeof(cpl_error_code));

  return results;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Initialize parameter mf_cacltrans_convolution_results
 *
 * @param params     Parameter structure
 *
 * @return mf_calctrans_convolution_results     Initializate structure
 *
 */
/* ---------------------------------------------------------------------------*/
static mf_calctrans_convolution_results * mf_calctrans_convolution_results_create(void)
{
  mf_calctrans_convolution_results *results = cpl_malloc(sizeof(mf_calctrans_convolution_results));

  results->spec_telluriccorr_format        = NULL;
  results->kernel_resampled_normalized = NULL;

  return results;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Set array with the best-fit parameters to mf_fit.
 *
 * @param fitpar             out: cpl_array with the fit_parameters
 * @param params             Parameter structure
 * @param res_table          Cpl_table with the best_fit parameters, result of mf_model
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_calctrans_best_fit_parameters(
    cpl_array                *fitpar,
    const mf_parameters      *params,
    cpl_table                *best_fit_parameters)
{
    int nmolec =     params->config->internal.molecules.n_molec;
    int nchip  =     params->config->internal.nchip;
    int ncont  = 1 + params->config->fitting.fit_continuum.n;
    int nwlc   = 1 + params->config->fitting.fit_wavelenght.n;

    /* Get number of fit parameters */
    int npar = nmolec + ((nwlc + ncont) * nchip) + 3;

    /* Set size of output CPL array */
    cpl_array_set_size(fitpar, npar);
    cpl_array_fill_window_double(fitpar, 0, nmolec, 1.);

    /* Set constant of continuum correction function to 1 and all other coefficients to 0 to obtain non-scaled transmission curve */
    int contparmin = nmolec     + (nchip * nwlc );
    int contparmax = contparmin + (nchip * ncont) - 1;
    for (int j = -1, i = contparmin; i <= contparmax; i++) {

        j++;

        if (j == ncont) j = 0;

        if (j == 0    ) cpl_array_set(fitpar, i, 1);
        else            cpl_array_set(fitpar, i, 0);
    }

    int    wlcparmin    = nmolec;
    int    wlcparmax    = wlcparmin + (nchip * nwlc) - 1;
    int    fitCoefIndex = wlcparmin;

    int    status       = -9999;
    double boxfwhm      = -1;
    double gaussfwhm    = -1;
    double lorentzfwhm  = -1;

    cpl_msg_info(cpl_func, "(mf_calctrans ) wlcparmin = %d", wlcparmin);
    cpl_msg_info(cpl_func, "(mf_calctrans ) wlcparmax = %d", wlcparmax);
    cpl_msg_info(cpl_func, "(mf_calctrans ) nwlc      = %d", nwlc);

    for (unsigned index = 0; index < cpl_table_get_nrow(best_fit_parameters); ++index) {

        char const *par_name = cpl_table_get_string(best_fit_parameters, MF_COL_PARAMETER, index);
        double     value     = cpl_table_get_double(best_fit_parameters, MF_COL_VALUE, index, NULL);

        if (!strcmp(MF_MODEL_FIT_STATUS, par_name)) {

            status = (int)value;
            if (status < 1) {
                cpl_array_set_size(fitpar, 0);
                cpl_msg_info(cpl_func, "(mf_calctrans ) No fit results available -> Return");
                return CPL_ERROR_NONE;
            }
            cpl_msg_info(cpl_func, "(mf_calctrans ) Read mf_model(...) fit parameters[index=%02d], %15s --> %15s : %d",                   index, "STATUS",       par_name, status);

        } else if (!strcmp(MF_MODEL_FIT_BOX_FWHM, par_name)) {

            boxfwhm = value;
            cpl_array_set(fitpar, npar - 3, boxfwhm);
            cpl_msg_info(cpl_func, "(mf_calctrans ) Read mf_model(...) fit parameters[index=%02d], %15s --> %15s : %g",                   index, "BOX_FWHM",     par_name, boxfwhm);

        } else if (!strcmp(MF_MODEL_FIT_GAUSS_FWHM, par_name)) {

            gaussfwhm = value;
            cpl_array_set(fitpar, npar - 2, gaussfwhm);
            cpl_msg_info(cpl_func, "(mf_calctrans ) Read mf_model(...) fit parameters[index=%02d], %15s --> %15s : %g",                   index, "GAUSS_FWHM",   par_name, gaussfwhm);

        } else if (!strcmp(MF_MODEL_FIT_LORENTZFWHM, par_name)) {

            lorentzfwhm = value;
            cpl_array_set(fitpar, npar - 1, lorentzfwhm);
            cpl_msg_info(cpl_func, "(mf_calctrans ) Read mf_model(...) fit parameters[index=%02d], %15s --> %15s : %g",                   index, "LORENTZ_FWHM", par_name, lorentzfwhm);

        } else if (!strncmp(MF_COL_CHIP, par_name, strlen(MF_COL_CHIP))) {

            cpl_msg_info(cpl_func, "(mf_calctrans ) Read mf_model(...) fit parameters[index=%02d], %15s --> %15s : %g [fitCoefIndex=%d]", index, "FIT_CHIP",     par_name, value, fitCoefIndex);

            cpl_array_set(fitpar, fitCoefIndex++, value);
            if(fitCoefIndex > wlcparmax + 1) {
                return cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                                             "Unexpected best-fit structure: Too many coef. for wavelength solution");
            }

        } else {
            //cpl_msg_info(cpl_func, "(mf_calctrans ) Read mf_model(...) fit parameters[index=%02d, not used] --> %s : %g", index, par_name, value);
        }
    }

    if (status == -9999 || boxfwhm < 0 || lorentzfwhm < 0 || fitCoefIndex == wlcparmin) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                                     "Unexpected best-fit structure: Coef. for wavelength solution not found");
    }

    return CPL_ERROR_NONE;
}

/** @endcond */


/**@}*/
