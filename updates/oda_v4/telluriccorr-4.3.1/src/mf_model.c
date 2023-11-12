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
#include "mf_molecules.h"
#include "mf_spectrum.h"
#include "mf_gdas.h"
#include "mf_atm_combined_gdas.h"
#include "mf_kernel_user.h"
#include "mf_fit.h"
#include "mf_lnfl.h"
#include "mf_lblrtm.h"

#include "mf_model.h"

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

/*  */
static mf_model_results * mf_model_results_create(
    const mf_parameters      *params);

/*  */
static mf_fit_results * mf_model_batch(
    mf_io_lblrtm_config      *lblrtm_config,
    mf_io_lnfl_config        *lnfl_config,
    mf_parameters            *params,
    mf_model_results         *results);

/*  */
static cpl_error_code mf_model_parameters_set_fit_flags(
    mf_parameters            *params,
    const cpl_array          *fit_molec,
    const cpl_array          *fit_cont,
    const cpl_array          *fit_wlc,
    const cpl_array          *fit_res);

/*  */
static cpl_error_code mf_model_calculate_func_telluric_absortion_correction(
    cpl_table                *spec,
    const mf_parameters      *params);

/*  */
static cpl_table * mf_model_res_create(
    mf_parameters            *params,
    cpl_table                *spec,
    cpl_table                *atm_profile,
    const mf_fit_results     *fit_results);

/*  */
static void mf_model_fill_parameter(
    cpl_table                *outtable,
    cpl_size                 index,
    const char               *parameter,
    const double             val,
    const double             uncertainty);

/* Calculates water vapor column in mm from profile in ppmv */
static cpl_error_code mf_model_calculate_atm_colums(
    cpl_table                *atm_profile,
    double                   *h2ocol,
    mf_parameters            *params);

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_model       Calculate the atmospheric profile with the best-fit parameters and model.
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
 * @brief Run telluriccorr in temporary folder and return the structure results.
 *
 * @param config                .
 * @param gdas_user             .
 * @param atm_profile_standard  .
 * @param atm_profile_combined  .
 * @param header_spec           property list containing values that telluriccorr reads (ESO TEL/INS)
 * @param spec_telluriccorr     array of input spectra, each entry is considered a chip
 * @param header_kernel         Property list containing the keywords of kernel
 * @param kernel                kernel, one row per total input spectra pixel
 * @param molecules             table defining molecules to fit
 * @param inc_wranges           wavelength inclusion table
 * @param exc_wranges           wavelength exclusion table
 * @param exc_pranges           pixel exclusion table
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
mf_model_results * mf_model(
    mf_configuration         *config,
    const cpl_table          *molecules,
    const cpl_propertylist   *header_spec,
    const cpl_table          *spec_telluriccorr,
    const cpl_table          *inc_wranges,
    const cpl_table          *exc_wranges,
    const cpl_table          *exc_pranges,
    const cpl_propertylist   *header_kernel,
    const cpl_matrix         *kernel,
    const cpl_table          *gdas_user,
    const cpl_table          *atm_profile_standard,
    const cpl_table          *atm_profile_combined)
{
    /* Check mandatory inputs in the external call */
    cpl_error_ensure(config,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Configuration parameters == NULL");
    cpl_error_ensure(molecules,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Input molecule/s (molectab:*cpl_table) == NULL");
    cpl_error_ensure(header_spec && spec_telluriccorr,
                     CPL_ERROR_NULL_INPUT, return NULL,
                     "Input spectrum (header_spec:*cpl_propertylist || spec:*cpl_table) == NULL");

    /* Check configuration structures */
    cpl_error_code err = CPL_ERROR_NONE;

    config->lnfl->use_ODA = CPL_TRUE;
    char* oda_option=getenv("MF_ODA_DISABLE");
    if (oda_option!=NULL) {
        if (strncasecmp( oda_option, "true", 4 ) == 0) {
            config->lnfl->use_ODA = CPL_FALSE;
        }
    }

    err += mf_lblrtm_config_check(    config->lblrtm    );
    err += mf_lnfl_config_check(      config->lnfl      );
    err += mf_parameters_config_check(config->parameters, config->lnfl);
    cpl_error_ensure(!err,
                     cpl_error_get_code(), return NULL,
                     "Configuration structures failed: %s", cpl_error_get_message());

    /* Check molecules cpl_table */
    err = mf_molecules_check(molecules);
    cpl_error_ensure(!err,
                     cpl_error_get_code(), return NULL,
                     "Configuration structures failed: %s", cpl_error_get_message());

    /* Check wave ranges cpl_tables */
    if (inc_wranges) err += mf_spectrum_ranges_check(spec_telluriccorr, inc_wranges, MF_INPUT_WAVE_INCLUDE );
    if (exc_wranges) err += mf_spectrum_ranges_check(spec_telluriccorr, exc_wranges, MF_INPUT_WAVE_EXCLUDE );
    if (exc_pranges) err += mf_spectrum_ranges_check(spec_telluriccorr, exc_pranges, MF_INPUT_PIXEL_EXCLUDE);

    cpl_error_ensure(!err,
                     cpl_error_get_code(), return NULL,
                     "Configuration structures failed: %s", cpl_error_get_message());

    /* Check user GDAS */
    if (gdas_user) err = mf_gdas_check(gdas_user);
    cpl_error_ensure(!err,
                     cpl_error_get_code(), return NULL,
                     "Configuration structures failed: %s", cpl_error_get_message());

    /* Check user ATM_PROFILES */
    if (atm_profile_standard) err += mf_atm_profile_standard_check(atm_profile_standard);
    if (atm_profile_combined) err += mf_atm_profile_combined_check(atm_profile_combined);
    cpl_error_ensure(!err,
                     cpl_error_get_code(), return NULL,
                     "Configuration structures failed: %s", cpl_error_get_message());


#ifdef _USE_OPENMP
    /* Set OMP_NUM_THREADS */
    /*cpl_msg_info(cpl_func,    "(mf_model     ) Setting telluriccorr OMP_NUM_THREADS = %i", config->parameters->inputs.omp_num_threads);
    omp_set_num_threads(config->parameters->inputs.omp_num_threads);*/
#else
    cpl_msg_warning(cpl_func, "(mf_model     ) OPENMP unsupported by the compiler!");
#endif


    /* Initialize parameters in a temporary directory */
    mf_parameters *params = mf_parameters_create(config->parameters, molecules, NULL);
    cpl_error_ensure(params,
                     cpl_error_get_code(), return NULL,
                     "mf_parameters structure non-create, configuration problems: %s", cpl_error_get_message());

    /* Get initial time */
    params->timers.time_start = cpl_test_get_walltime();

    /* Create result parameters */
    mf_model_results *results = mf_model_results_create(params);

    /* Make output cpl_table telluriccorr spectrum from header and spec and update params */
    results->spec = mf_spectrum_ranges_apply(params, spec_telluriccorr, inc_wranges, exc_wranges, exc_pranges);
    if (!(results->spec)) err = CPL_ERROR_NULL_INPUT;
    //cpl_table_dump( results->spec, 0,cpl_table_get_nrow(results->spec), NULL );

    /* Create internal kernel */
    if (!err && header_kernel && kernel) {

        cpl_size spec_n_lambdas = cpl_table_get_nrow(spec_telluriccorr);
        results->kernel = mf_kernel_user_create(spec_n_lambdas, header_spec, header_kernel, kernel);
        if (!(results->kernel)) err = CPL_ERROR_ILLEGAL_INPUT;
    }

    /* Create atmospheric profiles from reference profiles and GDAS data */
    if (!err) {
        
        /* Check if atm_profile_combined provided and keep time */
        if (atm_profile_combined) {

            cpl_msg_info(cpl_func, "(mf_atm       ) ATM_PROFILE_COMBINED user provided (ATM_PROFILE_STANDARD merged with GDAS_USER)");
            results->atm_profile_combined = cpl_table_duplicate(atm_profile_combined);
            results->gdas_src = "USER";

        } else {

            results->atm_profile_combined = mf_atm_combined_gdas(params,
                                                                 gdas_user,
                                                                 atm_profile_standard,
                                                                 &(results->gdas_before), &(results->gdas_after),
                                                                 &(results->gdas_interpolate),
                                                                 &(results->atm_profile_standard),
                                                                 &(results->gdas_src));
        }

        err = cpl_error_get_code();
        if (!(results->atm_profile_combined)) err = CPL_ERROR_UNSPECIFIED;
    }

    /* Run LNFL */
    if (!err) err = mf_lnfl(config->lnfl, params);

    if (config->lnfl->use_ODA) {
        /* Run ODA LBLRTM */
        if (!err) err = mf_io_lblrtm_oda(config->lnfl,
                                        config->lblrtm,
                                        params,
                                        results->atm_profile_combined);
    }

    /* Perform fitting procedure */
    if (!err) {

        /* Kernel spectrum convolution */
        if (results->kernel) cpl_msg_info(cpl_func, "(mf_model     ) Using user-defined kernel");
        else                 cpl_msg_info(cpl_func, "(mf_model     ) Using synthetic kernel"   );


        /*** CALL FIT SEQUENCE ***/
        cpl_msg_info(cpl_func, "(mf_model     ) Number of ranges = %d", params->config->internal.n_range);


        /* Allocate array with the error code status in each range */
        results->n_range         = params->config->internal.n_range;
        results->spec_out        = cpl_calloc(results->n_range, sizeof(cpl_table *));
        for (cpl_size i = 0; i < results->n_range; i++) {
            results->spec_out[i] = NULL;
        }


        results->range_status    = cpl_calloc(results->n_range, sizeof(cpl_error_code));

	cpl_table *tmp_table = params->rangetab;
        cpl_size tmp_size= cpl_table_get_nrow(tmp_table);
        cpl_msg_info(cpl_func,"<---------------------DUMPING RANGE TABLE HERE size %lld------------------->",tmp_size);
        //cpl_table_dump(tmp_table,0,tmp_size,NULL);
        
        for (cpl_size i = 0; i<tmp_size; i++) {
            const cpl_array* tmp_array;
            tmp_array=cpl_table_get_array (tmp_table,"cont_coef",i);
            cpl_size asize=cpl_array_get_size(tmp_array);
            cpl_array_dump(tmp_array,0,asize,NULL);
        }    
	 tmp_table = params->chiptab;
         tmp_size= cpl_table_get_nrow(tmp_table);
         //cpl_table_dump(tmp_table,0,tmp_size,NULL);
         
         for (cpl_size i = 0; i<tmp_size; i++) {
            const cpl_array* tmp_array;
            tmp_array=cpl_table_get_array (tmp_table,"wlc_coef",i);
            cpl_size asize=cpl_array_get_size(tmp_array);
            cpl_array_dump(tmp_array,0,asize,NULL);
        }    
        
        cpl_msg_info(cpl_func,"<---------------------DUMPING RANGE TABLE TO HERE------------------->");



        double         ts           = cpl_test_get_walltime();
        mf_fit_results *fit_results = mf_model_batch(config->lblrtm, config->lnfl, params, results);
        double         te           = cpl_test_get_walltime();
        double         runtime      = (te - ts) / 60;

        cpl_msg_info(cpl_func, "(mf_model     ) TIME FIT SEQUENCE  : %.2f  min", runtime);

        /* Write a summary of the fit results into an ASCII file in output folder */
        err = cpl_error_get_code();
	
	
	
        if (fit_results) {

            if (!err) {

                double ratio_lblrtm_call = params->timers.time_bin_lblrtm / (double)(fit_results->lblrtm_calls);

                cpl_msg_info(cpl_func, "(mf_model     ) time_fit           : %.2lf min", params->timers.time_fit / 60.);

                cpl_msg_info(cpl_func, "(mf_model     ) lblrtm_calls       : %d",        fit_results->lblrtm_calls    );
                cpl_msg_info(cpl_func, "(mf_model     ) time_lblrtm/calls  : %.2lf s",   ratio_lblrtm_call            );

                cpl_msg_info(cpl_func, "(mf_model     ) orignorm           : %.3e",      fit_results->orignorm        );
                cpl_msg_info(cpl_func, "(mf_model     ) bestnorm           : %.3e",      params->config->internal.chi2);
                cpl_msg_info(cpl_func, "(mf_model     ) npar               : %d",        fit_results->npar            );
                cpl_msg_info(cpl_func, "(mf_model     ) npix               : %d",        fit_results->nfunc           );
                cpl_msg_info(cpl_func, "(mf_model     ) niter              : %d",        fit_results->niter           );
                cpl_msg_info(cpl_func, "(mf_model     ) nfev               : %d",        fit_results->mpfit_calls     );

                results->res = mf_model_res_create(params, results->spec, results->atm_profile_fitted, fit_results);
                if (!results->res) err = CPL_ERROR_NULL_INPUT;

            } else if (fit_results->status >= 0) {

                fit_results->status = -99;
            }

            cpl_msg_info(cpl_func,     "(mf_model     ) status             : %d", fit_results->status);

            /* Cleanup */
            mf_fit_results_delete(fit_results);
        }
    }

    /* Calculate function for telluric absorption correction if a transmission spectrum was fitted */
    if (!err) err = mf_model_calculate_func_telluric_absortion_correction(results->spec, params);

    /* Set model flux to 0 for non-fitted pixels if desired */
    if (!err && params->config->inputs.clean_flux) {

        /* Get pointers to CPL table columns mflux and weight */
        double *mflux  = cpl_table_get_data_double(results->spec, MF_COL_MOD_FLUX);
        double *weight = cpl_table_get_data_double(results->spec, MF_COL_WEIGHT  );

        /* Set mflux to 0 for pixels with weight = 0 */
        int nrow = cpl_table_get_nrow(results->spec);
        for (cpl_size i = 0; i < nrow; i++) {
            if (weight[i] == 0.) mflux[i] = 0.;
        }

        err = cpl_error_get_code();
    }

    /* Show info review */
    if (!err) {

        params->timers.time_end = cpl_test_get_walltime();
        double runtime = params->timers.time_end - params->timers.time_start;

        double time_external = params->timers.time_bin_lnfl + params->timers.time_bin_lblrtm;

        cpl_msg_info(cpl_func, "(mf_model     ) time_bin_lnfl      : %.2lf min",  params->timers.time_bin_lnfl   / 60.);
        cpl_msg_info(cpl_func, "(mf_model     ) time_bin_lblrtm    : %.2lf min",  params->timers.time_bin_lblrtm / 60.);
        cpl_msg_info(cpl_func, "(mf_model     ) time_bin_calls     : %.2lf min",            time_external        / 60.);
        cpl_msg_info(cpl_func, "(mf_model     ) time_telluriccorr  : %.2lf min", (runtime - time_external)       / 60.);
        cpl_msg_info(cpl_func, "(mf_model     ) time_total         : %.2lf min",  runtime                        / 60.);
        //cpl_msg_info(cpl_func, "(mf_model     ) gdas_src           : %s",  results->gdas_src);
    }

    /* Cleanup parameters */
    mf_parameters_delete(params);
    mf_io_oda_delete_tableDB();
    cpl_msg_info(cpl_func,"FINISHED FREEING ODA TABLE");

    /* Check error and return */
    if (!err) {
        return results;
    } else {
        mf_model_results_delete(results);
        return NULL;
    }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Deallocate results structure.
 *
 * @param results            .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 */
/* ---------------------------------------------------------------------------*/
void mf_model_results_delete(
    mf_model_results         *results)
{
  if (results) {

      /* Remove telluriccorr temporary files */
      if (results->tmp_folder) {
        cpl_boolean rm_tmp_dir = CPL_TRUE;
        char* oda_option=getenv("MF_ODA_RM_TMP_DIR");
        if (oda_option!=NULL) {
            if (strncasecmp( oda_option, "false", 5 ) == 0) {
                rm_tmp_dir=CPL_FALSE;
            }
            if (rm_tmp_dir) {
                mf_io_rm_rf(results->tmp_folder, 10);
            }
        }
        cpl_free(results->tmp_folder);
      }

      if (results->gdas_before         ) cpl_table_delete( results->gdas_before         );
      if (results->gdas_after          ) cpl_table_delete( results->gdas_after          );
      if (results->gdas_interpolate    ) cpl_table_delete( results->gdas_interpolate    );
      if (results->atm_profile_combined) cpl_table_delete( results->atm_profile_combined);
      if (results->atm_profile_fitted  ) cpl_table_delete( results->atm_profile_fitted  );
      if (results->res                 ) cpl_table_delete( results->res                 );
      if (results->spec                ) cpl_table_delete( results->spec                );
      if (results->kernel              ) cpl_matrix_delete(results->kernel              );

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


/** @cond PRIVATE */

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param params             .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static mf_model_results * mf_model_results_create(
    const mf_parameters      *params)
{
  mf_model_results *results     = cpl_malloc(sizeof(mf_model_results));

  results->tmp_folder           = cpl_sprintf("%s", params->config->internal.tmp_folder);

  results->gdas_before          = NULL;
  results->gdas_after           = NULL;
  results->gdas_interpolate     = NULL;
  results->atm_profile_standard = NULL;
  results->atm_profile_combined = NULL;
  results->atm_profile_fitted   = NULL;
  results->res                  = NULL;
  results->spec                 = NULL;
  results->kernel               = NULL;

  results->n_range              = -1;
  results->spec_out             = NULL;
  results->range_status         = NULL;
  results->gdas_src             = "NONE";

  return results;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static mf_fit_results * mf_model_batch(
    mf_io_lblrtm_config      *lblrtm_config,
    mf_io_lnfl_config        *lnfl_config,
    mf_parameters            *params,
    mf_model_results         *results)
{
  /*!
   * \callgraph
   *
   * This routine runs the fitting procedure several times. For each fit
   * the fit flags are changed.
   *
   * \b INPUT:
   * \param spec     CPL table with observed spectrum
   * \param prof     CPL table with atmospheric profiles
   * \param params   mf_parameters parameter structure
   *
   * \b OUTPUT:
   * \param spec     CPL table with observed and best-fit model spectrum
   * \param result   CMPFIT structure for fit results
   *
   * \b ERRORS:
   * - No data
   * - Insufficient memory
   * - see mf_fit
   */

  cpl_error_code status = CPL_ERROR_NONE;
  cpl_array *fit_molec_0, *fit_molec_1, *fit_cont_0, *fit_cont_1;
  cpl_array *fit_wlc_0, *fit_wlc_1, *fit_res_0, *fit_res_1;
  int i = 0, runcheck[5] = {1, 1, 1, 1, 1}, niter = 0;
  double orignorm = 0.;

  /* Get size of fit flag arrays */
  cpl_size nmolec = cpl_table_get_nrow(params->molectab);
  if (nmolec < 1) {
      cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                            "No data: cpl_table params->molectab");
      return NULL;
  }

  cpl_size nrange = cpl_table_get_nrow(params->rangetab);
  if (nrange < 1) {
      cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                            "No data: cpl_table params->rangetab");
      return NULL;
  }

  /* This ncont is the dimensions of the tables/array that hold the 
   * continuum order values - not the continuum orders themselves!
   * It looked like this was truncating the fitting order artificially,
   * but it doesn't as we get 'Input data do not match: Invalid object 
   * structure: cpl_array *fit_cont (size != nrange)' (13) at 
   * mf_model_parameters_set_fit_flags:mf_model.c:890 
   */
  cpl_size ncont = nrange;
  /*cpl_size ncont  = 1 + params->config->fitting.fit_continuum.n;*/
  
  if (params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE) ncont++;

  cpl_size nchip = cpl_table_get_nrow(params->chiptab);
  if (nchip < 1) {
      cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                            "No data: cpl_table params->chiptab");
      return NULL;
  }

  cpl_size nwlc = nchip;
  cpl_size nres = 3;

  /* Initialize fit flag arrays */

  fit_molec_0 = cpl_array_new(nmolec, CPL_TYPE_INT);
  fit_molec_1 = cpl_array_new(nmolec, CPL_TYPE_INT);
  cpl_array_fill_window_int(fit_molec_0, 0, nmolec, 0);
  cpl_array_fill_window_int(fit_molec_1, 0, nmolec, 0);

  fit_cont_0  = cpl_array_new(ncont,  CPL_TYPE_INT);
  fit_cont_1  = cpl_array_new(ncont,  CPL_TYPE_INT);
  cpl_array_fill_window_int(fit_cont_0,  0, ncont,  0);
  cpl_array_fill_window_int(fit_cont_1,  0, ncont,  0);

  fit_wlc_0   = cpl_array_new(nwlc,   CPL_TYPE_INT);
  fit_wlc_1   = cpl_array_new(nwlc,   CPL_TYPE_INT);
  cpl_array_fill_window_int(fit_wlc_0,   0, nwlc,   0);
  cpl_array_fill_window_int(fit_wlc_1,   0, nwlc,   0);

  fit_res_0   = cpl_array_new(nres,   CPL_TYPE_INT);
  fit_res_1   = cpl_array_new(nres,   CPL_TYPE_INT);
  cpl_array_fill_window_int(fit_res_0,   0, nres,   0);
  cpl_array_fill_window_int(fit_res_1,   0, nres,   0);


  /* Set fit_molec_1 fit flags */
  for (i = 0; i < nmolec; i++) {
       cpl_msg_info(cpl_func,"fit_molec_1 for molec %d = %f",i,cpl_table_get(params->molectab, MF_COL_FIT_MOLECULES, i, NULL));
     cpl_array_set_int(fit_molec_1, i, cpl_table_get(params->molectab, MF_COL_FIT_MOLECULES, i, NULL));
  }

  /* Set fit_cont_1 fit flags */
  for (i = 0; i < nrange; i++) {
      cpl_msg_info(cpl_func,"fit_cont_1 for range %d = %f",i,cpl_table_get(params->rangetab, MF_COL_FIT_RANGE, i, NULL));
      cpl_array_set_int(fit_cont_1, i, cpl_table_get(params->rangetab, MF_COL_FIT_RANGE, i, NULL));
  }

  if (params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE) {
      /* Fit telescope background if relevant */
      cpl_array_set_int(fit_cont_1, nrange, params->config->fitting.fit_telescope_background.fit);
  }

  /* Set fit_wlc_1 fit flag */
  for (i = 0; i < nchip; i++) {
      cpl_msg_info(cpl_func,"fit_wlc_1 for chip %d = %f",i,cpl_table_get(params->chiptab, MF_COL_FIT_CHIP, i, NULL));
      cpl_array_set_int(fit_wlc_1, i, cpl_table_get(params->chiptab, MF_COL_FIT_CHIP, i, NULL));
  }
  
 
  /* Set fit_res_1 fit flags */

  cpl_array_set_int(fit_res_1, 0, params->config->fitting.fit_res_box.fit);
  cpl_array_set_int(fit_res_1, 1, params->config->fitting.fit_gauss.fit  );
  cpl_array_set_int(fit_res_1, 2, params->config->fitting.fit_lorentz.fit);

  /* Find "non-zero" fit flag arrays with all elements = 0 */
  if (cpl_array_get_max(fit_molec_1) == 0) runcheck[0] = 0;
  if (cpl_array_get_max(fit_cont_1 ) == 0) runcheck[1] = 0;
  if (cpl_array_get_max(fit_wlc_1  ) == 0) runcheck[2] = 0;
  if (cpl_array_get_max(fit_res_1  ) == 0) runcheck[3] = 0;


  /*** Run CMPFIT ***/
  cpl_size       run_execution = 0;
  mf_fit_results *fit_results  = mf_fit_results_create();

  /* 1st run: get the continuum first approximation */
  if (runcheck[1] == 0 || status != CPL_ERROR_NONE) {
      cpl_msg_info(cpl_func, "(mf_model     ) RUN skipped - 1st: Get the continuum first approximation!");
  } else {

      status = mf_model_parameters_set_fit_flags(params, fit_molec_0, fit_cont_1, fit_wlc_0, fit_res_0);
      if (!status) {

          /* Fitting */
          cpl_msg_info(cpl_func, "(mf_model     ) RUN INI - 1st: Get the continuum first approximation");

          double ts = cpl_test_get_walltime();
          status = mf_fit(++run_execution, lblrtm_config, lnfl_config,
                          results->spec, results->kernel, results->atm_profile_combined, &(results->atm_profile_fitted),
                          params, fit_results, results->spec_out, results->range_status);
          double te = cpl_test_get_walltime();
          double runtime = te - ts;
          cpl_msg_info(cpl_func, "(mf_model     ) RUN END - 1st: mf_fit(...) --> Time execution = %.2f s. (status=%d)", runtime, status);

          niter    += fit_results->niter;
          orignorm  = fit_results->orignorm;
      }
  }

  /* 2nd run: get the wavelength calibration and resolution first approximation */
  if ((runcheck[2] == 0 && runcheck[3] == 0) || status != CPL_ERROR_NONE) {
      cpl_msg_info(cpl_func, "(mf_model     ) RUN skipped - 2nd: Get the wavelength calibration and resolution first approximation!");
  } else {

      status = mf_model_parameters_set_fit_flags(params, fit_molec_0, fit_cont_0, fit_wlc_1, fit_res_1);
      if (!status) {

          /* Fitting */
          cpl_msg_info(cpl_func, "(mf_model     ) RUN INI - 2nd: Get the wavelength calibration and resolution first approximation");

          double ts = cpl_test_get_walltime();
          status = mf_fit(++run_execution, lblrtm_config, lnfl_config,
                          results->spec, results->kernel, results->atm_profile_combined, &(results->atm_profile_fitted),
                          params, fit_results, results->spec_out, results->range_status);
          double te = cpl_test_get_walltime();
          double runtime = te - ts;
          cpl_msg_info(cpl_func, "(mf_model     ) RUN END - 2nd: mf_fit(...) --> Time execution = %.2f s. (status=%d)", runtime, status);

          if (runcheck[1] == 0) {
              orignorm = fit_results->orignorm;
          }

          niter += fit_results->niter;
      }
  }

  /* 3rd run: get the continuum second approximation */
  if (runcheck[1] == 0 || status != CPL_ERROR_NONE) {
      cpl_msg_info(cpl_func, "(mf_model     ) RUN skipped - 3rd: Get the continuum second approximation!");
  } else {

      status = mf_model_parameters_set_fit_flags(params, fit_molec_0, fit_cont_1, fit_wlc_0, fit_res_0);
      if (!status) {

          /* Fitting */
          cpl_msg_info(cpl_func, "(mf_model     ) RUN INI - 3rd: Get the continuum second approximation");

          double ts = cpl_test_get_walltime();
          status = mf_fit(++run_execution, lblrtm_config, lnfl_config,
                          results->spec, results->kernel, results->atm_profile_combined, &(results->atm_profile_fitted),
                          params, fit_results, results->spec_out, results->range_status);
          double te = cpl_test_get_walltime();
          double runtime = te - ts;
          cpl_msg_info(cpl_func, "(mf_model     ) RUN END - 3rd: mf_fit(...) --> Time execution = %.2f s. (status=%d)", runtime, status);

          if (runcheck[1] == 0) {
              orignorm = fit_results->orignorm;
          }

          niter += fit_results->niter;
      }
  }

  /* 4th run: get the column densities first approximation */
  if (runcheck[0] == 0 || status != CPL_ERROR_NONE) {
      cpl_msg_info(cpl_func, "(mf_model     ) RUN skipped - 4th: Get the column densities first approximation!");
  } else {

      status = mf_model_parameters_set_fit_flags(params, fit_molec_1, fit_cont_0, fit_wlc_0, fit_res_0);
      if (!status) {

          /* Fitting */
          cpl_msg_info(cpl_func, "(mf_model     ) RUN INI - 4th: Get the column densities first approximation");

          double ts = cpl_test_get_walltime();
          status = mf_fit(++run_execution, lblrtm_config, lnfl_config,
                          results->spec, results->kernel, results->atm_profile_combined, &(results->atm_profile_fitted),
                          params, fit_results, results->spec_out, results->range_status);
          double te = cpl_test_get_walltime();
          double runtime = te - ts;
          cpl_msg_info(cpl_func, "(mf_model     ) RUN END - 4th: mf_fit(...) --> Time execution = %.2f s. (status=%d)", runtime, status);

          if (runcheck[1] == 0 && runcheck[2] == 0 && runcheck[3] == 0) {
              orignorm = fit_results->orignorm;
          }

          niter += fit_results->niter;
      }
  }

  /* 5th run: redo continuum, wavelength, resolution */
  if ((runcheck[1] == 0 && runcheck[2] == 0 && runcheck[3] == 0) || status != CPL_ERROR_NONE) {
      cpl_msg_info(cpl_func, "(mf_model     ) RUN skipped - 5th: Redo continuum, wavelength, resolution!");
  } else {

      status = mf_model_parameters_set_fit_flags(params, fit_molec_0, fit_cont_1, fit_wlc_1, fit_res_1);
      if (!status) {

          /* Fitting */
          cpl_msg_info(cpl_func, "(mf_model     ) RUN INI - 5th: Redo continuum, wavelength, resolution");

          double ts = cpl_test_get_walltime();
          status = mf_fit(++run_execution, lblrtm_config, lnfl_config,
                          results->spec, results->kernel, results->atm_profile_combined, &(results->atm_profile_fitted),
                          params, fit_results, results->spec_out, results->range_status);
          double te = cpl_test_get_walltime();
          double runtime = te - ts;
          cpl_msg_info(cpl_func, "(mf_model     ) RUN END - 5th: mf_fit(...) --> Time execution = %.2f s. (status=%d)", runtime, status);

          niter += fit_results->niter;
      }
  }

  /* 6th run: final adjustment */
  if (runcheck[0] == 0 || status != CPL_ERROR_NONE) {
      cpl_msg_info(cpl_func, "(mf_model     ) RUN skipped - 6th: Final adjustment!");
  } else {

      status = mf_model_parameters_set_fit_flags(params, fit_molec_1, fit_cont_1, fit_wlc_1, fit_res_1);
      if (!status) {

          /* Fitting */
          cpl_msg_info(cpl_func, "(mf_model     ) RUN INI - 6th: Final adjustment");

          double ts = cpl_test_get_walltime();
          status = mf_fit(++run_execution, lblrtm_config, lnfl_config,
                          results->spec, results->kernel, results->atm_profile_combined, &(results->atm_profile_fitted),
                          params, fit_results, results->spec_out, results->range_status);
          double te = cpl_test_get_walltime();
          double runtime = te - ts;
          cpl_msg_info(cpl_func, "(mf_model     ) RUN END - 6th: mf_fit(...) --> Time execution = %.2f s. (status=%d)", runtime, status);

          niter += fit_results->niter;
      }
  }

  /* Allocate memory for results structure if all fit flags = 0 */
  if (   runcheck[0] == 0
      && runcheck[1] == 0
      && runcheck[2] == 0
      && runcheck[3] == 0) {
      fit_results->status = -19;
      status = CPL_ERROR_ILLEGAL_OUTPUT;
  }

  /* Correct run-dependent parameters of results structure */
  fit_results->niter    = niter;
  fit_results->orignorm = orignorm;

  /* Cleanup */
  cpl_array_delete(fit_molec_0);
  cpl_array_delete(fit_molec_1);
  cpl_array_delete(fit_cont_0);
  cpl_array_delete(fit_cont_1);
  cpl_array_delete(fit_wlc_0);
  cpl_array_delete(fit_wlc_1);
  cpl_array_delete(fit_res_0);
  cpl_array_delete(fit_res_1);

  /* Check errors */
  if (status != CPL_ERROR_NONE) {
      cpl_msg_error(cpl_func, "(mf_model     ) RUN END with status[cpl_error_code] = %d", status);
      cpl_error_set_message(cpl_func, status, "Fit error!");
  }

  return fit_results;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_model_parameters_set_fit_flags(
    mf_parameters            *params,
    const cpl_array          *fit_molec,
    const cpl_array          *fit_cont,
    const cpl_array          *fit_wlc,
    const cpl_array          *fit_res)
{
  /*!
   * Modifies the fit flags (0 or 1) in the telluriccorr driver parameter
   * structure mf_parameters. The input CPL arrays have to consist of nmolec
   * (fit_molec), nrange (fit_cont), nchip (fit_wlc), and 3 (fit_res)
   * elements, respectively. In the case of a radiance spectrum a fit flag
   * for the telescope background has to be added to fit_cont.
   *
   * \b INPUT:
   * \param params     mf_parameters parameter structure
   * \param fit_molec  CPL array for flags for fitting of molecules
   * \param fit_cont   CPL array for flags for continuum fit
   * \param fit_wlc    CPL array for flags for wavelength fit
   * \param fit_res    CPL array for flags for resolution fit
   *
   * \b OUTPUT:
   * \param params     mf_parameters parameter structure with modified fit flags
   *
   * \b ERRORS:
   * - Invalid object value(s)
   * - Invalid object structure
   */

  /* Get number of molecules and compare with format of parameter table */
  int nmolec = params->config->internal.molecules.n_molec;
  int nrange = params->config->internal.n_range;
  int nchip  = params->config->internal.nchip;

  if (nmolec != cpl_table_get_nrow(params->molectab)) {
      return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                   "Invalid object value(s): nmolec of mf_parameters *params (number of molecules)");
  }

  /* Check number of fit_molec array elements */
  if (nmolec != cpl_array_get_size(fit_molec)) {
      return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                   "Invalid object structure: cpl_array *fit_molec (size != nmolec)");
  }

  /* Get number of ranges and compare with format of parameter table */
  if (nrange != cpl_table_get_nrow(params->rangetab)) {
      return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                   "Invalid object value(s): nrange of mf_parameters *params (number of ranges)");
  }

  /* Check number of fit_cont array elements */
  if (params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE) {
      if ((nrange + 1) != cpl_array_get_size(fit_cont)) {
          return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                       "Invalid object structure: cpl_array *fit_cont (size != nrange + 1)");
      }
  } else {
      if (nrange != cpl_array_get_size(fit_cont)) {
          return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                       "Invalid object structure: cpl_array *fit_cont (size != nrange)");
      }
  }

  /* Get number of ranges and compare with format of parameter table */
  if (cpl_table_get_nrow(params->chiptab) != nchip) {
      return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                   "Invalid object value(s): nchip of mf_parameters *params (number of chips)");
  }

  /* Check number of fit_wlc array elements */
  if (cpl_array_get_size(fit_wlc) != nchip) {
      return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                   "Invalid object structure: cpl_array *fit_wlc (size != nchip)");
  }

  /* Check number of fit_res array elements */
  if (cpl_array_get_size(fit_res) != 3) {
      return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                   "Invalid object structure: cpl_array *fit_res (size != 3)");
  }

  /* Check validity of fit_molec fit flags */
  for (cpl_size i = 0; i < cpl_array_get_size(fit_molec); i++) {
      int flag = cpl_array_get(fit_molec, i, NULL);
      if (flag < 0 || flag > 1) {
          return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                       "Invalid object value(s): cpl_array *fit_molec (fit flag != 0 or 1)");
      }
  }

  /* Check validity of fit_cont fit flags */
  for (cpl_size i = 0; i < cpl_array_get_size(fit_cont); i++) {
      int flag = cpl_array_get(fit_cont, i, NULL);
      if (flag < 0 || flag > 1) {
          return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                       "Invalid object value(s): cpl_array *fit_cont (fit flag != 0 or 1)");
      }
  }

  /* Check validity of fit_wlc fit flags */
  for (cpl_size i = 0; i < cpl_array_get_size(fit_wlc); i++) {
      int flag = cpl_array_get(fit_wlc, i, NULL);
      if (flag < 0 || flag > 1) {
          return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                       "Invalid object value(s): cpl_array *fit_wlc (fit flag != 0 or 1)");
      }
  }

  /* Check validity of fit_res fit flags */
  for (cpl_size i = 0; i < cpl_array_get_size(fit_res); i++) {
      int flag = cpl_array_get(fit_res, i, NULL);
      if (flag < 0 || flag > 1) {
          return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                       "Invalid object value(s): cpl_array *fit_res (fit flag != 0 or 1)");
      }
  }

  /* Set fit_molec fit flags */
  for (cpl_size i = 0; i < nmolec; i++) {
      cpl_table_set_int(params->molectab, MF_COL_FIT_MOLECULES, i, cpl_array_get(fit_molec, i, NULL));
  }

  /* Set fit_cont fit flags */
  for (cpl_size range_index = 0; range_index < nrange; range_index++) {
      cpl_table_set_int(params->rangetab, MF_COL_FIT_RANGE, range_index, cpl_array_get(fit_cont, range_index, NULL));
  }

  /* Fit telescope background if relevant */
  if (params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE) {
      /* Fit telescope background if relevant */
      params->config->fitting.fit_telescope_background.fit = cpl_array_get(fit_cont, nrange, NULL);
  }

  /* Set fit_wlc fit flag */
  for (cpl_size chip_index = 0; chip_index < nchip; chip_index++) {
      cpl_table_set_int(params->chiptab, MF_COL_FIT_CHIP, chip_index, cpl_array_get(fit_wlc, chip_index, NULL));
  }


  /*** Set fit_res fit flags ***/

  params->config->fitting.fit_res_box.fit = cpl_array_get(fit_res, 0, NULL);
  params->config->fitting.fit_gauss.fit   = cpl_array_get(fit_res, 1, NULL);
  params->config->fitting.fit_lorentz.fit = cpl_array_get(fit_res, 2, NULL);

  return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_model_calculate_func_telluric_absortion_correction(
    cpl_table                *spec,
    const mf_parameters      *params)
{
  /*!
   * Adds the column MF_COL_MOD_TRANS, which includes a transmission curve that can
   * be used for telluric absorption correction. It differs from the
   * best-fit model by the neglection of the continuum fit, i.e. the
   * unabsorbed continuum is characterised by a value of 1. No column is
   * created if the fitted spectrum is an telluric emission spectrum.
   *
   * \b INPUT:
   * \param spec    CPL table with observed and modelled spectrum
   * \param params  mf_parameters parameter structure
   *
   * \b OUTPUT:
   * \param spec    table with additional column for model transmission
   *                curve (only for telluric features in absorption)
   *
   * \b ERRORS:
   * - none
   */

  /* Return if plain sky spectrum was fitted */
  if (params->config->inputs.transmission != MF_PARAMETERS_TRANSMISSION_FALSE) {

      /* Create new column and set values to model flux */
      cpl_table_duplicate_column(spec, MF_COL_OUT_TELLURIC_CORR, spec, MF_COL_MOD_FLUX);

      /* Divide model spectrum by continuum scaling function */
      cpl_table_divide_columns(  spec, MF_COL_OUT_TELLURIC_CORR, MF_COL_MOD_SCALE);
  }

  return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_table * mf_model_res_create(
    mf_parameters            *params,
    cpl_table                *spec,
    cpl_table                *atm_profile,
    const mf_fit_results     *fit_results)
{
    /*!
     * Writes a summary of the CMPFIT results into the outtable.
     *
     * \b INPUT:
     * \param params   mf_parameters parameter structure
     * \param spec     CPL table with observed and best-fit spectrum
     * \param result   CMPFIT structure for fit results
     * \param fittime  fit run time in minutes
     *
     * \b OUTPUT:
     * \param params   mf_parameters parameter structure with \f${\chi^2}\f$
     *
     *
     * \b ERRORS:
     * - File opening failed
     * - Insufficient data points
     * - Invalid object value(s)
     */


    /* Get number of parameters */
    int ncont  = 1 + params->config->fitting.fit_continuum.n;
    int nwlc   = 1 + params->config->fitting.fit_wavelenght.n;
    int nrange =     params->config->internal.n_range;
    int nchip  =     params->config->internal.nchip;
    int nmolec =     params->config->internal.molecules.n_molec;

    int npar   = nmolec + (nchip * nwlc) + (nrange * ncont) + 3;

    /* For the table */
    cpl_size nRows =   1               /* H2O MM                         */
                     + nmolec          /* RELATIVE MOLECULAR GAS COLUMNS */
                     + nmolec          /* MOLECULAR GAS COLUMNS IN PPMV  */
                     + 1               /* TELESCOPE BACKGROUND           */
                     + nrange * ncont  /* CONTINUUM CORRECTION           */
                     + nchip  * nwlc   /* WAVELENGTH SOLUTION            */
                     + 1               /* FWHM Lorentzian                */
                     + 1               /* FWHM Gaussian                  */
                     + 2               /* FWHM boxcar + in pixels        */
                     + 15;             /* MPFIT                          */


    /* If not the first call, rewrite it. Otherwise, we end up with a memory leak. */
    cpl_table *res = cpl_table_new(nRows);
    cpl_table_new_column(res, MF_COL_PARAMETER,   CPL_TYPE_STRING);
    cpl_table_new_column(res, MF_COL_VALUE,       CPL_TYPE_DOUBLE);
    cpl_table_new_column(res, MF_COL_UNCERTAINTY, CPL_TYPE_DOUBLE);

    /* Write status ID and return in the case of errors */
    cpl_size index = 0;
    if (fit_results->status <= 0) {
        return res;                     /* FIXME ? */
    } else {
        mf_model_fill_parameter(res, index++, MF_MODEL_FIT_STATUS, fit_results->status, -1.);
    }


    /* Get number of data points with valid model fluxes and non-zero weight (= DOF+1) and calculate chi^2 */
    const double *weight  = cpl_table_get_data_double_const(spec, MF_COL_WEIGHT    );
    const int    *mrange  = cpl_table_get_data_int_const(   spec, MF_COL_MOD_RANGE );
    const double *mweight = cpl_table_get_data_double_const(spec, MF_COL_MOD_WEIGHT);
    const double *dev     = cpl_table_get_data_double_const(spec, MF_COL_DEV       );

    cpl_size     nmr      = 0;
    cpl_size     nmw      = 0;
    cpl_size     nw       = 0;
    double       chi2     = 0.;

    cpl_size nrows = cpl_table_get_nrow(spec);
    for (cpl_size l = 0; l < fit_results->nfunc && l < nrows; l++) {

        if (mrange[l]  > 0) nmr++;

        if (mweight[l] > 0) {
            nmw++;
            if (weight[l] > 0) {
                nw++;
                chi2 += dev[l] * dev[l];
            }
        }
    }

    if (nw == 1) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                              "Insufficient data points: cpl_table *spec (only 1 data point with weight > 0)");
        cpl_table_delete(res);
        return NULL;
    }

    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_N_PAR,        fit_results->npar,         -1.);
    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_DATA_POINTS,  fit_results->nfunc,        -1.);
    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_POS_WEIGHTS,  nw,                        -1.);
    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_PIX_FRAC,     (double)nmw / (double)nmr, -1.);
    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_N_ITER,       fit_results->niter,        -1.);
    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_FUNC_VAL,     fit_results->mpfit_calls,  -1.);
    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_LBLRTM_CALLS, fit_results->lblrtm_calls, -1.);
    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_INIT_CHI2,    fit_results->orignorm,     -1.);

    /* Write chi^2 and reduced chi2 */
    params->config->internal.chi2 = chi2;
    mf_model_fill_parameter(     res, index++, MF_MODEL_FIT_BEST_CHI2,    chi2,                      -1.);
    if (nw > 1) {
        double chi2red = chi2 / (nw - 1);
        mf_model_fill_parameter( res, index++, MF_MODEL_FIT_REDUCED_CHI2, chi2red,                   -1.);
        mf_model_fill_parameter( res, index++, MF_MODEL_FIT_RMS_REL_ERR,  sqrt(chi2red),             -1.);
    } else {
        mf_model_fill_parameter( res, index++, MF_MODEL_FIT_REDUCED_CHI2, -1.,                       -1.);
        mf_model_fill_parameter( res, index++, MF_MODEL_FIT_RMS_REL_ERR,  -1.,                       -1.);
    }

    /* Compute RMS relative to weighted mean and write it to file */
    const double *flux = cpl_table_get_data_double_const(spec, MF_COL_IN_FLUX);
    double       wsum  = 0.;
    double       wfsum = 0.;
    double       w2sum = 0.;
    for (cpl_size l = 0; l < fit_results->nfunc && l < nrows; l++) {
        if (mweight[l] > 0) {
            wsum  += weight[l];
            wfsum += weight[l] * flux[l];
            w2sum += weight[l] * weight[l];
        }
    }


    /* Check RMS errors and compute final value */
    double rms = 0.;
    if (wsum  == 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Invalid object value(s): cpl_table *spec (all weights = 0)");
        cpl_table_delete(res);
        return NULL;
    } else if (wfsum == 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Invalid object value(s): cpl_table *spec (all fluxes  = 0)");
        cpl_table_delete(res);
        return NULL;
    } else {
        rms = sqrt(chi2 / w2sum) * (wsum / wfsum);
    }

    if (wsum == 0 || wfsum == 0) mf_model_fill_parameter(res, index++, MF_MODEL_FIT_RMS_REL_MEAN, -1., -1.);
    else                         mf_model_fill_parameter(res, index++, MF_MODEL_FIT_RMS_REL_MEAN, rms, -1.);

    /* Boxcar: Spectral resolution */
    double slit_width  = params->config->instrumental.slit_width.value;
    double pixel_scale = params->config->instrumental.pixel_scale.value;

    cpl_boolean fit_res_box = params->config->fitting.fit_res_box.fit;
    double      rel_res_box = params->config->fitting.fit_res_box.const_val;
    if (fit_res_box) {
        mf_model_fill_parameter(               res, index++, MF_MODEL_FIT_BOX_FWHM,     rel_res_box,  fit_results->xerror[npar-3]                );
        mf_model_fill_parameter(               res, index++, MF_MODEL_FIT_BOX_FWHM_PIX, rel_res_box,  fit_results->xerror[npar-3] * slit_width / pixel_scale);
    } else {
        mf_model_fill_parameter(               res, index++, MF_MODEL_FIT_BOX_FWHM,     rel_res_box,  -1.                                   );
        mf_model_fill_parameter(               res, index++, MF_MODEL_FIT_BOX_FWHM_PIX, rel_res_box,  -1.                                   );
    }


    /* Gauss */
    cpl_boolean fit_gauss = params->config->fitting.fit_gauss.fit;
    double      res_gauss = params->config->fitting.fit_gauss.const_val;
    if (fit_gauss) mf_model_fill_parameter(    res, index++, MF_MODEL_FIT_GAUSS_FWHM,   res_gauss,   fit_results->xerror[npar-2]                );
    else           mf_model_fill_parameter(    res, index++, MF_MODEL_FIT_GAUSS_FWHM,   res_gauss,   -1.                                   );


    /* Lorentz */
    cpl_boolean fit_lorentz = params->config->fitting.fit_lorentz.fit;;
    double      res_lorentz = params->config->fitting.fit_lorentz.const_val;
    if (fit_lorentz) mf_model_fill_parameter(  res, index++, MF_MODEL_FIT_LORENTZFWHM,  res_lorentz, fit_results->xerror[npar-1]                );
    else             mf_model_fill_parameter(  res, index++, MF_MODEL_FIT_LORENTZFWHM,  res_lorentz, -1.                                   );


    /* Wavelength */
    int ipar = nmolec;
    for (cpl_size chip_index = 0; chip_index < nchip; chip_index++) {

        cpl_size chip = chip_index + 1;

        int fit_wavelength = cpl_table_get(params->chiptab, MF_COL_FIT_CHIP, chip_index, NULL);

        cpl_array *wlc_coef = cpl_array_duplicate(cpl_table_get_array(params->chiptab, MF_COL_WLC_COEF, chip_index));

        for (cpl_size coef = 0; coef < nwlc; coef++) {

            double value = cpl_array_get_double(wlc_coef, coef, NULL);

            char *pname = cpl_sprintf("%s %lld, coef %lld", MF_COL_CHIP, chip, coef);
            if (fit_wavelength == 1) mf_model_fill_parameter(res, index++, pname, value, fit_results->xerror[ipar++]);
            else                     mf_model_fill_parameter(res, index++, pname, value, -1.);
            cpl_free(pname);
        }

        cpl_array_delete(wlc_coef);
    }

    /* Continuum correction */
    for (cpl_size j = 0; j < nrange; j++) {

        int       chip          = cpl_table_get(params->rangetab, MF_COL_CHIP,      j, NULL);
        int       fit_continuum = cpl_table_get(params->rangetab, MF_COL_FIT_RANGE, j, NULL);

        cpl_array *cont_coef    = cpl_array_duplicate(cpl_table_get_array(params->rangetab, MF_COL_CONT_COEF, j));
        for (cpl_size i = 0; i < ncont; i++) {

            double value = cpl_array_get_double(cont_coef, i, NULL);

            char *pname = cpl_sprintf("Range %lld, chip %d, coef %lld", j + 1, chip, i);
            if (fit_continuum == 1) mf_model_fill_parameter(res, index++, pname, value, fit_results->xerror[ipar++]);
            else                    mf_model_fill_parameter(res, index++, pname, value, -1);
            cpl_free(pname);
        }
        cpl_array_delete(cont_coef);
    }

    /* Molecular spectrum in emission? (for radiance spectrum only) */
    if (params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE) {

        npar++;
        cpl_boolean fit_emission = params->config->fitting.fit_telescope_background.fit;
        double      telback      = params->config->fitting.fit_telescope_background.const_val;
        if (fit_emission) mf_model_fill_parameter( res, index++, MF_MODEL_FIT_TELBACK, telback, fit_results->xerror[npar-1]);
        else              mf_model_fill_parameter( res, index++, MF_MODEL_FIT_TELBACK, telback, -1.);

    } else {

        mf_model_fill_parameter(res, index++, MF_MODEL_FIT_TELBACK, -1., -1.);
    }


    /* Get molecules and properties */
    char   **mol   = cpl_table_get_data_string(params->molectab, MF_COL_LIST_MOLECULES);
    int    *fitmol = cpl_table_get_data_int(   params->molectab, MF_COL_FIT_MOLECULES );
    double *relcol = cpl_table_get_data_double(params->molectab, MF_COL_REL_COL   );

    /* Relative molecular gas columns */
    for (cpl_size i = 0; i < nmolec; i++) {
        char *pname = cpl_sprintf(MF_MODEL_FIT_REL_MOL_COL_"%s", mol[i]);
        if (fitmol[i] == 1) mf_model_fill_parameter(res, index++, pname, relcol[i], fit_results->xerror[i]);
        else                mf_model_fill_parameter(res, index++, pname, relcol[i], -1.);
        cpl_free(pname);
    }

    /* Calcule atm */
    double col_H2O = 0.;
    mf_model_calculate_atm_colums(atm_profile, &col_H2O, params);

    /* Molecular gas columns in ppmv */
    double   *ppmv          = cpl_table_get_data_double( params->molectab, MF_COL_PPMV);
    cpl_size nrows_molectab = cpl_table_get_nrow(        params->molectab);
    int      fit_H2O        = 0;
    double   reldel_H2O     = 0.;
    for (cpl_size i = 0; i < nmolec && i < nrows_molectab && ppmv; i++) {

        double reldel = fit_results->xerror[i] / relcol[i];

        if (strcmp(mol[i], MF_MOLECULES_H2O) == 0) {
            fit_H2O    = fitmol[i];
            reldel_H2O = reldel;
        }

        char *pname = cpl_sprintf(MF_MODEL_FIT_PPMV_"%s", mol[i]);
        if (fitmol[i] == 1) mf_model_fill_parameter(res, index++, pname, ppmv[i], reldel * ppmv[i]);
        else                mf_model_fill_parameter(res, index++, pname, ppmv[i], -1.);
        cpl_free(pname);
    }

    /* H2O column in mm */
    if (fit_H2O == 1) mf_model_fill_parameter(res, index++, MF_MODEL_FIT_H20_COL_MM, col_H2O, reldel_H2O * col_H2O);
    else              mf_model_fill_parameter(res, index++, MF_MODEL_FIT_H20_COL_MM, col_H2O, -1.);

    /* Return */
    return res;
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief .
 *
 * @param .                  .
 * @param                    .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @description .
 *
 * @note .
 *
 */
/* ---------------------------------------------------------------------------*/
static void mf_model_fill_parameter(
    cpl_table                *outtable,
    cpl_size                 index,
    const char               *parameter,
    const double             val,
    const double             uncertainty)
{
  cpl_table_set_string(outtable, MF_COL_PARAMETER,   index, parameter  );
  cpl_table_set_double(outtable, MF_COL_VALUE,       index, val        );
  cpl_table_set_double(outtable, MF_COL_UNCERTAINTY, index, uncertainty);
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Calculates water vapor column in mm from profile in ppmv.
 *
 * @param h2ocol             output: water vapor column in mm
 * @param params             input/output: mf_parameters parameter structure (supplemented by ppmv of provided molecules)
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - File opening failed.
 *                           - Invalid object structure
 *
 * @description The input profiles are interpreted as functions.
 *              The starting height is the altitude of the observing site.
 *              The profile information is read from the output profiles file.
 *              Moreover, the ppmv of the entire atmospheric column is provided for
 *              all molecules provided by the mf_parameters parameter structure.
 *
 * @note This information is written into a special ppmv column of the table of molecules of mf_parameters.
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_model_calculate_atm_colums(
    cpl_table                *atm_profile,
    double                   *h2ocol,
    mf_parameters            *params)
{
    /* Check input atm profile */
    if(!atm_profile) {
        return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                     "NULL input atm profile");
    }

    /* Default water vapor column */
    *h2ocol = 0.;

    /* Create column for ppmv of molecules in mf_parameters structure if not present */
    if (cpl_table_has_column(params->molectab, MF_COL_PPMV) != 1) {
        cpl_table_new_column(params->molectab, MF_COL_PPMV, CPL_TYPE_DOUBLE);
    }

    /* Set initial ppmv values to 0 */
    cpl_size nmolec = cpl_table_get_nrow(params->molectab);
    cpl_table_fill_column_window(params->molectab, MF_COL_PPMV, 0, nmolec, 0.);

    /* Get pointer to names of molecules in driver parameter structure */
    char **mol = cpl_table_get_data_string(params->molectab, MF_COL_LIST_MOLECULES);

    /* Get pointer to ppmv of molecules in driver parameter structure */
    double *ppmv = cpl_table_get_data_double(params->molectab, MF_COL_PPMV);

    /* Check existence of selected molecular columns in profile table */
    for (cpl_size j = 0; j < nmolec; j++) {
        if (cpl_table_has_column(atm_profile, mol[j]) != 1) {
            return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                         "Invalid object structure: cpl_table *prof (no %s column)",
                                         mol[j]);
        }
    }

    /* Altitude of observing site [in km] */
    double elevation = params->config->ambient.elevation.value * MF_CONV_M_2_KM;

    /* Number of layers, initialisation of upper layer height, and total column height */
    cpl_size nlayer = cpl_table_get_nrow(atm_profile);
    double   uhgt   = cpl_table_get(atm_profile, MF_COL_ATM_HGT, 0, NULL);

    /* Sum up molecular columns of all layers */
    double nmol0   = 0.;
    for (cpl_size i = 0; i < nlayer - 1; i++) {

        /* Lower and upper limit of layer */
        double lhgt = uhgt;

        /* Skip layers below height of observing site */
        uhgt = cpl_table_get(atm_profile, MF_COL_ATM_HGT, i+1, NULL);
        if (uhgt > elevation) {

          /* Thickness of layer in m */
          double dlayer;
          if (*h2ocol == 0. && uhgt > elevation) dlayer = (uhgt - elevation) * MF_CONV_KM_2_M;
          else                                   dlayer = (uhgt - lhgt   ) * MF_CONV_KM_2_M;

          /* Average pressure and temperature for layer */
          double pressure = (cpl_table_get(atm_profile, MF_COL_ATM_PRE, i, NULL) + cpl_table_get(atm_profile, MF_COL_ATM_PRE, i+1, NULL)) / 2 ;
          double temp     = (cpl_table_get(atm_profile, MF_COL_ATM_TEM, i, NULL) + cpl_table_get(atm_profile, MF_COL_ATM_TEM, i+1, NULL)) / 2;
          if (temp <= 0.) temp = MF_TOL;

          /* Number of mols per unit area for each molecule and PWV */
          for (cpl_size j = 0; j < nmolec; j++) {

              /* Average ppmv of molecule for layer */
              double molcol = (cpl_table_get(atm_profile, mol[j], i, NULL) + cpl_table_get(atm_profile, mol[j], i+1, NULL)) / 2;

              /* Column height [m] of molecule for layer */
              double ch = 1e-6 * molcol * dlayer;

              /* Number of mols per unit area [mol m^-2] for layer */
              double nmol = ch * pressure * MF_CONV_MBAR_2_PA / (MF_R * temp);

              /* Number of mols per unit area for atmosphere */
              ppmv[j] += nmol;

              /* Number of mols per unit area for air */
              if (j == 0) nmol0 += dlayer * pressure * MF_CONV_MBAR_2_PA / (MF_R * temp);

              /* H2O column in mm (PWV) : Mass per unit area [kg m^-2] */
              if (strcmp(mol[j], MF_MOLECULES_H2O) == 0) {
                  *h2ocol += nmol * MF_MOL_MASS_H2O;
              }
          }
        }
    }

    /* Volume mixing ratio for each molecule */
    for (cpl_size j = 0; j < nmolec; j++) {
        if (nmol0 <= 0.) ppmv[j]  = 0.;
        else             ppmv[j] *= 1e6 / nmol0;
    }

    return CPL_ERROR_NONE;
}

/** @endcond */


/**@}*/
