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

#include "cpl_mpfit.h"

#include "mf_lblrtm.h"
#include "mf_convolution.h"
#include "mf_spectrum.h"


#include "mf_fit.h"

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

/* Parameter constraint structure */
typedef struct {
  cpl_boolean           fit;             /* 1 = fixed; 0 = free                               */
  int                   limited[2];      /* 1 = low/upper limit; 0 = no limit                 */
  double                limits[2];       /* lower/upper limit boundary value                  */
  char                  *parname;        /* Name of parameter, or 0 for none                  */
  double                relstep;         /* Relative step size for finite difference          */
} mf_fit_parameter;

/* Structure for holding the fit parameters and their constraints */
typedef struct {
  int                   n;               /* Number of parameters                              */
  double                *p;              /* Parameter vector for fitting                      */
  mf_fit_parameter      *pars;           /* Structure for parameter constraints               */
} mf_fit_parameters;

/* Structure for holding the spectrum data and mf_parameters */
typedef struct {
  cpl_size              run_execution;   /* Number of run execution                           */
  cpl_table             *spec;           /* CPL table with wavelengths, fluxes, and weights   */
  cpl_matrix            *kernel;         /* CPL matrix with kernel to convolve                */
  cpl_table             *atm_prof;       /* CPL table with atmospheric profiles               */
  mf_parameters         *params;         /* Telluric Correction parameter structure           */
  cpl_table             **prof_out;      /* Output CPL table with atmospheric profiles        */
  cpl_array             *reffitpar;      /* Reference parameter fit array                     */
  int                   *mpfit_calls;    /* Number of lblrtm_calls                            */
  int                   *lblrtm_calls;   /* Number of lblrtm_calls                            */
  mf_io_lblrtm_config   *lblrtm_config;  /* Configuration of LBLRTM                           */
  mf_io_lnfl_config     *lnfl_config;    /* Configuration of LNFL                             */
  cpl_table             **spec_out;      /* Array with the cpl_table calculated in each range */
  cpl_error_code        *range_status;   /* Array with the error code status in each range    */
} mf_fit_calcdev_vars;

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/*  */
static cpl_error_code mf_fit_set_mf_fit_parameters(
    mf_fit_parameters        *fitpar,
    const mf_parameters      *params);

/*  */
static cpl_error_code mf_fit_set_mf_parameters(
    mf_parameters            *params,
    const mf_fit_parameters  *fitpar);

/* */
static int mf_fit_get_npar(
    const mf_parameters      *params,
    int                      *nmolec,
    int                      *nchip,
    int                      *nwlc,
    int                      *nrange,
    int                      *ncont,
    int                      *trans);

/*  */
static cpl_error_code mf_fit_parameters_create(
    mf_fit_parameters        *fitpar,
    const int                npar);

/*  */
static cpl_error_code mf_fit_parameters_delete(
    mf_fit_parameters        *fitpar);

/*  */
static void mf_fit_parameter_init(
    int                      *index,
    mf_fit_parameters        *fitpar,
    const double             value,
    const cpl_boolean        fit,
    const int                limited0,
    const int                limited1,
    const double             limits0,
    const double             limits1,
    const double             relstep,
    const char               *parname);

/*  */
static int mf_fit_execute_calcdev(
    int                      m,
    int                      n,
    double                   *p,
    double                   *dy,
    double                   **dvec,
    void                     *generic_vars);

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_fit     .
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
 * @brief Execute fit.
 *
 * @param run_execution      .
 * @param lblrtm_config       .
 * @param spec               .
 * @param kernel             CPL table with observed spectrum.
 * @param atm_profile        CPL table with atmospheric profiles.
 * @param atm_profile_out    .
 * @param params             mf_parameters parameter structure.
 * @param results            .
 * @param spec_out           CPL table with observed and best-fit model spectrum.
 * @param range_status       CMPFIT structure for fit results.
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - No data.
 *                           - Insufficient memory
 *                           - Error in subroutine (see subroutines).
 *
 * @note Handles the fitting routine CMPFIT (see mpfit.c for details on the
 *         fitting algorithm). Needs the observed spectrum plus weights,
 *         atmospheric profiles for pressure, temperature, and relevant molecules,
 *         and fit-related parameters from the mf_parameters driver parameter structure.
 *         The input data are used to compute a model spectrum which is compared
 *         to the observed spectrum by deriving a vector of weighted deviations.
 *         CMPFIT optimises the model spectrum by manipulating a vector of the
 *         fit parameters which is used as input for the model calculation.
 *         Information on the fit quality is written into a special results
 *         structure. The CPL table with the input spectrum is supplemented by the
 *         best-fit model and the weighted deviations taken for the \f${\chi^2}\f$
 *         computation. Moreover, an updated driver file, an ASCII file which
 *         summarises the fit results (including the ppmv of the different
 *         molecules and especially the water-related PWV in mm), and an ASCII
 *         file with the best-fit atmospheric profiles are written.
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_fit(
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
    cpl_error_code           *range_status)
{
  /* Set parameter vector and parameter constraints structure */
  mf_fit_parameters fitpar;
  mf_fit_set_mf_fit_parameters(&fitpar, params);
  if (fitpar.p == NULL || fitpar.pars == NULL) {
      results->status = -99;
      return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                   "No data: mf_fit_parameters fitpar");
  }

  /* Memory allocation for parameter errors */
  if (results->xerror) cpl_free(results->xerror);
  results->xerror = cpl_calloc(fitpar.n, sizeof(double));

  /* Reset number of fitting calls */
  cpl_msg_info(cpl_func, "(mf_fit       ) Fitting function calls and fit parameter changes -> Call par newval reldev");

  /* Start timer */
  double ts = cpl_test_get_walltime();

  /*** MINPACK-1 : Least Square Minimization ***/

  /* Memory reservation */
  mp_par *mpfit_constrains = cpl_calloc(fitpar.n, sizeof(mp_par));

  /* Structure for constrains */
  for (cpl_size i = 0; i < fitpar.n; i++) {
      mpfit_constrains[i].fixed      = (int)(fitpar.pars[i].fit); /* Casting need for cpl_mpfit call */
      mpfit_constrains[i].limited[0] = fitpar.pars[i].limited[0];
      mpfit_constrains[i].limited[1] = fitpar.pars[i].limited[1];
      mpfit_constrains[i].limits[0]  = fitpar.pars[i].limits[0];
      mpfit_constrains[i].limits[1]  = fitpar.pars[i].limits[1];
      mpfit_constrains[i].parname    = fitpar.pars[i].parname;
      mpfit_constrains[i].relstep    = fitpar.pars[i].relstep;
  }

  /* Structure for config */
  mp_config mpfit_config = { params->config->fitting.ftol,
                             params->config->fitting.xtol,
                             0, 0, 0, 0,
                             MF_FIT_ITER_MAX,
                             MF_FIT_ITER_MAX * fitpar.n,
                             0, 0, 0, 0};

  /* Structure for parameters for used in the user mf_fit_calcdev function */
  mf_fit_calcdev_vars vars;
  vars.run_execution = run_execution;
  vars.spec          = spec;
  vars.kernel        = kernel;
  vars.atm_prof      = atm_profile;
  vars.params        = params;
  vars.prof_out      = atm_profile_out;
  vars.reffitpar     = cpl_array_wrap_double(fitpar.p, fitpar.n);
  vars.mpfit_calls   = &(results->mpfit_calls);
  vars.lblrtm_calls  = &(results->lblrtm_calls);
  vars.lblrtm_config = lblrtm_config ;
  vars.lnfl_config   = lnfl_config;
  vars.spec_out      = spec_out;
  vars.range_status  = range_status;

  /* Memory reservation for results */
  int    nPoints = cpl_table_get_nrow(spec);
  double *resid  = cpl_calloc(nPoints, sizeof(double));
  mp_result mpfit_results = {0., 0., 0, 0, 0, 0, 0, 0, 0, resid, results->xerror, NULL, "" };

  /* Call fitting function for nPoints data points and n parameters and save results */
  int status = cpl_mpfit(mf_fit_execute_calcdev, nPoints, fitpar.n, fitpar.p, mpfit_constrains, &mpfit_config, (void *)&vars, &mpfit_results);

  /* Save results */
  results->orignorm    = mpfit_results.orignorm;
  results->niter       = mpfit_results.niter;
  results->mpfit_calls = mpfit_results.nfev;
  results->status      = mpfit_results.status;
  results->npar        = mpfit_results.npar;
  results->nfunc       = mpfit_results.nfunc;

  /* Cleanup */
  cpl_free(mpfit_constrains);
  if(resid) cpl_free(resid);

  /* Stop timer: FIT Time execution */
  double te        = cpl_test_get_walltime();
  double runtime   = te - ts;
  params->timers.time_fit += runtime;
  cpl_msg_info(cpl_func, "(mf_fit       ) RUN MF_FIT --> Time execution = %.2f s. (TOTAL_TIME_FIT = %.2f min., TIME_EXE = %.2f min., nfev = %3d)",
               runtime, params->timers.time_fit / 60., (te - params->timers.time_start) / 60., mpfit_results.nfev);

  /* Check status */
  cpl_error_code err = CPL_ERROR_NONE;
  if (status < 0) {
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                   "Error in subroutine: cpl_mpfit(...)");
  } else {

      /* Write best fit parameters to the telluriccorr parameter structure */
      mf_fit_set_mf_parameters(params, &fitpar);

      /* RUN LBLRTM */
      cpl_boolean last_call = CPL_TRUE;

      if (!err) {
          err = mf_lblrtm(     run_execution,
                               params, last_call, spec_out, range_status,
                               vars.reffitpar, results->mpfit_calls, vars.lblrtm_calls, lblrtm_config , lnfl_config,
                               atm_profile, atm_profile_out,
                               vars.reffitpar);
          if (err != CPL_ERROR_NONE) cpl_msg_error(cpl_func, "(mf_fit       ) Call to mf_lblrtm(...) fails, cpl_error_code=%d",        err);
      }


      if (!err) {
          err = mf_convolution(params, CPL_TRUE, last_call, spec_out, range_status,
                               vars.reffitpar, results->mpfit_calls,
                               spec,
                               kernel);
          if (err != CPL_ERROR_NONE) cpl_msg_error(cpl_func, "(mf_fit       ) Call to mf_convolution(...) fails, cpl_error_code=%d",        err);
      }
  }

  /* Unwrap reference parameter vector */
  cpl_array_unwrap(vars.reffitpar);

  /* Cleanup */
  mf_fit_parameters_delete(&fitpar);

  return err;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Run fit of telluriccorr in a temporary folder.
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
mf_fit_results * mf_fit_results_create(void)
{
  mf_fit_results *results = cpl_malloc(sizeof(mf_fit_results));

  results->mpfit_calls  = 0;
  results->lblrtm_calls = 0;
  results->xerror       = NULL;

  return results;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Delete fit results structure
 *
 * @param results            .
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
void mf_fit_results_delete(
    mf_fit_results           *results)
{
  if (results) {

      if (results->xerror) cpl_free(results->xerror);

      cpl_free(results);
  }
}

/** @cond PRIVATE */

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
static cpl_error_code mf_fit_set_mf_fit_parameters(
    mf_fit_parameters        *fitpar,
    const mf_parameters      *params)
{
  /*!
   * Provides a vector of parameters which are variables of the CMPFIT
   * fitting process. Moreover, constraints for these parameters are
   * delivered by the CMPFIT structure mf_par. Both objects are returned by
   * a container structure.
   *
   * \b INPUT:
   * \param params  mf_parameters parameter structure
   *
   * \b OUTPUT:
   * \param fitpar  structure containing the fit parameters
   *
   * \b ERRORS:
   * - see mf_fit_allocmempar
   */


  /* Allocate memory for mf_fit_parameters structure */
  int nmolec, nchip, nwlc, nrange, ncont, trans;
  int npar = mf_fit_get_npar(params, &nmolec, &nchip, &nwlc, &nrange, &ncont, &trans);
  if (mf_fit_parameters_create(fitpar, npar) != CPL_ERROR_NONE) {
      return cpl_error_get_code();
  }

  /* Write fit parameter values into mf_fit_parameters structure */
  int index = -1;

  /* Amount of molecular gas */
  for (cpl_size i = 0; i < nmolec; i++) {

      cpl_boolean fit   = 1 - cpl_table_get(params->molectab, MF_COL_FIT_MOLECULES, i, NULL);
      double      value =     cpl_table_get(params->molectab, MF_COL_REL_COL,   i, NULL);

      mf_fit_parameter_init(&index, fitpar, value, fit, 1, 1, 1e-5, 1e2, 0.01, cpl_table_get_string(params->molectab, MF_COL_LIST_MOLECULES, i));
  }

  /* Wavelength solution */
  cpl_boolean fitwlc1   = params->config->internal.wlc_lin;
  for (cpl_size chip_index = 0; chip_index < nchip; chip_index++) {

      cpl_size chip = chip_index + 1;

      cpl_boolean fitwlc    = cpl_table_get(                          params->chiptab, MF_COL_FIT_CHIP, chip_index, NULL);
      cpl_array   *wlc_coef = cpl_array_duplicate(cpl_table_get_array(params->chiptab, MF_COL_WLC_COEF, chip_index)     );

      for (cpl_size i = 0; i < nwlc; i++) {

          cpl_boolean fit    = (i == 1 && fitwlc) ? !fitwlc1 : !fitwlc;
          double      value  = cpl_array_get_double(wlc_coef, i, NULL);
          char        *pname = cpl_sprintf("c%lld_wlc_coef_%lld", chip, i);

          if (i == 1) mf_fit_parameter_init(&index, fitpar, value, fit, 0, 0, 0., 0., 1e-4, pname);
          else        mf_fit_parameter_init(&index, fitpar, value, fit, 0, 0, 0., 0., 0.1,  pname);

          cpl_free(pname);
      }
      cpl_array_delete(wlc_coef);
  }

  /* Continuum */
  for (cpl_size j = 0; j < nrange; j++) {

      cpl_boolean fitcont    = cpl_table_get(      params->rangetab,                     MF_COL_FIT_RANGE, j, NULL);
      cpl_array   *cont_coef = cpl_array_duplicate(cpl_table_get_array(params->rangetab, MF_COL_CONT_COEF, j));
      int range_poly_order   = params->config->fitting.cont_poly_order[j];

      for (cpl_size i = 0; i < ncont; i++) {

          cpl_boolean fit    = !fitcont;
          double      value  = cpl_array_get_double(cont_coef, i, NULL);
          char        *pname = cpl_sprintf("c%lld_cont_coef_%lld", j + 1, i);
	  cpl_boolean fit_coeff=CPL_TRUE;
	  if (i>range_poly_order) fit_coeff=CPL_FALSE;;
	  
	  fit = !(fitcont && fit_coeff);
	  

          if (i == 0) mf_fit_parameter_init(&index, fitpar, value, fit, 1, 0, MF_TOL * value, 0., 0.01, pname);
          else        mf_fit_parameter_init(&index, fitpar, value, fit, 0, 0, MF_TOL * value, 0., 0.01, pname);

          cpl_free(pname);
      }
      cpl_array_delete(cont_coef);
  }

  /* Boxcar */
  cpl_boolean fit_res_box = !(params->config->fitting.fit_res_box.fit);
  double      rel_res_box = params->config->fitting.fit_res_box.const_val;
  mf_fit_parameter_init(&index, fitpar, rel_res_box,   fit_res_box,   1, 1, 0.,   2.,       0.1,    MF_PARAMETERS_RES_BOX );

  /* Gaussian */
  cpl_boolean fit_gauss = !(params->config->fitting.fit_gauss.fit);
  double      res_gauss = params->config->fitting.fit_gauss.const_val;
  mf_fit_parameter_init(&index, fitpar, res_gauss,    fit_gauss,    1, 1, 0., 100.,       0.1,    MF_PARAMETERS_RES_GAUSS   );

  /* Lorentzian */
  cpl_boolean fit_lorentz = !(params->config->fitting.fit_lorentz.fit);
  double      res_lorentz     =  params->config->fitting.fit_lorentz.const_val;;
  mf_fit_parameter_init(&index, fitpar, res_lorentz,  fit_lorentz,  1, 1, 0., 100.,       0.1,    MF_PARAMETERS_RES_LORENTZ );

  /* Telescope background scaling */
  if (trans == MF_PARAMETERS_TRANSMISSION_FALSE) {
      cpl_boolean fit_telback = !(params->config->fitting.fit_telescope_background.fit);
      double      telback     = params->config->fitting.fit_telescope_background.const_val;
      mf_fit_parameter_init(&index, fitpar, telback, fit_telback, 1, 1, 0., 1 - MF_TOL, 0.01,   MF_PARAMETERS_TELESCOPE_BACK_CONST     );
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
static cpl_error_code mf_fit_set_mf_parameters(
    mf_parameters            *params,
    const mf_fit_parameters  *fitpar)
{
  /*!
   * Writes fit parameter values as provided by an mf_fit_parameters structure into
   * the mf_parameters parameter structure.
   *
   * \b INPUT:
   * \param params  mf_parameters parameter structure
   * \param fitpar  structure containing the fit parameters
   *
   * \b OUTPUT:
   * \param params  updated mf_parameters parameter structure
   *
   * \b ERRORS:
   * - Invalid object structure
   */

  /* Get parameters */
  int nmolec, nchip, nwlc, nrange, ncont, trans;
  int npar = mf_fit_get_npar(params, &nmolec, &nchip, &nwlc, &nrange, &ncont, &trans);
  if (fitpar->n != npar) {
      return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                   "Invalid object structure: fitpar->n != npar derived from mf_parameters *params");
  }

  /* Write fit parameter values into mf_parameters structure */
  int idx = -1;

  /* Amount of molecular gas */
  for (cpl_size i = 0; i < nmolec; i++) {
      cpl_table_set_double(params->molectab, MF_COL_REL_COL, i, fitpar->p[++idx]);
  }

  /* Wavelength solution */
  cpl_array *wlc_coef = cpl_array_new(nwlc, CPL_TYPE_DOUBLE);
  for (cpl_size chip_index = 0; chip_index < nchip; chip_index++) {
      for (cpl_size coef = 0; coef < nwlc; coef++) {
          cpl_array_set_double(wlc_coef, coef, fitpar->p[++idx]);
      }
      cpl_table_set_array(params->chiptab, MF_COL_WLC_COEF, chip_index, wlc_coef);
  }
  cpl_array_delete(wlc_coef);

  /* Continuum */
  cpl_array *cont_coef = cpl_array_new(ncont, CPL_TYPE_DOUBLE);
  for (cpl_size j = 0; j < nrange; j++) {
      for (cpl_size i = 0; i < ncont; i++) {
          cpl_array_set_double(cont_coef, i, fitpar->p[++idx]);
      }
      cpl_table_set_array(params->rangetab, MF_COL_CONT_COEF, j, cont_coef);
  }
  cpl_array_delete(cont_coef);

  /* Resolution */
  params->config->fitting.fit_res_box.const_val = fitpar->p[++idx];
  params->config->fitting.fit_gauss.const_val   = fitpar->p[++idx];
  params->config->fitting.fit_lorentz.const_val = fitpar->p[++idx];

  /* Telescope background scaling */
  if (trans == MF_PARAMETERS_TRANSMISSION_FALSE) {
      params->config->fitting.fit_telescope_background.const_val = fitpar->p[++idx];
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
static int mf_fit_get_npar(
    const mf_parameters      *params,
    int                      *nmolec,
    int                      *nchip,
    int                      *nwlc,
    int                      *nrange,
    int                      *ncont,
    int                      *trans)
{
  *nmolec =     params->config->internal.molecules.n_molec;
  *nrange =     params->config->internal.n_range;
  *nchip  =     params->config->internal.nchip;
  *trans  =     params->config->inputs.transmission;
  *ncont  = 1 + params->config->fitting.fit_continuum.n;
  *nwlc   = 1 + params->config->fitting.fit_wavelenght.n;

  int npar = *nmolec + ((*nchip) * (*nwlc)) + ((*nrange) * (*ncont)) + 3;

  if (*trans == MF_PARAMETERS_TRANSMISSION_FALSE) return (npar + 1);
  else                                            return  npar;
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
static cpl_error_code mf_fit_parameters_create(
    mf_fit_parameters        *fitpar,
    const int                npar)
{
  /*!
   * Allocates memory for an mf_fit_parameters structure.
   *
   * \b INPUT:
   * \param npar    number of fit parameters
   *
   * \b OUTPUT:
   * \param fitpar  mf_fit_parameters structure containing npar fit parameters
   *
   * \b ERRORS:
   * - Invalid input parameter(s)
   * - Insufficient memory
   */

  fitpar->n    = npar;
  fitpar->p    = NULL;
  fitpar->pars = NULL;

  /* Check number of parameters */
  if (fitpar->n < 0) {
      fitpar->n = 0;
      return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                   "Invalid input parameter(s): npar < 0");
  } else if (fitpar->n == 0) {
      return CPL_ERROR_NONE;
  }

  /* Allocate memory for parameter vector */
  fitpar->p    = cpl_calloc(fitpar->n, sizeof(double));

  /* Allocate memory for parameter constraints structure */
  fitpar->pars = cpl_calloc(fitpar->n, sizeof(mf_fit_parameter));

  /* Allocate memory for parameter names in parameter constraints structure */
  int nchar = 1024;
  cpl_boolean fl_mem = CPL_TRUE;
  for (cpl_size it = 0; it < 2; it++) {

      for (cpl_size i = 0; i < fitpar->n; i++) {

          fitpar->pars[i].parname = cpl_calloc(nchar, sizeof(char));

          if (it == 0 && fitpar->pars[i].parname == NULL) {
              nchar  = 0;
              fl_mem = CPL_FALSE;
              continue;
          }
      }

      if (       it == 0 && fl_mem == CPL_TRUE) {
          break;
      } else if (it == 1 && fl_mem == CPL_FALSE) {

          cpl_free(fitpar->p);
          cpl_free(fitpar->pars);

          fitpar->n    = 0;
          fitpar->p    = NULL;
          fitpar->pars = NULL;

          return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                       "Insufficient memory: mf_fit_parameters *fitpar");
      }
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
static cpl_error_code mf_fit_parameters_delete(
    mf_fit_parameters        *fitpar)
{
  /*!
   * Frees memory occupied by an mf_fit_parameters structure.
   *
   * \b INPUT:
   * \param fitpar  mf_fit_parameters structure containing npar fit parameters
   *
   * \b OUTPUT:
   * \param fitpar  mf_fit_parameters structure without allocated memory
   *
   * \b ERRORS:
   * - none
   */

  if (fitpar) {

    /* Free memory occupied by parameter vector */
    if (fitpar->p != NULL) {
        cpl_free(fitpar->p);
        fitpar->p = NULL;
    }

    /* Free memory occupied by parameter constraints structure */
    if (fitpar->pars != NULL) {
        for (cpl_size i = 0; i < fitpar->n; i++) {
            cpl_free(fitpar->pars[i].parname);
        }
        cpl_free(fitpar->pars);
        fitpar->pars = NULL;
    }

    /* Set number of parameters to 0 */
    fitpar->n = 0;
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
static void mf_fit_parameter_init(
    int                      *index,
    mf_fit_parameters        *fitpar,
    const double             value,
    const cpl_boolean        fit,
    const int                limited0,
    const int                limited1,
    const double             limits0,
    const double             limits1,
    const double             relstep,
    const char               *parname)
{
  (*index)++;

  fitpar->p[   *index]            = value;
  fitpar->pars[*index].fit        = fit;
  fitpar->pars[*index].limited[0] = limited0;
  fitpar->pars[*index].limited[1] = limited1;
  fitpar->pars[*index].limits[0]  = limits0;
  fitpar->pars[*index].limits[1]  = limits1;
  fitpar->pars[*index].relstep    = relstep;
  fitpar->pars[*index].parname    = cpl_sprintf("%s", parname);
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
static int mf_fit_execute_calcdev(
    int                      m,
    int                      n,
    double                   *p,
    double                   *dy,
    double                   **dvec,
    void                     *generic_vars)
{
  /*!
   * \callgraph
   *
   * User function for CMPFIT (The syntax of the function is predefined).
   * Returns weighted deviations between the model calculated by mf_lblrtm(...) and the observed spectrum.
   *
   * \b INPUT:
   * \param m     number of data points
   * \param n     number of parameters
   * \param p     array of fit parameters
   * \param dvec  derivatives (not used)
   * \param vars  private data -> observed spectrum and driver file
   *                              parameters
   *
   * \b OUTPUT:
   * \param dy    array of residuals ([model - obs. spectrum] * weight)
   *
   * \b ERRORS:
   * - Invalid input parameter(s)
   */

  /* Check dvec variable is NULL */
  cpl_ensure(!dvec, CPL_ERROR_ILLEGAL_INPUT, CPL_ERROR_ILLEGAL_INPUT);


  /* Unpack observed spectral data and input parameters */
  mf_fit_calcdev_vars *vars             = (mf_fit_calcdev_vars *)generic_vars;
  cpl_size            run_execution     = vars->run_execution;
  cpl_table           *spec             = vars->spec;
  cpl_matrix          *kernel           = vars->kernel;
  cpl_table           *atm_profile      = vars->atm_prof;
  mf_parameters       *params           = vars->params;
  cpl_table           **atm_profile_out = vars->prof_out;
  cpl_array           *reffitpar        = vars->reffitpar;
  int                 *mpfit_calls      = vars->mpfit_calls;
  int                 *lblrtm_calls     = vars->lblrtm_calls;
  mf_io_lblrtm_config *lblrtm_config     = vars->lblrtm_config;
  mf_io_lnfl_config   *lnfl_config      = vars->lnfl_config;
  cpl_table           **spec_out        = vars->spec_out;
  cpl_error_code      *range_status     = vars->range_status;

  /* Update number of fitting function call */
  (*mpfit_calls)++;

  /* Check for parameters with "nan" as value */
  cpl_size index = 0; /* This parameter is user after */
  for (index = 0; index < n; index++) {
      if (isnan(p[index]) != 0) {
          cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                "Invalid input parameter(s): double *p has nan values");
          break;
      }
  }

  /* Put fit parameters in CPL array */
  cpl_array *fitpar = cpl_array_wrap_double(p, n);

  /* Take reference parameter values in the case of errors */
  if (isnan(p[index - 1]) != 0) {
      for (cpl_size i = 0; i < n; i++) {
          double value = cpl_array_get(reffitpar, i, NULL);
          cpl_array_set(fitpar, i, value);
      }
  }

  /* RUN LBLRTM : Model calculation. It isn't the last call, it's a intermediate call in the fitting */
  cpl_boolean last_call = CPL_FALSE;

  cpl_error_code err = mf_lblrtm(run_execution,
                                 params, last_call, spec_out, range_status,
                                 reffitpar, *mpfit_calls, lblrtm_calls, lblrtm_config , lnfl_config,
                                 atm_profile, atm_profile_out,
                                 fitpar);

  if (!err) {
      err = mf_convolution(params, CPL_TRUE, last_call, spec_out, range_status,
                           fitpar,    *mpfit_calls,
                           spec,
                           kernel);
  }

  cpl_array_unwrap(fitpar);

  /* Fill array of residuals */
  if (!err) {
      for (cpl_size i = 0; i < m; i++) {
          dy[i] = cpl_table_get(spec, MF_COL_DEV, i, NULL);
      }
  } else {
      for (cpl_size i = 0; i < m; i++) {
          dy[i] = 0.;
      }
  }

  return err;
}

/** @endcond */


/**@}*/
