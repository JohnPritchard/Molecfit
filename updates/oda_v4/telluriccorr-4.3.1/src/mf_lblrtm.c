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

#include "mf_lblrtm.h"

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

/* Execute the lblrtm binary in all of the wavenumbers in a concrete range */
static cpl_error_code mf_lblrtm_range_execution(
    mf_io_lblrtm_config      *lblrtm_config ,
    mf_io_lnfl_config        *lnfl_config,
    mf_parameters            *params,
    const int                ismolcol,
    const cpl_size           run_execution,
    int                      *lblrtm_calls,
    cpl_table                *atm_profile,
    cpl_table                **spec_out,
    cpl_vector               *mol_abuns,
    cpl_error_code           *range_status,
    mf_io_oda_parameters     *oda_parameters);

/*  */
static cpl_error_code mf_lblrtm_molecular_correct_atm(
    cpl_table                *atm_profile,
    cpl_table                **atm_profile_out,
    const mf_parameters      *params,
    const cpl_array          *fitpar);

/* Make output spectrum : Combined and convolved */
static cpl_table * mf_lblrtm_range_combined_and_convolved(
    const char               *lblrtm_out_filename,
    const double             pixel_res,
    const double             wn1,
    const double             wn2,
    const char               *w_dir_range, const int range, cpl_vector* mol_abuns, cpl_boolean USE_HYBRID,
    mf_io_lnfl_config        *lnfl_config);

/*  */
static cpl_error_code mf_lblrtm_create_wavelength_grid(
    cpl_table                *outspec,
    const double             limlam[2],
    const double             resol);

/* Spectrum linear interpolation */
static cpl_error_code mf_lblrtm_linear_interpolate_spectrum(
    double                   *lambda,
    double                   *flux,
    const cpl_size           lim_min,
    const cpl_size           lim_max);

/* Break up a single range into a LBLRTM valid set of wavenumber subranges */
static cpl_vector* wavenumber_subranges(double nu1,double nu2);

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_lblrtm       .
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
mf_io_lblrtm_config * mf_lblrtm_config_create(void)
{
    mf_io_lblrtm_config *config = cpl_malloc(sizeof(mf_io_lblrtm_config));

    config->icntnm    = MF_LBLRTM_ICNTNM_INIT;
    config->iaersl    = MF_LBLRTM_IAERSL_INIT;
    config->mpts      = MF_LBLRTM_MPTS_INIT;
    config->npts      = MF_LBLRTM_NPTS_INIT;
    config->v[0]      = MF_LBLRTM_V1_INIT;
    config->v[1]      = MF_LBLRTM_V2_INIT;
    config->sample    = MF_LBLRTM_SAMPLE_INIT;
    config->alfal0    = MF_LBLRTM_ALFAL0_INIT;
    config->avmass    = MF_LBLRTM_AVMASS_INIT;
    config->dptmin    = MF_LBLRTM_DPTMIN_INIT;
    config->dptfac    = MF_LBLRTM_DPTFAC_INIT;
    config->tbound    = MF_LBLRTM_TBOUND_INIT;
    config->sremis[0] = MF_LBLRTM_SREMIS1_INIT;
    config->sremis[1] = MF_LBLRTM_SREMIS2_INIT;
    config->sremis[2] = MF_LBLRTM_SREMIS3_INIT;
    config->srrefl[0] = MF_LBLRTM_SRREFL1_INIT;
    config->srrefl[1] = MF_LBLRTM_SRREFL2_INIT;
    config->srrefl[2] = MF_LBLRTM_SRREFL3_INIT;
    config->model     = MF_LBLRTM_MODEL_INIT;
    config->itype     = MF_LBLRTM_ITYPE_INIT;
    config->nozero    = MF_LBLRTM_NOZERO_INIT;
    config->noprnt    = MF_LBLRTM_NOPRNT_INIT;
    config->ipunch    = MF_LBLRTM_IPUNCH_INIT;
    config->re        = MF_LBLRTM_RE_INIT;
    config->hspace    = MF_LBLRTM_HSPACE_INIT;
    config->ref_lat   = MF_PARAMETERS_LATITUDE_VALUE_INIT; //MF_LBLRTM_REF_LAT_INIT;
    config->h[0]      = MF_PARAMETERS_ELEVATION_VALUE_INIT/1000.0; // MF_LBLRTM_H1_INIT;
    config->h[1]      = MF_LBLRTM_H2_INIT;
    config->range     = MF_LBLRTM_RANGE_INIT;
    config->beta      = MF_LBLRTM_BETA_INIT;
    config->len       = MF_LBLRTM_LEN_INIT;
    config->hobs      = MF_LBLRTM_HOBS_INIT;
    config->avtrat    = MF_LBLRTM_AVTRAT_INIT;
    config->tdiff[0]  = MF_LBLRTM_TDIFF1_INIT;
    config->tdiff[1]  = MF_LBLRTM_TDIFF2_INIT;
    config->altd[0]   = MF_LBLRTM_ALTD1_INIT;
    config->altd[1]   = MF_LBLRTM_ALTD2_INIT;
    config->delv      = MF_LBLRTM_DELV_INIT;

    return config;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Deallocate the user parameter structure.
 *
 * @param config             .
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
void mf_lblrtm_config_delete(
    mf_io_lblrtm_config      *config)
{
  if (config) {

      cpl_free(config);
  }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Check the user parameter structure. For check the user modifications.
 *
 * @param config             .
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
cpl_error_code mf_lblrtm_config_check(
    mf_io_lblrtm_config      *config)
{
  cpl_ensure(config, CPL_ERROR_NULL_INPUT, CPL_ERROR_NULL_INPUT);

  // TODO: Add check for all variables and values
  cpl_error_code err = CPL_ERROR_NONE;

  if (!err && (   config->icntnm < MF_LBLRTM_ICNTNM_MIN
               || config->icntnm > MF_LBLRTM_ICNTNM_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_ICNTNM);
  }

  if (!err && (   config->iaersl < MF_LBLRTM_IAERSL_MIN
               || config->iaersl > MF_LBLRTM_IAERSL_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_IAERSL);
  }

  if (!err && (   config->sample < MF_LBLRTM_SAMPLE_MIN
               || config->sample > MF_LBLRTM_SAMPLE_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_SAMPLE);
  }

  if (!err && (   config->model < MF_LBLRTM_MODEL_MIN
               || config->model > MF_LBLRTM_MODEL_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_MODEL);
  }

  if (!err && (   config->itype < MF_LBLRTM_ITYPE_MIN
               || config->itype > MF_LBLRTM_ITYPE_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_ITYPE);
  }

  if (!err && (   config->nozero < MF_LBLRTM_NOZERO_MIN
               || config->nozero > MF_LBLRTM_NOZERO_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_NOZERO);
  }

  if (!err && (   config->noprnt < MF_LBLRTM_NOPRNT_MIN
               || config->noprnt > MF_LBLRTM_NOPRNT_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_NOPRNT);
  }

  if (!err && (   config->ipunch < MF_LBLRTM_IPUNCH_MIN
               || config->ipunch > MF_LBLRTM_IPUNCH_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_IPUNCH);
  }

  if (!err && (   config->ref_lat < MF_PARAMETERS_LATITUDE_VALUE_MIN // MF_LBLRTM_REF_LAT_MIN
               || config->ref_lat > MF_PARAMETERS_LATITUDE_VALUE_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_PARAMETERS_LATITUDE_VALUE);
  }

  if (!err && (   config->len < MF_LBLRTM_LEN_MIN
               || config->len > MF_LBLRTM_LEN_MAX) ){
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LBLRTM_LEN);
  }


  return err;
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief Execute LBLRTM
 *
 * @param run_execution      .
 * @param params             .
 * @param last_call          .
 * @param spec_out           .
 * @param range_status       .
 * @param reffitpar          .
 * @param mpfit_calls        .
 * @param lblrtm_calls       .
 * @param lblrtm_config       .
 * @param atm_profile        .
 * @param atm_profile_out    .
 * @param fitpar             .
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
cpl_error_code mf_lblrtm(
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
    cpl_array                *fitpar)
{
    /*** First Part of the call to LBLRTM ***/

    /* Check whether a relative molecular column density was changed (only for mf_model(...) not mf_calctrans...(...) */
    int ismolcol = 1;
    if (mpfit_calls != 1 || !last_call) {

        /* Check number of array values */
         cpl_size n    = cpl_array_get_size(fitpar);
         cpl_size nref = cpl_array_get_size(reffitpar);
         if (n == 0 && n == nref) {
             ismolcol = 1;
             return CPL_ERROR_NONE;
         } else if (n != nref) {
             ismolcol = -1;
             return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                          "Invalid object structure: cpl_array *fitpar != cpl_array *reffitpar (size)");
         }

         /* Get number of molecules in driver parameter structure */
         int nmolec = params->config->internal.molecules.n_molec;

         /* Check modification of molecular column density parameter and copy fitpar to reffitpar */
         if (mpfit_calls == 1) cpl_msg_debug(cpl_func, "(mf_lblrtm    ) FIT ITERATION: %4d ... Getting first approximation",                mpfit_calls        );
         else                  cpl_msg_debug(cpl_func, "(mf_lblrtm    ) FIT ITERATION: %4d ... Checking tolerance [MF_TOLERANCE =  %5.2e]", mpfit_calls, MF_TOL);

         cpl_boolean all_in_tolerance = CPL_TRUE;

         ismolcol = 0;
         for (cpl_size i = 0; i < n; i++) {

             double par    = cpl_array_get(fitpar,    i, NULL);
             double refpar = cpl_array_get(reffitpar, i, NULL);

             double reldev = fabs(par - refpar);

             if (     refpar == 0. && par == 0.) reldev  = 0.;
             else if (refpar == 0.             ) reldev /= par;
             else                                reldev /= refpar;

             if (mpfit_calls != 1 && reldev > MF_TOL) {
                 cpl_msg_debug(cpl_func, "(mf_lblrtm    ) [ELEM = %3lld, VAL = %14.5f, REF = %14.5f] => |VAL-REF|/VAL = %5.2e > %5.2e (NOT FITTING)",
                               i, par, refpar, reldev, MF_TOL);
                 if (all_in_tolerance) all_in_tolerance = CPL_FALSE;
             }

             if (i < nmolec && reldev > MF_TOL) ismolcol = 1;

             /* Check if reffitpar have enough elements */
             if (cpl_array_get_size(reffitpar) > i) {
                 cpl_array_set(reffitpar, i, par);
             }
         }

         if (mpfit_calls != 1 && all_in_tolerance) cpl_msg_info(cpl_func, "(mf_lblrtm    ) ALL FITTING !");
    }

    /* Run LBLRTM in any case if first or last function call */
    if (mpfit_calls == 1 || last_call) ismolcol = 1;

    /* Modify profiles of atmospheric molecules if necessary */
    cpl_table *atm_profile_local = NULL;
    if (ismolcol == 1) {
        atm_profile_local = cpl_table_duplicate(atm_profile);
        mf_lblrtm_molecular_correct_atm(atm_profile_local, atm_profile_out, params, fitpar);
    }

    /* Execute the LBLRTM calls */
    /*
    if (lnfl_config->use_ODA) {
        cpl_vector *mol_abuns = mf_io_molecule_abundancies(params,fitpar);
    }
    */
    cpl_vector *mol_abuns = mf_io_molecule_abundancies(params,fitpar);
    cpl_error_code err = mf_lblrtm_range_execution(lblrtm_config , lnfl_config, params,
                                                   ismolcol,
                                                   run_execution, lblrtm_calls,
                                                   atm_profile_local,
                                                   spec_out, mol_abuns, range_status,NULL);

    cpl_vector_delete(mol_abuns);
    /* Cleanup */
    if (atm_profile_local) cpl_table_delete(atm_profile_local);

    return err;
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
static cpl_error_code mf_lblrtm_molecular_correct_atm(
    cpl_table                *atm_profile,
    cpl_table                **atm_profile_out,
    const mf_parameters      *params,
    const cpl_array          *fitpar)
{
    /*!
     * \callgraph
     *
     * Multiplies profiles of molecules listed in the mf_parameters structure by
     * correction factors given by the CPL array of fit parameters.
     * The resulting CPL table of atmospheric profiles is written into an
     * ASCII output file (parameter MF_PARAMETERS_OUTPUT_NAME in mf_parameters structure +
     * MF_FILE_FIT_ATM_ASCII) which is required as input file for the LBLRTM code.
     *
     * \b INPUT:
     * \param modprof  CPL table with atmospheric profiles
     * \param params   mf_parameters parameter structure
     * \param fitpar   CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param modprof  atmospheric profiles with corrected molecular
     *                 frequencies
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    /* Initialise the return error code */
    cpl_error_code err = CPL_ERROR_NONE;

    /* Get number of columns with molecular profiles */
    cpl_size nmolcol = cpl_table_get_ncol(atm_profile) - 3;

    /* Get pointer to column names */
    cpl_array *colarray = cpl_table_get_column_names(atm_profile);
    char      **colname = cpl_array_get_data_string(colarray);

    /* Get number of molecules in driver parameter structure */
    int nmolec = params->config->internal.molecules.n_molec;

    /* Get pointer to names of molecules in driver parameter structure */
    char **mol = cpl_table_get_data_string(params->molectab, MF_COL_LIST_MOLECULES);

    /* Get pointer to CPL array with fit parameters */
    const double *par = cpl_array_get_data_double_const(fitpar);

    /* Modify column densities of molecules */
    for (cpl_size i = 0; i < nmolec; i++) {
        cpl_boolean ext = CPL_FALSE;
        for (cpl_size j = 0; j < nmolcol; j++) {
           if (strcmp(colname[j + 3], mol[i]) == 0) {
                ext = CPL_TRUE;
                cpl_table_multiply_scalar(atm_profile, colname[j + 3], par[i]);
            }
        }

        /* Molecule of the driver file is not present in list of CPL table column names */
        if (!ext) err = cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                              "Invalid object structure: cpl_table *modprof (column %s not found)",
                                              mol[i]);
    }

    /* Create output table. If not the first call, rewrite it. Otherwise, we end up with a memory leak. */
    if (*atm_profile_out) cpl_table_delete(*atm_profile_out);
    *atm_profile_out = cpl_table_duplicate(atm_profile);

    /* Hack for the first value of the pressure */
    cpl_table_set_double(*atm_profile_out, MF_COL_ATM_PRE, 0, cpl_table_get_double(atm_profile, MF_COL_ATM_PRE, 0, NULL) + 0.1);

    /* Cleanup */
    cpl_array_delete(colarray);

    return err;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Execute the LBLRTM binary in all of the wavenumbers in a concrete range.
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
static cpl_error_code mf_lblrtm_range_execution(
    mf_io_lblrtm_config      *lblrtm_config ,
    mf_io_lnfl_config        *lnfl_config,
    mf_parameters            *params,
    const int                ismolcol,
    const cpl_size           run_execution,
    int                      *lblrtm_calls,
    cpl_table                *atm_profile,
    cpl_table                **spec_out,
    cpl_vector               *mol_abuns,
    cpl_error_code           *range_status,
    mf_io_oda_parameters     *oda_parameters)
{

    /* Local parameters */
    cpl_error_code err               = CPL_ERROR_NONE;
    cpl_boolean    no_single_spec    = !params->config->internal.single_spectrum;
    int            nrange            =  params->config->internal.n_range;

    char           **run_w_dir_range = cpl_calloc(nrange, sizeof(char *));
    double         *wn_end           = cpl_calloc(nrange, sizeof(double));
    double         *pixel_res        = cpl_calloc(nrange, sizeof(double));
    cpl_boolean    *execute_range    = cpl_calloc(nrange, sizeof(cpl_boolean));

    /* MNB-HACK Define a flag to use the ODA Tables method from an env var*/
    /* MNB-HACK Define a flag to output the std LBLRTM  from an env var*/
    cpl_boolean USE_ODATABLE =mf_io_use_odatable(lnfl_config->use_ODA);
    cpl_boolean USE_STDLBLRTM=mf_io_use_stdlblrtm(lnfl_config->use_ODA);
    cpl_boolean WRITE_OUTPUTS=CPL_FALSE;
    if (oda_parameters!=NULL || USE_STDLBLRTM ) WRITE_OUTPUTS=CPL_TRUE;

    if (params->config->internal.single_spectrum) {
        cpl_msg_info(cpl_func, "(mf_lblrtm    ) Compute single spectrum with mf_lblrtm(...)");
    }

    /* Create folders */
    for (cpl_size range = 0; range < nrange; range++) {

        /* Compute LBLRTM spectrum. Get range (number of LBLRTM spectra) */
        if (oda_parameters==NULL) {
            run_w_dir_range[range] = cpl_sprintf("%s/run_%lld_%s_%lld_lblrtm_call_%lld",
                                             params->config->internal.tmp_folder,
                                             run_execution,
                                             MF_AER_WDIR_LBLRTM_RANGE_PATH,
                                             range + 1,
                                             *lblrtm_calls + range + 1);
        } else {
            run_w_dir_range[range] = cpl_sprintf("%s/range_%lld",
                                             oda_parameters->lblrtm_wdir,range + 1);
        }

        /* Get wavenumber range in [cm-1] */
        wn_end[   range] = cpl_table_get(params->rangetab, MF_COL_WN_END,  range, NULL)             ;
        pixel_res[range] = cpl_table_get(params->rangetab, MF_COL_PIX_RES, range, NULL) * MF_OVERSAMPLING_FACTOR;

        if (range == 0 || no_single_spec) {

            /* Skip empty ranges */
            if (wn_end[range] == 0.) {

                spec_out[range]     = NULL;
                range_status[range] = CPL_ERROR_NONE;

            } else if (ismolcol == 1) {

                if (pixel_res[range] == 0.) {

                    spec_out[range]     = NULL;
                    range_status[range] = CPL_ERROR_NONE;

                } else if(WRITE_OUTPUTS)  {

                    /* Create working folders */
                    if (mf_io_mkdir(run_w_dir_range[range]) != CPL_ERROR_NONE) {
                        err += cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                                        "Could not create working directory");
                    }
                }
            }
        }
    }


    if (!err) {

        /* Get name LBLRTM output : Emission/Transmission */
        char *lblrtm_out_filename;
        cpl_boolean spec_emission = params->config->inputs.transmission == MF_PARAMETERS_TRANSMISSION_FALSE;
        if (spec_emission) lblrtm_out_filename = cpl_sprintf("%s", MF_AER_TAPE27_FILE);
        else               lblrtm_out_filename = cpl_sprintf("%s", MF_AER_TAPE28_FILE);

        double   *wn1                 = cpl_calloc(nrange, sizeof(double    ));
        double   *wn2                 = cpl_calloc(nrange, sizeof(double    ));
        cpl_size *num_wavenumber      = cpl_calloc(nrange, sizeof(cpl_size  ));
        cpl_size **name_wavenumber    = cpl_calloc(nrange, sizeof(cpl_size *));
        char     ***folder_wavenumber = cpl_calloc(nrange, sizeof(char *    ));

        /* Execute ranges and update the number of lblrtm function calls */
        for (cpl_size range = 0; range < nrange; range++) {

            /* Get wavenumber range */
            wn1[range] = cpl_table_get(params->rangetab, MF_COL_WN_START, range, NULL);
            wn2[range] = cpl_table_get(params->rangetab, MF_COL_WN_END,   range, NULL);

            /* Prepare the execution files */
            if ( (range == 0 || no_single_spec)  && wn_end[range] != 0. && ismolcol == 1 && pixel_res[range] != 0. ) {

                cpl_size wn_jmax = 0;
                double   min_wn  = wn2[range];
                double   max_wn  = min_wn;
                for (cpl_size j = 0; min_wn > wn1[range]; j++) {

                    double delta = MF_LBLRTM_DELTA_FACTOR * pow(max_wn, MF_LBLRTM_DELTA_TOL);

                    if (delta > MF_LBLRTM_DELTA_MAX) delta = MF_LBLRTM_DELTA_MAX;
                    if (delta < MF_LBLRTM_DELTA_MIN) delta = MF_LBLRTM_DELTA_MIN;
                    //MNB-delta=delta/40;
                    min_wn = CPL_MAX(max_wn - delta, wn1[range]);
                    max_wn = min_wn;

                    wn_jmax = j;
                }

                /* Calculation of frequency for refractive geometry. It need to convert wavelength to wavenumber and calculate */
                double minc_config = MF_CONV_K_LAM / lblrtm_config ->v[0];
                double maxc_config = MF_CONV_K_LAM / lblrtm_config ->v[1];

                /* Check wavelength range */
                if (minc_config < 0 || maxc_config < 0 || fabs(maxc_config - minc_config) > MF_LBLRTM_DELTA_ABS) {
                    err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                                "LBLRTM variable delta too big (|%g - %g| = %g) => |maxc - minc| > %g cm-1",
                                                maxc_config, minc_config, fabs(maxc_config - minc_config), MF_LBLRTM_DELTA_ABS);
                }

                /* In case the wavenumber range is too large for LBLRTM, break it up
                 * into smaller subranges */
                cpl_msg_info(cpl_func,"MNBXWaveNumber nu1=%f nu2=%f",wn1[range],wn2[range]);
                cpl_vector* end_points=wavenumber_subranges(wn1[range],wn2[range]);
                int npts=cpl_vector_get_size(end_points);
                cpl_msg_info(cpl_func,"MNBXWaveNumber Subranges = %d",npts);
                cpl_vector_delete(end_points);

                double     vbar        = (minc_config + maxc_config) / 2.;
                double     angle       = 90. - params->config->ambient.telescope_angle.value;
                const char *lbl_molecs;
                if (oda_parameters==NULL) {
                    lbl_molecs=params->config->internal.molecules.lbl_molecs;
                } else {
                    lbl_molecs=oda_parameters->lbl_molecs;
                }
                /* Ranges */
                double         minc;
                double         maxc;

                num_wavenumber[   range] = 0;
                name_wavenumber[  range] = NULL;
                folder_wavenumber[range] = NULL;

                for (cpl_size i = 0; i < 2; i++) {

                    /* 1er run : Count the number of valid wavenumbers */
                    /* 2nd run : Execute the configuration             */

                    if (i != 0) {
                        name_wavenumber[  range] = cpl_calloc(num_wavenumber[range], sizeof(cpl_size));
                        folder_wavenumber[range] = cpl_calloc(num_wavenumber[range], sizeof(char *  ));
                    }

                    cpl_size j            = -1;
                    double   min_wn_local = wn2[range];
                    double   max_wn_local = min_wn_local;

                    for (cpl_size wavenumber = wn_jmax; min_wn_local > wn1[range] && !err; wavenumber--) {

                        double delta = MF_LBLRTM_DELTA_FACTOR * pow(max_wn_local, MF_LBLRTM_DELTA_TOL);
                        if (delta > MF_LBLRTM_DELTA_MAX) delta = MF_LBLRTM_DELTA_MAX;
                        if (delta < MF_LBLRTM_DELTA_MIN) delta = MF_LBLRTM_DELTA_MIN;
                        //delta=delta/40;
                        cpl_msg_info(cpl_func,"MARBEL-> [%f,%f],[%f,%f],wn_jmax=%lld,delat=%f",wn1[range],wn2[range],min_wn_local,max_wn_local,wn_jmax,delta);

                        min_wn_local = CPL_MAX(max_wn_local - delta, wn1[range]);

                        if (lnfl_config->use_ODA) {
                            if (wavenumber == wn_jmax) maxc = max_wn_local;
                            else                       maxc = max_wn_local + 1; // MF_EXTRA_WN_COVERAGE;

                            if (wavenumber == 0      ) minc = min_wn_local;
                            else                       minc = min_wn_local - 1; //- MF_EXTRA_WN_COVERAGE;
                        } else {
                            if (wavenumber == wn_jmax) maxc = max_wn_local;
                            else                       maxc = max_wn_local + MF_EXTRA_WN_COVERAGE;

                            if (wavenumber == 0      ) minc = min_wn_local;
                            else                       minc = min_wn_local - MF_EXTRA_WN_COVERAGE;
                        }

                        /* Update for the next loop */
                        max_wn_local = min_wn_local;

                        if (i == 0) {

                            /* Increase number of waves */
                            num_wavenumber[range]++;

                        } else {

                            /* Save the concrete number of wave */
                            j++;
                            name_wavenumber[  range][j] = wavenumber;
                            folder_wavenumber[range][j] = cpl_sprintf("%s/wv_number_%lld", run_w_dir_range[range], wavenumber);

                            /* Check inputs parameters */
                            if (minc < 0 || minc >= maxc || fabs(maxc - minc) > MF_LBLRTM_DELTA_ABS) {

                                err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                                             "(|%g - %g| = %g) => |maxc - minc| > %g cm-1",
                                                             maxc, minc, fabs(maxc - minc), MF_LBLRTM_DELTA_ABS);
                            } else {

                                if (WRITE_OUTPUTS) {
                                    err = mf_io_mkdir(folder_wavenumber[range][j]);
                                    if (!err) {

                                      /* Create a symbolic link to TAPE3 (output of LNFL) and
                                       *write the TAPE5 input file by LNFL execution */
                                      char* tape3;
                                      if (oda_parameters==NULL) {
                                          tape3 = cpl_sprintf("%s",
                                            cpl_table_get_string(params->rangetab,MF_COL_LNFL,range));
                                      } else {
                                          tape3=cpl_sprintf("%s/range_%lld/%s",
                                                            oda_parameters->lnfl_wdir,range + 1,
                                                            MF_AER_TAPE3_FILE);
                                      }

                                      err = mf_io_write_lblrtm_configuration(folder_wavenumber[range][j], tape3,
                                                                             minc, maxc,
                                                                             vbar, angle, spec_emission, lbl_molecs,
                                                                             lblrtm_config ,
                                                                             atm_profile);
                                      cpl_free(tape3);

                                    } /*end if (!err) */

                                } /*end if (WRITE_OUTPUTS) */

                            } /*end if (minc < 0 || minc >= maxc || fabs(maxc - minc) > MF_LBLRTM_DELTA_ABS)*/

                        } /* end if (i==0) */

                    }/* end wavenumber loop */

                }/* end i loop */

                execute_range[range] = (!err) ? CPL_TRUE : CPL_FALSE;

            } else {

                execute_range[range] = CPL_FALSE;

            }/* end if (range == 0 || no_single_spec)  && wn_end[range] != 0. && ismolcol == 1 && pixel_res[range] != 0. )*/

        } /* end range loop */

        if (!err) {

            cpl_msg_info(cpl_func, "(mf_lblrtm    ) Executing external calls to the LBLRTM binary ...");

            char *lblrtm_syscall;
            if (params->config->inputs.silent_external_bins) {
                lblrtm_syscall = cpl_sprintf("%s/%s", params->config->directories.telluric_path, MF_BIN_LBLRTM_SILENT);
            } else {
                lblrtm_syscall = cpl_sprintf("%s/%s", params->config->directories.telluric_path, MF_BIN_LBLRTM);
            }

            int    num_calls      = 0;
            double runtime_lblrtm = 0.;
            #ifdef _USE_OPENMP
                /* OMP parallel version */
                /*if (params->config->inputs.omp_num_threads > 1) {
                    cpl_msg_info(cpl_func, "(mf_lblrtm   ) OPENMP in parallel by wavelength ranges");
                }*/
                #pragma omp parallel for shared(lblrtm_syscall, execute_range, num_wavenumber, folder_wavenumber) reduction(+ : err, num_calls, runtime_lblrtm)
            #endif
            for (cpl_size range = 0; range < nrange; range++) {
                if (execute_range[range]) {
                    for (cpl_size wavenumber = 0; wavenumber < num_wavenumber[range]; wavenumber++) {
                        double runtime;
                        cpl_matrix* odTable;
                        odTable = mf_io_oda_get_tableDB(range);
                        /*if (!USE_ODATABLE || (USE_ODATABLE && odTable==NULL)) {*/
                        if ((!lnfl_config->use_ODA) || (lnfl_config->use_ODA && odTable==NULL)){
                            #ifdef _USE_OPENMP
                                err += mf_io_systemOpenMP(lblrtm_syscall, folder_wavenumber[range][wavenumber], &runtime);
                            #else
                                err += mf_io_system      (lblrtm_syscall, folder_wavenumber[range][wavenumber], &runtime);
                            #endif
                            runtime_lblrtm += runtime;
                        }
                    }
                    num_calls++;
                }
            }
            *lblrtm_calls                  += num_calls;
            params->timers.time_bin_lblrtm += runtime_lblrtm;

            mf_io_sync();

            cpl_free(lblrtm_syscall);

            if (!err) {

                for (cpl_size range = 0; range < nrange && !err; range++) {

                    if (execute_range[range]) {

                        /* Rename and move output file */
                        for (cpl_size wavenumber = 0; wavenumber < num_wavenumber[range]; wavenumber++) {

                                char *old_output_name = cpl_sprintf("%s/%s",      folder_wavenumber[range][wavenumber], lblrtm_out_filename                                        );
                                char *new_output_name = cpl_sprintf("%s/%s_%lld", run_w_dir_range[range],               lblrtm_out_filename, name_wavenumber[range][wavenumber] + 1);
                                cpl_matrix* odTable;
                                odTable = mf_io_oda_get_tableDB(range);
                                /*if (!USE_ODATABLE || (USE_ODATABLE && odTable==NULL)) {*/
                                if ((!lnfl_config->use_ODA) || (lnfl_config->use_ODA && odTable==NULL)){
                                    err = mf_io_mv(old_output_name, new_output_name);
                                }

                                if (err!=0) exit(1);

                                cpl_free(old_output_name);
                                cpl_free(new_output_name);
                        }


                        if (err != CPL_ERROR_NONE) {

                            cpl_error_set_message(cpl_func, err,
                                                    "LBLRTM failed for wavelength interval: %5.2f-%5.2f [µm]",
                                                    MF_CONV_K_LAM / wn1[range], MF_CONV_K_LAM / wn2[range]);
                        } else if (spec_out!=NULL) {

                            /* Get pixel resolution */
                            double local_pixel_res = CPL_MIN(1e6, pixel_res[range]);

                            /* Create spectrums combine the outputs for all the wavenumbers */
                            spec_out[range] = mf_lblrtm_range_combined_and_convolved(
                                lblrtm_out_filename, local_pixel_res,
                                wn1[range], wn2[range],
                                run_w_dir_range[range],range,mol_abuns,
                                USE_ODATABLE,
                                lnfl_config
                            );

                            /* Save output error */
                            err = (spec_out[range]) ? CPL_ERROR_NONE : CPL_ERROR_ILLEGAL_OUTPUT;
                            range_status[range] = err;
                        }
                    }
                }
            }
        }

        /* Cleanup */
        for (cpl_size range = 0; range < nrange; range++) {

            if (name_wavenumber[range]) cpl_free(name_wavenumber[range]);

            if (folder_wavenumber[range]) {
                for (cpl_size wavenumber = 0; wavenumber < num_wavenumber[range]; wavenumber++){
                    cpl_free(folder_wavenumber[range][wavenumber]);
                }
                cpl_free(folder_wavenumber[range]);
            }
        }
        cpl_free(num_wavenumber);
        cpl_free(name_wavenumber);
        cpl_free(folder_wavenumber);
        cpl_free(wn1);
        cpl_free(wn2);
        cpl_free(lblrtm_out_filename);
    }

    /* Cleanup */
    for (cpl_size range = 0; range < nrange; range++) {
        cpl_free(run_w_dir_range[range]);
    }
    cpl_free(run_w_dir_range);
    cpl_free(wn_end         );
    cpl_free(pixel_res      );
    cpl_free(execute_range  );

    return err;
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
static cpl_table * mf_lblrtm_range_combined_and_convolved(
    const char               *lblrtm_out_filename,
    const double             pixel_res,
    const double             wn1,
    const double             wn2,
    const char               *w_dir_range,const int range, cpl_vector* mol_abuns, cpl_boolean USE_ODATABLE,
    mf_io_lnfl_config        *lnfl_config)
{
    /* Merges and rebins a set of LBLRTM radiance ("*R.*") or transmission ("*T.*") spectra */
    double limlam[2] = { MF_CONV_K_LAM / wn2,
                         MF_CONV_K_LAM / wn1 };

    /* Initialize wavelength grid */
    cpl_table *spec = cpl_table_new(0);
    cpl_table_new_column(spec, MF_COL_IN_LAMBDA, CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, MF_COL_IN_FLUX,   CPL_TYPE_DOUBLE);

    /* Create grid */
    cpl_error_code err = mf_lblrtm_create_wavelength_grid(spec, limlam, pixel_res);
    if (err != CPL_ERROR_NONE) {
        cpl_error_set_message(cpl_func, err,
                              "Error in mf_lib_createlibspec");
        cpl_table_delete(spec);
        return NULL;
    }

    /*** Rebin set of LBLRTM spectra (wrapper output) in wavelength units [mu m] (variable step size possible) **/

    /* Find number of files and wavenumber ranges. */
    cpl_array *klim_all;
    if (USE_ODATABLE) {
        klim_all=mf_io_klim_from_odatable(range);
        if (klim_all==NULL) klim_all = mf_io_find_klim(w_dir_range, lblrtm_out_filename);
    } else {
        klim_all = mf_io_find_klim(w_dir_range, lblrtm_out_filename);
    }
    if (!klim_all) {
        err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                    "Invalid input parameter(s): klim_all NULL!");
    } else {

        /* Rebin spectra */
        cpl_size       size = cpl_array_get_size(klim_all);
        for (cpl_size wavenumber = 0; wavenumber < size - 1 && !err; wavenumber++) {

            double klim[2] = { cpl_array_get(klim_all, wavenumber,     NULL),
                               cpl_array_get(klim_all, wavenumber + 1, NULL)  };

            /* Check wavenumber range */
            if (klim[1] <= klim[0]) {
                err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                             "Invalid input parameter(s): klim[1] <= klim[0] (wavenumber limits)");
            } else {

                /* Convert wavenumbers to wavelengths */
                double llim[2] = { MF_CONV_K_LAM / klim[1],
                                   MF_CONV_K_LAM / klim[0] };

                /* Check number of data points in spectrum */
                int nrow = cpl_table_get_nrow(spec);
                if (nrow == 0) {
                    err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                                "No data: cpl_table *spec");
                }

                /* Get pointers to CPL table columns */
                double *lamv  = cpl_table_get_data_double(spec, MF_COL_IN_LAMBDA);
                double *fluxv = cpl_table_get_data_double(spec, MF_COL_IN_FLUX  );

                /* Compare wavelength ranges */
                if (llim[1] < lamv[0] || llim[0] > lamv[nrow - 1]) {

                    /* Desired wavelength range not covered by the input model spectrum */
                    err = CPL_ERROR_NONE;

                } else {

                    /* Output variable of Rebin -> mf_io_read_lblrtm_and_update_spec(...) */
                    cpl_boolean usampl =  CPL_FALSE;
                    int         jmin   = -1;
                    int         jmax   =  0;

                    /* Rebin spectrum */
                    char *spectrum_filename = cpl_sprintf("%s/%s_%lld", w_dir_range, lblrtm_out_filename, wavenumber + 1);
                    cpl_bivector* bvec;
                    bvec=NULL;
                    if (lnfl_config->use_ODA) {
                        bvec=mf_io_mergeODTables(range,mol_abuns,spectrum_filename);
                        if (bvec!=NULL && USE_ODATABLE==CPL_FALSE) {
                            cpl_bivector_delete(bvec);
                            bvec=NULL;
                        }
                        if (bvec==NULL) bvec = mf_io_read_lblrtm_spec(spectrum_filename);
                    }
                    err = mf_io_read_lblrtm_and_update_spec(
                        nrow, lamv, fluxv,
                        spectrum_filename,
                        bvec,
                        llim, &usampl, &jmin, &jmax,
                        lnfl_config
                    );
                    cpl_free(spectrum_filename);
                    if (lnfl_config->use_ODA) {
                        cpl_bivector_delete(bvec);
                    }

                    /* Interpolate "empty" bins by means of the closest bins that contain data points (valid flux values) */
                    if (usampl) mf_lblrtm_linear_interpolate_spectrum(lamv, fluxv, jmin, jmax);
                }
            }
        }
    }

    /* Cleanup */
    cpl_array_delete(klim_all);

    /* Check error and return */
    if (!err) {
        return spec;
    } else {
      cpl_table_delete(spec);
      return NULL;
    }
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
static cpl_error_code mf_lblrtm_create_wavelength_grid(
    cpl_table                *outspec,
    const double             limlam[2],
    const double             resol)
{
    /*!
     * Initialisation of a wavelength grid that is characterised by constant
     * resolution lambda / delta lambda
     *
     * \b INPUT:
     * \param limlam  lower and upper limit of wavelength range
     * \param resol   resolution lambda / delta lambda
     *
     * \b OUTPUT:
     * \param outspec  CPL table that provides wavelengths and fluxes (= 0.)
     *                 of a spectrum
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    /* Test input parameters */
    if (limlam[1] <= limlam[0] || resol <= 0.) {
        /* Return spectrum with zero data points */
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Invalid input parameter(s): limlam[1] <= limlam[0] || resol <= 0. (wavelength grid)");
    }

    /* Get number of data points and resize CPL table */
    cpl_size nrow = ceil(log(limlam[1] / limlam[0]) / log1p(1 / resol));
    cpl_table_set_size(outspec, nrow);

    /*** Create wavelength grid ***/

    /* Set first position */
    double lam = limlam[0];
    cpl_table_set(outspec, MF_COL_IN_LAMBDA, 0, lam);

    /* Set next positions */
    for (cpl_size i = 1; i < nrow; i++) {
        lam *= (1 + 1 / resol);
        cpl_table_set(outspec, MF_COL_IN_LAMBDA, i, lam);
    }

    /* Set flux = 0 */
    cpl_table_fill_column_window(outspec, MF_COL_IN_FLUX, 0, nrow, 0.);

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
static cpl_error_code mf_lblrtm_linear_interpolate_spectrum(
    double                   *lambda,
    double                   *flux,
    const cpl_size           lim_min,
    const cpl_size           lim_max)
{
  /*!
   * Linear interpolation for "empty" bins in spectra (marked by -HUGE_VAL).
   *
   * \b INPUT:
   * \param lim_min  pixel min limit of spectrum
   * \param lim_max  pixel max limit of spectrum
   *
   * \b OUTPUT:
   * \param spec    spectrum corrected by interpolation
   *
   * \b ERRORS:
   * - No data
   */

  /* Perform interpolation/extrapolation procedure */
  int            jlim[2] = {-1, -1};
  cpl_boolean    ok      = CPL_TRUE;
  cpl_error_code status  = CPL_ERROR_NONE;

  for (cpl_size j = lim_min; j <= lim_max; j++) {

      if (flux[j] != -HUGE_VAL) {

          if (ok == CPL_TRUE) {
              jlim[0] = j;
          } else {
              jlim[1] = j;
              if (jlim[0] == -1) {
                  /* Constant value for extrapolation of lower margin */
                  for (cpl_size i = lim_min; i <= jlim[1] - 1; i++) {
                      flux[i] = flux[jlim[1]];
                  }
              } else {
                  /* Linear interpolation */
                  double m = (flux[jlim[1]] - flux[jlim[0]]) / (lambda[jlim[1]] - lambda[jlim[0]]);
                  for (cpl_size i = jlim[0] + 1; i <= jlim[1] - 1; i++) {
                      flux[i] = m * (lambda[i] - lambda[jlim[0]]) + flux[jlim[0]];
                  }
              }
              jlim[0] = jlim[1];
              ok = CPL_TRUE;
          }

      } else {

          ok = CPL_FALSE;

          if (j == lim_max) {
              if (jlim[0] == -1) {
                  status = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                                 "No data: No valid flux value");
              } else {
                  /* Constant value for extrapolation of upper margin */
                  for (cpl_size i = jlim[0] + 1; i <= lim_max; i++) {
                      flux[i] = flux[jlim[0]];
                  }
              }
          }

      }

  }

  return status;
}


cpl_error_code mf_io_lblrtm_oda(mf_io_lnfl_config   *lnfl_config,
                                mf_io_lblrtm_config *lblrtm_config,
                                mf_parameters       *params,
                                cpl_table           *atm) {

    /*
     * --------
     * PURPOSE:
     * --------
     * Generating and loading into memory the ODa tables needed for the ODa process.
     *
     * LBLRTM has to be executed for individual molecule profiles, i.e. TAPE5 has to be
     * generated for each molecule; the relevant TAPE3 associated with that molecule has
     * to be soft linked; and the trasmission data has to be obtained from the output
     * TAPE28 files for each molecule.
     * As much usage as possible is made from mf_lblrtm_range_execution to prep and call
     * LBLRTM.
     * Note:
     * mf_lblrtm_range_execution may decide that wavenumber range is too long and needs
     * breaking down into several subregions: wavenumber_0,wavenumber_1 etc and call LBLRTM
     * executions in each of these subdirectores and place the outputs (TAPE28 files) in
     * the specified working directory with wavenumber range suffixes eg TAPE28_1, TAPE28_2
     * etc.
     *
     */

    cpl_msg_info(cpl_func,"Calling ODA LBLRTM");
    cpl_boolean debug = mf_io_use_debug();
    cpl_error_code err=CPL_ERROR_NONE;


    int nrows=cpl_table_get_nrow(atm);
    int ncols=cpl_table_get_ncol(atm);
    cpl_msg_info(cpl_func,"Table Details nrows=%d, ncols=%d",nrows,ncols);

    if (lnfl_config==NULL)   cpl_msg_info(cpl_func,"lnfl_config is NULL");
    if (lblrtm_config==NULL) cpl_msg_info(cpl_func,"lblrtm_config is NULL");
    if (params==NULL)        cpl_msg_info(cpl_func,"params is NULL");

    /* Initialise the molecule order permutation for the oda tables       */
    /* Note: this is needed because the oda table works with molecules in */
    /* the LBLRTM order, but the parameters listed in the recipie may be  */
    /* in a different order                                               */
    char **mol = cpl_table_get_data_string(params->molectab, MF_COL_LIST_MOLECULES);
    err=mf_oda_mol_idx(mol,params->config->internal.molecules.n_molec);


    /* Molecule string */
    /* We need to parse the molecule string to determine which molecules we
     * are using and to be able to create strings relevent for simulatint
     * single molecules only
     */
    char *lbl_molecs = params->config->internal.molecules.lbl_molecs;
    int lbl_size=strlen(lbl_molecs);

    /* Create a default all '0' molecule string */
    char zero_str[MF_MOLECULES_NUMBER_MAX+1];
    for (int i=0;i<MF_MOLECULES_NUMBER_MAX;i++) zero_str[i]='0';
    zero_str[MF_MOLECULES_NUMBER_MAX]='\0';
    cpl_msg_info(cpl_func,"molecular string = %s",lbl_molecs);

    /* Convert the molecule string into an array of molecule names*/
    cpl_array* mol_array=mf_io_molecstring2Names(lbl_molecs);
    int nmols=cpl_array_get_size(mol_array);

    /* Get number of ranges that we are dealing with */
    int nrange = params->config->internal.n_range;
    cpl_msg_info(cpl_func,"Nranges=%d", nrange);
    if (params->config->internal.single_spectrum) {
        cpl_msg_info(cpl_func,"=========<SINGLE SPECTRUM IS TRUE>=======");
        cpl_msg_info(cpl_func,"Nranges=%d", nrange);
        nrange = 1;
    }

    /* Create an nrange by nmols array placeholder for bivectors
     * in which we will store the optical depth values.
     */
    cpl_bivector* optical_depths[nrange][nmols];
    for (int range=0; range <nrange; range++) {
        for (int mol_idx=0;mol_idx<nmols;mol_idx++) {
            optical_depths[range][mol_idx]=NULL;
        }
    }

    /*Create Optical Depth working directory with molecule subdirectories*/
    char* od_dir = cpl_sprintf("%s/%s",
                               params->config->internal.tmp_folder,
                               "OpticalDepthsWorkDir");
    mf_io_mkdir(od_dir);
    for (int mol_idx=0;mol_idx<nmols;mol_idx++) {

        /* Generate a lbl string for this molecule only */
        for (int i=0, cnt=0;i<lbl_size;i++) {
            zero_str[i]='0';
            if (lbl_molecs[i]=='1') if (cnt++==mol_idx) zero_str[i]='1';
        }

        /* Create a molecule specific sub directory*/
        const char* mol_name=cpl_array_get_string(mol_array,mol_idx);
        char* mol_dir=cpl_sprintf("%s/%s", od_dir,mol_name);
        mf_io_mkdir(mol_dir);

        /* Create an lblrtm sub work directory for this molecule */
        char* lblrtm_wdir=cpl_sprintf("%s/%s", mol_dir,"lblrtm");
        mf_io_mkdir(lblrtm_wdir);

        /* Create a lnfl pseudo sub work directory for this molecule
           (pseudo because we are using it to softlink TAPE3s from the true
           workng dir to make it easier for lblrtm_range_execution to find it
        */
        char*   lnfl_wdir=cpl_sprintf("%s/%s", mol_dir,"lnfl"  );
        mf_io_mkdir(  lnfl_wdir);

        /* To leverage usage of mf_lblrtm_range_execution, we added an extra parameter
         * oda_parameter which is a pointer to a structure that specifies the molecule
         * string to use and the string names for the working directories
         */
        mf_io_oda_parameters oda_parameter;
        oda_parameter.lbl_molecs  = zero_str;
        oda_parameter.lnfl_wdir   = lnfl_wdir;
        oda_parameter.lblrtm_wdir = lblrtm_wdir;

        cpl_msg_info(cpl_func,"molecule %s wdir=%s",mol_name, mol_dir);
        cpl_msg_info(cpl_func,"lnfl___wdir=%s",oda_parameter.lnfl_wdir);
        cpl_msg_info(cpl_func,"lblrtm_wdir=%s",oda_parameter.lblrtm_wdir);

        /* Create range subdirectories for the lnfl wdir*/
        for (int range = 0; range < nrange; range++) {

            /* Create the range dir*/
            char* w_dir_range = cpl_sprintf("%s/%s_%d",lnfl_wdir,"range", range + 1);
            mf_io_mkdir(w_dir_range);

            /* Softlink TAPE3 from the true working directory to this one
             * (we do this to make it easier to find the correct TAPE3 in
             * lblrtm_range_execution */
            char* tape3_dest=cpl_sprintf("%s/%s",w_dir_range,"TAPE3");
            cpl_msg_info(cpl_func,"TAPE3 dest=%s",tape3_dest);
            char* assoc_lnfl_dir=cpl_sprintf("%s/%s_%d/%s",
                                             params->config->internal.tmp_folder,
                                             MF_AER_WDIR_LNFL_RANGE_PATH,
                                             range+1,mol_name);

            cpl_msg_info(cpl_func,"lnfl working dir for range %d = %s",
                                                     range,w_dir_range);
            cpl_msg_info(cpl_func,"lnfl associa dir for range %d = %s",
                                                  range,assoc_lnfl_dir);
            char* tape3_targ=cpl_sprintf("%s/%s",assoc_lnfl_dir,
                                                               "TAPE3");
            mf_io_oda_symlink(tape3_targ,tape3_dest);
            cpl_free(tape3_targ);
            cpl_free(tape3_dest);
            cpl_free(assoc_lnfl_dir);
            cpl_free(w_dir_range);
        }/* end of range loop*/

        /* Call LBLRTM executions for ths range definition
            * Note: there may be more than one execution as mf_lblrtm_range_execution
            * and the outputs will be in one or more TAPE28 file that need to be merged.
            */
        cpl_msg_info(cpl_func,"Attempt to call lblrtm from here");
        int lblrtm_calls=0;
        cpl_size run_execution=0;
        const int ismolcol=1;
        err=mf_lblrtm_range_execution(
                                            lblrtm_config,
                                            lnfl_config,
                                            params,
                                            ismolcol,
                                            run_execution,
                                            &lblrtm_calls,
                                            atm,
                                            NULL,
                                            NULL,
                                            NULL,
                                            &oda_parameter);
        cpl_msg_info(cpl_func,"Finished lblrtm with error %d",err);

        /* Load the transmission data from LBLRTM's output TAPE28 files */
        for (int range = 0; range < nrange; range++) {
            char* range_dir=cpl_sprintf("%s/range_%d",lblrtm_wdir,range+1);
            cpl_bivector* bvec = mf_io_merge_wavefiles(range_dir,"TAPE28");
            cpl_free(range_dir);
            if (bvec==NULL) {
                cpl_msg_info(cpl_func,"Failed to load bvec");
            } else {
                optical_depths[range][mol_idx]=bvec;
                int n=cpl_bivector_get_size(bvec);
                cpl_vector* x=cpl_bivector_get_x(bvec);
                cpl_msg_info(cpl_func,"Loaded %d numbers [%f,%f]",n,
                             cpl_vector_get(x,0),cpl_vector_get(x,n-1));
            }

            char * tape28=cpl_sprintf("%s/range_%d/%s",lblrtm_wdir,
                                      range+1,"TAPE28_1");
            cpl_msg_info(cpl_func,"Attempting to read %s",tape28);
            cpl_bivector* bvec2 = mf_io_read_lblrtm_spec(tape28);
            cpl_free(tape28);
            if (bvec2==NULL) {
                cpl_msg_info(cpl_func,"Failed to load bvec2");
            }

            FILE *stream;
            stream=fopen("bvec0.dat","w");
            cpl_bivector_dump(bvec,stream);
            fclose(stream);
            stream=fopen("bvec2.dat","w");
            cpl_bivector_dump(bvec2,stream);
            fclose(stream);
            cpl_bivector_delete(bvec2);

            /* Pause for debug*/
            if (debug) {
                printf("Pausing for debug purposes\n");
                getchar();
            }

        }/* end of range loop*/
        cpl_free(oda_parameter.lnfl_wdir);
        cpl_free(oda_parameter.lblrtm_wdir);

    }/* end for mol_idx*/

    /* Define a optical depth table for each range. Note this means rebinning the bivectors
     in each range into a single wavenumber axis*/
    for (int range=0; range <nrange; range++) {

        /* Determine which molecule specific bivector in this range has the most
         * data points. Note the wavenumber steps will be uniform between two fixed
         * points so the bivector with the most points has the smallest stepsize */
        int max_n=0;
        int ref_idx=0;
        for (int mol_idx=0; mol_idx <nmols; mol_idx++) {
            int n=cpl_bivector_get_size(optical_depths[range][mol_idx]);
            if (n>max_n) {
                max_n=n;
                ref_idx=mol_idx;
            }
        }

        cpl_msg_info(cpl_func,"Range %d Max bivector size=%d",range,max_n);

        /* Rebin via interpolation all other bivectors to the axis of max_idx */
        cpl_vector* ref_axis = cpl_bivector_get_x(optical_depths[range][ref_idx]);
        double ref_axis_min=cpl_vector_get(ref_axis,0);
        double ref_axis_max=cpl_vector_get(ref_axis,max_n-1);
        cpl_msg_info(cpl_func,"Got XMNBX chosen axis %d Range=[%f,%f]",ref_idx,ref_axis_min,ref_axis_max);
        for (int mol_idx=0; mol_idx <nmols; mol_idx++) {

            if (mol_idx==ref_idx) continue;

            /* Get the old bivector */
            cpl_bivector* old_bvec=optical_depths[range][mol_idx];
            int old_size=cpl_bivector_get_size(old_bvec);

            /* Define a new bivector as with the same axis as that of max_i*/
            cpl_vector* new_x=cpl_vector_new(max_n);
            cpl_vector* new_y=cpl_vector_new(max_n);
            for (int i=0;i<max_n;i++) cpl_vector_set(new_y,i,0.0);
            cpl_vector_copy(new_x,ref_axis);
            cpl_msg_info(cpl_func,"Creates newxy vector for %d",mol_idx);
            cpl_bivector* new_bvec=cpl_bivector_wrap_vectors(new_x,new_y);
            cpl_msg_info(cpl_func,"Created new bivector for %d",mol_idx);

            /* Now rebin/interpolate the old bivec to the new axis */
            cpl_msg_info(cpl_func,"About to Rebin new bivector for %d s1=%d, s2=%d",mol_idx,old_size,max_n);
            cpl_msg_info(cpl_func,"old bvec size =%d",old_size);
            int old_n= cpl_bivector_get_size(old_bvec);
            cpl_vector* old_x=cpl_bivector_get_x(old_bvec);
            cpl_vector* old_y=cpl_bivector_get_y(old_bvec);
            double old_axis_min=cpl_vector_get(old_x,0);
            double old_axis_max=cpl_vector_get(old_x,old_n-1);
            cpl_msg_info(cpl_func,"compare min axes %f %f",ref_axis_min,old_axis_min);
            if (ref_axis_min<old_axis_min) {
                cpl_msg_info(cpl_func,"Expect error because ref_axis_min<old_axis_min ie %f < %f",ref_axis_min,old_axis_min);
                cpl_msg_info(cpl_func,"Changeing old axis min");
                cpl_vector_set(old_x,0,ref_axis_min);
            }
            cpl_msg_info(cpl_func,"compare max axes %f %f",ref_axis_max,old_axis_max);
            if (ref_axis_max>old_axis_max) {
                cpl_msg_info(cpl_func,"Expect error because ref_axis_max>old_axis_max ie %f > %f",ref_axis_max,old_axis_max);
                cpl_msg_info(cpl_func,"Changeing old axis max");
                cpl_vector_set(old_x,old_n-1,ref_axis_max);
            }

            for (int i=0;i<5;i++) {
                double xval=cpl_vector_get(new_x,i);
                double yval=cpl_vector_get(new_y,i);
                double xold=cpl_vector_get(old_x,i);
                double yold=cpl_vector_get(old_y,i);
                cpl_msg_info(cpl_func,"New Axis i=%d x=%f y=%f cf old x=%f y=%f",i,xval,yval,xold,yold);
            }

            err=cpl_bivector_interpolate_linear(new_bvec,old_bvec);
            new_x=cpl_bivector_get_x(new_bvec);
            new_y=cpl_bivector_get_y(new_bvec);
            cpl_msg_info(cpl_func,"Rebinned new bivector for %d err=%d",mol_idx,err);
            if (err) {
                cpl_msg_info(cpl_func,"Rebinned Error err=%d",err);
                cpl_msg_info(cpl_func,"CPL_ERROR_NULL_INPUT=%d",CPL_ERROR_NULL_INPUT);
                cpl_msg_info(cpl_func,"CPL_ERROR_DATA_NOT_FOUND=%d",CPL_ERROR_DATA_NOT_FOUND);
                cpl_msg_info(cpl_func,"CPL_ERROR_ILLEGAL_INPUT=%d",CPL_ERROR_ILLEGAL_INPUT);

               // cpl_error_reset();
                if (old_bvec==NULL) {
                    cpl_msg_info(cpl_func,"old bvec=NULL");
                } else {
                    cpl_msg_info(cpl_func,"old bvec size =%d",old_size);
                }
                return err;
            }

            for (int i=0;i<5;i++) {
                double xval=cpl_vector_get(old_x,i);
                double yval=cpl_vector_get(old_y,i);
                double xval2=cpl_vector_get(new_x,i);
                double yval2=cpl_vector_get(new_y,i);
                cpl_msg_info(cpl_func,"i=%d x=%f y=%f ->  x=%f y=%f",i,xval,yval,xval2,yval2);
            }

            /* Replace the old bivec in the optical_depths table with the new*/
            optical_depths[range][mol_idx]=new_bvec;
            cpl_msg_info(cpl_func,"Replaced old bivector for new for %d",mol_idx);

            /* Dont forget to delete the old bvec*/
            cpl_bivector_delete(old_bvec);
            cpl_msg_info(cpl_func,"Deleted old bivector for %d",mol_idx);

        }/* end mol_idx loop */

        /* Now we can create the OD Table for this range*/
        for (int mol_idx=0; mol_idx<nmols;mol_idx++) {

            cpl_bivector* bvec=optical_depths[range][mol_idx];

            /* Get the size and the y data of this bivector */
            int m=cpl_bivector_get_size(bvec);
            double* y=cpl_bivector_get_y_data(bvec);

            /* Optical depth = -log(transission value) */
            for (int i=0; i<m;i++) {
                double trans=y[i];
                trans=CPL_MAX(trans,1.0e-08);
                double odval=-1.0*log(trans);
                y[i]=odval;
            }

            if (mol_idx==0) {

                /* Allocate the in memory Optical Depth Table size for this range*/
                int NCOLS=nmols+1;
                mf_io_oda_init_tableDB(range,NCOLS,m);

                /* Store the wavenumber values as the 0th column*/
                double* x=cpl_bivector_get_x_data(bvec);
                int WAVNUM_COLUMN_idx=0;
                mf_io_oda_set_tableDB(range,WAVNUM_COLUMN_idx,x, m);


            } /* endif mol_idx==0 */

            /*Load up the optical depth values for the column assigned to this molecule*/
            int column_idx=mol_idx+1;
            mf_io_oda_set_tableDB(range,column_idx,y,m);

        }/* end mol_idx loop */

    }/* end range loop */


    /* Cleanup*/
    cpl_array_delete(mol_array);
    cpl_free(od_dir);

    for (int range=0; range <nrange; range++) {
        for (int mol_idx=0; mol_idx <nmols; mol_idx++) {
            cpl_bivector* bvec=optical_depths[range][mol_idx];
            if (bvec!=NULL) cpl_bivector_delete(bvec);
        }
    }

    return CPL_ERROR_NONE;
}

cpl_vector* wavenumber_subranges(double nu1,double nu2) {

    /* COMMENTS:
    LBLRTM uses a single wave number value for refractive calculation.
    This value is usually chosen to be the mid point of the range.
    However, this single value approximation becomes poor if the the
    wave number range for simulatiom is too large. That is, for a wave
    number range, [ν1,ν2] , the range is too large if:

        v2 – v1 > delta

    where

        delta = 0.14 x (v2)^1.1

    and delta is clamped between [125.0,1750.0] cm-1;

    (NOTE: these numbers are stored as macros)

    Under these circumstances, the wave range has to be broken up into
    n sub-ranges that do not violate this rule, i.e:

    [vb0,vb1] , [vb1,vb2], [vb2,vb3], ..., [vbi-1,vbi], ... , [vbn-1,vbn]

    The puprose of this routine is to break down a supplied range
    [nu1,nu2] as specified above and to return a vector of calculated
    end points such that the first value is nu1 and the last is nu2.
    */

    /* Derive an minumum estimate of the DELTA factor */
    double min_delta = MF_LBLRTM_DELTA_FACTOR *
                       pow(nu1, MF_LBLRTM_DELTA_TOL);
    if (min_delta>MF_LBLRTM_DELTA_MAX) min_delta=MF_LBLRTM_DELTA_MAX;
    if (min_delta<MF_LBLRTM_DELTA_MIN) min_delta=MF_LBLRTM_DELTA_MIN;


    /* Derive a maximum number for subranges required within [nu1,nu2]
       i.e. how many times will min_delta fit in the range [nu1,nu2]
    */
    int max_n = ceil((nu2-nu1)/min_delta);

    /* Allocate a buffer vector to store new end points */
    cpl_vector* buffer = cpl_vector_new(max_n);

    /* Define the first end point in the buffer as nu2 */
    cpl_size stack_idx=0;
    cpl_vector_set(buffer,stack_idx,nu2);

    /* Iterate calculating subranges [mu1,mu2] that satisfy the range
     * criteria.
    */
    for (cpl_size i=1;i<max_n; i++) {

        /* Get upper range mu2 from the last value pushed in the
         * buffer stack
        */
        double mu2=cpl_vector_get(buffer,i-1);

        /* Define the maximum valid range for [mu1,mu2] */
        double delta = MF_LBLRTM_DELTA_FACTOR *
                       pow(mu2, MF_LBLRTM_DELTA_TOL);
        if (delta > MF_LBLRTM_DELTA_MAX) delta = MF_LBLRTM_DELTA_MAX;
        if (delta < MF_LBLRTM_DELTA_MIN) delta = MF_LBLRTM_DELTA_MIN;

        /* Use delta to estimate mu1 */
        double mu1 = mu2-delta;

        /* If mu1 is greater than nu1 then add mu1 to the stack otherwise
           break from loop as we are outside of the range [nu1,nu2] */
        if (mu1>nu1) {
            /* Push this value to the stack */
            cpl_vector_set(buffer,++stack_idx, mu1);
        } else {
            break;
        }
    } /* end for i */

    /*
     Required return vector is the sequence of values: nu1 and then
     the values popped from the buffer.
    */

    /* Allocate the size for the return vector */
    cpl_vector *retv = cpl_vector_new(stack_idx+1);

    /* Put nu1 as the first value */
    cpl_size i=0;
    cpl_vector_set(retv,i,nu1);

    /* Iterate through the stack elements */
    while (stack_idx>0) {

        /* Pop value from stack */
        double nu=cpl_vector_get(buffer,stack_idx--);

        /* Add value to the return vector */
        cpl_vector_set(retv,++i,nu);
    }

    /* Cleanup */
    cpl_vector_delete(buffer);

    return retv;

}

/** @endcond */


/**@}*/
