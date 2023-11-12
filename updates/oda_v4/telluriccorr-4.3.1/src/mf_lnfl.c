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

#include "mf_lnfl.h"

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

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup mf_lnfl       .
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create LNFL configuration.
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
mf_io_lnfl_config * mf_lnfl_config_create(void)
{
    mf_io_lnfl_config *config = cpl_malloc(sizeof(mf_io_lnfl_config));

    config->line_db     = cpl_sprintf("%s", MF_LNFL_LINE_DB_INIT);
    config->line_db_fmt = MF_LNFL_LINE_DB_FMT_INIT;
    config->use_ODA     = CPL_FALSE;

    return config;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Delete LNFL configuration.
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
void mf_lnfl_config_delete(
    mf_io_lnfl_config        *config)
{
  if (config) {

      if (config->line_db) cpl_free(config->line_db);

      cpl_free(config);
  }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Check LNFL configuration.
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
cpl_error_code mf_lnfl_config_check(
    mf_io_lnfl_config        *config)
{
  cpl_ensure(config, CPL_ERROR_NULL_INPUT, CPL_ERROR_NULL_INPUT);

  // TODO: Add check for all variables and values
  cpl_error_code err = CPL_ERROR_NONE;

  if (!err && !(strcmp(config->line_db, MF_LNFL_LINE_DB_AER_VERSION) == 0)) {
      //err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LNFL_LINE_DB);       /* Only allow the last version, included in the third party */
      cpl_msg_warning(cpl_func,"LNFL_LINE_DB: Using AER version %s instead of default %s",config->line_db,MF_LNFL_LINE_DB_AER_VERSION);
  }

  if (!err && config->line_db_fmt <= 0) {
      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, MF_LNFL_LINE_DB_FMT);   /* FIXME: Do I need check MF_LNFL_LINE_DB_FMT_OLD_HITRAN and MF_LNFL_LINE_DB_FMT_NEW_HITRAN? */
  }

  return err;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Calls LNFL for the different wavelength ranges saved in the mf_parameters parameter structure.
 *
 * @param config             .
 * @param params             mf_parameters parameter structure
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - Invalid object structure
 *                           - Error in subroutine (see subroutines).
 *
 * @note AVAILABLE MOLECULAR SPECIES :
 *  1) H2O   2) CO2   3)   O3   4)  N2O   5)  CO   6) CH4   7)    O2
 *  8)  NO   9) SO2  10)  NO2  11)  NH3  12)HNO3  13)  OH  14)    HF
 * 15) HCL  16) HBR  17)   HI  18)  CLO  19) OCS  20)H2CO  21)  HOCL
 * 22)  N2  23) HCN  24)CH3CL  25) H2O2  26)C2H2  27)C2H6  28)   PH3
 * 29)COF2  30) SF6  31)  H2S  32)HCOOH  33) HO2  34)   O  35)CLONO2
 * 36) NO+  37)HOBR  38) C2H4  39)CH3OH  40)CH3Br 41)CH3CN 42)   CF4
 * 43) C4H2 44) HC3N 45)    H2 46)    CS 47)   SO3
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_lnfl(
    mf_io_lnfl_config        *config,
    mf_parameters            *params)
{
    /* Get number of ranges */
    int nrange = params->config->internal.n_range;

    /* Set and check size of range table */
    if (params->config->internal.single_spectrum) {
        nrange = 1;
    } else if (cpl_table_get_nrow(params->rangetab) != nrange) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                     "Invalid object structure: mf_parameters *params (unexpected size of rangetab)");
    }

    cpl_error_code err = CPL_ERROR_NONE;

    /* Create LNFL configuration files */
    char   **w_dir_range  = cpl_calloc(nrange, sizeof(char *));
    double *wn_start      = cpl_calloc(nrange, sizeof(double));
    double *wn_end        = cpl_calloc(nrange, sizeof(double));
    for (int range = 0; range < nrange; range++) {

        w_dir_range[range] = cpl_sprintf("%s/%s_%d", params->config->internal.tmp_folder, MF_AER_WDIR_LNFL_RANGE_PATH, range + 1);

        if (mf_io_mkdir(w_dir_range[range]) != CPL_ERROR_NONE) {
            err += cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                            "Could not create working directory");
        }

        /* Get wavenumber range in [cm-1] */
        wn_start[range] = cpl_table_get(params->rangetab, MF_COL_WN_START, range, NULL);
        wn_end[range]   = cpl_table_get(params->rangetab, MF_COL_WN_END,   range, NULL);
        if (wn_start[range] != 0. || wn_end[range] != 0.) {

            /* Update LNFL global time execution and show results info */
            cpl_msg_info(cpl_func, "(mf_lnfl      ) Configuring LNFL   for %.3f - %.3f micron (Range: %3d)", 1.e4 / wn_end[range], 1.e4 / wn_start[range], range + 1);

            /* Create a symbolic link to TAPE1 (initial AER file) and write the TAPE5 input file by LNFL execution */
            err += mf_io_write_lnfl_configuration(params->config->directories.telluriccorr_data_path, w_dir_range[range],
                                                  wn_start[range], wn_end[range],
                                                  params->config->internal.molecules.lbl_molecs, config);
        }
    }

    if (!err) {

        cpl_msg_info(cpl_func, "(mf_lnfl      ) Executing external calls to the LNFL binary ...");

        char *lnfl_syscall;
        if (params->config->inputs.silent_external_bins) {
            if (config->use_ODA) {
                lnfl_syscall = cpl_sprintf("%s/%s", params->config->directories.telluric_path, MF_BIN_LNFLODA_SILENT);
            } else {
                lnfl_syscall = cpl_sprintf("%s/%s", params->config->directories.telluric_path, MF_BIN_LNFL_SILENT);
            }
        } else {
            if (config->use_ODA) {
                lnfl_syscall = cpl_sprintf("%s/%s", params->config->directories.telluric_path, MF_BIN_LNFLODA);
            } else {
                lnfl_syscall = cpl_sprintf("%s/%s", params->config->directories.telluric_path, MF_BIN_LNFL);
            }
        }

        /* EXECUTE LNFL for EVERY different range numbers */
        double lnfl_runtime = 0.;
#ifdef _USE_OPENMP
        /* OMP parallel version */
        /*if (params->config->inputs.omp_num_threads > 1) {
            cpl_msg_info(cpl_func, "(mf_lnfl      ) OPENMP in parallel by wavelength ranges");
        } */
        #pragma omp parallel for shared(lnfl_syscall, wn_start, wn_end, w_dir_range) reduction(+ : err, lnfl_runtime)
#endif
        for (int range = 0; range < nrange; range++) {
            if (wn_start[range] != 0. || wn_end[range] != 0.) {
                double runtime;
                #ifdef _USE_OPENMP
                    err      += mf_io_systemOpenMP(lnfl_syscall, w_dir_range[range], &runtime);
                #else
                    err      += mf_io_system      (lnfl_syscall, w_dir_range[range], &runtime);
		#endif
		lnfl_runtime += runtime;
            }
        }

        mf_io_sync();

        cpl_free(lnfl_syscall);

        if (!err) {

            /* Update global time */
            params->timers.time_bin_lnfl += lnfl_runtime;
            cpl_msg_info(cpl_func, "(mf_lnfl      ) User-time execution of LNFL : %g min (TOTAL : %g min)", lnfl_runtime / 60., params->timers.time_bin_lnfl / 60.);

            /* Storage the TAPE3 filename */
            for (int range = 0; range < nrange; range++) {
                if (wn_start[range] != 0. || wn_end[range] != 0.) {
                    char *tape3 = cpl_sprintf("%s/%s", w_dir_range[range], MF_AER_TAPE3_FILE);
                    err += cpl_table_set_string(params->rangetab, MF_COL_LNFL, range, tape3);
                    cpl_free(tape3);
                }
            }
        }
    }

    /* Cleanup */
    for (int range = 0; range < nrange; range++) {
        cpl_free(w_dir_range[ range]);
    }
    cpl_free(w_dir_range);
    cpl_free(wn_start);
    cpl_free(wn_end);

    /* Check and return */
    if (err != CPL_ERROR_NONE) {
       cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                             "Error in subroutine: LNFL failed!");
    }
    return err;
}


/** @cond PRIVATE */


/** @endcond */


/**@}*/
