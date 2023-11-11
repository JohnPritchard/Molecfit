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

#ifndef MF_CONSTANTS_H
#define MF_CONSTANTS_H

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

#if __GNUC__ >= 4 || __clang__
     #define MF_INTERNAL __attribute__((visibility("hidden" )))
     #define MF_EXPORT   __attribute__((visibility("default")))
#endif

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#if defined(_OPENMP) && !defined(__APPLE__)
  #define _USE_OPENMP
  #include <omp.h>
#else
  #undef _USE_OPENMP
#endif

#include <string.h>
#include <math.h>

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

#define MF_WAVELENGTH_MIN_MICRONS      0.28                                     /* Minimum expected wavelength in 'microns - mu m', warning emitted for smaller wavelength */
#define MF_WAVELENGTH_MAX_MICRONS      8600.                                    /* Maximum expected wavelength in 'microns - mu m', warning emitted for larger  wavelength */
#define MF_EXTRA_WN_COVERAGE           5                                        /* Extra wavenumber coverage at both sides of the model spectrum in 'cm^-1'                */
#define MF_OVERSAMPLING_FACTOR         5                                        /* Oversampling factor for model spectrum                                                  */

#define MF_FIT_N_POLYNOME_MIN          0                                        /* Minimum degree of the fitting polynome : Continuum / Wavelength */
#define MF_FIT_N_POLYNOME_MAX          10                                       /* Maximum degree of the fitting polynome : Continuum / Wavelength */

/*** LOCAL DIRECTORIES ***/

#define MF_BIN_PATH                    "bin"                                    /* molecif_third_party binaries folder                  */
#define MF_BINARIES_FLAGS              ">/dev/null 2>&1"                        /* Flags to SILENT the ouptut                           */

#define MF_LNFL_NAME                   "lnfl"                                   /* Binary executable of LNFL                            */
#define MF_BIN_LNFL                    MF_BIN_PATH"/"MF_LNFL_NAME               /* LNFL third_party binary                              */
#define MF_BIN_LNFL_SILENT             MF_BIN_LNFL" "MF_BINARIES_FLAGS          /* Command for exec LNFL third_party binary             */

#define MF_LNFLODA_NAME                "lnflOda"                                /* Binary executable of LNFL                            */
#define MF_BIN_LNFLODA                 MF_BIN_PATH"/"MF_LNFLODA_NAME            /* LNFL third_party binary                              */
#define MF_BIN_LNFLODA_SILENT          MF_BIN_LNFLODA" "MF_BINARIES_FLAGS       /* Command for exec LNFL third_party binary             */

#define MF_LBLRTM_NAME                 "lblrtm"                                 /* Binary executable of LBLRTM                          */
#define MF_BIN_LBLRTM                  MF_BIN_PATH"/"MF_LBLRTM_NAME             /* LBLRTM third_party binary                            */
#define MF_BIN_LBLRTM_SILENT           MF_BIN_LBLRTM" "MF_BINARIES_FLAGS        /* Command for exec LBLRTM third_party binary           */

#define MF_HITRAN_PATH                 "hitran"                                 /* Folder for HITRAN                                    */

#define MF_PROFILES_LIB_PATH           "profiles/lib"                           /*                                                      */
#define MF_PROFILES_MIPAS_PATH         "profiles/mipas"                         /*                                                      */
#define MF_PROFILES_GDAS_PATH          "profiles/gdas"                          /*                                                      */

#define MF_AER_WDIR_LNFL_RANGE_PATH    "wdir_"MF_LNFL_NAME"_range"              /* Working AER directory for LBLRTM and LNFL            */
#define MF_AER_WDIR_LBLRTM_RANGE_PATH  "wdir_"MF_LBLRTM_NAME"_range"            /* Working AER directory for LBLRTM and LNFL            */

#define MF_AER_TAPE1_FILE              "TAPE1"                                  /* AER TAPE file  1:                                    */
#define MF_AER_TAPE3_FILE              "TAPE3"                                  /* AER TAPE file  3:                                    */
#define MF_AER_TAPE5_FILE              "TAPE5"                                  /* AER TAPE file  5:                                    */
#define MF_AER_TAPE6_FILE              "TAPE6"                                  /* AER TAPE file  6:                                    */
#define MF_AER_TAPE7_FILE              "TAPE7"                                  /* AER TAPE file  6:                                    */
#define MF_AER_TAPE9_FILE              "TAPE9"                                  /* AER TAPE file  9:                                    */
#define MF_AER_TAPE10_FILE             "TAPE10"                                 /* AER TAPE file 10:                                    */
#define MF_AER_TAPE11_FILE             "TAPE11"                                 /* AER TAPE file 11:                                    */
#define MF_AER_TAPE12_FILE             "TAPE12"                                 /* AER TAPE file 12:                                    */
#define MF_AER_TAPE27_FILE             "TAPE27"                                 /* AER TAPE file 27:                                    */
#define MF_AER_TAPE28_FILE             "TAPE28"                                 /* AER TAPE file 28:                                    */


/*** GDAS ***/

#define MF_GDAS_FTP                    "ftp.eso.org"
#define MF_GDAS_FTP_ABSOLUTE           "eso-ftp02.hq.eso.org"


#define MF_GDAS_URL                    "/pub/dfs/pipelines/skytools/molecfit/gdas"
#define MF_GDAS_DOWNLOAD_TIMEOUT       10                                       /* Minutes                                              */
#define MF_GDAS_FILE_INVERVAL          3                                        /* Every GDAS file cover 3h                             */
#define MF_GDAS_FILE_SECURITY_SEARCH   2                                        /* Number of GDAS files search closest to the image     */

/**********************************/

#define MF_TOL                         1.e-7                                    /* Required relative accuracy for data comparisons      */

#define MF_LEN_MAX                     1024                                     /* Maximum number of string characters                  */

/* Conversions */
#define MF_CONV_T_2_K                  273.15                                   /* Convert from Celsius to Kelvin                       */
#define MF_CONV_K_2_T                  -MF_CONV_T_2_K                           /* Convert from Kelvin  to Celsius                      */
#define MF_CONV_K_LAM                  1e4                                      /* Convert wavelength to wavenumber: cm^{-1} to \mu m   */
#define MF_CONV_KM_2_M                 1e3                                      /* Convert kilometer to meter                           */
#define MF_CONV_M_2_KM                 1./MF_CONV_KM_2_M                        /* Convert kilometer to meter                           */
#define MF_CONV_MBAR_2_PA              1e2                                      /* Convert mili-bar to pascal                           */

/* Units */
#define MF_UNIT_VOL                    "ppmv"                                   /*  */
#define MF_UNIT_DIST                   "km"                                     /*  */
#define MF_UNIT_PRESS                  "mb"                                     /*  */
#define MF_UNIT_TEMP                   "K"                                      /*  */

/* Output table columns */
#define MF_COL_PARAMETER               "parameter"                              /*  */
#define MF_COL_VALUE                   "value"                                  /*  */
#define MF_COL_UNCERTAINTY             "uncertainty"                            /*  */


/*** Telluric Correction CPL table INPUTS ***/
#define MF_INPUT_MOLECULES             "MOLECULES"                              /*  */
#define MF_INPUT_WAVE_INCLUDE          "WAVE_INCLUDE"                           /*  */
#define MF_INPUT_WAVE_EXCLUDE          "WAVE_EXCLUDE"                           /*  */
#define MF_INPUT_PIXEL_EXCLUDE         "PIXEL_EXCLUDE"                          /*  */
#define MF_INPUT_GDAS                  "GDAS"                                   /*  */
#define MF_INPUT_ATM_PROFILE_STANDARD  "ATM_PROFILE_STANDARD"                   /*  */
#define MF_INPUT_ATM_PROFILE_COMBINED  "ATM_PROFILE_COMBINED"                   /*  */


/*** Telluric Correction CPL table Columns ***/

/* GDAS standard */
#define MF_COL_GDAS_RELHUM             "relhum"                                 /*  */
#define MF_COL_GDAS_HEIGHT             "height"                                 /*  */
#define MF_COL_GDAS_PRESS              "press"                                  /*  */
#define MF_COL_GDAS_TEMP               "temp"                                   /*  */

/* Atmospheric profile */
#define MF_COL_ATM_HGT                 "HGT"                                    /*  */
#define MF_COL_ATM_PRE                 "PRE"                                    /*  */
#define MF_COL_ATM_TEM                 "TEM"                                    /*  */

/* Input/Output columns */
#define MF_COL_IN_LAMBDA               "lambda"                                 /*  */
#define MF_COL_IN_FLUX                 "flux"                                   /*  */
#define MF_COL_IN_DFLUX                "dflux"                                  /*  */
#define MF_COL_IN_MASK                 "mask"                                   /*  */
#define MF_COL_OUT_FLUX                "cflux"                                  /*  */
#define MF_COL_OUT_DFLUX               "cdflux"                                 /*  */
#define MF_COL_OUT_MASK                "qual"                                   /*  */
#define MF_COL_OUT_TELLURIC_CORR       "mtrans"                                 /*  */

/* TAC columns */
#define MF_COL_TAC_FLUX                "tacflux"                                /*  */
#define MF_COL_TAC_DFLUX               "tacdflux"                               /*  */
#define MF_COL_TAC_MASK                "tacqual"                                /*  */
#define MF_COL_TAC_TELLURIC_CORR       "tacmtrans"                              /*  */

#define MF_COL_RANGE_FLUX              "oflux"                                  /*  */
#define MF_COL_RANGE_SCALE             "scal"                                   /*  */

#define MF_COL_MOD_RANGE               "mrange"                                 /*  */
#define MF_COL_MOD_LAMBDA              "mlambda"                                /*  */
#define MF_COL_MOD_FLUX                "mflux"                                  /*  */
#define MF_COL_MOD_SCALE               "mscal"                                  /*  */
#define MF_COL_MOD_WEIGHT              "mweight"                                /*  */

#define MF_COL_LAMBDA_0                "lambda0"                                /*  */
#define MF_COL_SLAMBDA                 "slambda"                                /*  */
#define MF_COL_WEIGHT                  "weight"                                 /*  */
#define MF_COL_DEV                     "dev"                                    /*  */

/* MOLEC_TAB */
#define MF_COL_LIST_MOLECULES          "LIST_MOLEC"                             /*  */
#define MF_COL_FIT_MOLECULES           "FIT_MOLEC"                              /*  */
#define MF_COL_REL_COL                 "REL_COL"                                /*  */
#define MF_COL_PPMV                    MF_UNIT_VOL                              /*  */

/* RANGE_TAB */
#define MF_COL_WAVE_RANGE_LOWER        "LOWER_LIMIT"                            /*  */
#define MF_COL_WAVE_RANGE_UPPER        "UPPER_LIMIT"                            /*  */
#define MF_COL_WAVE_RANGE_CONT_FIT     "CONT_FIT_FLAG"                          /*  */
#define MF_COL_WAVE_RANGE_CONT_ORDER   "CONT_POLY_ORDER"                        /*  */
#define MF_COL_WAVE_RANGE_MAP2CHIP     "MAPPED_TO_CHIP"                         /*  */
#define MF_COL_WAVE_RANGE_WLC_FIT      "WLC_FIT_FLAG"                           /*  */

/* RANGE_TAB */
#define MF_COL_CHIP                    "chip"                                   /*  */
#define MF_COL_FIT_RANGE               "fit_range"                              /*  */
#define MF_COL_CONT_RANGE              "cont_range"                             /*  */
#define MF_COL_CONT_COEF               "cont_coef"                              /*  */
#define MF_COL_PIX_RES                 "pixres"                                 /*  */
#define MF_COL_WN_START                "wn_start"                               /*  */
#define MF_COL_WN_END                  "wn_end"                                 /*  */
#define MF_COL_WN_STEP                 "wn_step"                                /*  */
#define MF_COL_LNFL                    "lnfl"                                   /*  */

/* CHIP_TAB */
#define MF_COL_FIT_CHIP                "fit_chip"                               /*  */
#define MF_COL_WLC_CHIP                "wlc_chip"                               /*  */
#define MF_COL_WLC_COEF                "wlc_coef"                               /*  */
#define MF_COL_WL_MIN                  "wl_min"                                 /*  */
#define MF_COL_WL_MAX                  "wl_max"                                 /*  */

/* KERNEL_TAB */
#define MF_COL_PIX_MIN                 "pixmin"                                 /*  */
#define MF_COL_PIX_MAX                 "pixmax"                                 /*  */
#define MF_COL_PIX_0                   "pix0"                                   /*  */
#define MF_COL_N_PIX                   "npix"                                   /*  */
#define MF_COL_T_N_PIX                 "tnpix"                                  /*  */
#define MF_COL_KERNEL                  "kernel"                                 /*  */
#define MF_COL_T_KERNEL                "tkernel"                                /*  */
#define MF_COL_KERNCEN                 "kerncen"                                /*  */

/* mf_model Output cpl_table parameters */
#define MF_MODEL_FIT_COEFFICIENT       "Chip 1, coef "                          /*  */
#define MF_MODEL_FIT_STATUS            "status"                                 /*  */
#define MF_MODEL_FIT_N_PAR             "fit_params"                             /*  */
#define MF_MODEL_FIT_DATA_POINTS       "data_points"                            /*  */
#define MF_MODEL_FIT_POS_WEIGHTS       "positive_weights"                       /*  */
#define MF_MODEL_FIT_PIX_FRAC          "valid_pix_frac"                         /*  */
#define MF_MODEL_FIT_N_ITER            "iterations"                             /*  */
#define MF_MODEL_FIT_FUNC_VAL          "func_eval"                              /*  */
#define MF_MODEL_FIT_LBLRTM_CALLS      "lblrtm_calls"                           /*  */
#define MF_MODEL_FIT_INIT_CHI2         "initial_chi2"                           /*  */
#define MF_MODEL_FIT_BEST_CHI2         "best_chi2"                              /*  */
#define MF_MODEL_FIT_REDUCED_CHI2      "reduced_chi2"                           /*  */
#define MF_MODEL_FIT_RMS_REL_ERR       "rms_rel_to_err"                         /*  */
#define MF_MODEL_FIT_RMS_REL_MEAN      "rms_rel_to_mean"                        /*  */
#define MF_MODEL_FIT_BOX_FWHM          "boxfwhm"                                /*  */
#define MF_MODEL_FIT_BOX_FWHM_PIX      "boxfwhm_pix"                            /*  */
#define MF_MODEL_FIT_GAUSS_FWHM        "gaussfwhm"                              /*  */
#define MF_MODEL_FIT_LORENTZFWHM       "lorentzfwhm"                            /*  */
#define MF_MODEL_FIT_TELBACK           "telback"                                /*  */
#define MF_MODEL_FIT_REL_MOL_COL_      "rel_mol_col_"                           /*  */
#define MF_MODEL_FIT_PPMV_     		   "ppmv_"                      			/*  */
#define MF_MODEL_FIT_H20_COL_MM        "h2o_col_mm"                             /*  */

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


#endif /* MF_CONSTANTS_H */
