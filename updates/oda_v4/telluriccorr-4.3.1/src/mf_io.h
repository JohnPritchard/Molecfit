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

#ifndef MF_IO_H
#define MF_IO_H

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include <cpl.h>

#include "mf_constants.h"
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

#define MF_IO_MAX_NUMBER_FILES         4096
#define MF_IO_DEFAULT_CHUNK            262144          /* 256k */

#define MF_IO_BUFFER_SOCKET_LENGTH     256

#define MF_IO_TMP_FOLDER_ENV           "TMPDIR"                                 /* Generic environment variable                                    */
#define MF_IO_TMP_FOLDER_TEMPLATE      "telluriccorr_tmp_folder"                /* Generic name for the TMP folder                                 */
#define MF_IO_TMP_FOLDER_INIT          MF_IO_TMP_FOLDER_TEMPLATE"_XXXXXX"       /* The string "_XXXXXX" is part mandatory for mktemp(...)          */
#define MF_IO_TMP_FOLDER_MAX_ATTEMPTS  10                                       /* Number of attempts to create the TMP Telluric Correction folder */
#define MF_IO_ODA_PROFILE_DIR          "INDIVIDUAL_MOLECULE_PROFS"              /* Name of directory containing the optical depth profiles for each mol*/
#define ODA_MAX_MOL_IDX 50                                                      /* Fix the langth of the number of molecul indecies that have to be permutated*/

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

/******************************************************************************/
/* mf_lnfl: CONFIGURATION                                                     */
/******************************************************************************/

/* For overriding default values in oreder to run in ODA mode*/
typedef struct {
    char * work_dir;
    char * lnfl_wdir;
    char * lblrtm_wdir;
    char * tape3;
    char * lbl_molecs;
} mf_io_oda_parameters;

/******************************************************************************/
/* mf_lblrtm: CONFIGURATION                                                   */
/******************************************************************************/

/* Exposed structure : Initially default values for to be modify by the user */
typedef struct {
  int                        icntnm;                 /* Continua and Rayleigh extinction [0,1,2,3,4,5]              */
  int                        iaersl;                 /* Aerosols [0,1]                                              */
  int                        mpts;                   /* Number of optical depth values                              */
  int                        npts;                   /* Number of values for each panel                             */
  double                     v[2];                   /* Ending wavenumber value for the calculation                 */
  int                        sample;                 /* Number of sample points per mean halfwidth [between 1 to 4] */
  double                     alfal0;                 /* Average collision broadened halfwidth [cm-1/atm]            */
  double                     avmass;                 /* Average molecular mass [amu] for Doppler halfwidth          */
  double                     dptmin;                 /* Min molecular optical depth below lines will be rejected    */
  double                     dptfac;                 /* Factor multiplying molecular continuum optical depth        */
  double                     tbound;                 /* Temperature of boundary [K]                                 */
  double                     sremis[3];              /* Emissivity coefficients                                     */
  double                     srrefl[3];              /* Reflectivity coefficients                                   */
  int                        model;                  /* Atmospheric profile [0,1,2,3,4,5,6]                         */
  int                        itype;                  /* Type of path [1,2,3]                                        */
  int                        nozero;                 /* Zeroing of small amounts of absorbers [0,1]                 */
  int                        noprnt;                 /* Do not print output? [0,1]                                  */
  int                        ipunch;                 /* Write out layer data to TAPE7 [0,1]                         */
  double                     re;                     /* Radius of earth [km]                                        */
  double                     hspace;                 /* Altitude definition for space [km]                          */
  double                     ref_lat;                /* Latitude of location of calculation [degrees] [-90.-90]     */
  double                     h[2];                   /* h[0] = Observer altitude, h[1] = Upper height limit; [km]   */
  double                     range;                  /* Length of a straight path from H1 to H2 [km]                */
  double                     beta;                   /* Earth centered angle from H1 to H2 [degrees]                */
  int                        len;                    /* Path length [0,1]                                           */
  double                     hobs;                   /* Height of observer                                          */
  double                     avtrat;                 /* Maximum Voigt width ratio across a layer                    */
  double                     tdiff[2];               /* Maximum layer temperature difference at ALTD1, ALTD2 [K]    */
  double                     altd[2];                /* Altitude of TDIFF1, TDIFF2 [km]                             */
  double                     delv;                   /* Number of wavenumbers [cm-1] per major division             */
} mf_io_lblrtm_config;

/*----------------------------------------------------------------------------*/
/**
 *                 Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Get the absolute current working directory */
MF_EXPORT char * mf_io_pwd(void);

/* Get value from a environment variable, NULL if not exist */
MF_INTERNAL const char * mf_io_getenv(
    const char               *env);

/* Synchronize memory to disk */
MF_INTERNAL void mf_io_sync(void);

/* Check if the file exist in the disk */
MF_INTERNAL cpl_error_code mf_io_access(
    const char               *file);

/* Wrapper from system_calls */
MF_INTERNAL cpl_error_code mf_io_system(
    const char               *command,
    const char               *path,
    double                   *runtime);

/* Wrapper from system_calls used within OpenMP */
MF_INTERNAL cpl_error_code mf_io_systemOpenMP(
    const char               *command,
    const char               *path,
    double                   *runtime);

/* Remove one file from the disk */
MF_INTERNAL cpl_error_code mf_io_rm(
    const char               *file);

/* Move one file between paths */
MF_INTERNAL cpl_error_code mf_io_mv(
    const char               *source_path,
    const char               *dest_path);

/* Remove a path recursively without follow symbolic links */
MF_INTERNAL cpl_error_code mf_io_rm_rf(
    char                     *path,
    int                      max_dir_depth);

/* Create a temporary file in the disk with mkstemp */
MF_INTERNAL char * mf_io_mkstemp(
    const char               *config_tmp_path,
    const char               *default_tmp_path);

/* Create a user folder in disk, if not exist */
MF_INTERNAL cpl_error_code mf_io_mkdir(
    const char               *new_folder);

/* Substitute the call to the command system curl from download a file from the FTP and save in a local path */
MF_INTERNAL cpl_error_code mf_io_curl(
    const char               *ftp_host,
    const char               *url_path,
    const char               *local_dst);

/* Get number of entries in GDAS tarball */
MF_INTERNAL long mf_io_tarball_nfiles(
    const char               *path);

/* Read GDAS profile */
MF_INTERNAL cpl_table * mf_io_read_gdas_file_and_create_table(
    const char               *gdas_file_ASCII,
    const char               *hgt_units);

/* Create execute LNFL configuration */
MF_INTERNAL cpl_error_code mf_io_write_lnfl_configuration(
    const char               *data_path,
    const char               *w_dir,
    const double             wn_start,
    const double             wn_end,
    const char               *lbl_molecs,
    const mf_io_lnfl_config  *config);

/* Create execute LBLRTM configuration */
cpl_error_code mf_io_write_lblrtm_configuration(
    const char               *w_dir,
    const char               *tape3,
    const double             V1,
    const double             V2,
    const double             vbar,
    const double             angle,
    const cpl_boolean        emission_spec,
    const char               *lbl_molecs,
    const mf_io_lblrtm_config *config,
    const cpl_table          *prof);

/*  */
MF_INTERNAL cpl_array * mf_io_find_klim(
    const char               *w_dir_range,
    const char               *lblrtm_out_filename);

/* Rebins LBLRTM spectra (wrapper output) in wavelength units [mu m] (variable step size possible) */
cpl_error_code mf_io_read_lblrtm_and_update_spec(
    const cpl_size           nrow,
    double                   *lamv,
    double                   *fluxv,
    const char               *spectrum_filename,
    cpl_bivector             *bvec,
    const double             llim[2],
    cpl_boolean              *usampl,
    int                      *jmin,
    int                      *jmax,
    mf_io_lnfl_config        *lnfl_config);
void          mf_io_load_oda_table(int range, const char* lblrtm_out_filename);
double**      mf_io_read_oda_table(int range);
cpl_bivector* mf_io_read_lblrtm_spec(const char *spectrum_filename);
cpl_bivector* mf_io_mergeODTables(const int range, cpl_vector* mol_abuns,const char* lblrtm_out_filename);
cpl_vector*   mf_io_molecule_abundancies(mf_parameters* params, cpl_array* fitpar);
cpl_matrix*   mf_io_oda_tableDB(int range, int molecule, double* vec, int nrows, int nmols, int option);
void          mf_io_oda_init_tableDB(int range, int nmols, int nrows);
void          mf_io_oda_set_tableDB(int range, int molecule, double *vec, int nrows);
cpl_matrix*   mf_io_oda_get_tableDB(int range);
void          mf_io_oda_delete_tableDB(void);
cpl_array*    mf_io_klim_from_odatable(int range);
cpl_bivector* mf_io_merge_wavefiles(const char* w_dir,const char* filename);

void          mf_io_DeclareOpticalDepthTable(int nrange, char* molec_string);
cpl_array*    mf_io_molecstring2Names(char* molec_str);
cpl_boolean   mf_io_use_odatable(cpl_boolean use_ODA);
cpl_boolean   mf_io_use_stdlblrtm(cpl_boolean use_ODA);
cpl_boolean   mf_io_use_debug(void);
cpl_error_code  mf_oda_mol_idx(char** mol, int nmol);
void          mf_io_oda_symlink(char* target, char* destination);
CPL_END_DECLS


#endif /* MF_IO_H */
