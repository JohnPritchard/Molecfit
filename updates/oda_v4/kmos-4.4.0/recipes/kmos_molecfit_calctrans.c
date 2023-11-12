/* 
 * This file is part of the KMOS Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#include <regex.h>
#endif

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include "kmos_molecfit.h"

/*----------------------------------------------------------------------------*/
/**
 *                              Defines
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                 Typedefs: Structs and enum types
 */
/*----------------------------------------------------------------------------*/

typedef struct {
  kmos_molecfit_parameter mf;            /* Generic molecfit parameter                                                                           */
  cpl_propertylist                 *pl_atmos_params;     /* Primary property list header in the ATMOS_PARM file                                                  */
  cpl_table                        *atmprof[N_IFUS];     /* cpl_table in the ATMOS_PARM file    (one for each extension DATA executed in kmos_molecfit_model)    */
  cpl_propertylist                 *pl_best_fit_params;  /* Primary property list header in the BEST_FIT_PARM file                                               */
  cpl_table                        *res_table[N_IFUS];   /* cpl_table in the BEST_FIT_PARM file (one for each extension DATA executed in kmos_molecfit_model)    */
  mf_calctrans_convolution_results *results[N_IFUS];     /* Results of telluric corrections after execute mf_calctrans_convolution                               */
  const char					   *scale_pwv;
  double					  	   pwv_sci;
  double					  	   pwv_ratio;
  double                           h2o_col_mm;

} kmos_molecfit_calctrans_parameter;


/*----------------------------------------------------------------------------*/
/**
 *                              Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Fill the internal KMOS configuration file */
static kmos_molecfit_calctrans_parameter * kmos_molecfit_calctrans_conf(
    const cpl_parameterlist *list);

/* Use the kmos_molecfit_model output to configure kmos_molecfit_caltrans */
static cpl_error_code kmos_molecfit_calctrans_frames_conf(
    kmos_molecfit_calctrans_parameter *conf, cpl_frameset *frameset, const cpl_parameterlist *parlist);

/* Fill the molecfit specific recipe configuration file */
static cpl_error_code kmos_molecfit_calctrans_mf_conf(
    kmos_molecfit_calctrans_parameter *conf,
    double                            median,
    mf_parameters_config              *config_parameters);

/* Clean variables allocated in the recipe */
static void kmos_molecfit_calctrans_clean(
    kmos_molecfit_calctrans_parameter *conf);

/*----------------------------------------------------------------------------*/
/**
 *                          Static variables
 */
/*----------------------------------------------------------------------------*/

#define RECIPE_NAME      KMOS_MOLECFIT_CALCTRANS
#define CONTEXT          "kmos."RECIPE_NAME

static char kmos_molecfit_calctrans_description[] =
    "This recipe applies the results from kmos_molecfit_model and runs calctrans to calculate the \n"
    "telluric correction for scientific input data. Scientific input data can have category: \n"
    " - STAR_SPEC (24 DATA plus 24 NOISE extensions)\n"
    " - EXTRACT_SPEC (24 DATA extensions, additional 24 NOISE extensions are optional)\n"
    " - SCIENCE (24 DATA extensions, additional 24 NOISE extensions are optional)\n"
    " - SCI_RECONSTRUCTED (24 DATA extensions, additional 24 NOISE extensions are optional)\n"
    "It is not mandatory that all the DATA extension contains data.\n"
    "It accounts for the difference in airmass between the input model and the input scientific data. \n"
    "It accounts for different spectral resolutions between the various KMOS IFUs.\n"
    "\n"
    "Input files: (It's mandatory to provide 1 file only of type A) \n"
    "\n"
    "   DO                 KMOS                                                  \n"
    "   category           Type   required   Explanation                         \n"
    "   --------           -----  --------   -----------                         \n"
    "   STAR_SPEC          F1I       A       The science spectrum to compute the telluric correction for (1D spectrum, IMAGE format).\n"
    "   EXTRACT_SPEC       F1I       A       The science spectrum to compute the telluric correction for (1D spectrum, IMAGE format).\n"
    "   SCIENCE            F3I       A       The science spectrum to compute the telluric correction for (3D cube,     IMAGE format).\n"
    "   SCI_RECONSTRUCTED  F3I       A       The science spectrum to compute the telluric correction for (3D cube,     IMAGE format).\n"
    "   ATMOS_PARM         F1L       Y       Atmospheric model as computed by the recipe kmos_molecfit_model.\n"
    "   BEST_FIT_PARM      F1L       Y       Best fitting model and parameters as computed by the recipe kmos_molecfit_model.\n"
    "   KERNEL_LIBRARY     F2I       N       The kernel library; must be the same grating as the other inputs.\n"
    "\n"
    "Output files:                                                               \n"
    "\n"
    "   DO                 KMOS                                                  \n"
    "   category           Type              Explanation                         \n"
    "   --------           -----             -----------                         \n"
    "   TELLURIC_DATA      F1L               Telluric correction and intermediate products for each IFU. Output fits file with 48 ext. (24-data and 24-error, in TABLE format).\n"
    "                                           Each data extension contains the telluric correction for one IFU.\n"
    "   TELLURIC_CORR      F1I               Telluric correction for each IFU. Output fits file with 48 ext. (24-data and 24-error, in IMAGE format).\n"
    "                                           Each data extension contains the telluric correction for one IFU.\n"
    "\n"
    "----------------------------------------------------------------------------\n"
    "The input atmospheric model and best fitting parameters (as obtained in the kmos_molecfit_model) for a given IFU.X are used to generate a telluric correction for IFU.Y.\n"
    "The instrumental spectral resolution of IFU.Y is taken into account if a kernel library is provided.\n"
    "The mapping between IFU.X and IFU.Y is either specified by the user or automatic.\n"
    "The difference in airmass between the data used in kmos_molecfit_model and the input data of kmos_molecfit_calctrans are taken into account.\n"
    "\n";

/* Standard CPL recipe definition */
cpl_recipe_define(	kmos_molecfit_calctrans,
                  	KMOS_BINARY_VERSION,
                  	"Jose A. Escartin, Yves Jung",
                  	"https://support.eso.org/",
                  	"2017",
                  	"Read the results from kmos_molecfit_model and computes the telluric correction for a scientific input data.",
                  	kmos_molecfit_calctrans_description);


/*----------------------------------------------------------------------------*/
/**
 * @defgroup kmos_molecfit_calctrans  It runs molectift on KMOS standard star file to compute an atmospheric model.
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
 *                              Functions code
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief    Interpret the command line options and execute the data processing
 *
 * @param    frameset   the frames list
 * @param    parlist    the parameters list
 *
 * @return   CPL_ERROR_NONE if everything is OK or CPL_ERROR_CODE in other case
 *
 */
/*----------------------------------------------------------------------------*/
static int kmos_molecfit_calctrans(
    cpl_frameset *frameset, const cpl_parameterlist *parlist)
{
  /* Check initial Entries */
  if (kmos_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
      return cpl_error_get_code();
  }

  /* Get initial errorstate */
  cpl_errorstate initial_errorstate = cpl_errorstate_get();

  /* Extract and verify data in the parameters of the recipe */
  cpl_msg_info(cpl_func, "Configuring initial parameters in the recipe ...");
  kmos_molecfit_calctrans_parameter *conf = kmos_molecfit_calctrans_conf(parlist);
  cpl_error_ensure(conf, CPL_ERROR_ILLEGAL_INPUT,
                   return (int)CPL_ERROR_ILLEGAL_INPUT, "Problems with the configuration parameters");


  /* Complete configuration with the ATMOS_PARAM and BEST_FIT_PARM, outputs of the recipe kmos_molecfit_model */
  if (kmos_molecfit_calctrans_frames_conf(conf, frameset, parlist) != CPL_ERROR_NONE) {
      if (conf) kmos_molecfit_calctrans_clean(conf);
      return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                        "Configuration with kmos_calctrans_model output failed!");
  }

  if (!(conf->mf.fit_wavelenght.fit)) cpl_msg_info(cpl_func, "Not fit wavelenght!");
  if (!(conf->mf.fit_continuum.fit) ) cpl_msg_info(cpl_func, "Not fit continuum!" );

  /* Loading data spectrums */
  cpl_error_code err = kmos_molecfit_load_spectrums(frameset, &(conf->mf), RECIPE_NAME);
  if(err != CPL_ERROR_NONE) {
      kmos_molecfit_calctrans_clean(conf);
      return err;
  }

  /* Loading kernels */
  const cpl_frame *frmKernel = cpl_frameset_find(frameset, KERNEL_LIBRARY);
  cpl_error_code sKernels_e;
  if (frmKernel && conf->mf.use_input_kernel) {

      if (kmos_molecfit_load_kernels(frmKernel, &(conf->mf)) != CPL_ERROR_NONE) {
          kmos_molecfit_calctrans_clean(conf);
          return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                            "Loading convolution kernel data in the input frameset failed!");
      }

      sKernels_e = kmos_molecfit_save(frameset, frameset, parlist, RECIPE_NAME, conf->mf.parms, CALCTRANS_KERNEL_LIBRARY, conf->mf.grating.name, conf->mf.suppress_extension, NULL);
      if (sKernels_e != CPL_ERROR_NONE){
          kmos_molecfit_calctrans_clean(conf);
          return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                            "Saving generic multi-extension output fits file ('%s') failed!",
                                            CALCTRANS_KERNEL_LIBRARY);
      }

  } else {

      if (frmKernel) {
          cpl_msg_warning(cpl_func, "Kernel is provided, but use_input_kernel = false !");
      } else if (conf->mf.use_input_kernel){
          cpl_msg_warning(cpl_func, "Kernel isn't provided, but use_input_kernel = true !");
      }
      cpl_msg_info(cpl_func, "Using the default molecfit kernels -> With the values inside BEST_FIT_PARMS!");
  }

  cpl_msg_info(cpl_func, " +++ All input data loaded successfully! +++");


  /* Saving generic multi-extension output *.fits file */

  const char * telcorr_file = "telcorr.fits";

  cpl_msg_info(cpl_func, "Saving generic multi-extension output fits file ('%s','%s') ...", TELLURIC_DATA, TELLURIC);
  cpl_error_code sTel_data = kmos_molecfit_save(frameset, frameset, parlist, RECIPE_NAME, conf->mf.parms, TELLURIC_DATA, conf->mf.grating.name, conf->mf.suppress_extension, NULL);
  //cpl_error_code sTel_img  = kmos_molecfit_save(frameset, frameset, parlist, RECIPE_NAME, conf->mf.parms, TELLURIC_CORR, conf->mf.grating.name, conf->mf.suppress_extension, NULL);
  cpl_error_code sTel_img = cpl_propertylist_save(conf->mf.header_spectrums, telcorr_file , CPL_IO_CREATE);


  if (sTel_data != CPL_ERROR_NONE || sTel_img != CPL_ERROR_NONE) {
      kmos_molecfit_calctrans_clean(conf);
      return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                        "Saving generic multi-extension output fits files ('%s','%s') failed!", TELLURIC_DATA, TELLURIC_CORR);
  }


  /* Constant values in Molecfit executions */
  const double wlmin = -1.;
  const double wlmax = -1.;

  /* Get molecules in the grating */
  cpl_table *molecules = conf->mf.grating.molecules;

  /*** General configuration of Molecfit for LBLRTM and LNFL calls ***/

  /* Initialize mf_configuration structure */
  mf_configuration *mf_config = mf_configuration_create();
  if (!mf_config) err = CPL_ERROR_NULL_INPUT;

  /* Update the molecfit configuration with the data primary header */
  err = mf_parameters_config_update_with_header_keywords(mf_config->parameters, conf->mf.header_spectrums);
  if (err != CPL_ERROR_NONE) {
      cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "Error updating molecfit parameters with the raw primary header configuration");
  }

  /* Execute molecfit (lblrtm and convolution) */
  for (cpl_size n_ifuY = 0; n_ifuY < N_IFUS; n_ifuY++) {

      /* Get specific IFU */
      kmos_spectrum *ifuY = &(conf->mf.ifus[n_ifuY]);

      cpl_msg_info(cpl_func, "IFU NUMBER : %d", ifuY->num);

      /* Running Molecfit lblrtm only in IFU_Y (when ifuY->num is equal to ifuY->map) */
      if (ifuY->num == ifuY->map) {

          /* Configuration variables */
          cpl_msg_info(cpl_func, "Building molecfict configuration variable in %s ...", ifuY->name);

          /* Check previous executions in ifuY */
          if (conf->results[n_ifuY]) {
            err = cpl_error_set_message(cpl_func, CPL_ERROR_DUPLICATING_STREAM,
                                        "Unexpected error: Molecfit call mf_calctrans_lblrtm(...) already executed in this IFU");
          }

          /* Building molecfic configuration variable ifuY */
          if (!err) err =  kmos_molecfit_calctrans_mf_conf(conf, ifuY->median, mf_config->parameters);

          //mf_config->parameters->inputs.silent_external_bins = CPL_FALSE;
          //cpl_table_dump(conf->res_table[n_ifuY], 0, cpl_table_get_nrow(conf->res_table[n_ifuY]), stdout);

          /* CALL Molecfit MF_LBLRTM routine */
          if (!err) {

              /* Create input molecfit spec format */
              cpl_table *spec_telluriccorr_lblrtm = mf_spectrum_create(mf_config->parameters, ifuY->data);

              /* CALL MOLECFIT */
              cpl_msg_info(cpl_func, "Call mf_calctrans_lblrtm(...), in IFU_Y->map=%d),(num=%d,name=%s)", ifuY->map, ifuY->num, ifuY->name);
              mf_calctrans_lblrtm_results *lblrtm_results;
              lblrtm_results = mf_calctrans_lblrtm( mf_config,                  /* mf_configuration       *config            */
                                                    spec_telluriccorr_lblrtm,   /* const cpl_table        *spec_telluriccorr */
                                                    molecules,                  /* cpl_table              *molectab          */
                                                    wlmin,                      /* double                 wl_start           */
                                                    wlmax,                      /* double                 wl_end             */
                                                    conf->atmprof[n_ifuY],      /* cpl_table              *atmprof           */
                                                    conf->res_table[n_ifuY]);   /* cpl_table              *res_table         */

              cpl_table_delete(spec_telluriccorr_lblrtm);

              /* Check call results */
              err = cpl_error_get_code();
              if (err != CPL_ERROR_NONE) {
                  cpl_msg_info(cpl_func, "Molecfit call mf_calctrans_lblrtm(...) failed! in %s -> error: %s",  ifuY->name, cpl_error_get_message());
              } else if (!lblrtm_results) {
                  err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                              "Molecfit call (mf_run_calctrans_lblrtm) doesn't provide error but the results structure is NULL!");
              } else {
                  cpl_msg_info(cpl_func, "Molecfit call mf_calctrans_lblrtm(...) run successfully! in %s ...", ifuY->name);

                  /* Execute molecfit convolution: including ifuX == ifuY */
                  for (cpl_size n_ifuX = 0; n_ifuX < N_IFUS && !err; n_ifuX++) {

                      /* Get specifics ifus */
                      kmos_spectrum *ifuX = &(conf->mf.ifus[n_ifuX]);

                      /* Running Molecfit convolution in each IFU_X(ifu->num) match with the reference math with IFU_Y (when ifuX->map == ifu->num) */
                      if (ifuX->map == ifuY->num) {

                          /* Check previous executions in ifuX */
                          if (conf->results[n_ifuX]) {

                              err = cpl_error_set_message(cpl_func, CPL_ERROR_DUPLICATING_STREAM,
                                                          "Unexpected error: Molecfit call mf_calctrans_convolution(...) already executed in this IFU");
                          } else {

                              /* Building molecfic configuration variable */
                              cpl_msg_info(cpl_func, "Building molecfict configuration variable ... ");

                              /* Initialize config MF_PARAMETERS structure */
                              mf_parameters_config *config_parameters_ifuX = mf_parameters_config_create();
                              if (!config_parameters_ifuX) err = CPL_ERROR_NONE;

                              /* Update the molecfit configuration with the data primary header */
                              err = mf_parameters_config_update_with_header_keywords(config_parameters_ifuX, conf->mf.header_spectrums);
                              if (err != CPL_ERROR_NONE) {
                                  cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                                        "Error updating molecfit parameters with the raw primary header configuration");
                              }

                              /* Building molecfic configuration variable ifuX */
                              if (!err) err = kmos_molecfit_calctrans_mf_conf(conf, ifuY->median, config_parameters_ifuX);

                              if (!err) {

                                  /* Create input molecfit spec format */
                                  cpl_table *spec_telluriccorr_convolution = mf_spectrum_create(config_parameters_ifuX, ifuY->data);

                                  /* CALL MOLECFIT : Convolution of the model transmission spectrum with the kernel */
                                  cpl_msg_info(cpl_func, "Call mf_calctrans_convolution(...), in IFU_X->map=%d),(num=%d,name=%s) IFU_Y->map=%d),(num=%d,name=%s)",
                                               ifuX->map, ifuX->num, ifuX->name, ifuY->map, ifuY->num, ifuY->name);
                                  conf->results[n_ifuX] = mf_calctrans_convolution( config_parameters_ifuX,         /* mf_parameters_config        *config            */
                                                                                    mf_config->lnfl,
                                                                                    lblrtm_results,                 /* mf_calctrans_lblrtm_results *lblrtm_results    */
                                                                                    ifuY->header_ext_data,          /* const cpl_propertylist      *header_spec       */
                                                                                    spec_telluriccorr_convolution,  /* const cpl_table             *spec_telluriccorr */
                                                                                    ifuX->kernel.header_ext_data,   /* const cpl_propertylist      *header_kernel     */
                                                                                    ifuX->kernel.data,              /* const cpl_matrix            *kernel            */
                                                                                    wlmin,                          /* double                      wl_start           */
                                                                                    wlmax,                          /* double                      wl_end             */
                                                                                    conf->res_table[n_ifuY]);       /* cpl_table                   *res_table         */

                                  cpl_table_delete(spec_telluriccorr_convolution);

                                  /* Check call results */
                                  err = cpl_error_get_code();
                                  if (err != CPL_ERROR_NONE) {
                                      cpl_msg_info(cpl_func, "Molecfit call mf_calctrans_convolution(...) failed! in %s -> error: %s",  ifuX->name, cpl_error_get_message());
                                  } else if (!(conf->results[n_ifuX])) {
                                      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                                                  "Molecfit call (mf_run_calctrans_convolution) doesn't provide error but the results structure is NULL!");
                                  } else {
                                      cpl_msg_info(cpl_func, "Molecfit call mf_calctrans_convolution(...) run successfully! in %s ...", ifuX->name);
                                  }
                              }

                              /* Cleanup mf_calctrans_lnfl(...) call */
                              if (config_parameters_ifuX) mf_parameters_config_delete(config_parameters_ifuX);
                          }
                      }
                  }


                  /* Cleanup */
                  mf_calctrans_lblrtm_results_delete(lblrtm_results);
              }
          }
      }
  }

  /* Cleanup: General configuration of Molecfit for LBLRTM and LNFL calls */
  if (mf_config) mf_configuration_delete(mf_config);

  /* Check errors in Molecfit calls */
  if (err != CPL_ERROR_NONE) {
      kmos_molecfit_calctrans_clean(conf);
      return cpl_error_set_message(cpl_func, err,
                                   "Call to Molecfit (MF_LBLRTM/MF_LNFL) failed, err_str: %s",
                                   cpl_error_get_message());
  }

  /*QC Parameters to be used in kmos_combine*/
  cpl_propertylist * tel_primary_header = kmclipm_propertylist_load(telcorr_file, 0); //kmo_dfs_load_primary_header(frameset,TELLURIC_CORR);
  cpl_propertylist * best_fit_header = kmclipm_propertylist_load(
  		cpl_frame_get_filename(cpl_frameset_find(frameset, BEST_FIT_PARM)),0);

  if (cpl_propertylist_has(best_fit_header, "ESO TEL AIRM END")){

	  float airm_end = cpl_propertylist_get_double(best_fit_header, "ESO TEL AIRM END");
	  float airm_start = cpl_propertylist_get_double(best_fit_header, "ESO TEL AIRM START");
	  kmclipm_update_property_double(tel_primary_header, QC_PARAM_AIRM_STD,
			  (airm_end+airm_start)/2.0, QC_PARAM_AIRM_STD_TXT);
  }

  if (cpl_propertylist_has(best_fit_header, QC_PARAM_H2O_AVG)){

	  kmclipm_update_property_double(tel_primary_header, QC_PARAM_H2O_AVG,
			  cpl_propertylist_get_double(best_fit_header, QC_PARAM_H2O_AVG),
			  QC_PARAM_H2O_AVG_TXT);
  }

  if (cpl_propertylist_has(best_fit_header, QC_PARAM_H2O_RMS)){

	  kmclipm_update_property_double(tel_primary_header, QC_PARAM_H2O_RMS,
	      		  cpl_propertylist_get_double(best_fit_header, QC_PARAM_H2O_RMS),
	  			  QC_PARAM_H2O_RMS_TXT);
  }

  if (cpl_propertylist_has(best_fit_header, QC_PARAM_MEAS_IWV)){

	  float meas_iwv = cpl_propertylist_get_double(best_fit_header, QC_PARAM_MEAS_IWV);
	  kmclipm_update_property_double(tel_primary_header, QC_PARAM_MEAS_IWV, meas_iwv,
			  QC_PARAM_MEAS_IWV_TXT);
	  kmclipm_update_property_double(tel_primary_header, QC_PARAM_H2O_RATIO,
			  cpl_propertylist_get_double(best_fit_header, QC_PARAM_H2O_AVG)/meas_iwv,
			  QC_PARAM_H2O_RATIO_TXT);
  }

  if (cpl_propertylist_has(best_fit_header, QC_PARAM_MEAS_IWV30D)){
	  kmclipm_update_property_double(tel_primary_header, QC_PARAM_MEAS_IWV30D,
			  cpl_propertylist_get_double(best_fit_header, QC_PARAM_MEAS_IWV30D),
	  	  			  QC_PARAM_MEAS_IWV30D_TXT);
  }

  if (cpl_propertylist_has(best_fit_header, QC_PARAM_H2O_DELTA)){
	  kmclipm_update_property_double(tel_primary_header, QC_PARAM_H2O_DELTA,
			  cpl_propertylist_get_double(best_fit_header, QC_PARAM_H2O_DELTA),
	  			  QC_PARAM_H2O_DELTA_TXT);
  }

  /*QC Parameters to be used in kmos_combine --- MOVE TO NEW FUNCTION ---*/

  // CHANGE FOR SUPPRESS-EXTENSION
  const char * suffix= "";

  if (!conf->mf.suppress_extension) {
	  suffix = kmo_dfs_get_suffix(kmo_dfs_get_frame(frameset, TELLURIC_DATA), TRUE, FALSE);
  }

  cpl_frame *frame_telluric_corr = kmo_dfs_get_frame(frameset, TELLURIC_DATA);
  kmo_dfs_save_main_header(frameset, TELLURIC_CORR , suffix, frame_telluric_corr, tel_primary_header, parlist, cpl_func);


  /* Save data in each the IFU */
  for (cpl_size n_ifu = 0; n_ifu < N_IFUS; n_ifu++) {

      /* Get specifics ifus */
      kmos_spectrum *ifu = &(conf->mf.ifus[n_ifu]);

      cpl_table  *telluric_data = (conf->results[n_ifu])->spec_telluriccorr_format;
      double     *data          = NULL;
      cpl_vector *telluric_vec  = NULL;
      if (telluric_data) {

          /* Wrap the data */
          data = cpl_table_get_data_double(telluric_data, MF_COL_OUT_TELLURIC_CORR);
          cpl_vector *vAux = cpl_vector_wrap(cpl_table_get_nrow(telluric_data), data);
          telluric_vec = cpl_vector_duplicate(vAux);
          cpl_vector_unwrap(vAux);

          double mtrans_max = cpl_vector_get_max(telluric_vec);
          //cpl_vector_divide_scalar(telluric_vec, mtrans_max);
          cpl_msg_info(cpl_func, "Saving telluric data and image in IFU=%02d ... (Convolution executed) --> Max_value(%lf)", ifu->num, mtrans_max);

      } else {
          cpl_msg_info(cpl_func, "Saving data in IFU=%02d ...", ifu->num);
      }

      cpl_propertylist *header_data;
      cpl_propertylist *header_noise;
      if (ifu->num == ifu->map) {
          header_data  = cpl_propertylist_duplicate(ifu->header_ext_data);
          header_noise = cpl_propertylist_duplicate(ifu->header_ext_noise);
      } else {

          kmos_spectrum *ifu_ref = &(conf->mf.ifus[ifu->map - 1]);
          header_data  = cpl_propertylist_duplicate(ifu_ref->header_ext_data);
          header_noise = cpl_propertylist_duplicate(ifu_ref->header_ext_noise);

          char* name_data = cpl_sprintf("IFU.%d.DATA", ifu->num);
          cpl_propertylist_update_string(header_data, EXTNAME, name_data);
          cpl_free(name_data);

          char* name_noise = cpl_sprintf("IFU.%d.NOISE", ifu->num);
          cpl_propertylist_update_string(header_noise, EXTNAME, name_noise);
          cpl_free(name_noise);
      }

      if ((telluric_data== NULL) && (telluric_vec==NULL)) {
    	  kmos_all_clean_plist(header_data) ;
      }
	  kmos_all_clean_plist(header_noise) ;

      //cpl_propertylist_dump(header_data, stdout);

      /* Save cpl_table with the data of the mf_convolution execution */
      sTel_data = kmos_molecfit_save_mf_results(header_data, header_noise,
                                                TELLURIC_DATA, conf->mf.grating.name, conf->mf.suppress_extension, NULL, NULL, telluric_data, NULL);
      if (sTel_data != CPL_ERROR_NONE) {
          if (telluric_data) cpl_vector_delete(telluric_data);
          if (telluric_vec) cpl_vector_delete(telluric_vec);
          cpl_propertylist_delete(header_data);
          cpl_propertylist_delete(header_noise);
          kmos_molecfit_calctrans_clean(conf);
          return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                            "Saving Molecfit output fits file ('%s') failed!, ext=%lld", TELLURIC_DATA, n_ifu + 1);
      }


      /* Save data image with the data obtained in the mf_convolution execution */

      sTel_img = kmos_molecfit_save_mf_results(header_data, header_noise,
                                               TELLURIC_CORR, conf->mf.grating.name, conf->mf.suppress_extension, NULL, NULL, NULL, telluric_vec);
      if (sTel_img != CPL_ERROR_NONE) {
          if (telluric_vec) cpl_vector_delete(telluric_vec);
          cpl_propertylist_delete(header_data);
          cpl_propertylist_delete(header_noise);
          kmos_molecfit_calctrans_clean(conf);
          return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                            "Saving Molecfit output fits file ('%s') failed!, ext=%lld", TELLURIC_CORR, n_ifu + 1);
      }

      /* Save kernel used in Molecfit */
      if (frmKernel && conf->mf.use_input_kernel) {

          cpl_propertylist *header_kernel_data  = ifu->kernel.header_ext_data  ? ifu->kernel.header_ext_data  : ifu->header_ext_data;
          cpl_propertylist *header_kernel_noise = ifu->kernel.header_ext_noise ? ifu->kernel.header_ext_noise : ifu->header_ext_noise;

          cpl_matrix *kernel = NULL;
          if (conf->results[n_ifu]) {
              if ((conf->results[n_ifu])->kernel_resampled_normalized) kernel = (conf->results[n_ifu])->kernel_resampled_normalized;
          }

          sKernels_e = kmos_molecfit_save_mf_results(header_kernel_data, header_kernel_noise,
                                                     CALCTRANS_KERNEL_LIBRARY, conf->mf.grating.name, conf->mf.suppress_extension, NULL, kernel, NULL, NULL);

          if (sKernels_e != CPL_ERROR_NONE) {
              if (telluric_vec) cpl_vector_delete(telluric_vec);
              cpl_propertylist_delete(header_data);
              cpl_propertylist_delete(header_noise);
              kmos_molecfit_calctrans_clean(conf);
              err = cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                          "Saving Molecfit output fits file ('%s') failed!",
                                          CALCTRANS_KERNEL_LIBRARY);
          }
      }


      /* Cleanup */
      if (telluric_vec) cpl_vector_delete(telluric_vec);
      cpl_propertylist_delete(header_data);
      cpl_propertylist_delete(header_noise);
      //if (telluric_data) cpl_table_delete(telluric_data);
      unlink(telcorr_file);

  }


  if (conf->pwv_ratio!=1.0){

	  const cpl_frame *frmAtmosParm = cpl_frameset_find(frameset, ATMOS_PARM);
	  const char *fileAtmosParm = cpl_frame_get_filename(frmAtmosParm);

	  //cpl_msg_info(cpl_func, "We get to this at least. 0.25  - %s", cpl_error_get_message());
	  cpl_error_code sAtm_e =kmos_molecfit_save(frameset, frameset, parlist, RECIPE_NAME, conf->mf.parms, CALCTRANS_ATMOS_PARM,conf->mf.grating.name, conf->mf.suppress_extension, NULL);
	  //cpl_msg_info(cpl_func, "Error up to now -> %s ", cpl_error_get_message());

	  for (cpl_size n_ifu = 1; n_ifu <= N_IFUS; n_ifu++) {

    	  int j = 2*n_ifu-1;
          cpl_propertylist *header_atmos_data  = cpl_propertylist_load(fileAtmosParm,j);
          cpl_propertylist *header_atmos_noise = cpl_propertylist_load(fileAtmosParm, j+1);
          if (conf->atmprof[n_ifu-1]== NULL) {
        	  kmos_all_clean_plist(header_atmos_data) ;
          }
    	  kmos_all_clean_plist(header_atmos_noise) ;

          //cpl_table * atmos_tbl_noise = kmclipm_table_load( frmAtmosParm, j, 0 );

          //sKernels_e = kmos_molecfit_save_mf_results(header_atmos_data, header_atmos_noise,
          //	  CALCTRANS_ATMOS_PARM, conf->mf.grating.name, conf->mf.suppress_extension, NULL, conf->atmprof[n_ifu-1], NULL, NULL);
          sKernels_e = kmo_dfs_save_table(conf->atmprof[n_ifu-1], CALCTRANS_ATMOS_PARM, suffix ,header_atmos_data);
          //sKernels_e = kmo_dfs_save_table(atmos_tbl_noise, CALCTRANS_ATMOS_PARM, suffix ,header_atmos_noise);

          if (sKernels_e != CPL_ERROR_NONE) {
              cpl_propertylist_delete(header_atmos_data);
              cpl_propertylist_delete(header_atmos_noise);
              //cpl_table_delete(atmos_tbl_noise);
              kmos_molecfit_calctrans_clean(conf);
              err = cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                          "Saving Molecfit output fits file ('%s') failed!",
										  CALCTRANS_ATMOS_PARM);
          }

      }


  }



  /* Cleanup configuration */
  cpl_msg_info(cpl_func,"Cleaning variables ...");
  kmos_molecfit_calctrans_clean(conf);

  /* Check Recipe status and end */
  if (cpl_errorstate_is_equal(initial_errorstate) && cpl_error_get_code() == CPL_ERROR_NONE ) {
      cpl_msg_info(cpl_func,"Recipe successfully!");
  } else {
      /* Dump the error history since recipe execution start. At this point the recipe cannot recover from the error */
      cpl_errorstate_dump(initial_errorstate, CPL_FALSE, NULL);
      cpl_msg_info(cpl_func,"Recipe failed!, error(%d)=%s", cpl_error_get_code(), cpl_error_get_message());
  }

  return (int)cpl_error_get_code();
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Function needed by cpl_recipe_define to fill the input parameters
 *
 * @param  self   parameterlist where you need put parameters
 *
 * @return cpl_error_code
 *
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code kmos_molecfit_calctrans_fill_parameterlist(
    cpl_parameterlist *self)
{
  /* Add the different default parameters to the recipe */
  cpl_errorstate prestate = cpl_errorstate_get();

  /* Fill the parameters list */
  cpl_error_code e;
  cpl_boolean    range     = CPL_TRUE;
  const void     *dummyMin = NULL;
  const void     *dummyMax = NULL;


  /* --KMOS_MOLECFIT_PARAMETER_USE_INPUT_KERNEL */
  cpl_boolean use_input_kernel = CPL_TRUE;
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_USE_INPUT_KERNEL,
                                   !range, dummyMin, dummyMax, CPL_TYPE_BOOL, &use_input_kernel,
                                   "In order to use kernel library provide by the user.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;


  /*** Mapping IFUS ***/
  int automatic = -1.;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_1 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_1,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 1 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_2 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_2,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 2 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_3 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_3,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 3 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_4 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_4,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 4 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_5 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_5,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 5 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_6 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_6,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 6 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_7 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_7,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 7 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_8 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_8,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 8 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_9 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_9,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 9 .", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_10 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_10,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 10.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_11 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_11,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 11.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_12 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_12,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 12.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_13 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_13,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 13.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_14 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_14,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 14.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_15 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_15,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 15.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_16 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_16,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 16.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_17 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_17,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 17.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_18 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_18,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 18.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_19 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_19,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 19.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_20 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_20,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 20.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_21 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_21,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 21.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_22 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_22,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 22.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_23 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_23,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 23.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_IFU_24 */
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_IFU_24,
                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, &automatic,
                                   "IFU number in ATMOS_PARM and BEST_FIT_PARM to be used to compute the telluric correction for IFU 24.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_SUPPRESS_EXTENSION */
  cpl_boolean suppress_extension = CPL_FALSE;
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_SUPPRESS_EXTENSION,
                                   !range, dummyMin, dummyMax, CPL_TYPE_BOOL, &suppress_extension,
                                   "Suppress arbitrary filename extension.(TRUE (apply) or FALSE (don't apply).", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* --KMOS_MOLECFIT_PARAMETER_SCALE_PWV  - PIPE-10280 */
  const char *scale_pwv = "auto";
  e = kmos_molecfit_fill_parameter(RECIPE_NAME, self, KMOS_MOLECFIT_PARAMETER_SCALE_PWV,
                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, (const void *)scale_pwv,
                                   "Value to use when scaling the precipitable water vapor.", CPL_FALSE);
  if (e != CPL_ERROR_NONE) return (int)e;

  /* Check possible errors */
  if (!cpl_errorstate_is_equal(prestate)) {
      return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                   "kmos_molecfit_calctrans_fill_parameterlist failed!");
  }

  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief  Generate the internal configuration file for the recipe and check values
 *
 * @param  list   parameterlist with the parameters
 *
 * @return configuration file or NULL if exist an error
 *
 */
/*----------------------------------------------------------------------------*/
static kmos_molecfit_calctrans_parameter * kmos_molecfit_calctrans_conf(
    const cpl_parameterlist *list)
{
  /* Check input */
  cpl_error_ensure(list, CPL_ERROR_NULL_INPUT,
                   return NULL, "kmos_molecfit_calctrans_fill_conf input list NULL!");

  /* Get preState */
  cpl_errorstate preState = cpl_errorstate_get();
  const cpl_parameter *p;


  /* Create the configuration parameter */
  kmos_molecfit_calctrans_parameter *conf = (kmos_molecfit_calctrans_parameter *)cpl_malloc(sizeof(kmos_molecfit_calctrans_parameter));
  kmos_molecfit_nullify(&(conf->mf));
  conf->pl_atmos_params    = NULL;
  conf->pl_best_fit_params = NULL;
  for (cpl_size i = 0; i < N_IFUS; i ++) {
      conf->atmprof[i]     = NULL;
      conf->res_table[i]   = NULL;
      conf->results[i]     = NULL;
  }

  conf->pwv_ratio = 1.0;
  conf->pwv_sci = -99.9;
  conf->h2o_col_mm  = -99.9;

  /* Initialize input parameters propertylist */
  conf->mf.parms = cpl_propertylist_new();

  /* User input kernel ? */
  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_USE_INPUT_KERNEL);
  conf->mf.use_input_kernel = cpl_parameter_get_bool(p);

  /* Mapping IFUS in kmos_molecfit_calctrans recipe */
  int ifu = -1;

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_1);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_2);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_3);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_4);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_5);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_6);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_7);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_8);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_9);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_10);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_11);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_12);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_13);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_14);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_15);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_16);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_17);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_18);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_19);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_20);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_21);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_22);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_23);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_IFU_24);
  conf->mf.ifus[++ifu].map = cpl_parameter_get_int(p);

  /* suppress_extension */
  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_SUPPRESS_EXTENSION);
  conf->mf.suppress_extension = cpl_parameter_get_bool(p);

  /* pwv scaling */
  p = cpl_parameterlist_find_const(list, KMOS_MOLECFIT_PARAMETER_SCALE_PWV);
  conf->scale_pwv = cpl_parameter_get_string(p);

  /* Save parameter in the output propertylist */
  cpl_propertylist_update_bool(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_USE_INPUT_KERNEL,   conf->mf.use_input_kernel);

  ifu = -1;
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_1,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_2,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_3,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_4,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_5,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_6,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_7,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_8,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_9,              conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_10,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_11,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_12,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_13,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_14,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_15,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_16,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_17,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_18,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_19,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_20,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_21,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_22,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_23,             conf->mf.ifus[++ifu].map );
  cpl_propertylist_update_int( conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_IFU_24,             conf->mf.ifus[++ifu].map );

  cpl_propertylist_update_bool(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_SUPPRESS_EXTENSION, conf->mf.suppress_extension);
  cpl_propertylist_update_char(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_SCALE_PWV, conf->scale_pwv);

  cpl_msg_info(cpl_func, "SCALE_PWV VALUE-  '%s'.", conf->scale_pwv);

  // TODO: DOUBLE CHECK IF THESE KEYWORDS ARE REQUIRED OR IF OTHER QC PARAMS ARE NEEDED.
  /*cpl_propertylist_update_char(conf->pwv_ratio, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_SCALE_PWV, conf->pwv_ratio);
  cpl_propertylist_update_char(conf->pwv_sci, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_SCALE_PWV, conf->pwv_sci);*/




  /* Check status */
  if (!cpl_errorstate_is_equal(preState)) {
      /* Configuration failed */
      kmos_molecfit_calctrans_clean(conf);
      return NULL;
  } else {
      /* Configuration successfully */
      return conf;
  }
}

/*----------------------------------------------------------------------------*/
/**
 * @brief  Generate the internal configuration by the input frames
 *
 * @param  conf       configuration file in the recipe
 * @param  frameset   the frames list
 *
 * @return cpl_error_code
 *
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code kmos_molecfit_calctrans_frames_conf(
    kmos_molecfit_calctrans_parameter *conf, cpl_frameset *frameset, const cpl_parameterlist *parlist)
{
  /* Check input */
  cpl_error_ensure(conf && frameset, CPL_ERROR_NULL_INPUT,
                   return CPL_ERROR_NULL_INPUT, "kmos_molecfit_calctrans_recipe_model_conf inputs NULL!");

  /* Get preState */
  cpl_errorstate initialState = cpl_errorstate_get();

  regex_t regex;
  int reti;
  reti = regcomp(&regex,"[a-zA-Z0-9]",0);
  if(reti){
    cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,"molecfit_calctrans: could not compile regex required to check SCALE_PWV");
    return cpl_error_get_code();
  }
  reti = regexec(&regex, conf->scale_pwv,0,NULL,0);

  /*** Get frame, header and check: ATMOS_PARAM ***/
  cpl_msg_info(cpl_func, "Loading '%s' header, input from kmos_molecfit_modeL recipe ...", ATMOS_PARM);
  const cpl_frame *frmAtmosParm = cpl_frameset_find(frameset, ATMOS_PARM);
  cpl_error_ensure(frmAtmosParm, CPL_ERROR_DATA_NOT_FOUND,
                   return CPL_ERROR_DATA_NOT_FOUND, ATMOS_PARM" not found in input frameset!");
  const char *fileAtmosParm = cpl_frame_get_filename(frmAtmosParm);
  conf->pl_atmos_params = cpl_propertylist_load(fileAtmosParm, 0);
  cpl_error_ensure(conf->pl_atmos_params, cpl_error_get_code(),
                   return cpl_error_get_code(), "Cannot load ATMOS_PARM primary header propertylist from '%s'!", fileAtmosParm);


  /*** Get frame, header and check: BEST_FIT_PARM ***/
  cpl_msg_info(cpl_func, "Loading '%s' header, input from kmos_molecfit_model recipe ...", BEST_FIT_PARM);
  const cpl_frame *frmBestFitParm = cpl_frameset_find(frameset, BEST_FIT_PARM);
  cpl_error_ensure(frmBestFitParm, CPL_ERROR_DATA_NOT_FOUND,
                   return CPL_ERROR_DATA_NOT_FOUND, BEST_FIT_PARM" not found in input frameset!");
  const char *fileBestFitParm = cpl_frame_get_filename(frmBestFitParm);
  conf->pl_best_fit_params = cpl_propertylist_load(fileBestFitParm, 0);
  cpl_error_ensure(conf->pl_best_fit_params, cpl_error_get_code(),
                   return cpl_error_get_code(), "Cannot load BEST_FIT_PARM primary header propertylist from '%s'!", fileBestFitParm);


  /*** Get frame, header and check: SCIENCE ***/
  cpl_msg_info(cpl_func, "Loading the SCIENCE header for keywords ");
  const cpl_frame *frmSci = kmos_molecfit_get_frame_spectrums(frameset);
  cpl_error_ensure(frmSci, CPL_ERROR_DATA_NOT_FOUND,
                   return CPL_ERROR_DATA_NOT_FOUND, "Input Data not found in input frameset!");
  const char *fileScience = cpl_frame_get_filename(frmSci);
  const cpl_propertylist *pl_sci = cpl_propertylist_load(fileScience, 0);
  cpl_error_ensure(pl_sci , cpl_error_get_code(),
                   return cpl_error_get_code(), "Cannot load primary Data header propertylist from '%s'!", fileScience);


  /*Get MJD_OBS values from BEST_FIT_PARAMETERS and SCIENCE files for QC param ESO DRS PWV DELTA MJD */
  int raise_err = !strcmp(conf->scale_pwv,"none");

  /*if(!cpl_propertylist_has(conf->pl_best_fit_params, KEY_MJDOBS)){
      if(raise_err)
        cpl_msg_error(cpl_func,"Missing %s keyword from %s.\n\tCannot populate ESO DRS PWV DELTA MJD.\n", KEY_MJDOBS,BEST_FIT_PARM);
      else
          cpl_msg_info(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",KEY_MJDOBS);
  }

  if(!cpl_propertylist_has(pl_sci, KEY_MJDOBS)){
      if(raise_err)
        cpl_msg_error(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",KEY_MJDOBS);
      else
        cpl_msg_info(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",KEY_MJDOBS);
  }*/

  /*key_exp*/
  /*if(!cpl_propertylist_has(conf->pl_best_fit_params, KEY_EXPTIME)){
      if(raise_err)
        cpl_msg_error(cpl_func,"Missing %s keyword from %s.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",KEY_EXPTIME, BEST_FIT_PARM);
      else
        cpl_msg_info(cpl_func,"Missing %s keyword from %s.\n\tCannot populate ESO DRS PWV DELTA MJD.\n", KEY_EXPTIME, BEST_FIT_PARM);
  }
  if(!cpl_propertylist_has(pl_sci,KEY_EXPTIME)){
      if(raise_err)
        cpl_msg_error(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n", KEY_EXPTIME);
      else
        cpl_msg_info(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n", KEY_EXPTIME);
  } */

  /* This new QC Keyword is included in the Molecfit PWV calcultaions - check with Lodo what QC parameters are needed for this.
   *
     if(cpl_propertylist_has(bfp_hdr,KEY_EXPTIME) && cpl_propertylist_has(bfp_hdr,KEY_MJDOBS) &&
     cpl_propertylist_has(data_hdr,KEY_EXPTIME) && cpl_propertylist_has(data_hdr,KEY_MJDOBS)){
      double t0 = cpl_propertylist_get_double(bfp_hdr, KEY_MJDOBS);
      cpl_type exp_type = cpl_propertylist_get_type(bfp_hdr,KEY_EXPTIME);
      if(exp_type == CPL_TYPE_INT){
          t0 = t0 + 0.5*cpl_propertylist_get_int(bfp_hdr,KEY_EXPTIME);
      } else if (exp_type == CPL_TYPE_FLOAT){
          t0 = t0 + 0.5*cpl_propertylist_get_float(bfp_hdr,KEY_EXPTIME);
      } else if (exp_type == CPL_TYPE_DOUBLE){
          t0 = t0 + 0.5*cpl_propertylist_get_double(bfp_hdr,KEY_EXPTIME);
      }
      exp_type = cpl_propertylist_get_type(data_hdr,KEY_EXPTIME);
      double t1 = cpl_propertylist_get_double(data_hdr, KEY_MJDOBS);
      if(exp_type == CPL_TYPE_INT){
          t1 = t1 + 0.5*cpl_propertylist_get_int(data_hdr,KEY_EXPTIME);
      } else if (exp_type == CPL_TYPE_FLOAT){
          t1 = t1 + 0.5*cpl_propertylist_get_float(data_hdr,KEY_EXPTIME);
      } else if (exp_type == CPL_TYPE_DOUBLE){
          t1 = t1 + 0.5*cpl_propertylist_get_double(data_hdr,KEY_EXPTIME);
      }

      cpl_propertylist_append_double(qc_new,"ESO DRS PWV DELTA MJD",t1-t0);
      cpl_propertylist_set_comment(qc_new,"ESO DRS PWV DELTA MJD","Sci-Telluric Mean MJD");
  } */

  /*** Loading data form the input fits file ***/
  cpl_errorstate preState = cpl_errorstate_get();

  //PIPE-10280 -> correcting the h2o column in the atm profile can be done here.

  cpl_msg_info(cpl_func, "SCALE_PWV value -  %s.", conf->scale_pwv);

  /*  // END of getting scale_pwv/ pwv_sci */
  if (!strcmp(conf->scale_pwv,"auto")){
      if(cpl_propertylist_has(pl_sci,"ESO TEL AMBI IWV START") && cpl_propertylist_has(pl_sci,"ESO TEL AMBI IWV END")){
          double iwv_start = cpl_propertylist_get_double(pl_sci,"ESO TEL AMBI IWV START");
          double iwv_end = cpl_propertylist_get_double(pl_sci,"ESO TEL AMBI IWV END");
          if(iwv_start == 0.0 && iwv_end == 0.0){
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_calctrans: The ESO TEL AMBI IWV START/END science header keywords are both 0.0 which is not permitted by SCALE_PWV=auto");
            return cpl_error_get_code();
          }
          conf->pwv_sci = 0.5*(iwv_start+iwv_end);
      } else {
          /* IWV data not in headers, raise error */
          cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_calctrans: Missing ESO TEL AMBI IWV START/END science header keywords required by SCALE_PWV=auto. Deactivating transmission function scaling. ");
          conf->scale_pwv = "none";
          //return cpl_error_get_code();
      }
  }
  //if the supplied value is a header keyword (i.e. it has characters in it) and not equal to none
  else if (!reti && strcmp(conf->scale_pwv,"none")){// && parameters->scale_pwv has chars and != NONE) {
      if(cpl_propertylist_has(pl_sci,conf->scale_pwv)){
          conf->pwv_sci = cpl_propertylist_get_double(pl_sci,conf->scale_pwv);
      } else {
          cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_calctrans: Missing science header keyword %s required by SCALE_PWV",  conf->scale_pwv);
          return cpl_error_get_code();
      }
  }
  else if (reti && strcmp(conf->scale_pwv,"none")) {
      /* numerical value supplied, convert to float and use it */
       conf->pwv_sci = atof(conf->scale_pwv);
  } else if (strcmp(conf->scale_pwv,"none")){
          cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_calctrans: Error converting supplied value of SCALE_PWV: %s",conf->scale_pwv);
          return cpl_error_get_code();
  }

  // END of getting scale_pwv/ pwv_sci

  for (cpl_size n_ifu = 0; n_ifu < N_IFUS; n_ifu++) {

      /* Get number of extension */
      cpl_size ext = (n_ifu * 2) + 1;

      /* Get Best fit parameters */
      conf->res_table[n_ifu] = cpl_table_load(fileBestFitParm, ext, 0);
      if (!(conf->res_table[n_ifu])) {
          /* The extension doesn't have data spectrum */
          cpl_errorstate_set(preState);
      } else {
          cpl_msg_info(cpl_func, "%s, ext=%02lld: Loaded input data to kmos_molecfit_calctrans recipe.", BEST_FIT_PARM, ext);
      }

      /* Get Atmospheric profile */
      conf->atmprof[n_ifu] = cpl_table_load(fileAtmosParm, ext, 0);
      if (!(conf->atmprof[n_ifu])) {
          /* The extension doesn't have atmprof */
          cpl_errorstate_set(preState);
      } else {
          cpl_msg_info(cpl_func, "%s,    ext=%02lld: Loaded input data to kmos_molecfit_calctrans recipe.", ATMOS_PARM, ext);


          if (!(conf->res_table[n_ifu])) {
              /* The extension doesn't have data spectrum */
              cpl_errorstate_set(preState);
              cpl_msg_info(cpl_func, "%s, ext=%02lld Does not exist- cannot scale PWV transmission lines.", BEST_FIT_PARM, ext);
          } else {
              cpl_msg_info(cpl_func, "ext=%02lld Exists: Transmission Scaling starting - nScale Value -%d.", ext, (int)strcmp(conf->scale_pwv,"auto"));

              int ridx;
              for (ridx=0; ridx<cpl_table_get_nrow(conf->res_table[n_ifu]) ; ridx++){
            	  if(!strcmp(cpl_table_get_string(conf->res_table[n_ifu],"parameter",ridx),"h2o_col_mm")){
            		  conf->h2o_col_mm = cpl_table_get_double(conf->res_table[n_ifu],"value",ridx,NULL);
            	  }
              }


              /* ------------------  PWV SCALING APPLIED - START ------------------*/
              /*Calculate pwv_sci */
              if ((!strcmp(conf->scale_pwv,"auto")) || (!reti && strcmp(conf->scale_pwv,"none")) || (reti && strcmp(conf->scale_pwv,"none"))  ){
                  cpl_msg_info(cpl_func, "internal calc starting now, pwv_sci-> %f, h2o_col_value-> %f",conf->pwv_sci ,conf->h2o_col_mm );
                      if(conf->h2o_col_mm != -99.9 && conf->h2o_col_mm != 0.0){

                          conf->pwv_ratio = conf->pwv_sci / conf->h2o_col_mm ;
                          cpl_msg_info(cpl_func, "internal calc done");
                      }
              }

              cpl_msg_info(cpl_func,"Value of parameters->pwv_sci: %e; parameters->h2o_col_mm (pwv_tell): %e; parameters->pwv_ratio: %e; parameters->scale_pwv: %s\n",conf->pwv_sci, conf->h2o_col_mm, conf->pwv_ratio, conf->scale_pwv);
              /* Rescale values in parameters->atm_parameters and write out to CALCTRANS_ATM_PARAMETERS.fits */
              if(conf->pwv_ratio == 1.0){
                    cpl_msg_warning(cpl_func,"molecfit_calctrans: pwv_ratio is 1.0. No scaling based on PWV performed.\n");
              }

              cpl_table_multiply_scalar(conf->atmprof[n_ifu], MF_MOLECULES_H2O , conf->pwv_ratio);

              /* ------------------  PWV SCALING APPLIED - END ------------------*/

          }

      }
      //atmprof[n_ifu] updating H2O column for pwv scaling;

  }

  cpl_msg_info(cpl_func, "ENDING PWV SCALING SECTIONS");
  cpl_msg_info(cpl_func, "ERROR STATE -  %s.", cpl_error_get_message());

  //cpl_frame_delete(frmBestFitParm) ;
  //cpl_frame_delete(frmAtmosParm);
  //cpl_frame_delete(frmSci);
  cpl_propertylist_delete(pl_sci);

  regfree(&regex);

  /*** Check and assign mapping of IFUS ***/
  cpl_msg_info(cpl_func, "Mapping IFU's ... (Using BEST_FIT_PARM data in automatic assignment)");
  for (cpl_size n_ifu = 0; n_ifu < N_IFUS; n_ifu++) {

      kmos_spectrum *ifu    = &(conf->mf.ifus[n_ifu]);
      cpl_size      ifu_num = n_ifu + 1;

      /* Check if it's necessary to re-mapping automatically */
      if (ifu->map != -1) {

          if (!(conf->res_table[ifu->map - 1])) {
              return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                           "Mapping IFU_X.%02lld->IFU_Y.%02d by the user failed! It doesn't contain IFU_Y data in the file %s.fits",
                                           ifu_num, ifu->map, BEST_FIT_PARM);
          }

          cpl_msg_info(cpl_func, "IFU_X.%02lld mapping to IFU_Y.%02d --> "
                       "Defined by the user with the input parameter (IFU_%02lld=%02d)",
                       ifu_num, ifu->map, ifu_num, ifu->map);
      } else {

          /* First rule: ifuX_i (science) have ifuY_i (calibration) */
          if (conf->res_table[n_ifu]) {
              ifu->map = ifu_num;
              cpl_msg_info(cpl_func, "IFU_X.%02lld mapping to IFU_Y.%02d --> "
                           "Defined automatically (Same  IFU_Y [with calibration data] for                  IFU_X [science data] - 1st rule)",
                           ifu_num, ifu->map);
          }

          /* Second rule: ifuX_(1-8)-->ifu_Y(1-8), ifuX_(9-16)-->ifu_Y(9-16), ifuX_(17-24)-->ifu_Y(17-24) */
          if (ifu->map == -1) {
              cpl_size ifuY_detector      = (n_ifu / 8) + 1;           /* Get detector of Science (IFU_X) */
              cpl_size ifuY_detector_init = (ifuY_detector - 1) * 8;   /* Starting IFU_Y (calibration) in the detector of IFU_X number (science) [valid = 0, 8, 16] */
              for (cpl_size ifu_Y = ifuY_detector_init; ifu_Y < ifuY_detector_init + 8 && ifu->map == -1; ifu_Y++) {

                  /* If exist data --> Map */
                  if (conf->res_table[ifu_Y]) {
                      ifu->map = ifu_Y + 1;
                      cpl_msg_info(cpl_func, "IFU_X.%02lld mapping to IFU_Y.%02d --> "
                                   "Defined automatically (First IFU_Y [with calibration data] inside detector:%lld of IFU_X [science data] - 2nd rule)",
                                   ifu_num, ifu->map, ifuY_detector);
                  }
              }
          }

          /* Third rule: If not found it in the normal range, looking for in all the range */
          if (ifu->map == -1) {
              cpl_size ifuX_detector = (n_ifu / 8) + 1;                /* Get detector of Science (IFU_X) */
              for (cpl_size ifu_Y = 0; ifu_Y < N_IFUS && ifu->map == -1; ifu_Y++) {

                  /* If exist data --> Map */
                  if (conf->res_table[ifu_Y]) {
                      cpl_size ifuY_detector = (ifu_Y / 8) + 1;        /* Get detector of calibration (IFU_Y) */
                      ifu->map = ifu_Y + 1;
                      cpl_msg_info(cpl_func, "IFU_X.%02lld mapping to IFU_Y.%02d --> "
                                   "Defined automatically (First IFU_Y [with calibration date] in the instrument (detector=%lld) outside the IFU_X detector=%lld [science data] - 3rd rule)", ifu_num, ifu->map,
                                   ifuY_detector, ifuX_detector);
                  }
              }
          }

          cpl_error_ensure(ifu->map != -1, CPL_ERROR_ILLEGAL_INPUT,
                           return CPL_ERROR_ILLEGAL_INPUT, "Cannot search any IFU_Y data for the IFU_X");
      }
  }


  /*** Get the molecfit parameters from the BEST_FIL_PARAMS propertylist ***/
  cpl_msg_info(cpl_func, "Loading input parameters provided from kmos_molecfit_model recipe ...");

  conf->mf.grating.name		        = cpl_propertylist_get_string(conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_WRANGE        );
  conf->mf.grating.wave_range       = cpl_propertylist_get_string(conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_WRANGE        );

  conf->mf.grating.list_molec       = cpl_propertylist_get_string(conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_LIST          );
  conf->mf.grating.fit_molec        = cpl_propertylist_get_string(conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_FIT           );
  conf->mf.grating.rel_col          = cpl_propertylist_get_string(conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_RELATIVE_VALUE);

  conf->mf.kernmode                 = cpl_propertylist_get_bool(  conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_KERN_MODE             );
  conf->mf.kernfac                  = cpl_propertylist_get_double(conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_KERN_FAC              );
  conf->mf.varkern                  = cpl_propertylist_get_bool(  conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_VAR_KERN              );

  conf->mf.fit_continuum.fit        = cpl_propertylist_get_bool(  conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_FIT_CONTINUUM         );
  conf->mf.fit_continuum.n          = cpl_propertylist_get_int(   conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_CONTINUUM_N           );

  conf->mf.fit_wavelenght.fit       = cpl_propertylist_get_bool(  conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_FIT_WLC               );
  conf->mf.fit_wavelenght.n         = cpl_propertylist_get_int(   conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_WLC_N                 );
  conf->mf.fit_wavelenght.const_val = cpl_propertylist_get_double(conf->pl_best_fit_params, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_WLC_CONST             );

  cpl_msg_info(cpl_func, "Filling up the result headers  ...");

  /* Save parameter in the output propertylist */

  cpl_propertylist_update_string(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_WRANGE,         conf->mf.grating.wave_range      );

  cpl_propertylist_update_string(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_LIST,           conf->mf.grating.list_molec      );
  cpl_propertylist_update_string(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_FIT,            conf->mf.grating.fit_molec       );
  cpl_propertylist_update_string(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "KMOS_MOLECFIT_KEYWORD_RELATIVE_VALUE, conf->mf.grating.rel_col         );

  cpl_propertylist_update_bool(  conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_KERN_MODE,              conf->mf.kernmode                );
  cpl_propertylist_update_double(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_KERN_FAC,               conf->mf.kernfac                 );
  cpl_propertylist_update_bool(  conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_VAR_KERN,               conf->mf.varkern                 );

  cpl_propertylist_update_bool(  conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_FIT_CONTINUUM,          conf->mf.fit_continuum.fit       );
  cpl_propertylist_update_int(   conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_CONTINUUM_N,            conf->mf.fit_continuum.n         );

  cpl_propertylist_update_bool(  conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_FIT_WLC,                conf->mf.fit_wavelenght.fit      );
  cpl_propertylist_update_int(   conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_WLC_N,                  conf->mf.fit_wavelenght.n        );
  cpl_propertylist_update_double(conf->mf.parms, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_WLC_CONST,              conf->mf.fit_wavelenght.const_val);


  /* Check status */
  if (!cpl_errorstate_is_equal(initialState)) {
      /* Configuration failed */
      return cpl_error_get_code();
  }

  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Function needed to fill the molecfit configuration file
 *
 * @param  conf   Recipe configuration.
 *
 * @return parameterlist with contain the config to molecfit or NULL if error
 *
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code kmos_molecfit_calctrans_mf_conf(
    kmos_molecfit_calctrans_parameter *conf,
    double                            median,
    mf_parameters_config              *config_parameters)
{
  /* Check input */
  cpl_error_ensure(conf, CPL_ERROR_NULL_INPUT,
                   return CPL_ERROR_NULL_INPUT, "conf input is NULL!");

  /*** Building generic configuration molecfic file ***/
  cpl_error_code err = kmos_molecfit_conf_generic(&(conf->mf), config_parameters);
  if (err != CPL_ERROR_NONE) {
      return cpl_error_set_message(cpl_func, err,
                                   "Building molecfit configuration variable failed!");
  }

  /*** PARAMETERS NOT INCLUDED IN THE RECIPE: HARD-CODED ***/

  config_parameters->fitting.fit_continuum.fit        = conf->mf.fit_continuum.fit;
  config_parameters->fitting.fit_continuum.n          = conf->mf.fit_continuum.n;
  config_parameters->fitting.fit_continuum.const_val  = median;  /* Molecfit default: 1. */
  cpl_msg_info(cpl_func,"--.cont_const = %g", config_parameters->fitting.fit_continuum.const_val);

  config_parameters->fitting.fit_wavelenght.fit       = conf->mf.fit_wavelenght.fit;
  config_parameters->fitting.fit_wavelenght.n         = conf->mf.fit_wavelenght.n;
  config_parameters->fitting.fit_wavelenght.const_val = conf->mf.fit_wavelenght.const_val;

  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief    Deallocate the given parameter configuration object and its contents
 *
 * @param    conf       The parameter configuration variable in the recipe.
 */
/*----------------------------------------------------------------------------*/
static void kmos_molecfit_calctrans_clean(
    kmos_molecfit_calctrans_parameter *conf)
{
  if (conf) {

     if (&(conf->mf)) kmos_molecfit_clean(&(conf->mf));

      if (conf->pl_atmos_params)    cpl_propertylist_delete(                conf->pl_atmos_params   );
      if (conf->pl_best_fit_params) cpl_propertylist_delete(                conf->pl_best_fit_params);
      for (cpl_size i = 0; i < N_IFUS; i ++) {
          if (conf->atmprof[i]  )   cpl_table_delete(                       conf->atmprof[i]        );
          if (conf->res_table[i])   cpl_table_delete(                       conf->res_table[i]      );
          if (conf->results[i]  )   mf_calctrans_convolution_results_delete(conf->results[i]        );
      }

      cpl_free(conf);
  }
}
