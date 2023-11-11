/*
 * This file is part of the Molecfit Pipeline
 * Copyright (C) 2001-2019 European Southern Observatory
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

#ifdef HAVE_CONFIG_H
    #include <config.h>
#endif

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include <cpl.h>
#include <regex.h>

#include "molecfit.h"

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

  cpl_boolean                      use_only_input_pri_ext;       /* If the user want to use only the  SCIENCE_CALCTRANS/SCIENCE input FITS primary extension                                    */
  int                              dflux_extension_data;         /* If use_only_input_pri_ext == CPL_TRUE, you can provide a extension as DFLUX (error flux)                                    */
  int                              mask_extension_data;          /* If use_only_input_pri_ext == CPL_TRUE, you can provide a extension as MASK                                                  */

  cpl_boolean                      use_input_kernel;             /* If the user provided, use or not the user kernel                                                                            */

  const char                       *mapping_kernel;              /* Mapping extensions : SCIENCE-CALCTRANS_KERNEL_LIBRARY; i.e. {Science_ext1-Kernel_lib_ext, ..., Science_extN-Kernel_lib_ext} */
  cpl_table                        *mapping_kernel_table;        /* Mapping kernel extensions cpl_table. Contains 1 columns [KERNEL_LIBRARY_EXT]                                                */

  const char                       *mapping_atmospheric;         /* Mapping extensions : SCIENCE-ATM_PARAMETERS, i.e. {Science_ext1-Lblrtm_results_ext, ..., Science_extN-Lblrtm_results_ext}   */
  cpl_table                        *mapping_atmospheric_table;   /* Mapping lblrtm_results extensions cpl_table. Contains 1 columns [LBLRTM_RESULTS_EXT]                                        */

  const char                       *mapping_convolve;            /* Mapping extensions : SCIENCE-LBLRTM_RESULTS, i.e. {Science_ext1-Lblrtm_results_ext, ..., Science_extN-Lblrtm_results_ext}   */
  cpl_table                        *mapping_convolve_table;      /* Mapping lblrtm_results extensions cpl_table. Contains 1 columns [LBLRTM_RESULTS_EXT]                                        */

  /* Calibration */
  molecfit_fits                    *atm_parameters;              /* cpl_table's in the ATM_PARAMETERS      (output of molecfit_model)                                                           */
  molecfit_fits                    *best_fit_parameters;         /* cpl_table's in the BEST_FIT_PARAMETERS (output of molecfit_model) */
  const char					   *scale_pwv;
  const char					   *mjd_pwv;
  const char					   *exp_pwv;
  const char					   *air1_pwv;
  const char					   *air2_pwv;
  double					  	   pwv_sci;
  double					  	   pwv_ratio;
  double                           h2o_col_mm;                   /* Variable to store value read in from BEST_FIT_PARAMETERS table (output of molecfit_model) */

  int                              sg_window_length;             /* Savitzky-Golay filter smoothing window length in pixels */
  cpl_boolean                      sg_as_max_length;             /* Treat sg_window_length as a maximum length */


  /* Science */
  cpl_size                         n_ext;                        /* Number of extensions                                                                                                        */
  cpl_propertylist                 **telluriccorr_head;          /* Science spectrum header array                                                                                               */
  cpl_table                        **telluriccorr_data;          /* Science spectrum data   array                                                                                               */
  cpl_table                        **best_fit_parameters_table;  /* Best fit parameters array                                                                                                   */
  mf_calctrans_lblrtm_results      **results_lblrtm;             /* Results execute mf_calctrans_lblrtm for every input extension                                                               */
  mf_calctrans_convolution_results **results_convolution;        /* Results of telluric corrections after execute mf_calctrans_convolution                                                      */

  mf_configuration                 *mf_config;                   /* Telluric Correction configuration parameter                                                                                 */

  cpl_boolean                      chip_extensions;              /* If the user want to combine the input FITS extensions                                                                       */

  cpl_propertylist                 *pl;                          /* All input parameters to save in the output files                                                                            */

} molecfit_calctrans_parameter;

/*----------------------------------------------------------------------------*/
/**
 *                              Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Fill the internal Molecfit configuration file */
static molecfit_calctrans_parameter * molecfit_calctrans_parameters(
    const cpl_frameset              *frameset,
    const cpl_parameterlist         *list,
    const cpl_propertylist          *scientific_head_pri,
    const cpl_size                  n_ext);

/* Clean variables allocated in the recipe */
static void molecfit_calctrans_parameter_delete(
    molecfit_calctrans_parameter *parameters);

cpl_vector* savgol(int nl,int nr,int order,int degree);
cpl_vector* convol(cpl_vector* v, cpl_vector* kernel);
cpl_table* mean_absolute_difference(cpl_table* telluric_data,molecfit_calctrans_parameter* parameters);
cpl_vector* calc_rms(cpl_table* t,int wlength,int nrows_inc,double* w1_inc,double* w2_inc,int nrows_exc,double* w1_exc,double* w2_exc,double* fres,double* cres,double* sres);

/*----------------------------------------------------------------------------*/
/**
 *                          Static variables
 */
/*----------------------------------------------------------------------------*/

#define RECIPE_NAME      MOLECFIT_CALCTRANS
#define CONTEXT          "molecfit."RECIPE_NAME

static char molecfit_calctrans_description[] =
    "This recipe applies the results from molecfit_model and runs calctrans (lblrtm/convolution)\n"
    "   to calculate the telluric correction for scientific input data. "
    "If no kernel library is provided, molecfit_calctrans uses the kernel information provided in \n"
    "   the BEST_FIT_PARAMETERS.fits table (boxfwhm, gaussfwhm, lorentzfwhm) together with the \n"
    "   HIERARCH ESO keywords: DRS MOLECFIT PARAM KERNMODE, DRS MOLECFIT PARAM KERNFAC, and DRS MOLECFIT PARAM VARKERN.\n"
    "\n"
    "The input scientific input data can have category: \n"
    "  - BINTABLE  - Table  (single or multi-extension, optionally with some NULL extensions)\n"
    "  - IMAGE(1D) - Vector (single or multi-extension, optionally with some NULL extensions)\n"
    "  - IMAGE(2D) - Image  (single or multi-extension, optionally with some NULL extensions)\n"
    "It is not mandatory that all the extensions contain data.\n"
    "Molecfit will be run  on all the extensions that contain data.\n"
    "The recipe also accepts as input a kernel library.\n"
    "\n"
    "Input files:\n"
    "  - Mandatory,  the input "MOLECFIT_SCIENCE_CALCTRANS"/"MOLECFIT_SCIENCE" data (Only one).\n"
    "  - Mandatory,  the fit "MOLECFIT_MODEL_MOLECULES"              (Same output of molecfit_model).\n"
    "  - Optionally, the user kernel library.\n"
    "  - If the user kernel library is used, the parameter "MOLECFIT_CALCTRANS_MAPPING_KERNEL" or "MOLECFIT_MAPPING_KERNEL" is a list of the extensions to map from the kernel library file "MOLECFIT_CALCTRANS_KERNEL_LIBRARY" or "MOLECFIT_KERNEL_LIBRARY" to the extensions of the "MOLECFIT_SCIENCE_CALCTRANS" or "MOLECFIT_SCIENCE" input file.\n"
    "  - Mandatory,  "MOLECFIT_ATM_PARAMETERS"                       (Same output of molecfit_model).\n"
    "  - Mandatory,  "MOLECFIT_BEST_FIT_PARAMETERS"                  (Same output of molecfit_model).\n"
    "  - Mandatory,  "MOLECFIT_MAPPING_ATMOSPHERIC"                  (List of the extensions to map from the associated file "MOLECFIT_ATM_PARAMETERS" to the extensions of the "MOLECFIT_SCIENCE_CALCTRANS" or "MOLECFIT_SCIENCE" input file).\n"
    "  - Mandatory,  "MOLECFIT_MAPPING_CONVOLVE"                     (List of the extensions to map from the associated file "MOLECFIT_LBLRTM_RESULTS" to the extensions of the "MOLECFIT_TELLURIC_CORR" output).\n"
    "\n"
    "DO: Category                                 Required   Explanation\n"
    "-------------                                --------   -----------\n"
    "  "MOLECFIT_SCIENCE_CALCTRANS" or "MOLECFIT_SCIENCE"                 Y        The input spectrum        (Multi-extension FITS           format).\n"
    "  "MOLECFIT_MODEL_MOLECULES"                              Y        The molecules table       (         3-cols FITS BINTABLE  format).\n"
    "  "MOLECFIT_CALCTRANS_KERNEL_LIBRARY" or "MOLECFIT_KERNEL_LIBRARY"   N        The kernel library        (Multi-extension FITS IMAGE(2D) format).\n"
    "  "MOLECFIT_CALCTRANS_MAPPING_KERNEL" or "MOLECFIT_MAPPING_KERNEL"   N        Kernel mapping table      (         1-col  FITS BINTABLE  format) or the recipe parameter --"MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL".\n"
    "  "MOLECFIT_ATM_PARAMETERS"                               Y        Atmospheric  parameters   (Multi-extension FITS BINTABLE  format).\n"
    "  "MOLECFIT_BEST_FIT_PARAMETERS"                          Y        Best fitting parameters   (Multi-extension FITS BINTABLE  format).\n"
    "  "MOLECFIT_MAPPING_ATMOSPHERIC"                          N        Atmospheric mapping table (         1-col  FITS BINTABLE  format) or the recipe parameter --"MOLECFIT_PARAMETER_MAPPING_ATMOSPHERIC".\n"
    "  "MOLECFIT_MAPPING_CONVOLVE"                             N        Convolution mapping table (         1-col  FITS BINTABLE  format) or the recipe parameter --"MOLECFIT_PARAMETER_MAPPING_CONVOLVE".\n"
    "\n"
    "Notes:\n"
    "  - "MOLECFIT_CALCTRANS_MAPPING_KERNEL"  BINTABLE column  ["MOLECFIT_MAPPING_KERNEL_EXT"]\n"
    "  - "MOLECFIT_MAPPING_ATMOSPHERIC"       BINTABLE column  ["MOLECFIT_MAPPING_ATMOSPHERIC_EXT"]\n"
    "  - "MOLECFIT_MAPPING_CONVOLVE"          BINTABLE column  ["MOLECFIT_MAPPING_CONVOLVE_EXT"]\n"
    "\n"
    "\n"
    "Output files:\n"
    "\n"
    "DO: Category                      Explanation\n"
    "-------------                     -----------\n"
    "  "MOLECFIT_CALCTRANS_CHIPS_COMBINED"        Input file chips combined in molecfit format (FITS BINTABLE format).\n"
    "  "MOLECFIT_LBLRTM_RESULTS"                  For every input spectrum in the "MOLECFIT_MAPPING_ATMOSPHERIC" extension a new extension is generated (Multi-extension FITS BINTABLE format).\n"
    "  "MOLECFIT_TELLURIC_DATA"                   An intermediate data product BINTABLE that contains the input spectrum (lambda and flux columns), the telluric correction model (mlambda and mtrans columns) and the corrected input spectrum (lambda and cflux columns).\n"
    
    "  "MOLECFIT_TELLURIC_CORR"                   Telluric correction (transmission spectrum of the atmosphere) that needs to be applied (Multi-extension FITS IMAGE(1D) format).\n"
    "  "MOLECFIT_CALCTRANS_KERNEL_LIBRARY"        The kernel library specifies optional line spread functions that must be convolved with the atmospheric model.\n"



    "\n"
    "----------------------------------------------------------------------------\n"
    "The best fitting parameters as provided by the molecfit_model are used by molecfit_calctrans to generate a high-spectral resolution telluric transmission spectrum.\n"
    "Molecfit_calctrans takes into account the difference between the observed airmasses of the input data provided to molecfit_model and molecfit_calctrans.\n"
    "This spectrum is then convolved with either the kernel(s) with the best fitting parameters determined by molecfit_model or with the kernel library (if provided)\n"
    "to match the spectral resolution of the scientific spectrum.\n"
    "\n";

/* Standard CPL recipe definition */
cpl_recipe_define(	molecfit_calctrans,
                  	MOLECFIT_BINARY_VERSION,
                  	"Jose A. Escartin",
                  	"https://support.eso.org/",
                  	"2019",
                  	"Read the results from molecfit_model and computes the telluric correction for a scientific input data.",
                  	molecfit_calctrans_description);


/*----------------------------------------------------------------------------*/
/**
 * @defgroup molecfit_calctrans  It runs molectift on a generic input spectrum file to compute an atmospheric model.
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
static int molecfit_calctrans(
    cpl_frameset             *frameset,
    const cpl_parameterlist  *parlist)
{


  /* Get initial errorstate */
  cpl_error_ensure(frameset && parlist, CPL_ERROR_NULL_INPUT, return CPL_ERROR_NULL_INPUT, "NULL input : frameset and/or parlist");
  cpl_errorstate initial_errorstate = cpl_errorstate_get();

  /* Check frameset TAGS */
  const char* infomsg_str = \
  "SOF file must contain the following mandatory tags:\n\t\tSCIENCE;\n\t\tMODEL_MOLECULES;\n\t\tATM_PARAMETERS;\n\t\tBEST_FIT_PARAMETERS";
  cpl_msg_info(cpl_func, "Check frameset (SOF tags) ...");
  cpl_error_code err = molecfit_check_and_set_groups(frameset);

  /* Check mandatory TAGS/Parameters */
  cpl_frame *input_frame = NULL;
  if (!err) {

      /* SCIENCE_CALCTRANS/SCIENCE */
      cpl_errorstate pre_state = cpl_errorstate_get();
      input_frame = cpl_frameset_find(frameset, MOLECFIT_SCIENCE_CALCTRANS);
      if (!input_frame) {
          cpl_errorstate_set(pre_state);
          input_frame = cpl_frameset_find(frameset, MOLECFIT_SCIENCE);
          if (!input_frame) {
	      cpl_msg_info(cpl_func, "%s", infomsg_str);
              err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                          "SCIENCE_CALCTRANS/SCIENCE data not found in frameset!");
          }
      }

      /* MODEL_MOLECULES */
      if (!err) {
          const cpl_frame *input_model_molecules = cpl_frameset_find_const(frameset, MOLECFIT_MODEL_MOLECULES);
          if (!input_model_molecules) {
             cpl_msg_info(cpl_func, "%s", infomsg_str);
             err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                          "MODEL_MOLECULES data not found in frameset!");
	  }
      }

      /* ATM_PARAMETERS */
      if (!err) {
          const cpl_frame *input_atm_parameters = cpl_frameset_find_const(frameset, MOLECFIT_ATM_PARAMETERS);
          if (!input_atm_parameters) {
             cpl_msg_info(cpl_func, "%s", infomsg_str);
             err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                          "ATM_PARAMETERS data not found in frameset!");
	  }
      }

      /* BEST_FIT_PARAMETERS */
      if (!err) {
          const cpl_frame *input_best_fit_parameters = cpl_frameset_find_const(frameset, MOLECFIT_BEST_FIT_PARAMETERS);
          if (!input_best_fit_parameters) {
	     cpl_msg_info(cpl_func, "%s", infomsg_str);
             err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                          "BEST_FIT_PARAMETERS data not found in frameset!");
	  }
      }

      /* MAPPING_CONVOLVE */
      if (!err) {
          const cpl_frame *input_mapping_convolve = cpl_frameset_find_const(frameset, MOLECFIT_MAPPING_CONVOLVE);
          if (!input_mapping_convolve && (!strcmp(cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, MOLECFIT_PARAMETER_MAPPING_CONVOLVE)), MF_PARAMETERS_NULL) ) ){
              err = CPL_ERROR_ILLEGAL_INPUT;
          }
      }
  }

  /* TEST KERNEL MAPPING */
  if (!err) {
      cpl_error_code map_err=CPL_ERROR_NONE;
      map_err=molecfit_config_mappingchk(frameset,parlist,MOLECFIT_CALCTRANS_KERNEL_LIBRARY,
                                                          "MAPPING_KERNEL",
							  "KERNEL_LIBRARY_EXT",
							  MOLECFIT_CALCTRANS_MAPPING_KERNEL,
							  "NULL");
      if (map_err) cpl_msg_error(cpl_func,"Mapping Errors detected cannot proceed further");
      
  }

  /* TEST _ATMOSPHERIC: MAPPING */
  if (!err) {
      cpl_error_code map_err=CPL_ERROR_NONE;
      map_err=molecfit_config_mappingchk(frameset,parlist,"ATM_PARAMETERS",
                                                          "MAPPING_ATMOSPHERIC",
							  "ATM_PARAMETERS_EXT",
							  "MAPPING_ATMOSPHERIC",
							  "SCIENCE");
      if (map_err) cpl_msg_error(cpl_func,"Mapping Errors detected cannot proceed further");
      
  }
  /* TEST _CONVOLVE: MAPPING */
  if (!err) {
      cpl_error_code map_err=CPL_ERROR_NONE;
      map_err=molecfit_config_mappingchk(frameset,parlist, MOLECFIT_LBLRTM_RESULTS,
                                                           MOLECFIT_MAPPING_CONVOLVE,
							  "LBLRTM_RESULTS_EXT",
							  "MAPPING_CONVOLVE",
							   MOLECFIT_TELLURIC_CORR);
      if (map_err) cpl_msg_error(cpl_func,"Mapping Errors detected cannot proceed further");
      
  }


  /* Load TAG = SCIENCE_CALCTRANS/SCIENCE */
  molecfit_fits *data = NULL;
  if (!err) {
      const char *science_file = cpl_frame_get_filename(input_frame);

      data = molecfit_fits_load(science_file, CPL_FALSE);

      if (!data) err = CPL_ERROR_ILLEGAL_INPUT;
      else       err = cpl_error_get_code();
  }

  /* Recipe Parameters : Need scientific_header_primary */
  molecfit_calctrans_parameter *parameters = NULL;
  if (!err) {

      /* Get recipe parameters and update the molecfit default configuration */
      cpl_msg_info(cpl_func, "Load '%s' recipe parameters ...", MOLECFIT_CALCTRANS);
      parameters = molecfit_calctrans_parameters(frameset, parlist, data->v_ext[0].header, data->n_ext);

      if (!parameters) err = CPL_ERROR_ILLEGAL_INPUT;
      else             err = cpl_error_get_code();

      /* Check if the science_data is a FITS BINTABLE : If not, convert to input molecfit_data cpl_table */
      if (!err) {

          if (!(parameters->mf_config->parameters->fitting.fit_wavelenght.fit)) cpl_msg_info(cpl_func, "Not fit wavelength!");
          if (!(parameters->mf_config->parameters->fitting.fit_continuum.fit )) cpl_msg_info(cpl_func, "Not fit continuum!" );

          err = molecfit_data_convert_to_table( data,
                                                parameters->chip_extensions,
                                                parameters->use_only_input_pri_ext,
                                                parameters->dflux_extension_data,
                                                parameters->mask_extension_data,
                                                parameters->mf_config->parameters->inputs.column_lambda,
                                                parameters->mf_config->parameters->inputs.column_flux,
                                                parameters->mf_config->parameters->inputs.column_dflux,
                                                parameters->mf_config->parameters->inputs.column_mask);

          /* Save input DATA --> Combined (molecfit_spectrum in BINTABLE format) */
          if (parameters->chip_extensions) {

              cpl_size index_combined = data->v_ext[0].spectrum_data ? 0 : 1;
              cpl_msg_info(cpl_func, "Save combined multi-extension molecfit_spectrum input FITS file ('%s') in ext=%lld ...", MOLECFIT_CALCTRANS_CHIPS_COMBINED, index_combined);

              err += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_CALCTRANS_CHIPS_COMBINED, NULL);
              if (!err) err = molecfit_save_mf_results(data->v_ext[index_combined].spectrum_head, MOLECFIT_CALCTRANS_CHIPS_COMBINED, CPL_TRUE, NULL, data->v_ext[index_combined].spectrum_data, NULL);
          }
      }
  }

  /* Load TAG = MOLECULES */
  cpl_table *molecules = NULL;
  if (!err) {

      cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_MODEL_MOLECULES);
      molecules = molecfit_load_unique_table(frameset, MOLECFIT_MODEL_MOLECULES);

      if (!molecules) err = CPL_ERROR_ILLEGAL_INPUT;
      else            err = cpl_error_get_code();
  }


  /* Load TAG : CALCTRANS_KERNEL_LIBRARY */
  molecfit_fits *kernel_data    = NULL;
  cpl_table       *mapping_kernel = NULL;
  if (!err) {

      kernel_data = molecfit_load_kernel_tag(frameset, MOLECFIT_CALCTRANS_KERNEL_LIBRARY, parameters->use_input_kernel, parameters->mf_config->parameters,CPL_TRUE);
      err = cpl_error_get_code();

      if (kernel_data && !err) {

          if (parameters->mapping_kernel_table) {

              /* Mapping in the recipe parameters */
              mapping_kernel = cpl_table_duplicate(parameters->mapping_kernel_table);

          } else {

              /* Mapping in the static_calib input FITS files */
              cpl_errorstate pre_state = cpl_errorstate_get();
              const cpl_frame *input_mapping_kernel = cpl_frameset_find(frameset, MOLECFIT_CALCTRANS_MAPPING_KERNEL);
              if (input_mapping_kernel) {
                  cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_CALCTRANS_MAPPING_KERNEL);
                  mapping_kernel = molecfit_load_unique_table(frameset, MOLECFIT_CALCTRANS_MAPPING_KERNEL);
              } else {
                  cpl_errorstate_set(pre_state);
                  input_mapping_kernel = cpl_frameset_find(frameset, MOLECFIT_MAPPING_KERNEL);
                  if (input_mapping_kernel) {
                      cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_MAPPING_KERNEL);
                      mapping_kernel = molecfit_load_unique_table(frameset, MOLECFIT_MAPPING_KERNEL);
                  } else {
                      err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                                  "Mapping kernel data not found in frameset and neither in the recipe parameter!");
                  }
              }
          }

          if (!mapping_kernel) err = CPL_ERROR_ILLEGAL_INPUT;
          else                 err = cpl_error_get_code();
      }

      if (!err) cpl_msg_info(cpl_func, "Using the default molecfit kernels -> With the values inside BEST_FIT_PARMS!");
  }


  /* Load TAG = MAPPING_ATMOSPHERIC */
  cpl_table *mapping_atmospheric = NULL;
  if (!err) {

      if (parameters->mapping_atmospheric_table) {
          mapping_atmospheric = cpl_table_duplicate(parameters->mapping_atmospheric_table);
      } else {
          cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_MAPPING_ATMOSPHERIC);
          mapping_atmospheric = molecfit_load_unique_table(frameset, MOLECFIT_MAPPING_ATMOSPHERIC);
      }

      if (!mapping_atmospheric) err = CPL_ERROR_ILLEGAL_INPUT;
      else                      err = cpl_error_get_code();
  }


  /* Load TAG = MAPPING_CONVOLVE */
  cpl_table *mapping_convolve = NULL;
  if (!err) {

      if (parameters->mapping_convolve_table) {
          mapping_convolve = cpl_table_duplicate(parameters->mapping_convolve_table);
      } else {
          cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_MAPPING_CONVOLVE);
          mapping_convolve = molecfit_load_unique_table(frameset, MOLECFIT_MAPPING_CONVOLVE);
      }

      if (!mapping_convolve) err = CPL_ERROR_ILLEGAL_INPUT;
      else                   err = cpl_error_get_code();
  }


  /*** Save generic outut files */
  /*These molecfit_save functions just save the placeholder (primary header). The actual data extensions are written later */
  if (!err) {

      cpl_msg_info(cpl_func, "Save generic multi-extension output FITS file ('%s','%s') ...", MOLECFIT_TELLURIC_DATA, MOLECFIT_TELLURIC_CORR);
      /* Extract some additional QC info from the BEST_FIT_PARAMETERS input file */
      cpl_propertylist* qc_orig = cpl_propertylist_new();
      cpl_propertylist* qc_new = cpl_propertylist_new();
      /* This hdu has the full primary header copy in it */
      cpl_propertylist* bfp_hdr =  parameters->best_fit_parameters->v_ext[1].header;
      cpl_propertylist* data_hdr = data->v_ext[0].header;
      err     += cpl_propertylist_copy_property_regexp(qc_orig,bfp_hdr,"^ESO PRO REC2",0);
      /* adjust the keyword names so they do not conflict with the existing ESO PRO REC2 ones */
      for (int idx=0;idx < cpl_propertylist_get_size(qc_orig);idx++){
          cpl_property* p = cpl_property_duplicate(cpl_propertylist_get(qc_orig,idx));
          cpl_property_set_name(p,cpl_sprintf("%s BFP",cpl_property_get_name(p)));
          cpl_propertylist_append_property(qc_new,p);
      }
      /*Try and get the difference between airmass and observation time ; do each separately in case either are missing keywords*/
      /* Should probably use HIERARCH ESO DRS MF PARAM OBSERVING_DATE_VALUE but this seems broken at the moment? */
      const char* key_exp = parameters->exp_pwv;
      const char* key_mjd = parameters->mjd_pwv; 
      const char* key_a1 = parameters->air1_pwv;
      const char* key_a2 = parameters->air2_pwv;

      /*If not none, then raise an error */
      int raise_err = strcmp(parameters->scale_pwv,"none");
      /*key_mjd*/
      if(!cpl_propertylist_has(bfp_hdr,key_mjd)){
          if(raise_err) {
            cpl_msg_error(cpl_func,"Missing %s keyword from BEST_FIT_PARAMETERS.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",key_mjd);
            return CPL_ERROR_ILLEGAL_INPUT;
          }
          else
            cpl_msg_warning(cpl_func,"Missing %s keyword from BEST_FIT_PARAMETERS.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",key_mjd);
      }
      if(!cpl_propertylist_has(data_hdr,key_mjd)){
          if(raise_err) {
            cpl_msg_error(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",key_mjd);
            return CPL_ERROR_ILLEGAL_INPUT;
          }
          else
            cpl_msg_warning(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",key_mjd);
      }
      /*key_exp*/
      if(!cpl_propertylist_has(bfp_hdr,key_exp)){
          if(raise_err){
            cpl_msg_error(cpl_func,"Missing %s keyword from BEST_FIT_PARAMETERS.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",key_exp);
            return CPL_ERROR_ILLEGAL_INPUT;
          }
          else
            cpl_msg_warning(cpl_func,"Missing %s keyword from BEST_FIT_PARAMETERS.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",key_exp);
      }
      if(!cpl_propertylist_has(data_hdr,key_exp)){
          if(raise_err){
            cpl_msg_error(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",key_exp);
            return CPL_ERROR_ILLEGAL_INPUT;
          }
          else 
            cpl_msg_warning(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA MJD.\n",key_exp);
      }
      /*key_a1*/
      if(!cpl_propertylist_has(bfp_hdr,key_a1)){
          if(raise_err){
            cpl_msg_error(cpl_func,"Missing %s keyword from BEST_FIT_PARAMETERS.\n\tCannot populate ESO DRS PWV DELTA AIRM.\n",key_a1);
            return CPL_ERROR_ILLEGAL_INPUT;
          }
          else
            cpl_msg_warning(cpl_func,"Missing %s keyword from BEST_FIT_PARAMETERS.\n\tCannot populate ESO DRS PWV DELTA AIRM.\n",key_a1);
      }
      if(!cpl_propertylist_has(data_hdr,key_a1)){
          if(raise_err){
            cpl_msg_error(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA AIRM.\n",key_a1);
            return CPL_ERROR_ILLEGAL_INPUT;
          }
          else 
            cpl_msg_warning(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA AIRM.\n",key_a1);
      }
      /*key_a2*/
      if(!cpl_propertylist_has(bfp_hdr,key_a2)){
          if(raise_err){
            cpl_msg_error(cpl_func,"Missing %s keyword from BEST_FIT_PARAMETERS.\n\tCannot populate ESO DRS PWV DELTA AIRM.\n",key_a2);
            return CPL_ERROR_ILLEGAL_INPUT;
          }
          else
            cpl_msg_warning(cpl_func,"Missing %s keyword from BEST_FIT_PARAMETERS.\n\tCannot populate ESO DRS PWV DELTA AIRM.\n",key_a2);
      }
      if(!cpl_propertylist_has(data_hdr,key_a2)){
          if(raise_err){
            cpl_msg_error(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA AIRM.\n",key_a2);
            return CPL_ERROR_ILLEGAL_INPUT;
          }
          else 
            cpl_msg_warning(cpl_func,"Missing %s keyword from science spectrum.\n\tCannot populate ESO DRS PWV DELTA AIRM.\n",key_a2);
      }


      if(cpl_propertylist_has(bfp_hdr,key_exp) && cpl_propertylist_has(bfp_hdr,key_mjd) &&
         cpl_propertylist_has(data_hdr,key_exp) && cpl_propertylist_has(data_hdr,key_mjd)){
          double t0 = cpl_propertylist_get_double(bfp_hdr,key_mjd); 
          cpl_type exp_type = cpl_propertylist_get_type(bfp_hdr,key_exp);
          if(exp_type == CPL_TYPE_INT){
              t0 = t0 + 0.5*cpl_propertylist_get_int(bfp_hdr,key_exp);
          } else if (exp_type == CPL_TYPE_FLOAT){
              t0 = t0 + 0.5*cpl_propertylist_get_float(bfp_hdr,key_exp);
          } else if (exp_type == CPL_TYPE_DOUBLE){
              t0 = t0 + 0.5*cpl_propertylist_get_double(bfp_hdr,key_exp);
          } 
          exp_type = cpl_propertylist_get_type(data_hdr,key_exp);
          double t1 = cpl_propertylist_get_double(data_hdr,key_mjd); 
          if(exp_type == CPL_TYPE_INT){
              t1 = t1 + 0.5*cpl_propertylist_get_int(data_hdr,key_exp);
          } else if (exp_type == CPL_TYPE_FLOAT){
              t1 = t1 + 0.5*cpl_propertylist_get_float(data_hdr,key_exp);
          } else if (exp_type == CPL_TYPE_DOUBLE){
              t1 = t1 + 0.5*cpl_propertylist_get_double(data_hdr,key_exp);
          } 

          cpl_propertylist_append_double(qc_new,"ESO DRS PWV DELTA MJD",t1-t0);
          cpl_propertylist_set_comment(qc_new,"ESO DRS PWV DELTA MJD","Sci-Telluric Mean MJD");
      }

      if(cpl_propertylist_has(bfp_hdr,key_a1) && cpl_propertylist_has(bfp_hdr,key_a2) &&
         cpl_propertylist_has(data_hdr,key_a1) && cpl_propertylist_has(data_hdr,key_a2)){
          double a0 = 0.5*(cpl_propertylist_get_double(bfp_hdr,key_a1) +  cpl_propertylist_get_double(bfp_hdr,key_a2));
          double a1 = 0.5*(cpl_propertylist_get_double(data_hdr,key_a1) +  cpl_propertylist_get_double(data_hdr,key_a2));
          cpl_propertylist_append_double(qc_new,"ESO DRS PWV DELTA AIRM",a1-a0);
          cpl_propertylist_set_comment(qc_new,"ESO DRS PWV DELTA AIRM","Sci-Telluric Mean Airmass");
      }
      
      /*cpl_propertylist_dump(qc_new,stdout);*/
      /* add some other values */
      cpl_propertylist_append_double(qc_new,"ESO DRS PWV SCI",parameters->pwv_sci);
      cpl_propertylist_append_double(qc_new,"ESO DRS PWV TELL",parameters->h2o_col_mm);
      cpl_propertylist_append_double(qc_new,"ESO DRS PWV RATIO",parameters->pwv_ratio);

      cpl_propertylist_set_comment(qc_new,"ESO DRS PWV SCI","Value for the science spectrum");
      cpl_propertylist_set_comment(qc_new,"ESO DRS PWV TELL","Best fit from telluric spectrum");
      cpl_propertylist_set_comment(qc_new,"ESO DRS PWV RATIO","SCI/TELL");

      /*Append to the primary header so that all the output products have this info*/
      err     += cpl_propertylist_append(parameters->pl,qc_new);

      err     += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_LBLRTM_RESULTS, NULL);

      err     += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_TELLURIC_DATA,  NULL);
      err     += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_TELLURIC_CORR,  NULL);
      //err     += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, "TELLURIC_CORR_SMOOTH",  NULL);

      if (kernel_data) {
          err += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_CALCTRANS_KERNEL_LIBRARY, NULL);
      }
  }

  /*** Execute molecfit (lblrtm) ***/
  if (!err) {

      int null;

      /* Execution extensions */
      cpl_size n_ext;
      if (     parameters->use_only_input_pri_ext) n_ext = 1;
      else if (parameters->chip_extensions    ) n_ext = data->v_ext[0].spectrum_data ? 1 : 2;
      else n_ext = data->n_ext;

      for (cpl_size ext = 0; ext < n_ext && !err; ext++) {

          /* Create input molecfit spec format */
          if (data->v_ext[ext].spectrum_data) {

              /* Get extension header/table data */
              cpl_msg_info(cpl_func, "Load spectrum for execute mf_calctrans_lblrtm(...) in ext = %lld ...", ext);
              parameters->telluriccorr_head[ext] = cpl_propertylist_duplicate(data->v_ext[ext].spectrum_head);
              parameters->telluriccorr_data[ext] = mf_spectrum_create(parameters->mf_config->parameters, data->v_ext[ext].spectrum_data);

              /* Get ATM_PARAMETERS/BEST_FIT_PARAMETERS index */
              cpl_size index_atmospheric = cpl_table_get(mapping_atmospheric, MOLECFIT_MAPPING_ATMOSPHERIC_EXT, ext, &null);

              /* This adds a clarifying message to the 'Access beyond boundaries' error that pops up if the MAPPING is incorrect */
              if(parameters->atm_parameters->v_ext[index_atmospheric].table == NULL){
                  err = cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                                      "Cannot find atm_parameters. Please check the incorrectly set MAPPING_ATMOSPHERIC parameter.");
              }

              /* Save BEST_FIT_PARAMETERS cpl_table for CONVOLUTION */
              parameters->best_fit_parameters_table[ext] = cpl_table_duplicate(parameters->best_fit_parameters->v_ext[index_atmospheric].table);

              /* CALL CALCTRANS_LBLRTM : Select all wavelength range for the Molecfit executions : wl_start = -1 and wl_end = -1 */
              const double wl_start = -1.;
              const double wl_end   = -1.;
              parameters->results_lblrtm[ext] = mf_calctrans_lblrtm( parameters->mf_config,                                             /* mf_configuration *config              */
                                                                     parameters->telluriccorr_data[ext],                                    /* const cpl_table  *spec_telluriccorr   */
                                                                     molecules,                                                         /* cpl_table        *molecules           */
                                                                     wl_start,                                                          /* double           wl_start             */
                                                                     wl_end,                                                            /* double           wl_end               */
                                                                     parameters->atm_parameters->v_ext[index_atmospheric].table,        /* cpl_table        *atm_parameters      */
                                                                     parameters->best_fit_parameters->v_ext[index_atmospheric].table);  /* cpl_table        *best_fit_parameters */

              /* Check possible errors */
              if (!(parameters->results_lblrtm[ext])) {
                  err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                              "Unexpected error in the Molecfit call mf_calctrans_lblrtm(...)");
              } else {

                  /* Only one range allow in the mf_calctrans_lblrtm(...) execution --> Mandatory */
                  cpl_size n_range = parameters->mf_config->parameters->internal.n_range;
                  if (n_range != 1 && !parameters->chip_extensions) {

                      err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                                                   "Unexpected n_ranges in the mf_calctrans_lblrtm(...) Molecfit execution : n_ranges = %lld (Mandatory == 1)", n_range);

                  } else {

                      /* Check execution results for the only range (position == 0) */
                      if (parameters->results_lblrtm[ext]->range_status[0] != CPL_ERROR_NONE) {

                          err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                                      "mf_calctrans_lblrtm(...) Molecfit execution : n_ranges = %lld (Mandatory == 1)", n_range);

                      } else {

                          /* Check other errors */
                          err = cpl_error_get_code();
                      }
                  }
              }
          }
      }


      /*** Execute molecfit (convolution) ***/
      for (cpl_size ext = 0; ext < n_ext && !err; ext++) {

          /* Get ATM_PARAMETERS/BEST_FIT_PARAMETERS index */
          cpl_size index_convolve       = cpl_table_get(mapping_convolve, MOLECFIT_MAPPING_CONVOLVE_EXT, ext, &null);
          cpl_size index_lblrtm_results = index_convolve - parameters->use_only_input_pri_ext;

          /* This adds a clarifying message to the 'Access beyond boundaries' error that pops up if the MAPPING is incorrect */ 
/*
          if(parameters->results_lblrtm[index_lblrtm_results] == NULL){
              err = cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT, "Cannot find lblrtm_results. Please check the incorrectly set MAPPING_CONVOLVE parameter.");
          }
*/
          /* Create input molecfit spec format */
          if (parameters->results_lblrtm[index_lblrtm_results]) {

              /* Get extension header data */
              cpl_msg_info(cpl_func, "Convolve input spectrum (ext_orig[ATM_PARAMETERS/BEST_FIT_PARAMETERS] = %lld) for execute mf_calctrans_convolve(...) in ext = %lld ...", index_convolve, ext);

              cpl_propertylist *header_kernel = NULL;
              cpl_matrix       *kernel        = NULL;
              if (kernel_data) {

                  if (mapping_kernel) {
                      cpl_size index_kernel_ext = cpl_table_get(mapping_kernel, MOLECFIT_MAPPING_KERNEL_EXT, ext, &null);
                      header_kernel = kernel_data->v_ext[index_kernel_ext].header;
                      kernel        = kernel_data->v_ext[index_kernel_ext].matrix;
                  } else {
                      header_kernel = kernel_data->v_ext[ext].header;
                      kernel        = kernel_data->v_ext[ext].matrix;
                  }
              }

              /* CALL CALCTRANS_CONVOLUTION : Select all wavelength range for the Molecfit executions : wl_start = -1 and wl_end = -1 */
              const double wl_start = -1.;
              const double wl_end   = -1.;
              parameters->results_convolution[ext] = mf_calctrans_convolution( parameters->mf_config->parameters,                             /* mf_parameters_config        *config              */
                                                                               parameters->mf_config->lnfl,
                                                                               parameters->results_lblrtm[index_lblrtm_results],              /* mf_calctrans_lblrtm_results *lblrtm_results      */
                                                                               parameters->telluriccorr_head[index_lblrtm_results],           /* const cpl_propertylist      *header_spec         */
                                                                               parameters->telluriccorr_data[index_lblrtm_results],           /* const cpl_table             *spec_telluriccorr   */
                                                                               header_kernel,                                                 /* const cpl_propertylist      *header_kernel       */
                                                                               kernel,                                                        /* const cpl_matrix            *kernel              */
                                                                               wl_start,                                                      /* double                      wl_start             */
                                                                               wl_end,                                                        /* double                      wl_end               */
                                                                               parameters->best_fit_parameters_table[index_lblrtm_results]);  /* cpl_table                   *best_fit_parameters */

              /* Check possible errors */
              if (!(parameters->results_convolution[ext])) {
                  err = cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                              "Unexpected error in the Molecfit call mf_calctrans_convolution(...)");
              } else {
                  err = cpl_error_get_code();
              }
          }
      }

      /* Save results */
      for (cpl_size ext = 0; ext < n_ext && !err; ext++) {

          cpl_propertylist *lblrtm_header   = cpl_propertylist_duplicate(data->v_ext[ext].header);
          cpl_table        *lblrtm_spec_out          = NULL;
                                                  
          cpl_table        *telluric_data            = NULL;
          cpl_vector       *telluric_vec             = NULL;

          //cpl_table        *selected_data = NULL;
          cpl_table        *residual_data = NULL;
          cpl_size maxpos=-1;
          double cres_at_max_res_ratio=-99.99;
          int max_res_length=-1;

          //cpl_vector       *telluric_vec_smoothed    = NULL;
          //flux = data->v_ext[index_combined].spectrum_data
          //cflux = 

          cpl_matrix       *kernel_matrix            = NULL;

          if (parameters->results_lblrtm[ext]) {

              cpl_msg_info(cpl_func, "LBLRTM_RESULTS exist in this ext = %lld ... saving ...", ext);

              /* Save tmp_dir in the LBLRTM_RESULTS header */
              cpl_propertylist_update_string(lblrtm_header, MF_PARAMETERS_CONTEX_DEFAULT" "MF_PARAMETERS_TMP_PATH, parameters->results_lblrtm[ext]->tmp_folder);

              /* Get LBLRTM results spectrum */
              lblrtm_spec_out = parameters->results_lblrtm[ext]->spec_out[0];
          }

          if (parameters->results_convolution[ext]) {

              /* Get CONVOLUTION results cpl_table */
              telluric_data = (parameters->results_convolution[ext])->spec_telluriccorr_format;
              int tdata_rows = cpl_table_get_nrow(telluric_data);

              /* Wrap the data */
              double *telluric_corr_column = cpl_table_get_data_double(telluric_data, MF_COL_OUT_TELLURIC_CORR);

              /* Convert in a spectrum (cpl_vector) */
              cpl_vector *vAux = cpl_vector_wrap(tdata_rows, telluric_corr_column);
              telluric_vec = cpl_vector_duplicate(vAux);

              /* The inputs for the smoothing comparison are as follows:
                wave_in = wavelength or mlambda
                flux_in = flux
                sflux_in = weight
                cflux_in = cflux
                N.B. molecfit_calctrans does not currently know about the regions used for fitting in molecfit_model.
              */
              //cpl_table_dump(telluric_data,0,5,stdout);
              //cpl_msg_info(cpl_func,"pre mean_absolute_difference");
              residual_data = mean_absolute_difference(telluric_data,parameters);
              //cpl_msg_info(cpl_func,"after mean_absolute_difference");
              if(cpl_table_get_column_maxpos(residual_data,"res_ratio",&maxpos) == CPL_ERROR_NONE){
                cres_at_max_res_ratio = cpl_table_get_double(residual_data,"cflux_res",maxpos,NULL);
                //max_res_ratio = cpl_table_get_double(residual_data,"res_ratio",maxpos,NULL);
                max_res_length = cpl_table_get_int(residual_data,"window_length",maxpos,NULL);
                cpl_msg_info(cpl_func,"max_res_length cres_at_max_res_ratio: %d %e\n",max_res_length,cres_at_max_res_ratio);
              }
              cpl_vector_unwrap(vAux);
              /* Get Kernel results */
              if (parameters->results_convolution[ext]->kernel_resampled_normalized) {
                  kernel_matrix = parameters->results_convolution[ext]->kernel_resampled_normalized;
              }
          }

          /*** SAVE results : LBLRTM_RESULTS, TELLURIC_DATA, TELLURIC_CORR and KERNEL_MOLECFIT results ***/
          if (!err && (parameters->use_only_input_pri_ext || ext > 0)) {

              cpl_msg_info(cpl_func, "Saving %s, %s, %s, %s ... (ext =%lld : lblrtm_results ? %d, convolution_results ? %d)", MOLECFIT_LBLRTM_RESULTS, MOLECFIT_CALCTRANS_KERNEL_LIBRARY, MOLECFIT_TELLURIC_DATA, MOLECFIT_TELLURIC_CORR, ext, lblrtm_spec_out != NULL, telluric_data != NULL);
                
              cpl_propertylist* use_hdr = cpl_propertylist_duplicate(data->v_ext[ext].header);
              //Add Savitzky-Golay filtered QC params to header for TELLURIC_DATA and TELLURIC_CORR
              if(residual_data){
                cpl_propertylist_append_double(use_hdr,"ESO QC MEAN_ABS_DEV MAX",cres_at_max_res_ratio);
                cpl_propertylist_append_int(use_hdr,"ESO QC MEAN_ABS_DEV MAX WIDTH",max_res_length);
                //this only really makes sense for case where parameters->sg_as_max_length == CPL_FALSE
                cpl_propertylist_append_int(use_hdr,"ESO QC MEAN_ABS_DEV INPUT WIDTH",parameters->sg_window_length);
                //cpl_propertylist_update_string(use_hdr,"EXTNAME","RESIDUAL_DATA");
                //err     += molecfit_save_mf_results(sm_hdr,"TELLURIC_CORR_SMOOTH",            CPL_TRUE, NULL,          residual_data,     NULL);
              }

              err     += molecfit_save_mf_results( lblrtm_header,                   MOLECFIT_LBLRTM_RESULTS,           CPL_TRUE, NULL,          lblrtm_spec_out, NULL         );

              err     += molecfit_save_mf_results( use_hdr,         MOLECFIT_TELLURIC_DATA,            CPL_TRUE, NULL,          telluric_data,   NULL         );
              err     += molecfit_save_mf_results( use_hdr,         MOLECFIT_TELLURIC_CORR,            CPL_TRUE, NULL,          NULL,            telluric_vec );
              /*if(selected_data){ 
                cpl_propertylist* sm_hdr = cpl_propertylist_duplicate(data->v_ext[ext].header);
                cpl_propertylist_update_string(sm_hdr,"EXTNAME","SELECTED_DATA");
                err     += molecfit_save_mf_results(sm_hdr,"TELLURIC_CORR_SMOOTH",            CPL_TRUE, NULL,          selected_data,     NULL);
              }*/



              if (kernel_data) {
                  err += molecfit_save_mf_results( kernel_data->v_ext[ext].header,  MOLECFIT_CALCTRANS_KERNEL_LIBRARY, CPL_TRUE, kernel_matrix, NULL,            NULL         );
              }
          }

          /* Cleanup */
          if (lblrtm_header) cpl_propertylist_delete(lblrtm_header);
          if (telluric_vec ) cpl_vector_delete(telluric_vec);
      }
  }

  /* Cleanup */
  if (parameters         ) molecfit_calctrans_parameter_delete( parameters         );
  if (data               ) molecfit_fits_delete(                data               );
  if (molecules          ) cpl_table_delete(                    molecules          );
  if (kernel_data        ) molecfit_fits_delete(                kernel_data        );
  if (mapping_kernel     ) cpl_table_delete(                    mapping_kernel     );
  if (mapping_atmospheric) cpl_table_delete(                    mapping_atmospheric);
  if (mapping_convolve   ) cpl_table_delete(                    mapping_convolve   );


  /* Check Recipe status and end */
  if (!err && cpl_errorstate_is_equal(initial_errorstate)) {
      cpl_msg_info(cpl_func,"Recipe successfully!");
  } else {
      /* Dump the error history */
      cpl_errorstate_dump(initial_errorstate, CPL_FALSE, NULL);
      cpl_msg_error(cpl_func,"Recipe failed!, error(%d)=%s", err, cpl_error_get_message());
  }

  return err;
}

cpl_vector* calc_rms(cpl_table* t,int wlength,int nrows_inc,double* w1_inc,double* w2_inc,int nrows_exc,double* w1_exc,double* w2_exc,double* fres,double* cres,double* sres){
    //1. Duplicate tdata
    //2. Add columns for smoothed data using Savitzky-Golay filter with window_length=wlength
    //3. Filter table based on *inc/*exc data
    //4. Calculate statistics and populate the results 
    cpl_table* tdata = cpl_table_duplicate(t);
    int tdata_rows = cpl_table_get_nrow(tdata);
    //calculate smoothed versions of the columns flux and cflux
    cpl_msg_info(cpl_func,"calc_rms: Calculating Savitzky-Golay coefficients for window_length %d\n",wlength);
    cpl_vector* coeffs = savgol(floor(wlength/2),floor(wlength/2),0,3);

    /*before we do anything, we need to remove the bad quality (null rows), since these can introduce nan values */
    int *qual_col_orig = cpl_table_get_data_int(tdata, MF_COL_OUT_MASK);
    cpl_error_code err = cpl_table_unselect_all(tdata);
    for (int i=0;i<tdata_rows;i++){
        if(qual_col_orig[i] == 1){
            err = cpl_table_select_row(tdata,i);
        }
    }
            
    tdata=cpl_table_extract_selected(tdata);
    tdata_rows = cpl_table_get_nrow(tdata);
    //lambda == MF_COL_IN_LAMBDA, mlambda == MF_COL_MOD_LAMBDA 
    double *wave_col = cpl_table_get_data_double(tdata, MF_COL_MOD_LAMBDA);
    double *flux_col = cpl_table_get_data_double(tdata, MF_COL_IN_FLUX);
    double *sflux_col = cpl_table_get_data_double(tdata, MF_COL_MOD_WEIGHT);
    double *cflux_col = cpl_table_get_data_double(tdata, MF_COL_OUT_FLUX);
    cpl_vector *flux = cpl_vector_wrap(tdata_rows, flux_col);
    cpl_vector *sflux = cpl_vector_wrap(tdata_rows, sflux_col);
    cpl_vector *cflux = cpl_vector_wrap(tdata_rows, cflux_col);
    cpl_vector *flux_sm = convol(flux,coeffs);
    cpl_vector *cflux_sm = convol(cflux,coeffs);
    cpl_vector *flux_norm = cpl_vector_new(tdata_rows); 
    cpl_vector *cflux_norm = cpl_vector_new(tdata_rows); 
    cpl_vector *sflux_norm = cpl_vector_new(tdata_rows); 
    cpl_msg_info(cpl_func,"calc_rms: Calculating normalised flux columns");
    for(int i=0;i<tdata_rows;i++){
        cpl_vector_set(flux_norm,i,1.0-(cpl_vector_get(flux,i)/cpl_vector_get(flux_sm,i)));
        cpl_vector_set(cflux_norm,i,1.0-(cpl_vector_get(cflux,i)/cpl_vector_get(cflux_sm,i)));
        cpl_vector_set(sflux_norm,i,(cpl_vector_get(sflux,i)/cpl_vector_get(flux_sm,i)));
    }

    //calculate the normalised columns

    //populate the new columns
    //(add the new columns)
    cpl_table_new_column(tdata,"flux_sm",CPL_TYPE_DOUBLE);
    cpl_table_new_column(tdata,"cflux_sm",CPL_TYPE_DOUBLE);
    cpl_table_new_column(tdata,"flux_norm",CPL_TYPE_DOUBLE);
    cpl_table_new_column(tdata,"cflux_norm",CPL_TYPE_DOUBLE);
    cpl_table_new_column(tdata,"sflux_norm",CPL_TYPE_DOUBLE);
    cpl_table_copy_data_double(tdata,"flux_sm",cpl_vector_get_data(flux_sm));
    cpl_table_copy_data_double(tdata,"cflux_sm",cpl_vector_get_data(cflux_sm));
    cpl_table_copy_data_double(tdata,"flux_norm",cpl_vector_get_data(flux_norm));
    cpl_table_copy_data_double(tdata,"cflux_norm",cpl_vector_get_data(cflux_norm));
    cpl_table_copy_data_double(tdata,"sflux_norm",cpl_vector_get_data(sflux_norm));
    //(fill the new columns)

    //filter the table...
    cpl_msg_info(cpl_func,"calc_rms: Filtering results table");
    int *qual_col = cpl_table_get_data_int(tdata, MF_COL_OUT_MASK);
    err = cpl_table_unselect_all(tdata);
    for (int i=0;i<tdata_rows;i++){
        //cpl_msg_info(cpl_func,"qual_col[%d] = %d\n",i,qual_col[i]);
        if(qual_col[i] == 1){
            cpl_boolean sel_row = CPL_FALSE;
            //first select the row if it lies in an inclusion range
            if(nrows_inc != -1){
                for(int j=0;j<nrows_inc;j++){
                    if(wave_col[i] >= w1_inc[j] && wave_col[i] <= w2_inc[j])
                      sel_row = CPL_TRUE; 
                }
            }
            //then exclude any rows if they lie in an exclusion range
            if(nrows_exc != -1){
                for(int j=0;j<nrows_exc;j++){
                    if(wave_col[i] >= w1_exc[j] && wave_col[i] <= w2_exc[j])
                      sel_row = CPL_FALSE; 
                }
            }
            if(nrows_inc == -1 && nrows_exc == -1){
                sel_row = CPL_TRUE;
            }
            //if the row was selected, actually select it
            if(sel_row == CPL_TRUE){
              err = cpl_table_select_row(tdata,i);
            }
        }
    }
    tdata=cpl_table_extract_selected(tdata);
    tdata_rows = cpl_table_get_nrow(tdata);

    /* Calculate statistics from the cleaned up table */
    cpl_msg_info(cpl_func,"calc_rms: Calculating statistics on filtered results table");

    //get the table data
    double *fnorm_col = cpl_table_get_data_double(tdata, "flux_norm");
    double *cfnorm_col = cpl_table_get_data_double(tdata, "cflux_norm");
    double *sfnorm_col = cpl_table_get_data_double(tdata, "sflux_norm");
    cpl_vector* fnorm = cpl_vector_wrap(tdata_rows, fnorm_col);
    cpl_vector* cfnorm = cpl_vector_wrap(tdata_rows, cfnorm_col);
    cpl_vector* sfnorm = cpl_vector_wrap(tdata_rows, sfnorm_col);

    //now the calculations
    double flux_norm_mean = cpl_vector_get_mean(fnorm);
    double cflux_norm_mean = cpl_vector_get_mean(cfnorm);
    double sflux_norm_mean = cpl_vector_get_mean(sfnorm);

    cpl_vector* flux_norm_res  = cpl_vector_new(tdata_rows);
    cpl_vector* cflux_norm_res = cpl_vector_new(tdata_rows);
    cpl_vector* sflux_norm_res = cpl_vector_new(tdata_rows);
    for(int i=0;i<tdata_rows;i++){
        cpl_vector_set(flux_norm_res,i,fabs(cpl_vector_get(flux_norm,i) - flux_norm_mean));
        cpl_vector_set(cflux_norm_res,i,fabs(cpl_vector_get(cflux_norm,i) - cflux_norm_mean));
        cpl_vector_set(sflux_norm_res,i,fabs(cpl_vector_get(sflux_norm,i) - sflux_norm_mean));
    }
    *fres = cpl_vector_get_mean(flux_norm_res);
    *cres = cpl_vector_get_mean(cflux_norm_res);
    *sres = cpl_vector_get_mean(sflux_norm_res);
    
    return cpl_vector_duplicate(cflux_sm);

}
cpl_table* mean_absolute_difference(cpl_table* telluric_data,molecfit_calctrans_parameter* parameters){
    //TODO: very important, smooth the full spectrum and then filter wavelength ranges afterwards...ugh
    //this is necessary to avoid introducing bad artefacts into the smoothed functions...

    //TODO: ASSERT maxwl >= 3

    //First find the wavelength region information from molecfit_model via best_fit_parameters
    //If we have WAVE_INCLUDE or WAVE_EXCLUDE in BEST_FIT_PARAMETERS, use them, 
    //otherwise,  use entire wavelength range
    //extension indices
    //number of rows in the tables
    int nrows_inc = -1;
    int nrows_exc = -1;
    //wavelength ranges to include
    double* w1_inc = NULL; 
    double* w2_inc = NULL; 
    //wavelength ranges to exclude
    double* w1_exc = NULL; 
    double* w2_exc = NULL; 

    cpl_msg_info(cpl_func,"mean_absolute_difference: calculating wave ranges");
    for(int i=0;i<parameters->best_fit_parameters->n_ext;i++){
        if(cpl_propertylist_has(parameters->best_fit_parameters->v_ext[i].header,"EXTNAME")){
            if(!strcmp(cpl_propertylist_get_string(parameters->best_fit_parameters->v_ext[i].header,"EXTNAME"),"WAVE_INCLUDE")){
                if(cpl_table_has_column(parameters->best_fit_parameters->v_ext[i].table,MF_COL_WAVE_RANGE_LOWER) &&
                      cpl_table_has_column(parameters->best_fit_parameters->v_ext[i].table,MF_COL_WAVE_RANGE_UPPER)){
                    nrows_inc =cpl_table_get_nrow(parameters->best_fit_parameters->v_ext[i].table);
                    w1_inc = cpl_table_get_data_double(parameters->best_fit_parameters->v_ext[i].table,MF_COL_WAVE_RANGE_LOWER);
                    w2_inc = cpl_table_get_data_double(parameters->best_fit_parameters->v_ext[i].table,MF_COL_WAVE_RANGE_UPPER);
                }

            }
            if(!strcmp(cpl_propertylist_get_string(parameters->best_fit_parameters->v_ext[i].header,"EXTNAME"),"WAVE_EXCLUDE")){
                if(cpl_table_has_column(parameters->best_fit_parameters->v_ext[i].table,MF_COL_WAVE_RANGE_LOWER) &&
                      cpl_table_has_column(parameters->best_fit_parameters->v_ext[i].table,MF_COL_WAVE_RANGE_UPPER)){
                    nrows_exc =cpl_table_get_nrow(parameters->best_fit_parameters->v_ext[i].table);
                    w1_exc = cpl_table_get_data_double(parameters->best_fit_parameters->v_ext[i].table,MF_COL_WAVE_RANGE_LOWER);
                    w2_exc = cpl_table_get_data_double(parameters->best_fit_parameters->v_ext[i].table,MF_COL_WAVE_RANGE_UPPER);
                }
            }
        }
    }


    int wlength;
    if(parameters->sg_as_max_length == CPL_FALSE){
        double fres, cres, sres;
        cpl_vector* cflux_sm;
        wlength = parameters->sg_window_length;
        cpl_msg_info(cpl_func,"mean_absolute_difference: single window_length specified: %d",wlength);
        cflux_sm = calc_rms(telluric_data,wlength,nrows_inc,w1_inc,w2_inc,nrows_exc,w1_exc,w2_exc,&fres,&cres,&sres);

        /*First fixup the cflux_sm vector so that it includes any bad data from the original table */
        cpl_size nrows_orig = cpl_table_get_nrow(telluric_data);
        cpl_vector* cflux_sm_proper = cpl_vector_new(nrows_orig); 
        int *qual_col = cpl_table_get_data_int(telluric_data, MF_COL_OUT_MASK);
        int counter = 0;
        for (int i=0;i<nrows_orig; i++){
            if(qual_col[i] == 0){
                cpl_vector_set(cflux_sm_proper,i,0.0);
            } else {
                cpl_vector_set(cflux_sm_proper,i,cpl_vector_get(cflux_sm,counter));
                counter = counter + 1;
            }
        }
        /*Add column to telluric data with smoothed spectrum (cflux_sm) */
        cpl_table_new_column(telluric_data,"cflux_savgol",CPL_TYPE_DOUBLE);
        cpl_table_copy_data_double(telluric_data,"cflux_savgol",cpl_vector_get_data(cflux_sm_proper));

        //fres = cres = sres = 1.0;
        cpl_msg_info(cpl_func,"mean_absolute_difference: calc_rms called: %e %e %e",fres,cres,sres);
        cpl_table* res = cpl_table_new(1);
        cpl_table_new_column(res,"window_length",CPL_TYPE_INT);
        cpl_table_new_column(res,"flux_res",CPL_TYPE_DOUBLE);
        cpl_table_new_column(res,"cflux_res",CPL_TYPE_DOUBLE);
        cpl_table_new_column(res,"sflux_res",CPL_TYPE_DOUBLE);
        cpl_table_new_column(res,"res_ratio",CPL_TYPE_DOUBLE);
        cpl_table_set(res,"window_length",0,wlength);
        cpl_table_set(res,"flux_res",0,fres);
        cpl_table_set(res,"cflux_res",0,cres);
        cpl_table_set(res,"sflux_res",0,sres);
        cpl_table_set(res,"res_ratio",0,fres/cres);
        return cpl_table_duplicate(res);
    } else {
        int maxwl = parameters->sg_window_length;
        int minwl = 5;
        int count = 0;
        for(wlength=maxwl;wlength>=minwl;wlength-=2){
            count+=1;
        }
        cpl_table* res = cpl_table_new(count);
        cpl_table_new_column(res,"window_length",CPL_TYPE_INT);
        cpl_table_new_column(res,"flux_res",CPL_TYPE_DOUBLE);
        cpl_table_new_column(res,"cflux_res",CPL_TYPE_DOUBLE);
        cpl_table_new_column(res,"sflux_res",CPL_TYPE_DOUBLE);
        cpl_table_new_column(res,"res_ratio",CPL_TYPE_DOUBLE);
        int idx = 0;
        cpl_vector* cflux_sm = NULL;
        cpl_vector* cflux_sm_keep = NULL;
        double running_res_max;
        for(wlength=maxwl;wlength>=minwl;wlength-=2){
            double fres, cres, sres;
            cflux_sm = calc_rms(telluric_data,wlength,nrows_inc,w1_inc,w2_inc,nrows_exc,w1_exc,w2_exc,&fres,&cres,&sres);

            cpl_table_set(res,"window_length",idx,wlength);
            cpl_table_set(res,"flux_res",idx,fres);
            cpl_table_set(res,"cflux_res",idx,cres);
            cpl_table_set(res,"sflux_res",idx,sres);
            double res_ratio = fres/cres;
            cpl_table_set(res,"res_ratio",idx,res_ratio);
            //keep track of the max res ratio.
            if(wlength == maxwl){
                running_res_max = res_ratio;
                cflux_sm_keep = cpl_vector_duplicate(cflux_sm);
            }
            if(res_ratio > running_res_max){
                running_res_max = res_ratio;
                if(cflux_sm_keep) cpl_vector_delete(cflux_sm_keep);
                cflux_sm_keep = cpl_vector_duplicate(cflux_sm);
            }
            idx += 1;
        }
        /*First fixup the cflux_sm vector so that it includes any bad data from the original table */
        cpl_size nrows_orig = cpl_table_get_nrow(telluric_data);
        cpl_vector* cflux_sm_proper = cpl_vector_new(nrows_orig); 
        int *qual_col = cpl_table_get_data_int(telluric_data, MF_COL_OUT_MASK);
        int counter = 0;
        for (int i=0;i<nrows_orig; i++){
            if(qual_col[i] == 0){
                cpl_vector_set(cflux_sm_proper,i,0.0);
            } else {
                cpl_vector_set(cflux_sm_proper,i,cpl_vector_get(cflux_sm,counter));
                counter = counter + 1;
            }
        }
        /*Add column to telluric data with smoothed spectrum (cflux_sm) */
        cpl_table_new_column(telluric_data,"cflux_savgol",CPL_TYPE_DOUBLE);
        cpl_table_copy_data_double(telluric_data,"cflux_savgol",cpl_vector_get_data(cflux_sm_keep));
        
        return cpl_table_duplicate(res);
    }
}

/* Computes the Savitzky-Golay filter coefficients following the IDL function savgol.pro.
   Results match those from savgol_coeffs function from scipy.
*/
cpl_vector* savgol(int nl,int nr,int order,int degree){
    cpl_matrix* a = cpl_matrix_new(degree+1, degree+1);
    cpl_matrix* b = cpl_matrix_new(degree+1,1);
    //cpl_vector* b = cpl_vector_new(degree+1);
    cpl_vector* cr = cpl_vector_new(nr);
    cpl_vector* cl = cpl_vector_new(nl);
    cpl_vector* power = cpl_vector_new(degree+1);
    int np = nl + nr + 1;
    cpl_vector* c = cpl_vector_new(np);

    if(nl < 0 || nr < 0){
      cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,"savgol: Both nl (%d) and nr (%d) must be positive",nl,nr);
        return NULL;
    }
    if(order > degree){
      cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,"savgol: The order (%d) must be less than or equal to the degree (%d)",order,degree);
        return NULL;
    }
    if(degree >= np){
      cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,"savgol: The degree (%d) must be less than the filter width (nl+nr+1=%d).",degree,np);
        return NULL;
    }

    int i,j,k;
    for(i=0;i<nr;i++){
        cpl_vector_set(cr,i,1.0+i*1.0);
    }
    for(i=0;i<nl;i++){
        cpl_vector_set(cl,i,-(1.0+i*1.0));
    }
    for(i=0;i<degree+1;i++){
        cpl_vector_set(power,i,i*1.0);
        if(i == order){
            cpl_matrix_set(b,i,0,1.0);
            //cpl_vector_set(b,i,1.0);
        } else {
            cpl_matrix_set(b,i,0,0.0);
            //cpl_vector_set(b,i,0.0);
        }
    }

    for(i=0;i<=degree*2;i++){
        double sum = (i == 0) ? 1.0 : 0.0; 
        if(nr > 0){
            for(j=0;j<nr;j++){
                sum = sum + pow(cpl_vector_get(cr,j),i);
            }
        }
        if(nl > 0){
            for(j=0;j<nl;j++){
                sum = sum + pow(cpl_vector_get(cl,j),i);
            }
        }
        double mm = i < (2 * degree - i) ? i : (2 * degree - i);
        for(j=-mm;j <= mm; j+=2){
            cpl_matrix_set(a,(i+j)/2,(i-j)/2,sum);
        }
    }
    //gsl_matrix_fprintf(stdout,a,"%.2f");
    /*for(i=0;i<degree+1;i++){
        for(j=0;j<degree+1;j++){
            printf("a[%d][%d] = %.2f\n",i,j,cpl_matrix_get(a,i,j));
        }
    }*/
    //unsigned long N = cpl_matrix_get_nrow(a);//->size1;//, N = A->size2;
    //unsigned long N = a->size1;//, N = A->size2;

    //int signum;
    //int s =0;

    //gsl_permutation * perm = gsl_permutation_alloc(N);
    //gsl_matrix * lu  = cpl_matrix_new(N, N);
    //gsl_vector * res = cpl_vector_new(N);
  
    //printf("A size: %d %d\n",A->size1,A->size2);
    //printf("b size: %d\n",b->size);
    //gsl_matrix_memcpy(lu, a);
  
    //s += gsl_linalg_LU_decomp(lu, perm, &signum);
    //s += gsl_linalg_LU_solve(lu, perm, b, res);   
    cpl_matrix* res = cpl_matrix_solve(a,b);
    if(res == NULL){
      cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
                            "Could not solve matrices for Savitzky-Golay coefficients: %s",cpl_error_get_message());
        return NULL;
    }
    //cpl_msg_info(cpl_func,"Error: %s",cpl_error_get_message());
    //cpl_matrix_dump(res,stdout);
    for(k=-nl;k<=nr;k++){
        double sum = 0.0;
        for(i=0;i<degree+1;i++){
            sum = sum+(cpl_matrix_get(res,i,0)*pow(k,cpl_vector_get(power,i)));
            //sum = sum+(cpl_vector_get(res,i)*pow(k,cpl_vector_get(power,i)));
        }
        cpl_vector_set(c,k+nl,sum);
    }
    //these are the Savitzky-Golay coefficients
    //for(i=0;i<np;i++){
    //    printf("c[%d] = %.8f\n",i,cpl_vector_get(c,i));
    //}
    return c;

}
cpl_vector* convol(cpl_vector* v, cpl_vector* kernel){
    /* Code taken from cpl_vector_filter_lowpass_create */
    int nv = cpl_vector_get_size(v);
    int nk = cpl_vector_get_size(kernel);
    cpl_vector* filtered = cpl_vector_new(nv);
    int hw = floor(nk/2);
    double replace;
    cpl_size i,j;

    /* compute edge effects for the first hw elements */
    for (i = 0; i < hw; i++) {
        replace = 0.0;
        for (j = -hw; j <= hw; j++) {
            if (i + j < 0) {
                replace += cpl_vector_get(kernel,hw + j) * cpl_vector_get(v,0);
            }
            else {
                replace += cpl_vector_get(kernel,hw + j) * cpl_vector_get(v,i + j);
            }
        }
        cpl_vector_set(filtered,i,replace);
    }

    /* compute edge effects for the last hw elements */
    for (i = nv - hw; i < nv; i++) {
        replace = 0.0;
        for (j = -hw; j <= hw; j++) {
            if (i + j > nv - 1) {
                replace += cpl_vector_get(kernel,hw + j) * cpl_vector_get(v,nv - 1);
            }
            else {
                replace += cpl_vector_get(kernel,hw + j) * cpl_vector_get(v,i + j);
            }
        }
        cpl_vector_set(filtered,i,replace);
    }

    /* compute all other elements */
    for (i = hw; i < nv - hw; i++) {
        replace = 0.0;
        for (j = -hw; j <= hw; j++) {
            replace += cpl_vector_get(kernel,hw + j) * cpl_vector_get(v,i + j);
        }
        cpl_vector_set(filtered,i,replace);
    }
    //cpl_vector_delete(kernel);
    //cpl_tools_add_flops(4 * hw * nv);
    return filtered;
}   

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
static cpl_error_code molecfit_calctrans_fill_parameterlist(
    cpl_parameterlist *self)
{
  /* Add the different default parameters to the recipe */
  cpl_errorstate pre_state = cpl_errorstate_get();
  cpl_error_code e         = CPL_ERROR_NONE;

  /* Fill the parameters list */
  cpl_boolean    range     = CPL_TRUE;
  const void     *dummyMin = NULL;
  const void     *dummyMax = NULL;


  /* --MOLECFIT_PARAMETER_USE_ONLY_INPUT_PRI */
  cpl_boolean use_only_input_pri_ext = MOLECFIT_PARAMETER_USE_ONLY_INPUT_PRI_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_USE_ONLY_INPUT_PRI,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_BOOL, MF_STRLST_NULL,
						   &use_only_input_pri_ext,
                                                   MOLECFIT_PARAMETER_USE_ONLY_INPUT_PRI_DESC, CPL_FALSE);

  /* --MOLECFIT_PARAMETER_DFLUX_EXTENSION_DATA */
  int dflux_extension_data = MOLECFIT_PARAMETER_DFLUX_EXTENSION_DATA_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_DFLUX_EXTENSION_DATA,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, MF_STRLST_NULL,
						   &dflux_extension_data,
                                                   MOLECFIT_PARAMETER_DFLUX_EXTENSION_DATA_DESC, CPL_FALSE);

  /* --MOLECFIT_PARAMETER_MASK_EXTENSION_DATA */
  int mask_extension_data = MOLECFIT_PARAMETER_MASK_EXTENSION_DATA_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_MASK_EXTENSION_DATA,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, MF_STRLST_NULL,
						   &mask_extension_data,
                                                   MOLECFIT_PARAMETER_MASK_EXTENSION_DATA_DESC, CPL_FALSE);


  /* --MOLECFIT_PARAMETER_USE_INPUT_KERNEL */
  cpl_boolean use_input_kernel = MOLECFIT_PARAMETER_USE_INPUT_KERNEL_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_USE_INPUT_KERNEL,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_BOOL, MF_STRLST_NULL,
						   &use_input_kernel,
                                                   MOLECFIT_PARAMETER_USE_INPUT_KERNEL_DESC, CPL_FALSE);

    /* --MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL */
  const char *mapping_kernel = MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
						   (const void *)mapping_kernel,
                                                   MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL_DESC, CPL_FALSE);


  /* --MOLECFIT_PARAMETER_MAPPING_ATMOSPHERIC */
  const char *mapping_atmospheric = MOLECFIT_PARAMETER_MAPPING_ATMOSPHERIC_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_MAPPING_ATMOSPHERIC,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
						   (const void *)mapping_atmospheric,
                                                   MOLECFIT_PARAMETER_MAPPING_ATMOSPHERIC_DESC, CPL_FALSE);


  /* --MOLECFIT_PARAMETER_MAPPING_CONVOLVE */
  const char *mapping_convolution = MOLECFIT_PARAMETER_MAPPING_CONVOLVE_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_MAPPING_CONVOLVE,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
						   (const void *)mapping_convolution,
                                                   MOLECFIT_PARAMETER_MAPPING_CONVOLVE_DESC, CPL_FALSE);


  /* --MOLECFIT_PARAMETER_COMBINE_EXTENSIONS */
  cpl_boolean combine_extensions = MOLECFIT_PARAMETER_CHIP_EXTENSIONS_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_CHIP_EXTENSIONS,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_BOOL, MF_STRLST_NULL,
						   &combine_extensions,
                                                   MOLECFIT_PARAMETER_CHIP_EXTENSIONS_DESC, CPL_FALSE);

  /* --MOLECFIT_PARAMETER_CALCTRANS_SCALE_PWV */
  const char *scale_pwv = MOLECFIT_PARAMETER_CALCTRANS_SCALE_PWV_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_CALCTRANS_SCALE_PWV,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
						   (const void *)scale_pwv,
                                                   MOLECFIT_PARAMETER_CALCTRANS_SCALE_PWV_DESC, CPL_FALSE); 

  /* Some parameters to choose what fits header keywords to use */
  const char *mjd_pwv = MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_MJD_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_MJD,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
						   (const void *)mjd_pwv,
                                                  MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_MJD_DESC, CPL_FALSE); 

  const char *exp_pwv = MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_EXP_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_EXP,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
						   (const void *)exp_pwv,
                                                  MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_EXP_DESC, CPL_FALSE); 

  const char *air1_pwv = MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR1_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR1,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
						   (const void *)air1_pwv,
                                                  MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR1_DESC, CPL_FALSE); 

  const char *air2_pwv = MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR2_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR2,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
						   (const void *)air2_pwv,
                                                  MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR2_DESC, CPL_FALSE); 



  /* Savitzky-Golay parameters */
  int sgwl = MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self,MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_INT, MF_STRLST_NULL,&sgwl,
                                                   MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_DESC, CPL_FALSE); 

  cpl_boolean sgwl_asmax = MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_AS_MAX_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self,MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_AS_MAX,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_BOOL, MF_STRLST_NULL,&sgwl_asmax,
                                                   MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_AS_MAX_DESC, CPL_FALSE); 

  /* Check possible errors */
  if (!cpl_errorstate_is_equal(pre_state) || e != CPL_ERROR_NONE) {
      return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                   "molecfit_calctrans_fill_parameterlist failed!");
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
static molecfit_calctrans_parameter * molecfit_calctrans_parameters(
    const cpl_frameset              *frameset,
    const cpl_parameterlist         *list,
    const cpl_propertylist          *scientific_head_pri,
    const cpl_size                  n_ext)
{
  /* Check input */
  cpl_error_ensure(frameset && list, CPL_ERROR_NULL_INPUT,
                   return NULL, "list input is NULL!");

  /* Get preState */
  cpl_errorstate preState = cpl_errorstate_get();
  const cpl_parameter *p;

  /* Create the configuration parameter */
  molecfit_calctrans_parameter *parameters = cpl_malloc(sizeof(molecfit_calctrans_parameter));
  parameters->mapping_kernel_table               = NULL;
  parameters->mapping_atmospheric_table          = NULL;
  parameters->mapping_convolve_table             = NULL;
  parameters->atm_parameters                     = NULL;
  parameters->best_fit_parameters                = NULL;
  parameters->telluriccorr_head                  = NULL;
  parameters->telluriccorr_data                  = NULL;
  parameters->best_fit_parameters_table          = NULL;
  parameters->results_lblrtm                     = NULL;
  parameters->results_convolution                = NULL;
  parameters->mf_config                          = NULL;
  parameters->pl                                 = cpl_propertylist_new();

  parameters->scale_pwv  						 = MOLECFIT_PARAMETER_CALCTRANS_SCALE_PWV_INIT;
  parameters->mjd_pwv  						     = MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_MJD_INIT;
  parameters->exp_pwv  						     = MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_EXP_INIT;
  parameters->air1_pwv  						 = MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR1_INIT;
  parameters->air2_pwv  						 = MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR2_INIT;
  parameters->pwv_ratio  						 = 1.0;
  parameters->pwv_sci  							 = -99.9;
  parameters->h2o_col_mm                         = -99.9;
  parameters->sg_window_length                   = MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_INIT; //15
  parameters->sg_as_max_length                   = MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_AS_MAX_INIT; //CPL_FALSE


  /* Use only primary extension in the input FITS file ? */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_USE_ONLY_INPUT_PRI);
  parameters->use_only_input_pri_ext = cpl_parameter_get_bool(p);
  cpl_propertylist_update_bool(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_USE_ONLY_INPUT_PRI, parameters->use_only_input_pri_ext);

  /* If parameters->use_only_input_pri_ext == CPL_TRUE, you can provide a error flux extension */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_DFLUX_EXTENSION_DATA);
  parameters->dflux_extension_data = cpl_parameter_get_int(p);
  cpl_propertylist_update_int(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_DFLUX_EXTENSION_DATA, parameters->dflux_extension_data);

  /* If parameters->use_only_input_pri_ext == CPL_TRUE, you can provide a mask extension */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_MASK_EXTENSION_DATA);
  parameters->mask_extension_data = cpl_parameter_get_int(p);
  cpl_propertylist_update_int(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_MASK_EXTENSION_DATA, parameters->mask_extension_data);


  /* User input kernel ? */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_USE_INPUT_KERNEL);
  parameters->use_input_kernel = cpl_parameter_get_bool(p);
  cpl_propertylist_update_bool(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_USE_INPUT_KERNEL, parameters->use_input_kernel);

  /* Mapping kernel */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL);
  parameters->mapping_kernel = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL, parameters->mapping_kernel);
  /* Create mapping kernel cpl_table */
  if (strcmp(parameters->mapping_kernel, MF_PARAMETERS_NULL)) {
      parameters->mapping_kernel_table = molecfit_config_table_mapping(parameters->mapping_kernel, MOLECFIT_MAPPING_KERNEL_EXT);
  }


  /* Scale PWV  */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CALCTRANS_SCALE_PWV);
  if (p) {
	  parameters->scale_pwv = cpl_parameter_get_string(p);
	  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CALCTRANS_SCALE_PWV, parameters->scale_pwv);
  }
  
  /* Other header keywords used for airmass calculation. We do not ask for RA, DEC and LST as these are more standard. */
  
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_MJD);
  if (p) {
	  parameters->mjd_pwv = cpl_parameter_get_string(p);
	  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_MJD, parameters->mjd_pwv);
  }

  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_EXP);
  if (p) {
	  parameters->exp_pwv = cpl_parameter_get_string(p);
	  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_EXP, parameters->exp_pwv);
  }

  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR1);
  if (p) {
	  parameters->air1_pwv = cpl_parameter_get_string(p);
	  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR1, parameters->air1_pwv);
  }
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR2);
  if (p) {
	  parameters->air2_pwv = cpl_parameter_get_string(p);
	  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CALCTRANS_PWV_HDR_AIR2, parameters->air2_pwv);
  }

  /* Savitzky-Golay filter parameters */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH);
  if (p) {
	  parameters->sg_window_length= cpl_parameter_get_int(p);
	  cpl_propertylist_update_int(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH, parameters->sg_window_length);
  }

  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_AS_MAX);
  if (p) {
	  parameters->sg_as_max_length= cpl_parameter_get_bool(p);
	  cpl_propertylist_update_bool(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CALCTRANS_SG_WINDOW_LENGTH_AS_MAX, parameters->sg_as_max_length);
  }




  /* Mapping ATMOSPHERIC : SCIENCE_CALCTRANS/SCIENCE - ATM_PARAMETERS */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_MAPPING_ATMOSPHERIC);
  parameters->mapping_atmospheric = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_MAPPING_ATMOSPHERIC, parameters->mapping_atmospheric);
  /* Create mapping atmospheric cpl_table */
  if (strcmp(parameters->mapping_atmospheric, MF_PARAMETERS_NULL)) {
      parameters->mapping_atmospheric_table = molecfit_config_table_mapping(parameters->mapping_atmospheric, MOLECFIT_MAPPING_ATMOSPHERIC_EXT);
  }


  /* Mapping CONVOLVE : LBLRTM_RESULTS - TELLURIC_CORR */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_MAPPING_CONVOLVE);
  parameters->mapping_convolve = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_MAPPING_CONVOLVE, parameters->mapping_convolve);
  /* Create mapping convolve cpl_table */
  if (strcmp(parameters->mapping_convolve, MF_PARAMETERS_NULL)) {
      parameters->mapping_convolve_table = molecfit_config_table_mapping(parameters->mapping_convolve, MOLECFIT_MAPPING_CONVOLVE_EXT);
  }


  /* Load TAG = MOLECFIT_ATM_PARAMETERS */
  cpl_msg_info(cpl_func, "Load %s from molecfit_model(...)", MOLECFIT_ATM_PARAMETERS);
  const cpl_frame *frm_atm_parameters = cpl_frameset_find_const(frameset, MOLECFIT_ATM_PARAMETERS);
  if (frm_atm_parameters) {
      const char *filename = cpl_frame_get_filename(frm_atm_parameters);
      parameters->atm_parameters = molecfit_fits_load(filename, CPL_FALSE);
  }

  /* Load TAG = MOLECFIT_BEST_FIT_PARAMETERS */
  cpl_msg_info(cpl_func, "Load %s from molecfit_model(...)", MOLECFIT_BEST_FIT_PARAMETERS);
  const cpl_frame *frm_best_fit_parameters = cpl_frameset_find_const(frameset, MOLECFIT_BEST_FIT_PARAMETERS);
  if (frm_best_fit_parameters) {
      const char *filename = cpl_frame_get_filename(frm_best_fit_parameters);
      parameters->best_fit_parameters = molecfit_fits_load(filename, CPL_FALSE);
  }

  /* Check if MOLECFIT_ATM_PARAMETERS and MOLECFIT_BEST_FIT_PARAMETERS were loaded */
  if (   !(parameters->atm_parameters       ) || !(parameters->best_fit_parameters       )
      ){
      //||   parameters->atm_parameters->n_ext  !=   parameters->best_fit_parameters->n_ext) {
      molecfit_calctrans_parameter_delete(parameters);
      cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_model output files not correct or not provided!");
      return NULL;
  }

  /* Reserve memory for the output results */
  parameters->n_ext                              = n_ext;
  parameters->telluriccorr_head                  = cpl_calloc(parameters->n_ext, sizeof(cpl_table                        *));
  parameters->telluriccorr_data                  = cpl_calloc(parameters->n_ext, sizeof(cpl_table                        *));
  parameters->best_fit_parameters_table          = cpl_calloc(parameters->n_ext, sizeof(cpl_table                        *));
  parameters->results_lblrtm                     = cpl_calloc(parameters->n_ext, sizeof(mf_calctrans_lblrtm_results      *));
  parameters->results_convolution                = cpl_calloc(parameters->n_ext, sizeof(mf_calctrans_convolution_results *));
  for (cpl_size ext = 0; ext < n_ext; ext++) {
      parameters->telluriccorr_head[ext]         = NULL;
      parameters->telluriccorr_data[ext]         = NULL;
      parameters->best_fit_parameters_table[ext] = NULL;
      parameters->results_lblrtm[ext]            = NULL;
      parameters->results_convolution[ext]       = NULL;
  }


  /* Get recipe parameters in the recipe molecfit_model --> Best_fit_parameters */
  parameters->mf_config = molecfit_config_load_header_parameters(parameters->best_fit_parameters->v_ext[0].header, scientific_head_pri, parameters->pl);

  regex_t regex;
  int reti;
  reti = regcomp(&regex,"[a-zA-Z0-9]",0);
  if(reti){
    cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,"molecfit_calctrans: could not compile regex required to check SCALE_PWV");
          return NULL;
  }
  reti = regexec(&regex,parameters->scale_pwv,0,NULL,0);

  /* This needs to be updated to work for all n_ext */
  int ridx;
  int nrows;
  if(n_ext > 0 && parameters->atm_parameters->v_ext[1].table != NULL){
      nrows  = cpl_table_get_nrow(parameters->best_fit_parameters->v_ext[1].table);
      for (ridx=0;ridx<nrows;ridx++){
        if(!strcmp(cpl_table_get_string(parameters->best_fit_parameters->v_ext[1].table,"parameter",ridx),"h2o_col_mm")){
            parameters->h2o_col_mm = cpl_table_get_double(parameters->best_fit_parameters->v_ext[1].table,"value",ridx,NULL);
            break;
        }
      }
  }
  /*Do we need to throw an error here if we cannot read in h2o_col_mm ?*/

  /*Calculate pwv_sci */
  if (!strcmp(parameters->scale_pwv,"auto")){
      if(cpl_propertylist_has(scientific_head_pri,"ESO TEL AMBI IWV START") && cpl_propertylist_has(scientific_head_pri,"ESO TEL AMBI IWV END")){
          double iwv_start = cpl_propertylist_get_double(scientific_head_pri,"ESO TEL AMBI IWV START");
          double iwv_end = cpl_propertylist_get_double(scientific_head_pri,"ESO TEL AMBI IWV END");
          if(iwv_start == 0.0 && iwv_end == 0.0){
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_calctrans: The ESO TEL AMBI IWV START/END science header keywords are both 0.0 which is not permitted by SCALE_PWV=auto");
            return NULL;
          }
          parameters->pwv_sci = 0.5*(iwv_start+iwv_end);
          if(parameters->h2o_col_mm != -99.9 && parameters->h2o_col_mm != 0.0){
              parameters->pwv_ratio = parameters->pwv_sci / parameters->h2o_col_mm ;
          }
      } else {
          /* IWV data not in headers, raise error */
          cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_calctrans: Missing ESO TEL AMBI IWV START/END science header keywords required by SCALE_PWV=auto");
          return NULL;
      }
  }
  //if the supplied value is a header keyword (i.e. it has characters in it) and not equal to none
  else if (!reti && strcmp(parameters->scale_pwv,"none")){// && parameters->scale_pwv has chars and != NONE) {
      if(cpl_propertylist_has(scientific_head_pri,parameters->scale_pwv)){
          parameters->pwv_sci = cpl_propertylist_get_double(scientific_head_pri,parameters->scale_pwv);
          if(parameters->h2o_col_mm != -99.9 && parameters->h2o_col_mm != 0.0){
              parameters->pwv_ratio = parameters->pwv_sci / parameters->h2o_col_mm ;
          }
      } else {
          cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_calctrans: Missing science header keyword %s required by SCALE_PWV",parameters->scale_pwv);
          return NULL;
      }
  }
  else if (reti && strcmp(parameters->scale_pwv,"none")) {
      /* numerical value supplied, convert to float and use it */
       parameters->pwv_sci = atof(parameters->scale_pwv);
       if(parameters->h2o_col_mm != -99.9 && parameters->h2o_col_mm != 0.0){
           parameters->pwv_ratio = parameters->pwv_sci / parameters->h2o_col_mm ;
       }
  } else if (strcmp(parameters->scale_pwv,"none")){
          cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "molecfit_calctrans: Error converting supplied value of SCALE_PWV: %s",parameters->scale_pwv);
          return NULL;
  }
  regfree(&regex);
  cpl_msg_info(cpl_func,"Value of parameters->pwv_sci: %e; parameters->h2o_col_mm (pwv_tell): %e; parameters->pwv_ratio: %e; parameters->scale_pwv: %s\n",parameters->pwv_sci,parameters->h2o_col_mm,parameters->pwv_ratio,parameters->scale_pwv);
  /* Rescale values in parameters->atm_parameters and write out to CALCTRANS_ATM_PARAMETERS.fits */
  if(n_ext > 0 && parameters->atm_parameters->v_ext[1].table != NULL && cpl_table_has_column(parameters->atm_parameters->v_ext[1].table,"H2O")){
      if(parameters->pwv_ratio == 1.0){
        cpl_msg_warning(cpl_func,"molecfit_calctrans: pwv_ratio is 1.0. No scaling based on PWV performed.\n");
      }
      nrows  = cpl_table_get_nrow(parameters->atm_parameters->v_ext[1].table);
      for (ridx=0;ridx<nrows;ridx++){
          double rval = cpl_table_get_double(parameters->atm_parameters->v_ext[1].table,"H2O",ridx,NULL)*parameters->pwv_ratio;
          cpl_table_set_double(parameters->atm_parameters->v_ext[1].table,"H2O",ridx,rval);
      }
  } else {
      cpl_msg_warning(cpl_func,"molecfit_calctrans: ATM_PARAMETERS missing H2O column. No scaling based on PWV performed.\n");
  }
  /*write out the new table to disk */
  /* This does not work on complex datasets - probably have to use:
  molecfit_save / molecfit_save_mf_results (as in molecfit_model) 
  instead of molecfit_fits_write */
  cpl_error_code err;
  if(strcmp(parameters->scale_pwv,"none")) {
      cpl_frameset *frameset_output = cpl_frameset_new();
      cpl_size n_frames = cpl_frameset_get_size(frameset);
      for (cpl_size frame = 0; frame < n_frames; frame++) {
             /* Check all inputs in frameset to get only the ATM_PARAMETERS frame*/
             cpl_frame  *frame_data = cpl_frameset_get_position(frameset, frame);
             const char *tag        = cpl_frame_get_tag(frame_data);
             if(!strcmp(tag,"ATM_PARAMETERS")){
                cpl_frameset_insert(frameset_output, cpl_frame_duplicate(frame_data));
                err = molecfit_fits_write(frameset,frameset_output,list,RECIPE_NAME,"ATM_PARAMETERS",parameters->atm_parameters,"CALCTRANS_ATM_PARAMETERS.fits");
                break;
             }
      }
  }

  /* Combine extensions in the FITS file ? */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CHIP_EXTENSIONS);
  parameters->chip_extensions = cpl_parameter_get_bool(p);
  cpl_propertylist_update_bool(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CHIP_EXTENSIONS, parameters->chip_extensions);


  /* Check status */
  if (!cpl_errorstate_is_equal(preState) || !(parameters->mf_config)) {
      /* Configuration failed */
      molecfit_calctrans_parameter_delete(parameters);
      return NULL;
  }

  return parameters;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief    Deallocate the given parameter configuration object and its contents
 *
 * @param    parameters       The parameter configuration variable in the recipe.
 */
/*----------------------------------------------------------------------------*/
static void molecfit_calctrans_parameter_delete(
    molecfit_calctrans_parameter *parameters)
{
  if (parameters) {

      if (parameters->mapping_kernel_table)                   cpl_table_delete(                        parameters->mapping_kernel_table);

      if (parameters->mapping_atmospheric_table)              cpl_table_delete(                        parameters->mapping_atmospheric_table);

      if (parameters->mapping_convolve_table)                 cpl_table_delete(                        parameters->mapping_convolve_table);

      if (parameters->atm_parameters)                         molecfit_fits_delete(                  parameters->atm_parameters);
      if (parameters->best_fit_parameters)                    molecfit_fits_delete(                  parameters->best_fit_parameters);

      if (parameters->telluriccorr_head) {
          for (cpl_size ext = 0; ext < parameters->n_ext; ext++) {
              if (parameters->telluriccorr_head[ext])             cpl_propertylist_delete(                 parameters->telluriccorr_head[ext]);
          }
          cpl_free(parameters->telluriccorr_head);
      }

      if (parameters->telluriccorr_data) {
          for (cpl_size ext = 0; ext < parameters->n_ext; ext++) {
              if (parameters->telluriccorr_data[ext])             cpl_table_delete(                        parameters->telluriccorr_data[ext]);
          }
          cpl_free(parameters->telluriccorr_data);
      }

      if (parameters->best_fit_parameters_table) {
          for (cpl_size ext = 0; ext < parameters->n_ext; ext++) {
              if (parameters->best_fit_parameters_table[ext]) cpl_table_delete(                        parameters->best_fit_parameters_table[ext]);
          }
          cpl_free(parameters->best_fit_parameters_table);
      }

      if (parameters->results_lblrtm) {
          for (cpl_size ext = 0; ext < parameters->n_ext; ext++) {
              if (parameters->results_lblrtm[ext])            mf_calctrans_lblrtm_results_delete(      parameters->results_lblrtm[ext]);
          }
          cpl_free(parameters->results_lblrtm);
      }

      if (parameters->results_convolution) {
          for (cpl_size ext = 0; ext < parameters->n_ext; ext++) {
              if (parameters->results_convolution[ext])       mf_calctrans_convolution_results_delete( parameters->results_convolution[ext]);
          }
          cpl_free(parameters->results_convolution);
      }

      if (parameters->mf_config)                              mf_configuration_delete(                              parameters->mf_config);

      if (parameters->pl)                                     cpl_propertylist_delete(                 parameters->pl);

      cpl_free(parameters);
  }
}
