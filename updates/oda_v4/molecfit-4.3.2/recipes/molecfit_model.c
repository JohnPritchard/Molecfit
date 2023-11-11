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

  cpl_boolean                      use_only_input_pri_ext;       /* If the user want to use only the STD_MODEL/SCIENCE_CALCTRANS/SCIENCE input FITS primary extension                          */
  int                              dflux_extension_data;         /* If use_only_input_pri_ext == CPL_TRUE, you can provide a extension as DFLUX (error flux)                                   */
  int                              mask_extension_data;          /* If use_only_input_pri_ext == CPL_TRUE, you can provide a extension as MASK                                                 */

  cpl_boolean                      use_input_kernel;             /* If the user provided, use or not the user kernel                                                                           */

  const char                       *list_molec;                  /* list_molec,           i.e. {H2O,CO,CO2,CH4,O2} string                                                                      */
  const char                       *fit_molec;                   /* fit_molec,            i.e. {1,  0, 1,  1,  0 } flag                                                                        */
  const char                       *rel_col;                     /* relative density,     i.e. {1.,1.,1.06,1.,1. } double                                                                      */
  cpl_table                        *molecules_table;             /* Molecules to fit. Contains 3 columns [MF_COL_LIST_MOLECULES, MF_COL_FIT_MOLECULES, MF_COL_REL_COL]                         */

  const char                       *wave_range_include;          /* wave_ranges_include,  i.e. {ini1,end1,...,iniN,endN} double                                                                */
  cpl_table                        *wave_ranges_include_table;   /* Wavelength ranges include to fit. Contains 2 columns [MF_COL_WAVE_RANGE_LOWER, MF_COL_WAVE_RANGE_UPPER]                    */

  const char                       *wave_range_exclude;          /* wave_ranges_exclude,  i.e. {ini1,end1,...,iniN,endN} double                                                                */
  cpl_table                        *wave_ranges_exclude_table;   /* Wavelength ranges exclude to fit. Contains 2 columns [MF_COL_WAVE_RANGE_LOWER, MF_COL_WAVE_RANGE_UPPER]                    */

  const char                       *pixel_range_exclude;         /* pixel_ranges_exclude, i.e. {ini1,end1,...,iniN,endN} int                                                                   */
  cpl_table                        *pixel_ranges_exclude_table;  /* Pixel ranges exclude      to fit. Contains 2 columns [MF_COL_WAVE_RANGE_LOWER, MF_COL_WAVE_RANGE_UPPER]                    */

  const char                       *mapping_kernel;              /* Mapping extensions : STD_MODEL/SCIENCE-MODEL_KERNEL_LIBRARY; i.e. {Data_ext1-Kernel_lib_ext, ..., Data_extN-Kernel_lib_ext} */
  cpl_table                        *mapping_kernel_table;        /* Mapping kernel extensions cpl_table. Contains 1 columns [KERNEL_LIBRARY_EXT]                                               */

  mf_configuration                 *mf_config;                   /* Molecfit configuration parameter                                                                                           */

  cpl_boolean                      chip_extensions;              /* If the user want to combine the input FITS extensions                                                                      */

  cpl_propertylist                 *pl;                          /* All input parameters to save in the output files                                                                           */

} molecfit_model_parameter;

/*----------------------------------------------------------------------------*/
/**
 *                              Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Fill the internal Molecfit configuration file */
static molecfit_model_parameter * molecfit_model_parameters(
    const cpl_parameterlist            *list,
    const cpl_propertylist             *raw_head_pri);

/* Clean variables allocated in the recipe */
static void molecfit_model_parameter_delete(
    molecfit_model_parameter *parameters);

static cpl_error_code molecfit_model_check_extensions_and_ranges(cpl_size chip, double min_wav, double max_wav, cpl_table* range);
    
static cpl_table* molecfit_model_amalgamate_ranges_table(
     cpl_table* inc_wranges_ext_table,
     const char* wave_include_str, 
     const char* map_regions_str,  
     const char* fit_continuum_str,
     const char* continuum_n_str,  
     const char* fit_wlc_str,
     int         nchips,
     double      wlg2mn); 
         
static int parsMFparsReadIValsFromString (const char* str, int*    vec, int vsize,int max_strlen) ;
static int parsMFparsReadDValsFromString (const char* str, double* vec, int vsize,int max_strlen);
static cpl_error_code molecfit_model_expert_mode(const cpl_frame* inp_frame,cpl_table* molecules, mf_configuration *config);
static int parseStr4Coeffs(const char *str, int *r_idx, int* c_idx, int* coef_idx);

/*----------------------------------------------------------------------------*/
/**
 *                          Static variables
 */
/*----------------------------------------------------------------------------*/

#define RECIPE_NAME      MOLECFIT_MODEL
#define CONTEXT          "molecfit."RECIPE_NAME

/*    "  -            (If used, optionally "MOLECFIT_MODEL_MAPPING_KERNEL"/"MOLECFIT_MAPPING_KERNEL" between "MOLECFIT_MODEL_KERNEL_LIBRARY"/"MOLECFIT_KERNEL_LIBRARY" to "MOLECFIT_STD_MODEL"/"MOLECFIT_SCIENCE" input).\n"*/
static char molecfit_model_description[] =
    "This recipe runs Molecfit on a 1D spectrum (direct or extracted form the input).\n"
    "The input data "MOLECFIT_STD_MODEL"/"MOLECFIT_SCIENCE_CALCTRANS"/"MOLECFIT_SCIENCE" could be:\n"
    "  - BINTABLE  - Table  (single or multi-extension, optionally with some NULL extensions)\n"
    "  - IMAGE(1D) - Vector (single or multi-extension, optionally with some NULL extensions)\n"
    "  - IMAGE(2D) - Image  (single or multi-extension, optionally with some NULL extensions)\n"
    "It is not mandatory that all the extensions contain data.\n"
    "Molecfit will be run  on all the extensions that contain data.\n"
    "The recipe also accepts as input a kernel library.\n"
    "If no kernel is provided, this is considered as a free parameter and it will be determined by Molecfit itself.\n"
    "\n"
    "Input files:\n"
    "  - Mandatory,  the input "MOLECFIT_STD_MODEL"/"MOLECFIT_SCIENCE_CALCTRANS"/"MOLECFIT_SCIENCE" data (Only one).\n"
    "  - Mandatory,  the fit molecules.\n"
    "  - Optionally, the inc/exc ranges/pixel.\n"
    "  - Optionally, the user kernel library.\n"
    "  - If the user kernel library is used, the parameter "MOLECFIT_MODEL_MAPPING_KERNEL" or "MOLECFIT_MAPPING_KERNEL" is a list of the extensions to map from the kernel library file "MOLECFIT_MODEL_KERNEL_LIBRARY" or "MOLECFIT_KERNEL_LIBRARY" to the extensions of the "MOLECFIT_STD_MODEL" or "MOLECFIT_SCIENCE" input file.\n"
    "  - Optionally, the GDAS profile.\n"
    "  - Optionally, the atmospheric standard profile.\n"
    "\n"
    "DO: Category                                  Required    Explanation\n"
    "-------------                                 --------    -----------\n"
    "  "MOLECFIT_STD_MODEL" or "MOLECFIT_SCIENCE_CALCTRANS" or "MOLECFIT_SCIENCE"     Y         The input spectrum           (Multi-extension FITS           format).\n"
    "  "MOLECFIT_MOLECULES"                                     Y      The molecules          table (         3-cols FITS BINTABLE  format) or the recipe parameters --"MOLECFIT_PARAMETER_LIST", --"MOLECFIT_PARAMETER_FIT", --"MOLECFIT_PARAMETER_RELATIVE_VALUE".\n"
    "  "MOLECFIT_WAVE_INCLUDE"                                  N      The wavelength include table (         2-cols FITS BINTABLE  format) or the recipe parameter  --"MOLECFIT_PARAMETER_WAVE_RANGE_INCLUDE".\n"
    "  "MOLECFIT_WAVE_EXCLUDE"                                  N      The wavelength exclude table (         2-cols FITS BINTABLE  format) or the recipe parameter  --"MOLECFIT_PARAMETER_WAVE_RANGE_EXCLUDE".\n"
    "  "MOLECFIT_PIXEL_EXCLUDE"                                 N      The pixel      exclude table (         2-cols FITS BINTABLE  format) or the recipe parameter  --"MOLECFIT_PARAMETER_PIXEL_RANGE_EXCLUDE".\n"
    "  "MOLECFIT_MODEL_KERNEL_LIBRARY" or "MOLECFIT_KERNEL_LIBRARY"        N         The kernel library           (Multi-extension FITS IMAGE(2D) format).\n"
    "  "MOLECFIT_MODEL_MAPPING_KERNEL" or "MOLECFIT_MAPPING_KERNEL"        N      Kernel mapping table         (         1-col  FITS BINTABLE  format) or the recipe parameter  --"MOLECFIT_PARAMETER_MODEL_MAPPING_KERNEL".\n"
    "  "MOLECFIT_GDAS"                                          N         User define GDAS profile     (                FITS BINTABLE  format).\n"
    "  "MOLECFIT_ATM_PROFILE_STANDARD"                          N         User define GDAS profile     (                FITS BINTABLE  format).\n"
    "\n"
    "Notes:\n"
    "  - "MOLECFIT_MOLECULES"                            BINTABLE columns ["MF_COL_LIST_MOLECULES", "MF_COL_FIT_MOLECULES", "MF_COL_REL_COL"]\n"
    "  - "MOLECFIT_WAVE_INCLUDE"                         BINTABLE columns ["MF_COL_WAVE_RANGE_LOWER", "MF_COL_WAVE_RANGE_UPPER"]\n"
    "  - "MOLECFIT_WAVE_EXCLUDE"                         BINTABLE columns ["MF_COL_WAVE_RANGE_LOWER", "MF_COL_WAVE_RANGE_UPPER"]\n"
    "  - "MOLECFIT_PIXEL_EXCLUDE"                        BINTABLE columns ["MF_COL_WAVE_RANGE_LOWER", "MF_COL_WAVE_RANGE_UPPER"]\n"
    "  - "MOLECFIT_MODEL_MAPPING_KERNEL" or "MOLECFIT_MAPPING_KERNEL"  BINTABLE column  ["MOLECFIT_MAPPING_KERNEL_EXT"]\n"
    "\n"
    "\n"
    "Output files:\n"
    "\n"
    "DO: Category                      Explanation\n"
    "-------------                     -----------\n"
    "  "MOLECFIT_MODEL_CHIPS_COMBINED"            Input file chips combined in molecfit format    (                FITS BINTABLE format).\n"
    "  "MOLECFIT_GDAS"                            Used GDAS profile                               (                FITS BINTABLE format).\n"
    "  "MOLECFIT_GDAS_BEFORE"                     If ESO DB GDAS is used, file before the MJD-OBS (                FITS BINTABLE format).\n"
    "  "MOLECFIT_GDAS_AFTER"                      If ESO DB GDAS is used, file after  the MJD-OBS (                FITS BINTABLE format).\n"
    "  "MOLECFIT_MODEL_MOLECULES"                 The molecules          table                    (         3-cols FITS BINTABLE format).\n"
    "  "MOLECFIT_WAVE_INCLUDE"                    The wavelength include table                    (         2-cols FITS BINTABLE format).\n"
    "  "MOLECFIT_WAVE_EXCLUDE"                    The wavelength exclude table                    (         2-cols FITS BINTABLE format).\n"
    "  "MOLECFIT_PIXEL_EXCLUDE"                   The pixel      exclude table                    (         2-cols FITS BINTABLE format).\n"
    "  "MOLECFIT_ATM_PROFILE_STANDARD"            Atmospheric profile standard - Initial          (                FITS BINTABLE format).\n"
    "  "MOLECFIT_ATM_PROFILE_COMBINED"            Atmospheric profile combined with GDAS.         (                FITS BINTABLE format).\n"
    "  "MOLECFIT_ATM_PARAMETERS"                  Atmospheric parameters                          (Multi-extension FITS BINTABLE format).\n"
    "  "MOLECFIT_BEST_FIT_PARAMETERS"             Best fitting parameters                         (Multi-extension FITS BINTABLE format).\n"
    "  "MOLECFIT_BEST_FIT_MODEL"                  Best fit model and intermediate products        (Multi-extension FITS BINTABLE format).\n"
    "  "MOLECFIT_MODEL_KERNEL_LIBRARY"            The kernel library processed in mf_model(...) Molecfit function: Normalize and Resample.\n"
    "\n";

/* Standard CPL recipe definition */
cpl_recipe_define(        molecfit_model,
                          MOLECFIT_BINARY_VERSION,
                          "Jose A. Escartin",
                          "https://support.eso.org/",
                          "2019",
                          "Extract 1D spectrum's from a multi-extension input raw FITS files and compute an atmospheric model for each.",
                          molecfit_model_description);


/*----------------------------------------------------------------------------*/
/**
 * @defgroup molecfit_model  It runs Molecfit on a generic input spectrum file to compute an atmospheric model.
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
static int molecfit_model(
    cpl_frameset             *frameset,
    const cpl_parameterlist  *parlist)
{


  /* Check parameters and get initial errorstate */
  cpl_error_ensure(frameset && parlist, CPL_ERROR_NULL_INPUT, return CPL_ERROR_NULL_INPUT, "NULL input : frameset and/or parlist");
  cpl_errorstate initial_errorstate = cpl_errorstate_get();
  cpl_error_code err = CPL_ERROR_NONE;
  
  /* Check that all input string lists are in valid formats */  
  cpl_error_code list_err=CPL_ERROR_NONE;
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MOLECFIT_MODEL_MAPPING_KERNEL,	N0LIST );
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MOLECFIT_PARAMETER_LIST,           SLIST );
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MOLECFIT_PARAMETER_FIT,		 BLIST );
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MOLECFIT_PARAMETER_RELATIVE_VALUE, DLIST );
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MOLECFIT_WAVE_INCLUDE,		 DRANGE);
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MOLECFIT_WAVE_EXCLUDE,		 DRANGE);
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MOLECFIT_PIXEL_EXCLUDE,		 NRANGE);
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MF_PARAMETERS_FIT_CONTINUUM,	 BLIST_OLD );
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MF_PARAMETERS_CONTINUUM_N,	 N0LIST);
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MF_PARAMETERS_MAP_REGIONS_TO_CHIP, NLIST );
  if (!list_err) list_err = molecfit_config_strlstchk(parlist,MF_PARAMETERS_MAP_REGIONS_TO_CHIP, NLIST );
  if (!list_err) list_err = molecfit_config_stroptchk(parlist,MF_PARAMETERS_WAVELENGTH_FRAME, "VAC,AIR,VAC_RV,AIR_RV");
  if (list_err!=CPL_ERROR_NONE) {
      err=CPL_ERROR_ILLEGAL_INPUT;
      cpl_error_set_message(cpl_func,err,"Invalid string list parameters");
      return err;
  } 


  if (!err) {
      /* Check frameset TAGS */
      cpl_msg_info(cpl_func, "Check frameset (SOF tags) ...");
      err = molecfit_check_and_set_groups(frameset);
  }


  /* Test kernel mapping */
  if (!err) {
      cpl_error_code map_err=CPL_ERROR_NONE;
      map_err=molecfit_config_mappingchk(frameset,parlist,"KERNEL_LIBRARY",
                                                          "MAPPING_KERNEL",
							  "KERNEL_LIBRARY_EXT",
							  MOLECFIT_MODEL_MAPPING_KERNEL,
							  "NULL");
      if (map_err) cpl_msg_error(cpl_func,"Mapping Errors detected cannot proceed further");
      
  }
  
  /* Check mandatory TAGS/Parameters */
  const cpl_frame *input_frame = NULL;
  cpl_boolean science_data_is_table=CPL_FALSE;
  if (!err) {

      /* STD_MODEL/SCIENCE_CALCTRANS/SCIENCE */
      cpl_errorstate pre_state = cpl_errorstate_get();
      input_frame = cpl_frameset_find(frameset, MOLECFIT_STD_MODEL);
      if (!input_frame) {
          cpl_errorstate_set(pre_state);
          input_frame = cpl_frameset_find(frameset, MOLECFIT_SCIENCE_CALCTRANS);
          if (!input_frame) {
              cpl_errorstate_set(pre_state);
              input_frame = cpl_frameset_find(frameset, MOLECFIT_SCIENCE);
              if (!input_frame) {
                  err = cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                              "STD_MODEL/SCIENCE_CALCTRANS/SCIENCE data not found in frameset!");
              }
          }
      }
      
      /* Checking the SCIENCE file data type eg is it in table form */
      science_data_is_table=molecfit_config_chk_science_is_table(input_frame);
      if ( science_data_is_table) cpl_msg_info(cpl_func,"Science frame is in table data format");
      if (!science_data_is_table) cpl_msg_info(cpl_func,"Science frame is not in table data format");

      /* MOLECULES */
      if (!err) {
          const cpl_frame *input_molecules = cpl_frameset_find_const(frameset, MOLECFIT_MOLECULES);
          if (!input_molecules) {
	      /* There is no external molecule table so check on the string lists*/
	      cpl_error_code molstr_err=CPL_ERROR_NONE;
	      molstr_err=molecfit_config_molecule_strs_chk(
	         cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, MOLECFIT_PARAMETER_LIST          )),
	         cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, MOLECFIT_PARAMETER_FIT           )),
		 cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, MOLECFIT_PARAMETER_RELATIVE_VALUE)));

              if( molstr_err) {
	          err = CPL_ERROR_ILLEGAL_INPUT;
		  cpl_error_set_message(cpl_func, err,"Bad molecules definition!");
              }
	  }
      }

  }


  /* Load TAG = STD_MODEL/SCIENCE_CALCTRANS/SCIENCE */
  molecfit_fits *data = NULL;
  if (!err) {
      const char *data_file = cpl_frame_get_filename(input_frame);

      data = molecfit_fits_load(data_file, CPL_FALSE);

      if (!data) err = CPL_ERROR_ILLEGAL_INPUT;
      else           err = cpl_error_get_code();
  }

  /* If the data file has been loaded successfully check the header keywords */
  if (!err) {
  
       /* Check that the following keywords (if user specified) do exists in the fits data header */
	cpl_error_code kwd_err=CPL_ERROR_NONE;
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_OBSERVATORY_ERF_RV_KEYWORD  );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_OBSERVING_DATE_KEYWORD      );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_UTC_KEYWORD		        );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_TELESCOPE_ANGLE_KEYWORD     );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_RELATIVE_HUMIDITY_KEYWORD   );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_PRESSURE_KEYWORD	        );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_TEMPERATURE_KEYWORD	        );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_MIRROR_TEMPERATURE_KEYWORD  );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_ELEVATION_KEYWORD	        );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_LONGITUDE_KEYWORD	        );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_LATITUDE_KEYWORD	        );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_SLIT_WIDTH_KEYWORD	        );
	kwd_err+=mf_config_check_keyword(data,parlist, MF_PARAMETERS_PIXEL_SCALE_KEYWORD	        );

	if (kwd_err) {
            cpl_msg_error(cpl_func,"Header Keyword Names specifed in the .rc file not found in the target science frame!");
            err=CPL_ERROR_ILLEGAL_INPUT;
	}
  }

  /* Recipe Parameters : Need raw_header_primary */
  molecfit_model_parameter *parameters = NULL;
  
  /* Note n physical chips is defined via a combination of the science data extensions */
  /*  and paremeter flags                                                              */
  int n_chip_subranges=-1; /* -1 Implies that chip subranges are not being used */

  
  if (!err) {

      /* Get recipe parameters and update the molecfit default configuration */
      cpl_msg_info(cpl_func, "Load '%s' recipe parameters ...", MOLECFIT_MODEL);
      parameters = molecfit_model_parameters(parlist, data->v_ext[0].header);
      cpl_boolean use_only_input_pri_ext; 
      use_only_input_pri_ext = cpl_parameter_get_bool(cpl_parameterlist_find_const(parlist,MOLECFIT_PARAMETER_USE_ONLY_INPUT_PRI));  


      if (!parameters)  err = CPL_ERROR_ILLEGAL_INPUT;
      else              err = cpl_error_get_code();

      /* ---------------------------------------------------------------------*/
      /* Determine if using CHIPS mechanism and if so how many chips are there*/
      /* -------------------------------------------------------------------- */
      if (!err) {
                
          /* If USECHIPS flag is true then extensions are to be treated as "CHIP" sub ranges of 
              a full spectrum and chip specfic parameters must be checked.        */
          /* Use chips mechanism only if chip_extension flag is true and use_oly_pri is false*/
          if (!parameters->use_only_input_pri_ext && parameters->chip_extensions) {
              /* No of chip subregions is the number of extensions (not counting the primary) */
              n_chip_subranges= data->n_ext-1; 
          }
          
          if (n_chip_subranges>0) {
              cpl_msg_info(cpl_func,"Use CHIPS flag set to true");
              cpl_msg_info(cpl_func,"Use only primary data flag set to false.");
              cpl_msg_info(cpl_func,"Will treat the %d extensions of the data file as subregions of full spectrum",n_chip_subranges);
             
          }/* end if n_chip_subranges>0*/  
         
      }    

      /* ----------------------------------------------------------------------------------------------------*/
      /* Convert the imported fits file data which could be in one of many different formats into cpl tables */
      /* that molecfit/telluriccorr can process                                                              */
      /* ----------------------------------------------------------------------------------------------------*/
      /* Note: the imported fits data is already contained in the data structure *data    */ 
      /* Note: the new cpl_tables will be added to the structure *data                    */
      /* Note in the event of chip_extension=true, all extension data in the fits file is merged into a single cpl_table  */
      
      /* First check that all required column names exists (This test is only if USE_ONLY_PRIMARY_DATA=FALSE)*/
      if (!use_only_input_pri_ext && science_data_is_table) {
          cpl_error_code col_err=CPL_ERROR_NONE;

          col_err=mf_config_chk_column(data,parameters->mf_config->parameters->inputs.column_lambda);
          col_err=mf_config_chk_column(data,parameters->mf_config->parameters->inputs.column_flux  );
          col_err=mf_config_chk_column(data,parameters->mf_config->parameters->inputs.column_dflux );
          col_err=mf_config_chk_column(data,parameters->mf_config->parameters->inputs.column_mask  );
          if (col_err) err=CPL_ERROR_ILLEGAL_INPUT;
      }
    
 
      if (!err) {
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
              cpl_msg_info(cpl_func, "Save combined multi-extension molecfit_spectrum input FITS file ('%s') in ext=%lld ...", MOLECFIT_MODEL_CHIPS_COMBINED, index_combined);

              err += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_MODEL_CHIPS_COMBINED, NULL);
              if (!err) err = molecfit_save_mf_results(data->v_ext[index_combined].spectrum_head, MOLECFIT_MODEL_CHIPS_COMBINED, CPL_TRUE, NULL, data->v_ext[index_combined].spectrum_data, NULL);
          }
      }

  }

  /* Load TAG = MOLECULES */
  cpl_table *molecules = NULL;
  if (!err) {

      /* Get MOLECUTES table */
      cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_MOLECULES);
      if (parameters->molecules_table) {
          molecules = cpl_table_duplicate(parameters->molecules_table);
      } else {
          cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_MOLECULES);
          molecules = molecfit_load_unique_table(frameset, MOLECFIT_MOLECULES);
      }


      if (!molecules) err = CPL_ERROR_ILLEGAL_INPUT;
      else {
          if ((int)cpl_table_get_column_max(molecules, "FIT_MOLEC") < 1){
                  cpl_msg_warning(cpl_func, "No molecules to fit !");
          }
              err = cpl_error_get_code();
      }


     /* Now check that selected Reference Atmosphere file has the required molecules */
      cpl_error_code rat_err=CPL_ERROR_NONE;
      rat_err=mf_config_chk_ref_atmos(parameters->mf_config->parameters->directories.telluriccorr_data_path,
                                      parameters->mf_config->parameters->atmospheric.ref_atm,
				      molecules);        
      if (rat_err) {
          cpl_msg_info(cpl_func,"Bad Specifications of Molecules and Reference Atmosphere");
	  err=CPL_ERROR_ILLEGAL_INPUT;
      }
 
  }


  /* Load TAG = WAVE_INCLUDE, WAVE_EXCLUDE, PIXEL_EXCLUDE */
  cpl_table *inc_wranges = NULL;
  cpl_table *exc_wranges = NULL;
  cpl_table *exc_pranges = NULL;
  
  
  /* ------------------------------------------------------------*/
  /* Sort the wavlength ranges and any range specific parameters */
  /* ------------------------------------------------------------*/

  /* Obtain all relevant parameters as defined in the .rc file */
  const char * wave_include_str  = cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, "WAVE_INCLUDE"));
  const char * map_regions_str   = cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, "MAP_REGIONS_TO_CHIP"));
  const char * fit_continuum_str = cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, "FIT_CONTINUUM"      ));
  const char * continuum_n_str   = cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, "CONTINUUM_N"        ));
  const char * fit_wlc_str       = cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, "FIT_WLC"            ));
  double       wlg2mn            = cpl_parameter_get_double(cpl_parameterlist_find_const(parlist, "WLG_TO_MICRON"      ));
  cpl_table *inc_wranges_ext_table = NULL;

  if (!err) {
   
      /* Load up on any external fits table or previously defined table which may have some or even all required data*/
      if (parameters->wave_ranges_include_table ) {
          inc_wranges_ext_table = cpl_table_duplicate(parameters->wave_ranges_include_table);
     } else {
          cpl_msg_info (cpl_func, "Loading %s cpl_table as external", MOLECFIT_WAVE_INCLUDE);
          inc_wranges_ext_table = molecfit_load_unique_table(frameset, MOLECFIT_WAVE_INCLUDE);
     }
 
      /* Merge all the .rc string parameters and table data into a single inc_wranges table*/
      int    nchips     = (1>n_chip_subranges)?1:n_chip_subranges;
      inc_wranges = molecfit_model_amalgamate_ranges_table(inc_wranges_ext_table,
                                                                     wave_include_str, 
                                                                     map_regions_str,  
                                                                     fit_continuum_str,
                                                                     continuum_n_str,  
                                                                     fit_wlc_str,
								     nchips,
								     wlg2mn); 

       /* Free up the external table*/
      if (inc_wranges_ext_table!=NULL) cpl_free(inc_wranges_ext_table);

      if (inc_wranges==NULL) {
          cpl_msg_error(cpl_func,"Unable to construct full wavelength range table!");
          err=CPL_ERROR_NULL_INPUT;
      } else {
          err = cpl_error_get_code();
      }
      
  }
  

 
  
  /* -----------------------------------------------------------------------------*/
  /* Define the initial coefficient values to use for the WLC and Continuum models*/
  /* NOTE: May be overriden by EXPERT MODE lower down */
  /* -----------------------------------------------------------------------------*/

  if (!err) {
  
      /* The number of ranges is now defined in the inc_wranges table */
      cpl_size n_ranges=cpl_table_get_nrow(inc_wranges);
      
      /* Default initial shift value shuld be define in the parameters structure */
      double const_value = parameters->mf_config->parameters->fitting.fit_continuum.const_val;
      
      /* Go through each range and define the initial poly as a shift by the const_value*/
      int max_porder=0;
      for (int i=0;i<n_ranges;i++) {
      
          /* Get the polynomial order for this range and set all coeffs to zero*/
          int porder=cpl_table_get_int(inc_wranges,MF_COL_WAVE_RANGE_CONT_ORDER,i,NULL);
          if(porder>max_porder) max_porder=porder;
          for (int j=0;j<=porder;j++) {
              parameters->mf_config->parameters->fitting.cont_coeffs[i][j]=0.0;
          }
          
          /* Now set as y=x+shift*/
          parameters->mf_config->parameters->fitting.cont_coeffs[i][0]=const_value;
          parameters->mf_config->parameters->fitting.cont_coeffs[i][1]=1.0;
      }
      
      /* Originally there was only one poly order for all ranges and it was defined in: */
      /*       parameters->mf_config->parameters->fitting.fit_continuum.n               */
      /* This parameter is still used to allocate array sizes for the coeffs in the     */
      /* the range tables passed to molecfit execution. It would appear that array      */
      /* columns in tables must have the same dimension so we change the meaning of     */
      /* this parameter to be the maximum poly order of all the ranges                  */
      parameters->mf_config->parameters->fitting.fit_continuum.n=max_porder;
      
      
      /* Default initial shift value shuld be define in the parameters structure */
      int    n_chips     = (1>n_chip_subranges)?1:n_chip_subranges;
      int    porder      = parameters->mf_config->parameters->fitting.fit_wavelenght.n;
      double const_shift= parameters->mf_config->parameters->fitting.fit_wavelenght.const_val;
      
      cpl_msg_info(cpl_func,"N chips=%d, porder = %d const_shift=%f",n_chips,porder,const_shift);
       
      /* Go through each chip and define the initial poly as a shift by shft)value*/
      for (int i=0;i<n_chips;i++) {
      
          for (int j=0;j<=porder;j++) {
              parameters->mf_config->parameters->fitting.wlc_coeffs[i][j]=0.0;
          }
          
          /* Now set as y = x+shift */
          parameters->mf_config->parameters->fitting.wlc_coeffs[i][0]=const_shift;
          parameters->mf_config->parameters->fitting.wlc_coeffs[i][1]=1.0;
      }
             
  }
 
 
  
  
  /* ------------------------------------------------------------*/
  /* Temporaray solution to fit chips and fit ranges flag mechanism */
  /* ------------------------------------------------------------*/

  if (!err) {

      /* 
         We are still using the halfway house code which overrides wlc fit flags and 
         continuum fit flags in subroutine:
                                 mf_spectrum_replace_coefficientstine
         with parameters:
                          params->config->fitting.fit_ranges[i]
                          params->config->fitting.fit_chips[i]
      */   
      
      /* Sort out the fit_ranges parameter */
      int n_ranges  = cpl_table_get_nrow(inc_wranges);  
      for (int i=0;i<n_ranges;i++) {
          cpl_boolean flag=CPL_TRUE;
          int ival = cpl_table_get_int(inc_wranges,MF_COL_WAVE_RANGE_CONT_FIT,i,NULL);
          if (ival==0) flag=CPL_FALSE;
          parameters->mf_config->parameters->fitting.fit_ranges[i]=flag;
          int order = cpl_table_get_int(inc_wranges,MF_COL_WAVE_RANGE_CONT_ORDER,i,NULL);
          parameters->mf_config->parameters->fitting.cont_poly_order[i]=order;
      }
      
      /* Sort out the fit chips parameter */
      int n_chips=(1>n_chip_subranges)?1:n_chip_subranges;

      cpl_boolean fitv[n_chips]; 
      for (int i=0;i<n_chips;i++) fitv[i]=CPL_FALSE; 

      for (int i=0;i<n_ranges;i++) {
          int ival    = cpl_table_get_int(inc_wranges,MF_COL_WAVE_RANGE_WLC_FIT,i,NULL);
          if (ival==0) continue;
          int chip_no = cpl_table_get_int(inc_wranges,MF_COL_WAVE_RANGE_MAP2CHIP,i,NULL);
          int chip_idx=chip_no-1;
          if(chip_idx<0 || chip_idx>n_chips) continue;
          fitv[chip_idx]=CPL_TRUE;          
      }
      
      for (int i=0;i<n_chips;i++) {
          parameters->mf_config->parameters->fitting.fit_chips[i]=fitv[i];
          cpl_msg_info(cpl_func,"Sorting chip fits chip %d fit flag = %d",i+1,fitv[i]);
      }
   
   }
 
     
  /* ------------------------------------------------- ------------*/
  /* END Sort the wavlength ranges and any associated RANGE flagss */
  /* ------------------------------------------------- ------------*/

  /* ------------------------------------------------- ------------*/
  /* Sort any conflicts between the RANGES and the CHIPS           */
  /* ------------------------------------------------- ------------*/
  if (!err) {
      /* 
         Only need to check if flags:
               use_chips flag        = true
               use only primary flag = false 
      */
      cpl_msg_info(cpl_func,"Checking for any conflicts between ranges and chips");
      if (!parameters->use_only_input_pri_ext && parameters->chip_extensions) {
       /*
          err=molecfit_model_check_chips_and_ranges(parameters,data);
      */
      }
  }  

  
  /* ------------------------------------------------- ------------*/
  /* END Sort any conflicts between the RANGES and the CHIPS       */
  /* ------------------------------------------------- ------------*/


  if (!err) {
                     
      if (parameters->wave_ranges_exclude_table ) {
          exc_wranges = cpl_table_duplicate(parameters->wave_ranges_exclude_table);
      } else {
          cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_WAVE_EXCLUDE);
          exc_wranges = molecfit_load_unique_table(frameset, MOLECFIT_WAVE_EXCLUDE);
      }

      if (parameters->pixel_ranges_exclude_table) {
          exc_pranges = cpl_table_duplicate(parameters->pixel_ranges_exclude_table);
      } else {
          cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_PIXEL_EXCLUDE);
          exc_pranges = molecfit_load_unique_table(frameset, MOLECFIT_PIXEL_EXCLUDE);
      }

      err = cpl_error_get_code();
  }


  /* ------------*/ 
  /* EXPERT MODE */
  /* ------------*/ 
  if (!err) {
  
      
      /* Check if parameter EXPERT_MODE=TRUE*/
      cpl_boolean expert_mode_flag=cpl_parameter_get_bool(cpl_parameterlist_find_const(parlist, "EXPERT_MODE"));
      /*
      const char * expert_mode_str  = cpl_parameter_get_string(cpl_parameterlist_find_const(parlist, "EXPERT_MODE"))
      if (expert_mode_str) {
          if (expert_mode_str[0]=='T' || expert_mode_str[0]=='t') expert_mode_flag=CPL_TRUE;
      }
      */     
      parameters->mf_config->parameters->fitting.expert_mode=expert_mode_flag;    
      
      if (expert_mode_flag) {
          cpl_msg_info(cpl_func,"    ========================");
          cpl_msg_info(cpl_func,"     EXPERT MODE REQUESTED!");
          cpl_msg_info(cpl_func,"    ========================");
          

          input_frame = cpl_frameset_find_const(frameset,MOLECFIT_MODEL_INIT_FIT_PARAMETERS);
          if (!input_frame) {
              cpl_msg_info(cpl_func,"No Expert Mode Initial Fit Values File Found!");
          } else {
              cpl_msg_info(cpl_func,"Expert Mode Initial Fit Values File Located");
              err= molecfit_model_expert_mode(input_frame,molecules,parameters->mf_config);
          }
          cpl_msg_info(cpl_func,"    ===============================");
          cpl_msg_info(cpl_func,"     END OF EXPERT MODE PROCESSING ");
          cpl_msg_info(cpl_func,"    ===============================");
      } /* endif expert_mode_flag*/
  }



  /* Save MOLECULES table */
  if (!err) {
      err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_MODEL_MOLECULES, NULL);
      if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_MODEL_MOLECULES, CPL_TRUE, NULL, molecules, NULL);
  }


  /* Save WV_INCLUDE RANGE's table */
  if (!err && inc_wranges) {
      err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_WAVE_INCLUDE, NULL);

      if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_WAVE_INCLUDE, CPL_TRUE, NULL, inc_wranges, NULL);
/*
      if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_WAVE_INCLUDE, CPL_TRUE, NULL, inc_wranges_new_table, NULL);
*/ 
  }

  /* Save WV_EXCLUDE RANGE's table */
  if (!err && exc_wranges) {
      err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_WAVE_EXCLUDE, NULL);
      if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_WAVE_EXCLUDE, CPL_TRUE, NULL, exc_wranges, NULL);
  }

  /* Save PIX_INCLUDE RANGE's table */
  if (!err && exc_pranges) {
      err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_PIXEL_EXCLUDE, NULL);
      if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_PIXEL_EXCLUDE, CPL_TRUE, NULL, exc_pranges, NULL);
  }


  /* Load TAG : MODEL_KERNEL_LIBRARY/KERNEL_LIBRARY */
  molecfit_fits *kernel_data    = NULL;
  cpl_table     *mapping_kernel = NULL;
  if (!err) {

      kernel_data = molecfit_load_kernel_tag(frameset, MOLECFIT_MODEL_KERNEL_LIBRARY, parameters->use_input_kernel, parameters->mf_config->parameters,CPL_TRUE);
      err = cpl_error_get_code();

      if (kernel_data && !err) {

          if (parameters->mapping_kernel_table) {

              /* Mapping in the recipe parameters */
              mapping_kernel = cpl_table_duplicate(parameters->mapping_kernel_table);

          } else {

              /* Mapping in the static_calib input FITS files */
              cpl_errorstate pre_state = cpl_errorstate_get();
              const cpl_frame *input_mapping_kernel = cpl_frameset_find(frameset, MOLECFIT_MODEL_MAPPING_KERNEL);
              if (input_mapping_kernel) {
                  cpl_msg_info (cpl_func, "Loading %s cpl_table", MOLECFIT_MODEL_MAPPING_KERNEL);
                  mapping_kernel = molecfit_load_unique_table(frameset, MOLECFIT_MODEL_MAPPING_KERNEL);
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

          if (!mapping_kernel && !err) err = CPL_ERROR_ILLEGAL_INPUT;
          else                         err = cpl_error_get_code();
      }

      if (!err) cpl_msg_info(cpl_func, "Using the default molecfit kernels -> With the values inside BEST_FIT_PARMS!");
  }

  /* Load TAG = GDAS */
  cpl_boolean gdas_write = CPL_FALSE;
  cpl_table   *gdas_user = NULL;
  if (!err) {
      char * gdas_prof=parameters->mf_config->parameters->atmospheric.gdas_prof;
      cpl_error_code gdas_err;
      gdas_err = mf_config_chk_gdas_prof(gdas_prof);
      if (gdas_err) cpl_msg_info(cpl_func,"Invalid value for GDAS_PROF: %s",gdas_prof);

      //only search the SOF if GDAS_PROF = NULL (default is auto)
      if((strcmp(gdas_prof,"null") == 0 || strcmp(gdas_prof,"NULL") == 0)){
          if(cpl_frameset_find_const(frameset,MOLECFIT_GDAS) && !gdas_err){
            gdas_user = molecfit_load_unique_table(frameset, MOLECFIT_GDAS);
            if (gdas_user) cpl_msg_info (cpl_func, "Loaded %s cpl_table", MOLECFIT_GDAS);
            err = cpl_error_get_code();
          } else {
            /* Use the average if possible */
            cpl_msg_info(cpl_func,"GDAS profile missing from SOF: Reverting to GDAS_PROF=auto behaviour");
          //parameters->mf_config->parameters->atmospheric.gdas_prof = "auto";
          cpl_free(parameters->mf_config->parameters->atmospheric.gdas_prof);
          parameters->mf_config->parameters->atmospheric.gdas_prof = cpl_strdup("auto");
          }
      } 
      //if gdas_prof is a file...
      else if(!(!strcmp(gdas_prof,"auto") || !strcmp(gdas_prof,"AUTO") ||
           !strcmp(gdas_prof,"none") || !strcmp(gdas_prof,"NONE") ||
           !strcmp(gdas_prof,"null") || !strcmp(gdas_prof,"NULL"))){
          if(!gdas_err){
              cpl_msg_info(cpl_func,"Loading GDAS_PROF from file: %s",gdas_prof);
              gdas_user = cpl_table_load(gdas_prof,1,0);
              if(!gdas_user){
                /* Error loading file, even if gdas_prof suspects it is OK*/
                err = CPL_ERROR_ILLEGAL_INPUT;
                cpl_msg_error(cpl_func,"Error loading file specified by GDAS_PROF: %s",gdas_prof);
              }
          } else {
              /* Error parsing gdas_prof */
              err = CPL_ERROR_ILLEGAL_INPUT;
              cpl_msg_error(cpl_func,"Error loading file specified by GDAS_PROF: %s",gdas_prof);
          }
      } 
      /*this is no longer needed, since we force use of the average files in share/molecfit/data/profiles/lib/' 
        corresponding to the month of the observation in mf_gdas inside telluriccorr */
      /*else if(strcmp(gdas_prof,"none") == 0 || strcmp(gdas_prof,"NONE") == 0){
          //Otherwise do not check the SOF and instead try and use an auto detected one here 
          cpl_free(parameters->mf_config->parameters->atmospheric.gdas_prof);
          parameters->mf_config->parameters->atmospheric.gdas_prof = cpl_strdup("auto");
      }*/
  }

  /* Load TAG = ATM_PROFILE_STANDARD */
  cpl_boolean atm_profile_standard_write = CPL_FALSE;
  cpl_table   *atm_profile_standard      = NULL;
  cpl_table   *atm_profile_combined      = NULL;
  if (!err) {
      atm_profile_standard = molecfit_load_unique_table(frameset, MOLECFIT_ATM_PROFILE_STANDARD);
      if (atm_profile_standard) cpl_msg_info (cpl_func, "Loaded %s cpl_table", MOLECFIT_ATM_PROFILE_STANDARD);
      err = cpl_error_get_code();
  }
  
  /* Check Request LNFL DB is valid */
  if (!err) {
      cpl_error_code lnfl_db_err=CPL_ERROR_NONE;
      lnfl_db_err = mf_config_chk_lnfl_line_db(
                    parameters->mf_config->parameters->directories.telluriccorr_data_path,
		    parameters->mf_config->lnfl->line_db);
      if(lnfl_db_err) {
          cpl_msg_error(cpl_func,"Invalid parameter value for LNFL_LINE_DB: %s",parameters->mf_config->lnfl->line_db);
	  err = CPL_ERROR_ILLEGAL_INPUT;
	  cpl_error_set_message(cpl_func, err,"molecfit_model .rc parameter checks failed!");
      }
  }

  /*** Save generic output files */
  if (!err) {

      cpl_msg_info(cpl_func, "Save generic multi-extension output FITS file ('%s','%s','%s') ...", MOLECFIT_ATM_PARAMETERS, MOLECFIT_BEST_FIT_MODEL, MOLECFIT_BEST_FIT_PARAMETERS);

      err     += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_ATM_PARAMETERS,       NULL);
      err     += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_BEST_FIT_PARAMETERS,  NULL);
      err     += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_BEST_FIT_MODEL,       NULL);
      err     += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, "MOLECFIT_DATA",       NULL);


      if (kernel_data) {
          err += molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_MODEL_KERNEL_LIBRARY, NULL);
      }
  }


  int null;

  /*** Recipe execution ***/
  if (!err) {

      /* Execution extensions with spectral data in table data->v_ext[?].spectrum data  */
      /*  NOTE: If chip_extensions = TRUE then all spectra data from extensions with    */
      /*  index > 0 have been merged into the table data->v_ext[1].spectrum data        */

      cpl_size n_ext;
      if (     parameters->use_only_input_pri_ext) n_ext = 1;
      else if (parameters->chip_extensions    ) n_ext = data->v_ext[0].spectrum_data ? 1 : 2;
      else n_ext = data->n_ext;
      

      if (!parameters->use_only_input_pri_ext && n_ext == 1){
              cpl_msg_warning(cpl_func,"No data extensions found. USE_ONLY_INPUT_PRIMARY_DATA=FALSE and number of extensions ==  1");
      }

      
      /* --------------- */
      /* Execution Loops */
      /* --------------- */
      /* Iterate over all extensions in the data object that contain spectrum_data tables.*/
      for (cpl_size ext = 0; ext < n_ext && !err; ext++) {

          mf_model_results *results = NULL;

          /* Create input molecfit spec format */
          if (data->v_ext[ext].spectrum_data) {

              /* Get kernel */
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

              /* Get STD_MODEL/SCIENCE_CALCTRANS/SCIENCE extension header/table data */
              cpl_msg_info(cpl_func, "Load spectrum, ext = %lld ...", ext);
              cpl_table *molecfit_data = mf_spectrum_create(parameters->mf_config->parameters, data->v_ext[ext].spectrum_data);


              
              /* If the chip_extension flag is false then molecfit_data is taken from a single extension as we need to test */
              /* that all ranges assigned to this extension no have wavlength overlap with the data in this extension          */
              if (!parameters->chip_extensions) {
                  double max_wave = cpl_table_get_column_max(molecfit_data,MF_COL_IN_LAMBDA);
                  double min_wave = cpl_table_get_column_min(molecfit_data,MF_COL_IN_LAMBDA);
                  err= molecfit_model_check_extensions_and_ranges(ext,min_wave,max_wave,inc_wranges); 
                  if(!err){
                      inc_wranges = cpl_table_extract_selected(inc_wranges);
                  }
              }


 
              if (!err) {
 
                  /* CALL MOLECFIT */
                  results = mf_model( parameters->mf_config,               /* mf_configuration       *config               */
                                      molecules,                           /* cpl_table              *molecules            */
                                      data->v_ext[ext].spectrum_head,      /* const cpl_propertylist *header_spec          */
                                      molecfit_data,                       /* const cpl_table        *spec                 */
                                      inc_wranges,                         /* cpl_table              *inc_wranges          */
                                      exc_wranges,                         /* cpl_table              *exc_wranges          */
                                      exc_pranges,                         /* cpl_table              *exc_pranges          */
                                      header_kernel,                       /* const cpl_propertylist *header_kernel        */
                                      kernel,                              /* const cpl_matrix       *kernel               */
                                      gdas_user,                           /* const cpl_table        *gdas_user            */
                                      atm_profile_standard,                /* const cpl_table        *atm_profile_standard */
                                      atm_profile_combined);               /* const cpl_table        *atm_profile_combined */


                  /* Cleanup */
                  cpl_table_delete(molecfit_data);
                  err = cpl_error_get_code();
              
              }
              /*if(err){
                  cpl_msg_error(cpl_func,"ERROR AFTER MOLECFIT MODEL err=%d %d error=%s", err, cpl_error_get_code(),cpl_error_get_message());
              }*/
              cpl_propertylist_update_string(data->v_ext[0].header,"ESO DRS MOLECFIT GDAS_SOURCE",(results ? results->gdas_src : "NONE"));
              /*if(err){
                  cpl_msg_error(cpl_func,"ERROR BEFORE GDAS WRITE err=%d %d error=%s", err, cpl_error_get_code(),cpl_error_get_message());
              }*/

              //cpl_msg_error(cpl_func,"NROWS OF GDAS_USER table=%lld; gdas_write %d", cpl_table_get_nrow(gdas_user),gdas_write);
              

              if (!err) {

                  // Write GDAS files on disk 
                  if (results && !gdas_write) {

                      if (!(results->gdas_interpolate)) {
                          err = CPL_ERROR_ILLEGAL_OUTPUT;
                      } else {

                          if (gdas_user) {
                              cpl_msg_info(cpl_func, "Saving Molecfit output user provide fits files ('%s') ... [only first call to Molecfit!]", MOLECFIT_GDAS);
                          } else if (results->gdas_before && results->gdas_after && results->gdas_interpolate) {
                              cpl_msg_info(cpl_func, "Save Molecfit output GDAS automatic fits files ('%s','%s','%s') ... [only first call to Molecfit!]",
			                             MOLECFIT_GDAS_BEFORE, MOLECFIT_GDAS_AFTER, MOLECFIT_GDAS);

                              // Keep time, no access to ESO GDAS DB : Same gdas_interpolate for all the IFUs because all the arms have the same 'date' -> MDJ_OBS in the header 
                              gdas_user = cpl_table_duplicate(results->gdas_interpolate);
                          }

                          err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_GDAS, NULL);
                          if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_GDAS, CPL_TRUE, NULL, results->gdas_interpolate, NULL);

                          // If not default GDAS profiles, save the GDAS files (before/after) from the ESO GDAS DB 

                          if (results->gdas_before) {
                              err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_GDAS_BEFORE, NULL);
                              if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_GDAS_BEFORE, CPL_TRUE, NULL, results->gdas_before, NULL);
                          }

                          if (results->gdas_after) {
                              err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_GDAS_AFTER, NULL);
                              if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_GDAS_AFTER, CPL_TRUE, NULL, results->gdas_after, NULL);
                          }

                          gdas_write = CPL_TRUE;
                      }
                  }
              }
              /*if(err){
                  cpl_msg_error(cpl_func,"ERROR BEFORE ATM PROFILE err=%d %d error=%s", err, cpl_error_get_code(),cpl_error_get_message());
              }*/

              /* Write ATM_PROFILE files on disk */
              if (!err && results && !atm_profile_standard_write) {

                  if (!(results->atm_profile_standard) || !(results->atm_profile_combined)) {
                      err = CPL_ERROR_ILLEGAL_OUTPUT;
                  } else {

                      if (atm_profile_standard) {
                          cpl_msg_info(cpl_func, "Save Molecfit output user provide fits files ('%s') ... [only first call to Molecfit!]", MOLECFIT_ATM_PROFILE_STANDARD);
                      } else {
                          cpl_msg_info(cpl_func, "Saving Molecfit output automatic ATM_PROFILE fits files ('%s','%s') ... [only first call to Molecfit!]",
                                       MOLECFIT_ATM_PROFILE_STANDARD, MOLECFIT_ATM_PROFILE_COMBINED);

                          /* Keep time, no access to disk : Same atm_profile_standard for all the extensions because they have the same 'date' -> MDJ_OBS in the header */
                          atm_profile_standard = cpl_table_duplicate(results->atm_profile_standard);
                      }

                      /* Keep time, no access to disk : Same atm_profile_standard for all the extensions because they have the same 'date' -> MDJ_OBS in the header */
                      atm_profile_combined = cpl_table_duplicate(results->atm_profile_combined);

                      err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_ATM_PROFILE_STANDARD, NULL);
                      if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_ATM_PROFILE_STANDARD, CPL_TRUE, NULL, results->atm_profile_standard, NULL);

                      err = molecfit_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, MOLECFIT_ATM_PROFILE_COMBINED, NULL);
                      if (!err) err = molecfit_save_mf_results(data->v_ext[0].header, MOLECFIT_ATM_PROFILE_COMBINED, CPL_TRUE, NULL, results->atm_profile_combined, NULL);

                      atm_profile_standard_write = CPL_TRUE;
                  }
              }

              /* Save output kernel */
              if (!err && kernel_data) {

                  cpl_matrix *kernel_matrix = NULL;
                  if (results) {
                      if (results->kernel) kernel_matrix = results->kernel;
                  }

                  err += molecfit_save_mf_results(header_kernel, MOLECFIT_MODEL_KERNEL_LIBRARY, CPL_TRUE, kernel_matrix, NULL, NULL);
              }
	  }
          //if(err){
          //    cpl_msg_error(cpl_func,"ERROR %d error=%s", cpl_error_get_code(),cpl_error_get_message());
              //cpl_msg_error(cpl_func,"ERROR error(%d)=%s", err, cpl_error_get_message());
          //}

          /* Save outputs : extensions with/without execute molecfit */
          if (!err && (parameters->use_only_input_pri_ext || ext > 0) ) {
	          cpl_table *atm_profile_fitted = results ? results->atm_profile_fitted : NULL;
              cpl_table *best_fit_params    = results ? results->res                : NULL;
              cpl_table *best_fit_model     = results ? results->spec               : NULL;

              
              err += molecfit_save_mf_results(data->v_ext[ext].header, MOLECFIT_ATM_PARAMETERS,      CPL_TRUE, NULL, atm_profile_fitted, NULL);
              err += molecfit_save_mf_results(data->v_ext[ext].header, MOLECFIT_BEST_FIT_PARAMETERS, CPL_TRUE, NULL, best_fit_params,    NULL);
              /* Also add inc_wranges and exc_wranges to BEST_FIT_PARAMETERS so they be used by molecfit_calctrans */
              if(inc_wranges) {
                  cpl_propertylist* inc_hdr = cpl_propertylist_duplicate(data->v_ext[ext].header);
                  cpl_propertylist_update_string(inc_hdr,"EXTNAME","WAVE_INCLUDE");
                  err += molecfit_save_mf_results(inc_hdr, MOLECFIT_BEST_FIT_PARAMETERS, CPL_TRUE, NULL, inc_wranges,    NULL);
              }
              if(exc_wranges) {
                  cpl_propertylist* exc_hdr = cpl_propertylist_duplicate(data->v_ext[ext].header);
                  cpl_propertylist_update_string(exc_hdr,"EXTNAME","WAVE_EXCLUDE");
                  err += molecfit_save_mf_results(exc_hdr, MOLECFIT_BEST_FIT_PARAMETERS, CPL_TRUE, NULL, exc_wranges,    NULL);
              }
              if(exc_pranges){
                  cpl_propertylist* exc_hdr2 = cpl_propertylist_duplicate(data->v_ext[ext].header);
                  cpl_propertylist_update_string(exc_hdr2,"EXTNAME","PIXEL_EXCLUDE");
                  err += molecfit_save_mf_results(exc_hdr2, MOLECFIT_BEST_FIT_PARAMETERS, CPL_TRUE, NULL, exc_pranges,    NULL);
              }
              err += molecfit_save_mf_results(data->v_ext[ext].header, MOLECFIT_BEST_FIT_MODEL,      CPL_TRUE, NULL, best_fit_model,     NULL);
	      if (results) {
                  err += molecfit_save_mf_results(data->v_ext[ext].header, "MOLECFIT_DATA",              CPL_TRUE, NULL, results->spec,     NULL);
              }
      
 
          }

          /* Cleanup mf_model(...) results */
          if (results) mf_model_results_delete(results, NULL);
      
      }
      /* ----------------------- */
      /* END of Execution Loops */
      /* ---------------------- */

  }/* End if (!err) */

  /* Cleanup */
  if (parameters          ) molecfit_model_parameter_delete( parameters           );
  if (data                ) molecfit_fits_delete(            data                 );
  if (molecules           ) cpl_table_delete(                molecules            );
  if (inc_wranges         ) cpl_table_delete(                inc_wranges          );
  if (exc_wranges         ) cpl_table_delete(                exc_wranges          );
  if (exc_pranges         ) cpl_table_delete(                exc_pranges          );
  if (kernel_data         ) molecfit_fits_delete(            kernel_data          );
  if (mapping_kernel      ) cpl_table_delete(                mapping_kernel       );
  if (gdas_user           ) cpl_table_delete(                gdas_user            );
  if (atm_profile_standard) cpl_table_delete(                atm_profile_standard );
  if (atm_profile_combined) cpl_table_delete(                atm_profile_combined );


  /* Check Recipe status and end */
//if (!err && cpl_errorstate_is_equal(initial_errorstate)) {
//    cpl_msg_info(cpl_func,"Recipe successfully!");
//} else {
//    /* Dump the error history */
//    //only dump if errors
//    if(err){
//        cpl_errorstate_dump(initial_errorstate, CPL_FALSE, NULL);
//        cpl_msg_error(cpl_func,"Recipe failed!, error(%d)=%s", err, cpl_error_get_message());
//    }
//}

  return err;
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
static cpl_error_code molecfit_model_fill_parameterlist(
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

  /* --MOLECFIT_PARAMETER_MODEL_MAPPING_KERNEL */
  const char *mapping_kernel = MOLECFIT_PARAMETER_MODEL_MAPPING_KERNEL_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_MODEL_MAPPING_KERNEL,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
                                                   (const void *)mapping_kernel,
                                                   MOLECFIT_PARAMETER_MODEL_MAPPING_KERNEL_DESC, CPL_FALSE);


  /* --MOLECFIT_PARAMETER_MOLECULES_LIST */
  const char *molecules_list = MOLECFIT_PARAMETER_LIST_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_LIST,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
                                                   (const void *)molecules_list,
                                                   MOLECFIT_PARAMETER_LIST_DESC, CPL_FALSE);

  /* --MOLECFIT_PARAMETER_FIT */
  const char *molecules_fit = MOLECFIT_PARAMETER_FIT_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_FIT,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
                                                   (const void *)molecules_fit,
                                                   MOLECFIT_PARAMETER_FIT_DESC, CPL_FALSE);

  /* --MOLECFIT_PARAMETER_RELATIVE_VALUE */
  const char *molecules_relative_value = MOLECFIT_PARAMETER_RELATIVE_VALUE_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_RELATIVE_VALUE,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
                                                   (const void *)molecules_relative_value,
                                                   MOLECFIT_PARAMETER_RELATIVE_VALUE_DESC, CPL_FALSE);


  /* --MOLECFIT_PARAMETER_WAVE_RANGE_INCLUDE */
  const char *wave_range_include = MOLECFIT_PARAMETER_WAVE_RANGE_INCLUDE_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_WAVE_RANGE_INCLUDE,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
                                                   (const void *)wave_range_include,
                                                   MOLECFIT_PARAMETER_WAVE_RANGE_INCLUDE_DESC, CPL_FALSE);

  /* --MOLECFIT_PARAMETER_WAVE_RANGE_EXCLUDE */
  const char *wave_range_exclude = MOLECFIT_PARAMETER_WAVE_RANGE_EXCLUDE_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_WAVE_RANGE_EXCLUDE,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
                                                   (const void *)wave_range_exclude,
                                                   MOLECFIT_PARAMETER_WAVE_RANGE_EXCLUDE_DESC, CPL_FALSE);

  /* --MOLECFIT_PARAMETER_PIXEL_RANGE_EXCLUDE */
  const char *pixel_range_exclude = MOLECFIT_PARAMETER_PIXEL_RANGE_EXCLUDE_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_PIXEL_RANGE_EXCLUDE,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_STRING, MF_STRLST_NULL,
                                                   (const void *)pixel_range_exclude,
                                                   MOLECFIT_PARAMETER_PIXEL_RANGE_EXCLUDE_DESC, CPL_FALSE);


  /* Add Molecfit configure parameters to the recipe */
  if (!e) e = molecfit_config_fill_parameterlist(RECIPE_NAME, self);


  /* --MOLECFIT_PARAMETER_COMBINE_EXTENSIONS */
  cpl_boolean combine_extensions = MOLECFIT_PARAMETER_CHIP_EXTENSIONS_INIT;
  if (!e) e = molecfit_config_fill_parameter(RECIPE_NAME, self, MOLECFIT_PARAMETER_CHIP_EXTENSIONS,
                                                   !range, dummyMin, dummyMax, CPL_TYPE_BOOL, MF_STRLST_NULL,
                                                   &combine_extensions,
                                                   MOLECFIT_PARAMETER_CHIP_EXTENSIONS_DESC, CPL_FALSE);

  /* Check possible errors */
  if (!cpl_errorstate_is_equal(pre_state) || e != CPL_ERROR_NONE) {
      return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                   "molecfit_model_fill_parameterlist failed!");
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
static molecfit_model_parameter * molecfit_model_parameters(
    const cpl_parameterlist            *list,
    const cpl_propertylist             *raw_head_pri)
{
  /* Check input */
  cpl_error_ensure(list, CPL_ERROR_NULL_INPUT,
                   return NULL, "list input is NULL!");

  /* Get preState */
  cpl_errorstate preState = cpl_errorstate_get();

  /* Create the configuration parameter */
  molecfit_model_parameter *parameters = cpl_malloc(sizeof(molecfit_model_parameter));
  parameters->molecules_table                = NULL;
  parameters->wave_ranges_include_table      = NULL;
  parameters->wave_ranges_exclude_table      = NULL;
  parameters->pixel_ranges_exclude_table     = NULL;
  parameters->mapping_kernel_table           = NULL;
  parameters->mf_config                      = NULL;
  parameters->pl                             = cpl_propertylist_new();


  /*** Load recipe parameters ***/
  const cpl_parameter *p;


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
  p = cpl_parameterlist_find_const(list, MOLECFIT_MODEL_MAPPING_KERNEL);
  parameters->mapping_kernel = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_MODEL_MAPPING_KERNEL, parameters->mapping_kernel);
  /* Create mapping kernel cpl_table */
  if (strcmp(parameters->mapping_kernel, MF_PARAMETERS_NULL)) {
      parameters->mapping_kernel_table = molecfit_config_table_mapping(parameters->mapping_kernel, MOLECFIT_MAPPING_KERNEL_EXT);
  }


  /* Molecules list */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_LIST);
  parameters->list_molec = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_LIST, parameters->list_molec);

  /* Molecules fit */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_FIT);
  parameters->fit_molec = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_FIT, parameters->fit_molec);

  /* Molecules relative value */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_RELATIVE_VALUE);
  parameters->rel_col = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_RELATIVE_VALUE, parameters->rel_col);

  /* Create table */
  if (   strcmp(parameters->list_molec, MF_PARAMETERS_NULL)
      && strcmp(parameters->fit_molec,  MF_PARAMETERS_NULL)
      && strcmp(parameters->rel_col,    MF_PARAMETERS_NULL)) {
      parameters->molecules_table = molecfit_config_table_molecules(parameters->list_molec, parameters->fit_molec, parameters->rel_col);
  }


  /* Wavelength ranges included */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_WAVE_RANGE_INCLUDE);
  parameters->wave_range_include = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_WAVE_RANGE_INCLUDE, parameters->wave_range_include);
  /* Create wavelength ranges included cpl_table */
  if (strcmp(parameters->wave_range_include, MF_PARAMETERS_NULL)) {
      parameters->wave_ranges_include_table = molecfit_config_table_ranges(parameters->wave_range_include, CPL_TYPE_DOUBLE);
  }

  /* Wavelength ranges exclude */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_WAVE_RANGE_EXCLUDE);
  parameters->wave_range_exclude = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_WAVE_RANGE_EXCLUDE, parameters->wave_range_exclude);
  /* Create wavelength ranges exclude cpl_table */
  if (strcmp(parameters->wave_range_exclude, MF_PARAMETERS_NULL)) {
      parameters->wave_ranges_exclude_table = molecfit_config_table_ranges(parameters->wave_range_exclude, CPL_TYPE_DOUBLE);
  }

  /* Pixel ranges exclude */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_PIXEL_RANGE_EXCLUDE);
  parameters->pixel_range_exclude = cpl_parameter_get_string(p);
  cpl_propertylist_update_string(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_PIXEL_RANGE_EXCLUDE, parameters->pixel_range_exclude);
  /* Create pixel ranges exclude cpl_table */
  if (strcmp(parameters->pixel_range_exclude, MF_PARAMETERS_NULL)) {
      parameters->pixel_ranges_exclude_table = molecfit_config_table_ranges(parameters->pixel_range_exclude, CPL_TYPE_INT);
  }


  /*** Get Molecfit config parameters from the input recipe parameter list ***/
  parameters->mf_config = molecfit_config_get_parameters(list, raw_head_pri, parameters->pl);


  /* Combine extensions in the FITS file ? */
  p = cpl_parameterlist_find_const(list, MOLECFIT_PARAMETER_CHIP_EXTENSIONS);
  parameters->chip_extensions = cpl_parameter_get_bool(p);
  cpl_propertylist_update_bool(parameters->pl, MF_PARAMETERS_CONTEX_DEFAULT" "MOLECFIT_PARAMETER_CHIP_EXTENSIONS, parameters->chip_extensions);



  /* Check status */  
  if (!cpl_errorstate_is_equal(preState) || !(parameters->mf_config)) {
      /* Configuration failed */
      molecfit_model_parameter_delete(parameters);
      return NULL;
  } else {
      /* Configuration successfully */
      return parameters;
  }
}

/*----------------------------------------------------------------------------*/
/**
 * @brief    Deallocate the given parameter configuration object and its contents
 *
 * @param    parameters       The parameter configuration variable in the recipe.
 */
/*----------------------------------------------------------------------------*/
static void molecfit_model_parameter_delete(
    molecfit_model_parameter *parameters)
{
  if (parameters) {

      if (parameters->molecules_table)            cpl_table_delete(        parameters->molecules_table            );

      if (parameters->wave_ranges_include_table)  cpl_table_delete(        parameters->wave_ranges_include_table  );

      if (parameters->wave_ranges_exclude_table)  cpl_table_delete(        parameters->wave_ranges_exclude_table  );

      if (parameters->pixel_ranges_exclude_table) cpl_table_delete(        parameters->pixel_ranges_exclude_table );

      if (parameters->mapping_kernel_table)       cpl_table_delete(        parameters->mapping_kernel_table       );

      if (parameters->mf_config)                  mf_configuration_delete(              parameters->mf_config                  );

      if (parameters->pl)                         cpl_propertylist_delete( parameters->pl                         );

      cpl_free(parameters);
  }
}



/*----------------------------------------------------------------------------*/
/**
 * @brief    Checks the validity of the range wavlength coverage and that of the associated extension data
 *
 * @param    rangtab       The range table.
 * @param    data          The science data molecfit_fits object.
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code molecfit_model_check_extensions_and_ranges(cpl_size extension, double min_wav, double max_wav, cpl_table* range) {

  cpl_msg_info (cpl_func,"CHECK FOR extension %lld, min %f max %f",  extension, min_wav,max_wav);
  //cpl_table_dump(range,0,cpl_table_get_nrow(range),NULL);
  
  /* Check all ranges that are mapped to this extension */
  
  cpl_size nranges=cpl_table_get_nrow(range);
  
  cpl_boolean check_ok=CPL_TRUE;
  for (cpl_size i=0;i<nranges;i++) {
      
      //cpl_size range_idx=i+1;
     
      int    range_assigned_extension =  cpl_table_get_int   (range,MF_COL_WAVE_RANGE_MAP2CHIP,i,NULL);
      double range_low_limit          =  cpl_table_get_double(range,MF_COL_WAVE_RANGE_LOWER,        i, NULL);
      double range_upp_limit          =  cpl_table_get_double(range,MF_COL_WAVE_RANGE_UPPER,        i,NULL);
      
      /* Only check ranges that have been mapped to this extension data */
      if (range_assigned_extension!=extension) continue;
      
      /* Check that there is an overlap between this range and the wavlength coverage on this extension */
      //if (range_low_limit>=max_wav) check_ok=CPL_FALSE; 
      //if (range_upp_limit<=min_wav) check_ok=CPL_FALSE; 
      if(range_upp_limit<=min_wav || range_low_limit>=max_wav){
        cpl_table_unselect_row(range,i);
      }
  }
  if(cpl_table_count_selected(range) == 0){
      check_ok = CPL_FALSE;
  }
     
     
  /* If the check is not ok then set error message and return illegal input flag */
  if (!check_ok) {
         cpl_msg_error(cpl_func,"No valid ranges found");
         //cpl_msg_error(cpl_func,"Range %lld [%f , %f] Assigned to Extension %d [%f , %f]. Check Overlap NOT OK",
         //                    range_idx, range_low_limit, range_upp_limit, range_assigned_extension,min_wav,max_wav);
         return CPL_ERROR_ILLEGAL_INPUT;
  }
  
      /* Output that this range has been checked and show details */ 
     // cpl_msg_info(cpl_func,"Range %lld [%f , %f] Assigned to Extension %d [%f , %f]. Check Overlap OK", 
     //                          range_idx, range_low_limit, range_upp_limit, range_assigned_extension,min_wav,max_wav);


  return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief    Checks and merges all wavelength range specific data into a single
 *           table.
 *
 * @param  inc_wranges_new_table  The new wavelength range table to be built.
 * @param  inc_wranges_ext_table  The wavelength range table loaded from external fits file.
 * @param  wave_include_str       Parameter string from .rc file that defines wave_include
 * @param  map_regions_str        Parameter string from .rc file that defines map_regions
 * @param  fit_continuum_str      Parameter string from .rc file that defines fit_continuum
 * @param  continuum_n_str        Parameter string from .rc file that defines continuum_n
 * @param  fit_wlc_str            Parameter string from .rc file that defines fit_wlc
 * @param  nchips                 No of chips declared.
*/
/*----------------------------------------------------------------------------*/
static cpl_table* molecfit_model_amalgamate_ranges_table(
  cpl_table* inc_wranges_ext_table,
  const char* wave_include_str, 
  const char* map_regions_str,  
  const char* inp_fit_continuum_str,
  const char* continuum_n_str,  
  const char* inp_fit_wlc_str,
  int         nchips,
  double      wlg2mn)      


{

  /* Temporary Fix from old regression tests that sets continuum_str = TRUE or FALSE */
  const char * fit_continuum_str;
  if        (strcmp(inp_fit_continuum_str,"TRUE")==0) {
      fit_continuum_str="1";
  } else if (strcmp(inp_fit_continuum_str,"FALSE")==0) {
      fit_continuum_str="0";
  } else {
      fit_continuum_str=inp_fit_continuum_str;
  }
  
  const char * fit_wlc_str;
  if        (strcmp(inp_fit_wlc_str,"TRUE")==0) {
      fit_wlc_str="1";
  } else if (strcmp(inp_fit_wlc_str,"FALSE")==0) {
      fit_wlc_str="0";
  } else {
      fit_wlc_str=inp_fit_wlc_str;
  }


  cpl_table* inc_wranges_new_table = NULL; /* The return value if succesful */

  double rangev    [2*MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];  
  int    mapv      [  MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];
  int    fit_contv [  MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];
  int    cont_nv   [  MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];
  int    fit_wlcv  [  MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];
  cpl_boolean R2C_MAP;
  
 
  cpl_msg_info (cpl_func,"Amalgamating all region related parameters into single table");          
  

  /* ------------------------------------------------------------------------------
  If there is no external wave_include fits file then check that none of the string
     parameters is declared as "NULL"
     ------------------------------------------------------------------------------*/
  if (inc_wranges_ext_table==NULL) {    
  
      cpl_msg_info(cpl_func,"There is no external table");
                                  
      if (strcmp(wave_include_str,MF_PARAMETERS_NULL)==0) {
          cpl_msg_error(cpl_func,"No range values defined");
          /*return NULL;*/
      }

      if (strcmp(map_regions_str,MF_PARAMETERS_NULL)==0) {
          cpl_msg_info(cpl_func,"No map regions to chips defined. Assuming all are mapped to chip 1");
          R2C_MAP=CPL_FALSE;
      } else {
          R2C_MAP=CPL_TRUE;
      }

      if (strcmp(fit_continuum_str,MF_PARAMETERS_NULL)==0) {
          cpl_msg_error(cpl_func,"No range specific continuum fit flags defined");
          return NULL;
      }

      if (strcmp(continuum_n_str,MF_PARAMETERS_NULL)==0) {
          cpl_msg_error(cpl_func,"No range specific poly order of continuum modelling defined");
          return NULL;
      }

      if (strcmp(fit_wlc_str,MF_PARAMETERS_NULL)) {
          cpl_msg_error(cpl_func,"No range specific wavlength correction flags set");
          return NULL;
      }

  } else {
      cpl_msg_info(cpl_func,"There is an external table");

  }
                                
  /* At this stage we know that there is either an external table or range values defined in a string */

                               
  /* ------------------------------- */
  /* WAVE_INCLUDE RANGE LIMIT VALUES */
  /* ------------------------------- */

  int    nranges;
  double llimv[MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];
  double ulimv[MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];

  if (strcmp(wave_include_str,MF_PARAMETERS_NULL)==0) {
      
      /* No values defined in a string so get them from the external table */
      cpl_msg_info(cpl_func,"Wave Inlude string=NULL. Reading data from external table.");
      
      /* Check that we have the lower limit data column */
      if (!cpl_table_has_column(inc_wranges_ext_table,MF_COL_WAVE_RANGE_LOWER)) {
          cpl_msg_error(cpl_func,"External table does not have column %s",MF_COL_WAVE_RANGE_LOWER);
          return NULL;
      }

      /* Check that we have the upper limit data column */
      if (!cpl_table_has_column(inc_wranges_ext_table,MF_COL_WAVE_RANGE_UPPER)) {
          cpl_msg_error(cpl_func,"External table does not have column %s", MF_COL_WAVE_RANGE_UPPER);
          return NULL;
      }

      nranges=cpl_table_get_nrow(inc_wranges_ext_table); 
      for (int i=0;i<nranges;i++) {
          llimv[i] = cpl_table_get(inc_wranges_ext_table,MF_COL_WAVE_RANGE_LOWER,i,NULL);
          ulimv[i] = cpl_table_get(inc_wranges_ext_table,MF_COL_WAVE_RANGE_UPPER,i,NULL);
      }
 
 
  } else {
  
      /* Values are defined in a string as an even list of doubles so parse the string into a vector of doubles */
      cpl_msg_info(cpl_func,"Wave Include string=%s",wave_include_str);
      
      int full_string_size=strlen(wave_include_str);
      int ndvals=parsMFparsReadDValsFromString(wave_include_str,rangev, 2*MF_PARAMETERS_MAXIMUM_NO_OF_RANGES,full_string_size);

      /* Check that parsing string was successful */
      if (ndvals<1) {
          cpl_msg_error(cpl_func,"Invalid list of wavelength ranges : %s",wave_include_str);
          return NULL;
      }
      
      /* Check that we have an even no of doubles */
      if (ndvals % 2 !=0 ) {
          cpl_msg_error(cpl_func,"Odd number of range boundaries defined: %s",wave_include_str);
          return NULL;
      } else {
          cpl_msg_info(cpl_func,"Even number of range values defined");
      }
      
      nranges=ndvals/2;
      for (int i=0; i<nranges;i++) {
          llimv[i]=rangev[2*i  ];
          ulimv[i]=rangev[2*i+1];
      }
  }
  
  /* Now check that lower values < upper values */
  cpl_msg_info(cpl_func,"%d Ranges defined:",nranges);
  for (int i=0; i<nranges;i++) {
      cpl_msg_info(cpl_func,"Range %d: [%f, %f]um",i+1, llimv[i]*wlg2mn, ulimv[i]*wlg2mn);
      if (llimv[i]>=ulimv[i]) {
          cpl_msg_error(cpl_func,"Invalid range boundaries [%f, %f]um",llimv[i]*wlg2mn, ulimv[i]*wlg2mn);
          return NULL;
      }
  }
 

  /* ------------------------------- */
  /* RANGES MAPPED TO CHIP  DETAILS  */
  /* ------------------------------- */

  /* Assume there is no valid range to chip map column unless proven otherwise */
  R2C_MAP=CPL_FALSE;
      cpl_msg_info(cpl_func,"Sorting Range Mappings To %d Chips.",nchips);
  
  if (strcmp(map_regions_str,MF_PARAMETERS_NULL)==0) {
      
      /* No values defined in a string so get them from the external table */
      cpl_msg_info(cpl_func,"Map Regions string=NULL. Reading data from external table");
      
      /* Check that there is the required data column */
      if (!cpl_table_has_column(inc_wranges_ext_table,MF_COL_WAVE_RANGE_MAP2CHIP)) {
      
          /* No valid column so report and set R2C_MAP to false */
          cpl_msg_info(cpl_func,"External table does not have column %s. Assuming that all regions are mapped to chip 1",
          MF_COL_WAVE_RANGE_MAP2CHIP);
      
      } else {
      
          /* Load the column values into a vector */
          for (int i=0;i<nranges;i++) {
              mapv[i] = cpl_table_get_int(inc_wranges_ext_table,MF_COL_WAVE_RANGE_MAP2CHIP,i,NULL);
          }
          R2C_MAP=CPL_TRUE;
      }
      
  } else {
  
      /* Values are defined in a string as a list of integers so parse the string into a vector of integers */
      cpl_msg_info(cpl_func,"Map Regions string=%s",map_regions_str);
      
      int full_string_size=strlen(map_regions_str);
      int nivals=parsMFparsReadIValsFromString(map_regions_str,mapv, MF_PARAMETERS_MAXIMUM_NO_OF_RANGES,full_string_size);

      /* Check that we have the correct number of integers */
      if (nivals!=nranges && nivals!=1) {
          cpl_msg_error(cpl_func,"Range mapping parameter has incorrect list size");
          return NULL;
      } else if (nivals==1) {
          /* The exception to the rule is if there is only one value then that implies that this is the smae for all */
          cpl_msg_info(cpl_func,"Assuming that all regions are mapped to Chip %d",mapv[0]);
          for (int i=1;i<nranges;i++) mapv[i] = mapv[0];
      }
      /* The R2C MAP has now been derived from a string definition into the vector form */
      R2C_MAP=CPL_TRUE;
      
  }
  
  /* In the event that there are no maps defined, (i.e. the map string is set to "NULL" and either the external   */
  /* table does not exists or the required column does not exist in the external table) then the requirements are */
  /* that all ranges should be mapped to chip  1                                                                  */
  if (R2C_MAP==CPL_FALSE) {
      for (int i=0;i<nranges;i++) mapv[i] = 1;
  }



  /* ------------------------------- */
  /* RANGES FIT CONTINUUM FLAGS      */
  /* ------------------------------- */
  if (strcmp(fit_continuum_str,MF_PARAMETERS_NULL)==0) {
      
      /* No values defined in a string so get them from the external table */
      
      /* Check that there is the required data column */
      if (!cpl_table_has_column(inc_wranges_ext_table,MF_COL_WAVE_RANGE_CONT_FIT)) {
          cpl_msg_error(cpl_func,"External table does not have column %s",MF_COL_WAVE_RANGE_CONT_FIT);
          return NULL;
      }
      
      /* Load the values into a vector */
      for (int i=0;i<nranges;i++) {
          fit_contv[i] = cpl_table_get_int(inc_wranges_ext_table,MF_COL_WAVE_RANGE_CONT_FIT,i,NULL);
      }
      
  } else {
  
      /* Values are defined in a string as a list of integers so parse the string into a vector of integers */
      cpl_msg_info(cpl_func,"Fit Continuum string=%s",fit_continuum_str);
      
      int full_string_size=strlen(fit_continuum_str);
      int nivals=parsMFparsReadIValsFromString(fit_continuum_str,fit_contv, MF_PARAMETERS_MAXIMUM_NO_OF_RANGES,full_string_size);

      /* Check that we have the correct number of integers */
      if (nivals!=nranges && nivals!=1) {
          cpl_msg_error(cpl_func,"Range fit continuum parameter has incorrect list size");
          return NULL;
      } else if (nivals==1) {
          /* Special case for nivals=1 implies that this value should be used for all */
          cpl_msg_info(cpl_func,"Assuming that all ranges have fit continuum flag %d",fit_contv[0]);
          for (int i=1;i<nranges;i++) fit_contv[i] = fit_contv[0];
      }         
      
  }

  /* ------------------------------------------- */
  /* RANGES POLYNOMIAL ORDER FOR CONTINUUM MODEL */
  /* ------------------------------------------- */

  if (strcmp(continuum_n_str,MF_PARAMETERS_NULL)==0) {
      
      /* No values defined in a string so get them from the external table */
      
      /* Check that there is the required data column */
      if (!cpl_table_has_column(inc_wranges_ext_table,MF_COL_WAVE_RANGE_CONT_ORDER)) {
          cpl_msg_error(cpl_func,"External table does not have column %s",MF_COL_WAVE_RANGE_CONT_ORDER);
          return NULL;
      }
      
      /* Load the values into a vector */
      for (int i=0;i<nranges;i++) {
          cont_nv[i] = cpl_table_get_int(inc_wranges_ext_table,MF_COL_WAVE_RANGE_CONT_ORDER,i,NULL);
      }
      
  } else {
  
      /* Values are defined in a string as a list of integers so parse the string into a vector of integers */
      cpl_msg_info(cpl_func,"Continuum n string=%s",continuum_n_str);
      
      int full_string_size=strlen(continuum_n_str);
      int nivals=parsMFparsReadIValsFromString(continuum_n_str,cont_nv, MF_PARAMETERS_MAXIMUM_NO_OF_RANGES,full_string_size);
 
     /* Check that we have the correct number of integers */
      if (nivals!=nranges && nivals!=1) {
          cpl_msg_error(cpl_func,"Range continuum polynomial order parameter has incorrect list size");
          return NULL;
      } else if (nivals==1) {
          /* Special case for nivals=1 implies that this value should be used for all ranges */
          cpl_msg_info(cpl_func,"Assuming that all ranges have polynomial order %d to model continuum",cont_nv[0]);
          for (int i=1;i<nranges;i++) cont_nv[i] = cont_nv[0];
      }         
      
  }
  
   

  /* ------------------------------------------- */
  /* RANGES WAVLENGTH CORRECTION FIT FLAGS       */
  /* ------------------------------------------- */

  if (strcmp(fit_wlc_str,MF_PARAMETERS_NULL)==0) {
      cpl_msg_info(cpl_func,"Fit WLC string is NULL!");
      
      /* No values defined in a string so get them from the external table */
      
      /* Check that there is the required data column */
      if (!cpl_table_has_column(inc_wranges_ext_table,MF_COL_WAVE_RANGE_WLC_FIT)) {
          cpl_msg_error(cpl_func,"External table does not have column %s",MF_COL_WAVE_RANGE_WLC_FIT);
          return NULL;
      }
      
      /* Load the values into a vector */
      cpl_msg_info(cpl_func,"Loading WLC Fit flags from fits table");
      for (int i=0;i<nranges;i++) {
          fit_wlcv[i] = cpl_table_get_int(inc_wranges_ext_table,MF_COL_WAVE_RANGE_WLC_FIT,i,NULL);
       }
      
  } else {
  
      /* Values are defined in a string as a list of integers so parse the string into a vector of integers */
      cpl_msg_info(cpl_func,"Fit WLC string=%s",fit_wlc_str);
      
      int full_string_size=strlen(fit_wlc_str);
      int nivals=parsMFparsReadIValsFromString(fit_wlc_str,fit_wlcv, MF_PARAMETERS_MAXIMUM_NO_OF_RANGES,full_string_size);

      /* Check that we have the correct number of integers */
      if (nivals!=nranges && nivals!=1) {
          cpl_msg_error(cpl_func,"Range wavlength correction flags parameter has incorrect list size");
          return NULL;
      } else if (nivals==1){
          /* Special case for nivals=1 implies that this value should be used for all ranges */
          for (int i=1;i<nranges;i++) fit_wlcv[i] = fit_wlcv[0];
      }          
      
  }
  for (int i=0;i<nranges;i++) {
      cpl_msg_info(cpl_func,"WLC Fit flag for range %d =%d",i+1,fit_wlcv[i]);
  }


  /* ------------------------------------------- */
  /* ALL DATA OBTAINED SO NOW CREATE NEW TABLE   */
  /* ------------------------------------------- */
  
  cpl_error_code err=CPL_ERROR_NONE;
  
  cpl_msg_info(cpl_func,"About to create a new table with %d rows",nranges);
  
  /* Create a new table with same number of rows as there are ranges*/
  inc_wranges_new_table=cpl_table_new(nranges);
  
  if (inc_wranges_new_table==NULL) {
      err=CPL_ERROR_NULL_INPUT;
      cpl_msg_error(cpl_func,"Failed to create a new table");
  } else {  
      int tmpi = cpl_table_get_nrow(inc_wranges_new_table);
      cpl_msg_info(cpl_func,"Have created a new table with %d rows",tmpi);
  }
  
  /* Define the table columns  */  
  if (!err) err=cpl_table_new_column(inc_wranges_new_table,MF_COL_WAVE_RANGE_LOWER,     CPL_TYPE_DOUBLE);
  if (!err) err=cpl_table_new_column(inc_wranges_new_table,MF_COL_WAVE_RANGE_UPPER,     CPL_TYPE_DOUBLE);
  if (!err) err=cpl_table_new_column(inc_wranges_new_table,MF_COL_WAVE_RANGE_MAP2CHIP,  CPL_TYPE_INT);
  if (!err) err=cpl_table_new_column(inc_wranges_new_table,MF_COL_WAVE_RANGE_CONT_FIT,  CPL_TYPE_INT);
  if (!err) err=cpl_table_new_column(inc_wranges_new_table,MF_COL_WAVE_RANGE_CONT_ORDER,CPL_TYPE_INT);
  if (!err) err=cpl_table_new_column(inc_wranges_new_table,MF_COL_WAVE_RANGE_WLC_FIT,   CPL_TYPE_INT);
  
  /* Populate the table */
  for (int i=0;i<nranges;i++) {
      if (!err) err=cpl_table_set_double(inc_wranges_new_table, MF_COL_WAVE_RANGE_LOWER,     i,llimv[i]);
      if (!err) err=cpl_table_set_double(inc_wranges_new_table, MF_COL_WAVE_RANGE_UPPER,     i,ulimv[i]);
      if (!err) err=cpl_table_set_int    (inc_wranges_new_table,MF_COL_WAVE_RANGE_CONT_FIT,  i,fit_contv[i]);
      if (!err) err=cpl_table_set_int    (inc_wranges_new_table,MF_COL_WAVE_RANGE_CONT_ORDER,i,cont_nv[i]);
      if (!err) err=cpl_table_set_int    (inc_wranges_new_table,MF_COL_WAVE_RANGE_MAP2CHIP,  i,mapv[i]);
      if (!err) err=cpl_table_set_int    (inc_wranges_new_table,MF_COL_WAVE_RANGE_WLC_FIT,   i,fit_wlcv[i]);  
  }

  if (err) {
      cpl_msg_error(cpl_func,"Unable to create new wavelength range table");
      if (inc_wranges_new_table!=NULL) cpl_free(inc_wranges_new_table);
      return NULL;
  } 

  return inc_wranges_new_table;

}



/*----------------------------------------------------------------------------*/
/**
 * @brief    Parses a string defining a ',' deliminated list of integer values
 *           and populates values onto supplied vector
 *           table. Returns no of integer values parsed or -1 if failure.
 *
 * @param  str           The string with the list of values
 * @param  vec           The vector to populate
 * @param  vec           The size of the vector
 * @param  max_strlen    Maximum size of the string and deliminated substrings
*/
/*----------------------------------------------------------------------------*/

int parsMFparsReadIValsFromString (const char* str, int* vec, int vsize,int max_strlen) {

    char bufferA[vsize][max_strlen+1];
    
    /*Parse string lists into individual string elements and place in bufferA*/
    int str_size=strlen(str);        
    int ndelim=0;
    int i=0;
    for (int idx=0;idx<str_size;idx++) {
        if (str[idx]==',') {
           bufferA[ndelim][i]='\0';
           ndelim++;
           i=0;
        } else {
            bufferA[ndelim][i]=str[idx];
            i++;
	    bufferA[ndelim][i]='\0';
        }
    }

    /* Iterate each string stored in the buffer array 
       convert to an integer and store in the return array */
    int nelems=ndelim+1;
    int ival;
    
    for (i=0;i<nelems;i++) {
        int rn = sscanf(bufferA[i], "%d", &ival);
        if (rn==1) {
           vec[i]=ival;
        } else {
           return -1;
        }
    }
    
    return nelems;

}

/*----------------------------------------------------------------------------*/
/**
 * @brief    Parses a string defining a ',' deliminated list of doublevalues
 *           and populates values onto supplied vector
 *           table. Returns no of integer values parsed or -1 if failure.
 *
 * @param  str           The string with the list of values
 * @param  vec           The vector to populate
 * @param  vec           The size of the vector
 * @param  max_strlen    Maximum size of the string and deliminated substrings
*/
/*----------------------------------------------------------------------------*/


int parsMFparsReadDValsFromString (const char* str, double* vec, int vsize,int max_strlen) {

    char bufferA[vsize][max_strlen];
    
    /*Parse string lists into individual string elements and place in bufferA*/
    int str_size=strlen(str);        
    int ndelim=0;
    int i=0;
    for (int idx=0;idx<str_size;idx++) {
        if (str[idx]==',') {
           bufferA[ndelim][i]='\0';
           ndelim++;
           i=0;
        } else {
            bufferA[ndelim][i]=str[idx];
            i++;
	    bufferA[ndelim][i]='\0';
        }
    }

    /* Iterate each string stored in the buffer array 
       convert to an integer and store in the return array */
    int nelems=ndelim+1;
    double dval;
    
    for (i=0;i<nelems;i++) {
        
       
        int rn = sscanf(bufferA[i], "%lf", &dval);
        
        if (rn==1) {
           vec[i]=dval;
        } else {
           return -1;
        }
    }
    
    return nelems;

}

/*----------------------------------------------------------------------------*/
/**
 * @brief    FOR EXPERT MODE ONLY
 *           Will look for fits file that defines initial values and will 
 *           overide default or predeclared values
 *
 * @param  inp_frame              The fits frame conatning initial values to use.
 * @param  molecules_table        The current moleculs table with current values.
 */
/*----------------------------------------------------------------------------*/

static cpl_error_code molecfit_model_expert_mode(const cpl_frame* inp_frame, cpl_table* molecules_table, mf_configuration  *config) {
    
    cpl_msg_info(cpl_func,"=======================================");
    cpl_msg_info(cpl_func,"Parsing Expert Mode Initial Fit Values");
    cpl_msg_info(cpl_func,"(Warning still a work in progress!)");
    cpl_msg_info(cpl_func,"=======================================\n");
    
    
    cpl_msg_info(cpl_func,"===================MOLECULE TABLE DUMP========================");
    cpl_size n_molecs=cpl_table_get_nrow(molecules_table);
    cpl_table_dump(molecules_table,0,n_molecs,NULL);
 
    double  chipCoeffA[MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS][MF_FIT_N_POLYNOME_MAX];
    double rangeCoeffA[MF_PARAMETERS_MAXIMUM_NO_OF_RANGES][MF_FIT_N_POLYNOME_MAX];

    cpl_boolean chipCoeffADefined[MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS][MF_FIT_N_POLYNOME_MAX];
    for (int chip_idx=0;chip_idx<MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS;chip_idx++) {
        for (int coeff_idx=0;coeff_idx<MF_FIT_N_POLYNOME_MAX;coeff_idx++) {
            chipCoeffADefined[chip_idx][coeff_idx]=CPL_FALSE;
        }
    }
    
    const char *filename = cpl_frame_get_filename(inp_frame);
    
    cpl_msg_info(cpl_func,"Found file %s",filename);
    cpl_table* experttab=NULL;
    

    /* Iterate through all extensions in this fits file until a table
       with columns parameter and value is found                    */
    molecfit_fits* fits_obj = molecfit_fits_load(filename, CPL_FALSE);
    int n_exts= fits_obj->n_ext;
    cpl_table *tab;
    for (int i=0; i<n_exts;i++) {
    
        tab=fits_obj->v_ext[i].table;
        if (tab==NULL) continue;

        /* Check that this table has parameter and value columns */
        cpl_array* col_names= cpl_table_get_column_names (tab) ;
        if (col_names==NULL) continue;

        cpl_size n_cols=cpl_array_get_size(col_names);
        cpl_boolean found_parameter=CPL_FALSE;
        cpl_boolean found_value   =CPL_FALSE;
        for (cpl_size j=0; j<n_cols;j++) {
            const char * col_name=cpl_array_get_string(col_names,j);
            if (!strcmp(col_name,"parameter")) found_parameter=CPL_TRUE;
            if (!strcmp(col_name,"value"    )) found_value    =CPL_TRUE;
        }
        if (!found_parameter || !found_value) continue;
                        
        cpl_msg_info(cpl_func,"Extension %d has table with parameter and value columns", i);        
        experttab=tab;
    }

    /*  If we have expert mode table then parse for specific parameter and values ;*/
    if (!experttab) {
       cpl_msg_info(cpl_func,"NO TABLE LOADED!");
    } else {
    
       cpl_msg_info(cpl_func,"PARAMETER TABLE LOADED!");
    
       int n=cpl_table_get_nrow(experttab);
       cpl_msg_info(cpl_func,"TABLE SIZDE=%d",n);
       //cpl_table_dump(experttab,0,n,NULL);
       
       for (cpl_size i=0; i<n;i++) {
           const char* str=cpl_table_get_string(experttab,"parameter",i);
           double dval =cpl_table_get_double(experttab,"value",i,NULL);
           if (str==NULL) continue;


           /* ----------------------------- */
           /* Check for the FWHM parameters */
           /* ----------------------------- */
           if (!strcmp(str,"boxfwhm"    )) {
               config->parameters->fitting.fit_res_box.const_val=dval;
               cpl_msg_info(cpl_func,"BOXWIDTH FWHM  --->  %f",dval);
           }
           if (!strcmp(str,"gaussfwhm"  )) {
               config->parameters->fitting.fit_gauss.const_val=dval;
               cpl_msg_info(cpl_func,"GAUSSWIDTH   FWHM =  %f",dval);
           }
           if (!strcmp(str,"lorentzfwhm")) {
               config->parameters->fitting.fit_lorentz.const_val=dval;
               cpl_msg_info(cpl_func,"LORENTZWIDTH FWHM =  %f",dval);
           }
           

           /* ----------------------------------------------------------- */
           /* Check for the relative column values for specific molecules */
           /* ----------------------------------------------------------- */
           char substr[12];
           strncpy(substr,str,11);
           if (!strcmp(substr,"rel_mol_col")) {
               int str_size=strlen(str);
               char mol_str[str_size];
               for (int j=12;j<str_size;j++) {
                   mol_str[j-12]=str[j];
                   mol_str[j-11]='\0';
               }    
               cpl_error_code mol_err=CPL_ERROR_NONE;
               mol_err= mf_molecules_str_check(mol_str) ;
               if (mol_err==CPL_ERROR_NONE) {                    
                   cpl_msg_info(cpl_func,"REL_MOL_COL for %s =  %f",mol_str,dval);
                   cpl_boolean found_flag=CPL_FALSE;
                   for (cpl_size molec_idx=0;molec_idx<n_molecs;molec_idx++) {
                       const char* mol_str2=cpl_table_get_string(molecules_table,"LIST_MOLEC",molec_idx);
                       double dval2  =cpl_table_get_double(molecules_table,"REL_COL",molec_idx,NULL);
                       if (!strcmp(mol_str2,mol_str)) {
                           found_flag=CPL_TRUE;
                           cpl_msg_info(cpl_func,"Overriding REL_COL value for molecule %s from value %f to value %f",mol_str2,dval2,dval);
                           cpl_table_set_double(molecules_table,"REL_COL",molec_idx,dval);
                       }
                       if (!found_flag) cpl_msg_info(cpl_func,"Undeclared molecule %s. Ignoring", mol_str);
                   }
                                  
               }
           } /* endif rel_col_mol*/
           
           
           /* ----------------------------------------------------------- */
           /* Check for coeffs in the chip and range specific polynomials */
           /* ----------------------------------------------------------- */
           int ridx,cidx,coefidx;
           if (parseStr4Coeffs(str,&ridx,&cidx,&coefidx)>=0) {
           
               if (cidx   >=0 && cidx   < MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS &&
                   coefidx>=0 && coefidx< MF_FIT_N_POLYNOME_MAX) {
 
                   if(ridx==-1) {
                       /* A valid chip coeff indexed value */
                       cpl_msg_info(cpl_func,"CHIP %d,Coeff %d = %f",cidx,coefidx,dval);
                       chipCoeffA[cidx][coefidx]=dval;
                       chipCoeffADefined[cidx][coefidx]=CPL_TRUE;
                   
                   } else if (ridx>=0 && ridx< MF_PARAMETERS_MAXIMUM_NO_OF_RANGES) {
                   
                       /* A valid Range coeff indexed value */
                       cpl_msg_info(cpl_func,"RANGE %d,CHIP %d,Coeff %d = %f",ridx,cidx,coefidx,dval);
                       rangeCoeffA[ridx][coefidx]=dval;
                   
                   }/* endif ridx */
                    
               }

           }/*endif parseStr4Coeffs*/
           
       }/* end for i */
    
    }/* endifexperttab*/
    
    
    /* -----------------------------------------------*/
    /* MUST PUT HERE WHAT TO DO WITH THE FOUND VALUES */
    /*          (WORK IN PROGRESS)                    */  
    /* -----------------------------------------------*/
    
    /*---------------------------*/
    /* WAVLENGTH CORRECTION POLYS*/
    /*---------------------------*/
 
    /* Poly order is fixed for all chips*/
    int np = config->parameters->fitting.fit_wavelenght.n;
    cpl_msg_info(cpl_func,"wc poly order=%d",np);

    /* Initialising the chip coeff array to unit scale and zero else where*/
    for (int chip_idx=0;chip_idx<MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS;chip_idx++) {
        for (int coeff_idx=0;coeff_idx<MF_FIT_N_POLYNOME_MAX;coeff_idx++) {
            config->parameters->fitting.wlc_coeffs[chip_idx][coeff_idx]=0.0;
        }
            config->parameters->fitting.wlc_coeffs[chip_idx][1]=1.0;
    }

    /* Where defined in the valid range, set the coeffs values to those extracted*/
    for (int chip_idx=1;chip_idx<MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS;chip_idx++) {
        for (int coeff_idx=0;coeff_idx<=np;coeff_idx++) {
            double dval=0.0;
            if (chipCoeffADefined[chip_idx][coeff_idx]) dval=chipCoeffA[chip_idx][coeff_idx];
            config->parameters->fitting.wlc_coeffs[chip_idx-1][coeff_idx]=dval;
            
            if (chip_idx<=4) cpl_msg_info(cpl_func,"Initial Value for CHIP %d Coeff %d of %d set to %f",chip_idx,coeff_idx,np+1,dval);
        }
    }

 
    /*-----------------------*/
    /* CONTINUUM MODEL  POLYS*/
    /*-----------------------*/
    
 
    /* NOTE: Poly order is specific to each range */

    /* Initialising the range coeff array to unit scale and zero everywhere else */
    for (int range_idx=0;range_idx<MF_PARAMETERS_MAXIMUM_NO_OF_RANGES;range_idx++) {
        for (int coeff_idx=0;coeff_idx<MF_FIT_N_POLYNOME_MAX;coeff_idx++) {
            config->parameters->fitting.cont_coeffs[range_idx][coeff_idx]=0.0;
        }
        config->parameters->fitting.cont_coeffs[range_idx][1]=1.0;
    }

    /* Where defined in the valid range, set the coeffs values to those extracted*/
    for (int range_idx=1;range_idx<MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS;range_idx++) {
        int i = range_idx-1;
        int cont_porder = config->parameters->fitting.cont_poly_order[i];
        for (int coeff_idx=0;coeff_idx<=cont_porder;coeff_idx++) {
            double dval=rangeCoeffA[range_idx][coeff_idx];
            config->parameters->fitting.cont_coeffs[range_idx-1][coeff_idx]=dval;
            
            if (range_idx<=4) cpl_msg_info(cpl_func,
                  "Initial Value for RANGE %d Coeff %d (of %d coeffs) set to %f",
                  range_idx,coeff_idx,cont_porder+1,dval);
        }
    }



/*
  int                        cont_poly_order[MF_PARAMETERS_MAXIMUM_NO_OF_RANGES];
  double                     cont_coeffs    [MF_PARAMETERS_MAXIMUM_NO_OF_RANGES][MF_FIT_N_POLYNOME_MAX];
  int                        wlc_poly_order [MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS];
  double                     wlc_coeffs     [MF_PARAMETERS_MAXIMUM_NO_OF_CHIPS ][MF_FIT_N_POLYNOME_MAX];
*/

    cpl_msg_info(cpl_func,"-------------------------------------------------------------");
    
    cpl_free(fits_obj);
 
    return CPL_ERROR_NONE; 

}

/*----------------------------------------------------------------------------*/
/**
 * @brief  Parse a string for identification of RANGE or CHIP specific 
 *         polynomial coefficients.
  *
 * @param  str             The fits frame conatning initial values to use.
 * @param  r_idx           The current moleculs table with current values.
 * @param  c_idx           The current moleculs table with current values.
 * @param  coef_idx        The current moleculs table with current values.
 */
/*----------------------------------------------------------------------------*/
static int parseStr4Coeffs(const char *str, int *r_idx, int* c_idx, int* coef_idx) {
    
    int len=strlen(str);
    char onechar;
    int start_idx;
    char buffer[len];
    
    enum TagName {NONE, CHIP, RANGE, COEFF};
    enum TagName current=NONE;
    

    /* Default values */   
    *r_idx     = -1;
    *c_idx     = -1;
    *coef_idx  = -1;
    
    for (int i=0;i<=len;i++) {
        onechar=str[i];
        if (onechar==' ') continue;
        
        /* Look ahead for "RANGE " */
        if ((onechar   =='R' || onechar =='r') && i+5<len) {
           if ((str[i+1]=='A' || str[i+1]=='a') && 
               (str[i+2]=='N' || str[i+2]=='n') &&            
               (str[i+3]=='G' || str[i+3]=='g') &&            
               (str[i+4]=='E' || str[i+4]=='e') &&            
               str[i+5]==' '                     ) {
               current=RANGE;
               start_idx=i+5;
           }   
        }
        
        /* Look ahead for "CHIP " */
        if ((onechar=='C' || onechar=='c') && i+4<len) {
           if ((str[i+1]=='H' || str[i+1]=='h') && 
               (str[i+2]=='I' || str[i+2]=='i') &&            
               (str[i+3]=='P' || str[i+3]=='p') &&            
                str[i+4]==' '                     ) {           
               current=CHIP;
               start_idx=i+4;
           }   
        }
        
        
        /* Look ahead for "COEF " */
        if ((onechar=='C' || onechar=='c') && i+4<len) {
           if ((str[i+1]=='O' || str[i+1]=='o') && 
               (str[i+2]=='E' || str[i+2]=='e') &&            
               (str[i+3]=='F' || str[i+3]=='f') &&           
                str[i+4]==' '                     ) {           
               current=COEFF;
           }   
               start_idx=i+4;
        }
        
        if (current==NONE) continue;
        
        /* Look ahead for next whitespace or end of string */
        int k=0;                
        for (int j=start_idx;j<=len;j++) {
            if (str[j] ==' ' && k==0) continue;
            if (str[j] ==' ' || str[j]=='\0') break;
            buffer[k]=str[j];
            k++;
        }
        buffer[k]='\0';
        
        /* Now parse the buffer string into an integer*/
        int ival;
        int rn = sscanf(buffer, "%d", &ival);
        if (rn<0) return rn;
        if (current==RANGE) *r_idx   =ival;
        if (current==CHIP ) *c_idx   =ival;
        if (current==COEFF) *coef_idx=ival;
        
        current=NONE;
        
        
    }
    
    /* Return -1 if failed to find relevant compoents */
    if ( *coef_idx == -1 )          return -1;
    if ( *r_idx==-1 && *c_idx==-1 ) return -1;

    return 0;
}
