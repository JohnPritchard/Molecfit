/*
 * This file is part of the ESO X-Shooter Pipeline
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/* xshooter headers */
#include <xsh_error.h>
#include <xsh_utils.h>

/* Molecfit Model */
#include "xsh_molecfit_calctrans.h"
#include <mf_wrap_config.h>
#include <telluriccorr.h>
//#include <mf_spectrum.h>
//#include <mf_wrap.h>
/*----------------------------------------------------------------------------*/
/**
 *                 Typedefs: Enumeration types
 */
/*----------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------
                             Recipe Defines
 ---------------------------------------------------------------------------*/

#define RECIPE_ID "xsh_molecfit_calctrans"
#define RECIPE_AUTHOR "N. Fernando, B. Miszalski"
#define RECIPE_CONTACT "nuwanthika.fernando@partner.eso.org"

/*---------------------------------------------------------------------------
                            Functions prototypes
 ---------------------------------------------------------------------------*/

/*
 *   Plugin initalization, execute and cleanup handlers
 */

int xsh_molecfit_calctrans_create(cpl_plugin *);
int xsh_molecfit_calctrans_exec(cpl_plugin *);
int xsh_molecfit_calctrans_destroy(cpl_plugin *);

/* The actual executor function */
int xsh_molecfit_calctrans(cpl_frameset *frameset, const cpl_parameterlist  *parlist);

/*----------------------------------------------------------------------------*/
/**
 *                 static variables
 */
/*----------------------------------------------------------------------------*/

char xsh_molecfit_calctrans_description_short[] =
"Applies molecfit_calctrans to input spectra";

char xsh_molecfit_calctrans_description[] =
"Applies molecfit_calctrans to input spectra";

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
 * @defgroup molecfit_model  It runs Molecfit on a generic input spectrum file to compute an atmospheric model.
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
 *                              Functions code
 */
/*----------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/**
  @brief    Build the list of available plugins, for this module.
  @param    list    the plugin list
  @return   0 if everything is ok, -1 otherwise

  Create the recipe instance and make it available to the application using
  the interface. This function is exported.
 */
/*--------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist *list) {
  cpl_recipe *recipe = NULL;
  cpl_plugin *plugin = NULL;

  recipe = cpl_calloc(1, sizeof(*recipe));
  if ( recipe == NULL ){
    return -1;
  }

  plugin = &recipe->interface ;

  cpl_plugin_init(plugin,
                  CPL_PLUGIN_API,                   /* Plugin API */
                  XSH_BINARY_VERSION,            /* Plugin version */
                  CPL_PLUGIN_TYPE_RECIPE,           /* Plugin type */
                  RECIPE_ID,                        /* Plugin name */
                  xsh_molecfit_calctrans_description_short, /* Short help */
                  xsh_molecfit_calctrans_description,   /* Detailed help */
                  RECIPE_AUTHOR,                    /* Author name */
                  RECIPE_CONTACT,                   /* Contact address */
                  xsh_get_license(),                /* Copyright */
                  xsh_molecfit_calctrans_create,
                  xsh_molecfit_calctrans_exec,
                  xsh_molecfit_calctrans_destroy);

  cpl_pluginlist_append(list, plugin);

  return (cpl_error_get_code() != CPL_ERROR_NONE);
 }

/*--------------------------------------------------------------------------*/
/**
  @brief    Setup the recipe options
  @param    plugin  the plugin
  @return   0 if everything is ok

  Create the recipe instance and make it available to the application using
  the interface.

 */
/*--------------------------------------------------------------------------*/

int xsh_molecfit_calctrans_create(cpl_plugin *plugin){
  cpl_recipe *recipe = NULL;

  /* Reset library state */
  xsh_init();

  /* Check input */
  assure( plugin != NULL, CPL_ERROR_NULL_INPUT, "Null plugin");

  /* Get the recipe out of the plugin */
  assure( cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE,
          CPL_ERROR_TYPE_MISMATCH,
          "Plugin is not a recipe");

  recipe = (cpl_recipe *)plugin;

  /* Create the parameter list in the cpl_recipe object */
  recipe->parameters = cpl_parameterlist_new();
  assure( recipe->parameters != NULL,
          CPL_ERROR_ILLEGAL_OUTPUT,
          "Memory allocation failed!");

  //xsh_molecfit_calctrans parameters

  //MOLECFIT_PARAMETER_USE_INPUT_KERNEL
  //use_input_kernel
  check(xsh_parameters_new_boolean(recipe->parameters,RECIPE_ID,
  MOLECFIT_PARAMETER_USE_INPUT_KERNEL,CPL_TRUE,
  "If TRUE, then the input KERNEL_LIBRARY_XXX given in the SOF is used, where XXX is UVB, VIS or NIR. If FALSE, or if the KERNEL_LIBRARY_XXX is not given, then the information stored in BEST_FIT_PARAMETERS_YYY_XXX will be used to compute the line spread function, where YYY is SCI or STD, and XXX is UVB, VIS or NIR."));

  //apply_wlc_corr
  //check(xsh_parameters_new_boolean(recipe->parameters,RECIPE_ID,
  //"APPLY_WLC_CORR",CPL_FALSE,
  //""));

  //MF_PARAMETERS_SCALE_TO_PWV
  //scale_to_obs_pwv
//check(xsh_parameters_new_string(recipe->parameters,RECIPE_ID,
//MF_PARAMETERS_SCALE_TO_PWV,"TEL AMBI PWV START",
//" \
// "));
////MOLECFIT_PARAMETER_USE_MOLEC_DATABASE
////use_molec_database
//check(xsh_parameters_new_boolean(recipe->parameters,RECIPE_ID,
//MOLECFIT_PARAMETER_USE_MOLEC_DATABASE,false,
//" \
// "));

  cleanup:
    if ( cpl_error_get_code() != CPL_ERROR_NONE ){
      xsh_error_dump(CPL_MSG_ERROR);
      return 1;
    }
    else {
      return 0;
    }
}


/*--------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*--------------------------------------------------------------------------*/

int xsh_molecfit_calctrans_exec(cpl_plugin *plugin) {
  cpl_recipe *recipe = NULL;

  /* Check parameter */
  assure( plugin != NULL, CPL_ERROR_NULL_INPUT, "Null plugin" );

  /* Get the recipe out of the plugin */
  assure( cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE,
          CPL_ERROR_TYPE_MISMATCH, "Plugin is not a recipe");

  recipe = (cpl_recipe *)plugin;
  /* Check recipe */
  xsh_molecfit_calctrans( recipe->frames, recipe->parameters);

  cleanup:
    if ( cpl_error_get_code() != CPL_ERROR_NONE ) {
      xsh_error_dump(CPL_MSG_ERROR);
      cpl_error_reset();
      return 1;
    }
    else {
      return 0;
    }
}

/*--------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*--------------------------------------------------------------------------*/
int xsh_molecfit_calctrans_destroy(cpl_plugin *plugin)
{
    cpl_recipe *recipe = NULL;

    /* Check parameter */
    assure( plugin != NULL, CPL_ERROR_NULL_INPUT, "Null plugin" );

    /* Get the recipe out of the plugin */
    assure( cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE,
            CPL_ERROR_TYPE_MISMATCH, "Plugin is not a recipe");

    recipe = (cpl_recipe *)plugin;

    xsh_free_parameterlist(&recipe->parameters);

  cleanup:
    if (cpl_error_get_code() != CPL_ERROR_NONE)
        {
            return 1;
        }
    else
        {
            return 0;
        }
}


cpl_error_code xsh_molecfit_calc_setup_frameset(cpl_frameset* frameset,cpl_parameterlist* list,const char* arm,const char* input_name){

      const char* tag;
      //check the input frame (SCIENCE); we know the tag already thanks to the input_name parameter
      //so no need to search frameset for it
      //input_name may be one of the following:
      //SCI_SLIT_FLUX_IDP_XXX
      //SCI_SLIT_FLUX_MERGE1D_XXX
      //TELL_SLIT_MERGE1D_XXX
      //TELL_SLIT_FLUX_MERGE1D_XXX
      //STD_SLIT_FLUX_IDP_YYY_XXX
      //where XXX is the arm and YYY is NOD, STARE or OFFSET

      tag = input_name;
      cpl_frame* f = cpl_frameset_find(frameset,input_name);
      cpl_frame_set_group(f,CPL_FRAME_GROUP_RAW);

      //cpl_frame_set_type(f,CPL_FRAME_TYPE_NONE);
      //cpl_frame_set_tag(f,"SCIENCE");
      //cpl_frame_set_level(f,CPL_FRAME_LEVEL_NONE);

      //other inputs
/*
      //KERNEL_LIBRARY_XXX
      tag = cpl_sprintf("%s_%s",MOLECFIT_KERNEL_LIBRARY,arm);
      f = cpl_frameset_find(frameset,tag);
      if(f){
            cpl_frame_set_group(f, CPL_FRAME_GROUP_CALIB);
      }
      //ATM_PARAMETERS_YYY_XXX
      tag = cpl_sprintf("%s_%s",MOLECFIT_ATM_PARAMETERS,arm);
      f = cpl_frameset_find(frameset,tag);
      if(f){
            cpl_frame_set_group(f, CPL_FRAME_GROUP_CALIB);
      }
      //BEST_FIT_PARAMETERS_YYY_XXX
      tag = cpl_sprintf("%s_%s",MOLECFIT_BEST_FIT_PARAMETERS,arm);
      f = cpl_frameset_find(frameset,tag);
      if(f){
            cpl_frame_set_group(f, CPL_FRAME_GROUP_CALIB);
      }
      //BEST_FIT_MODEL_YYY_XXX
      tag = cpl_sprintf("%s_%s",MOLECFIT_BEST_FIT_MODEL,arm);
      f = cpl_frameset_find(frameset,tag);
      if(f){
            cpl_frame_set_group(f, CPL_FRAME_GROUP_CALIB);
      }
      //MODEL_MOLECULES_XXX
      tag = cpl_sprintf("%s_%s",MOLECFIT_MODEL_MOLECULES,arm);
      f = cpl_frameset_find(frameset,tag);
      if(f){
            cpl_frame_set_group(f, CPL_FRAME_GROUP_CALIB);
      }
      */

}

/*--------------------------------------------------------------------------*/
/**
  @brief    Build the list of available plugins, for this module.
  @param    list    the plugin list
  @return   0 if everything is ok, -1 otherwise

  Create the recipe instance and make it available to the application using
  the interface. This function is exported.
 */
/*--------------------------------------------------------------------------*/


    /* -----------------------------------------------------------------------------------------------------/
     *
     * OVERVIEW OF STEPS
     * 		1. DO SOME ERROR CHECKING ON FRAMESET? (may be done later in setup config, but just being thorough)
     * 			1.5 GET *TAG* OF INPUT SCIENCE from FRAMSET using either MOLECFIT_STD_MODEL, MOLECFIT_SCIENCE_CALCTRANS or MOLECFIT_SCIENCE
     * 				1.55 (Check for YYY in this FILENAME, in case we need it later...)
     * 		2. DETERMINE IF SUFFIX NEEDED FROM DATASET (e.g. UVB/VIS/NIR for xshooter)
     * 		3. SETUP INPUT FILE NAME TAGS
     * 		4. SETUP OUTPUT FILE NAME TAGS
     * 		5. SETUP PARAMETERS (e.g. --list_molec, etc.)
     * 			N.B. We expect to be provided a FITS binary table of WAVE_INCLUDE & MOLECULES params,
     * 			but FOR NOW, we setup these values internally via the string parameter versions.
     * 			6. SETUP *EXTRA* PARAMETERS
     *
     *
    -------------------------------------------------------------------------------------------------------*/

    cpl_error_code xsh_molecfit_calctrans_config(cpl_frameset *frameset, const cpl_parameterlist  *parlist,
        		cpl_parameterlist* ilist, cpl_parameterlist* iframelist){

            cpl_msg_info(cpl_func,"xsh_molecfit_calctrans_config");
            cpl_msg_info(cpl_func,"FRAMESET");
            //cpl_frameset_dump(frameset,stdout);
            cpl_msg_info(cpl_func,"PARLIST");
            cpl_parameterlist_dump(parlist,stdout);
            //cpl_msg_info(cpl_func,"");

            cpl_parameterlist* iframe = cpl_parameterlist_new();
            //INPUTNAME == input tag name
            //ARM == UVB, VIS, NIR
            //OBSMODE == NOD, STARE or OFFSET; If none of these, it is set to DEFAULT
            //IDP (bool) == TRUE or FALSE
            //1. Get all the input params into the iframe parameterlist (INPUTNAME,ARM,OBSMODE,IDP)
            cpl_error_code err= CPL_ERROR_NONE;
            err=xsh_molecfit_utils_find_input_frame(frameset, iframe);


            const char* input_name = cpl_parameter_get_string(cpl_parameterlist_find(iframe,"INPUTNAME"));
            const char* arm = cpl_parameter_get_string(cpl_parameterlist_find(iframe,"ARM"));
            const char* obsmode = cpl_parameter_get_string(cpl_parameterlist_find(iframe,"OBSMODE"));
            const char* is_idp = cpl_parameter_get_string(cpl_parameterlist_find(iframe,"IDP"));
            const char* fname = cpl_parameter_get_string(cpl_parameterlist_find(iframe,"INPUTFILENAME"));

            cpl_msg_info(cpl_func,"iframe details; INPUTNAME: %s; ARM: %s; IDP: %s; OBSMODE: %s; INPUTFILENAME: %s",input_name,arm,is_idp,obsmode,fname);

            //add iframe parameters (INPUTNAME,ARM,OBSMODE,IDP) to iframelist so that we can access them from xsh_molecfit_model
            //these are not added to ilist, as they are only meant to be instrument dependent molecfit parameters
            err = cpl_parameterlist_append(iframelist,cpl_parameterlist_find(iframe,"INPUTNAME"));
            err = cpl_parameterlist_append(iframelist,cpl_parameterlist_find(iframe,"ARM"));
            err = cpl_parameterlist_append(iframelist,cpl_parameterlist_find(iframe,"OBSMODE"));
            err = cpl_parameterlist_append(iframelist,cpl_parameterlist_find(iframe,"IDP"));
            err = cpl_parameterlist_append(iframelist,cpl_parameterlist_find(iframe,"INPUTFILENAME"));

            //cpl_msg_info(cpl_func,"finishing up xsh_molecfit_model list_molec: %s",val);


            //return err;

              //5. SETUP PARAMETERS (e.g. --list_molec, etc.)
              ///UVB - if we are using the UVB arm, just assume defaults for now.
              //(We may define specific values for its parameters later)
              /*By default, these are all NULL -- just to spell it out here, so we remind ourselves that we are satisfying the requirements*/


              cpl_boolean USE_INPUT_KERNEL;
              cpl_boolean APPLY_WLC_CORR;
              //const char* SCALE_TO_OBS_PWV;
              //cpl_boolean USE_MOLEC_DATABASE;

              const char* CALCTRANS_MAPPING_KERNEL;
              char * MAPPING_ATMOSPHERIC;
              char * MAPPING_CONVOLVE;
              cpl_boolean USE_ONLY_INPUT_PRIMARY_DATA;
              int USE_DATA_EXTENSION_AS_DFLUX;
              int USE_DATA_EXTENSION_AS_MASK;
              cpl_boolean CHIP_EXTENSIONS;

              //the input framset will have _YYY_XXX suffixes ; we have to be mindful of that!
            
              
              //this maps to parameters->mapping_kernel which is a const char* 
              //cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL,CPL_TYPE_STRING,NULL,NULL,"1"));
              cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL,CPL_TYPE_STRING,NULL,NULL,"0,1"));
              cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_PARAMETER_USE_INPUT_KERNEL,CPL_TYPE_BOOL,NULL,NULL,CPL_FALSE));
              //cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_PARAMETER_CALCTRANS_MAPPING_KERNEL,CPL_TYPE_INT,NULL,NULL,1));

              const char* input_tag;
              //const char* combined_suffix = cpl_sprintf("%s_%s","SCI",arm);
              //first try "SCI" as YYY 
              //if the sof does not have ATM_PARAMETER_YYY_XXX, add this default param
              /*input_tag = mf_wrap_tag_suffix(MOLECFIT_ATM_PARAMETERS,combined_suffix,CPL_FALSE);
              if(!cpl_frameset_find(frameset,input_tag)){
                cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_ATM_PARAMETERS,CPL_TYPE_STRING,NULL,NULL,"NULL"));
              }
              //if the sof does not have BEST_FIT_PARAMETERS_YYY_XXX, add this default param
              input_tag = mf_wrap_tag_suffix(MOLECFIT_BEST_FIT_PARAMETERS,combined_suffix,CPL_FALSE);
              if(!cpl_frameset_find(frameset,input_tag)){
                cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_BEST_FIT_PARAMETERS,CPL_TYPE_STRING,NULL,NULL,"NULL"));
              }
              //if the sof does not have BEST_FIT_MODEL_YYY_XXX, add this default param
              input_tag = mf_wrap_tag_suffix(MOLECFIT_BEST_FIT_MODEL,combined_suffix,CPL_FALSE);
              if(!cpl_frameset_find(frameset,input_tag)){
                cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_BEST_FIT_MODEL,CPL_TYPE_STRING,NULL,NULL,"auto"));
              }*/

/*              input_tag = mf_wrap_tag_suffix(MOLECFIT_KERNEL_LIBRARY,arm,CPL_FALSE);
              if(!cpl_frameset_find(frameset,input_tag)){
                //cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_MODEL_MAPPING_KERNEL,CPL_TYPE_STRING,NULL,NULL,"NULL"));
                cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_PARAMETER_USE_INPUT_KERNEL,CPL_TYPE_BOOL,NULL,NULL,false));
              }
              //if we have an actual input kernel, we will have to set MOLECFIT_MODEL_MAPPING_KERNEL!

*/



              //if the sof does not have MODEL_MOLECULES_XXX, add this default param
              //input_tag = mf_wrap_tag_suffix(MOLECFIT_MODEL_MOLECULES,arm,CPL_FALSE);
              //if(!cpl_frameset_find(frameset,input_tag)){
              //  cpl_parameterlist_append(ilist,cpl_parameter_new_value(MOLECFIT_MODEL_MOLECULES,CPL_TYPE_STRING,NULL,NULL,"auto"));
             // }

              //FOR ALL CASES, even if IDP or not

              //CALCTRANS_MAPPING_KERNEL = 1;
              //cpl_propertylist_update_string(ilist,"USE_DATA_EXTENSION_AS_DFLUX",USE_DATA_EXTENSION_AS_DFLUX);
              //cpl_parameterlist_append(ilist,cpl_parameter_new_value("CALCTRANS_MAPPING_KERNEL",CPL_TYPE_INT,NULL,NULL,CALCTRANS_MAPPING_KERNEL));

              //IDP format...
              if (!strcmp(is_idp,"TRUE")){
            	  USE_ONLY_INPUT_PRIMARY_DATA = CPL_FALSE;
            	  USE_DATA_EXTENSION_AS_DFLUX = 0;
            	  MAPPING_ATMOSPHERIC = "0,1";
            	  MAPPING_CONVOLVE = "0,1";

              }
              else {
            	  USE_ONLY_INPUT_PRIMARY_DATA = CPL_TRUE;
            	  USE_DATA_EXTENSION_AS_DFLUX = 1;
            	  MAPPING_ATMOSPHERIC = "1";
            	  MAPPING_CONVOLVE = "1";

              }
              cpl_parameterlist_append(ilist,cpl_parameter_new_value("USE_ONLY_INPUT_PRIMARY_DATA",CPL_TYPE_BOOL,NULL,NULL,USE_ONLY_INPUT_PRIMARY_DATA));
              cpl_parameterlist_append(ilist,cpl_parameter_new_value("USE_DATA_EXTENSION_AS_DFLUX",CPL_TYPE_INT,NULL,NULL,USE_DATA_EXTENSION_AS_DFLUX));
              cpl_parameterlist_append(ilist,cpl_parameter_new_value("MAPPING_ATMOSPHERIC",CPL_TYPE_STRING,NULL,NULL,MAPPING_ATMOSPHERIC));
              cpl_parameterlist_append(ilist,cpl_parameter_new_value("MAPPING_CONVOLVE",CPL_TYPE_STRING,NULL,NULL,MAPPING_CONVOLVE));



              //default is 0 - do not need to set?
              USE_DATA_EXTENSION_AS_MASK = 0;
              cpl_parameterlist_append(ilist,cpl_parameter_new_value("USE_DATA_EXTENSION_AS_MASK",CPL_TYPE_INT,NULL,NULL,USE_DATA_EXTENSION_AS_MASK));


              CHIP_EXTENSIONS = CPL_FALSE ;//e-3;//0.001;
              //cpl_propertylist_update_string(ilist,"WLG_TO_MICRON",WLG_TO_MICRON);
              cpl_parameterlist_append(ilist,cpl_parameter_new_value("CHIP_EXTENSIONS",CPL_TYPE_BOOL, NULL ,NULL, CHIP_EXTENSIONS));

              //cpl_msg_info(cpl_func,"calling xsh_molecfit_model_spec_header_calcs with fname = %s",fname);
              //err = xsh_molecfit_model_spec_header_calcs(fname,arm,ilist);

              //CONTINUUM_CONST
              //this is handled in xsh_molecfit_model_spec_data_calcs()

            //free up iframe
              //have to do this at the end, otherwise strings (e.g. fname) get deleted
              //cpl_parameterlist_delete(iframe);

              return err;
     }


/*Instrument sid- example instrument - XSHOOTER?*/
int xsh_molecfit_calctrans(cpl_frameset *frameset, const cpl_parameterlist  *parlist)
{
    //cpl_errorstate_ensure(frameset && parlist, CPL_ERROR_NULL_INPUT, return, CPL_ERROR_NULL_INPUT, 'NULL input: framset and/or parlist');
    cpl_error_code err= CPL_ERROR_NONE;

    //INPUT FRAMES WILL BE NAMED DIFFERENT TO MOLECFIT_SCIENCE, MOLECFIT_STD* etc.
    //err = xsh_check_and_set_input_tags(frameset); //INSTRUMENT_SPECIFIC

    //TO WRITE
    //err = mf_wrap_calc_preconfig_check(frameset, parlist);/*Check parameters and frameset for correct inputs - returns CPL_ERROR*/

    cpl_parameterlist* ilist = cpl_parameterlist_new();
    cpl_parameterlist* iframelist = cpl_parameterlist_new();
    cpl_errorstate initial_errorstate = cpl_errorstate_get();
    // get instrument specific parameter defaults
    err = xsh_molecfit_calctrans_config(frameset,parlist,ilist,iframelist);

    //these are in iframelist, so they don't end up as input parameters for molecfit via 
    const char* input_name = cpl_parameter_get_string(cpl_parameterlist_find(iframelist,"INPUTNAME"));
    const char* arm = cpl_parameter_get_string(cpl_parameterlist_find(iframelist,"ARM"));
    const char* obsmode = cpl_parameter_get_string(cpl_parameterlist_find(iframelist,"OBSMODE"));
    const char* is_idp = cpl_parameter_get_string(cpl_parameterlist_find(iframelist,"IDP"));
    const char* fname = cpl_parameter_get_string(cpl_parameterlist_find(iframelist,"INPUTFILENAME"));
    
    //if MOLECFIT_BEST_FIT_PARAMETERS_STD
    
    const char* bf_sci = cpl_sprintf("%s_%s_%s",MOLECFIT_BEST_FIT_PARAMETERS,"SCI",arm);
    const char* bf_std = cpl_sprintf("%s_%s_%s",MOLECFIT_BEST_FIT_PARAMETERS,"STD",arm);
    const char* combined_suffix = cpl_sprintf("%s_%s",(strstr(input_name,"SCI")) ? "SCI" : "STD",arm);

    //fix for PIPE-10006
    if(cpl_frameset_find(frameset, bf_sci)){
        combined_suffix = cpl_sprintf("%s_%s","SCI",arm);
    }
    else if(cpl_frameset_find(frameset, bf_std)){
        combined_suffix = cpl_sprintf("%s_%s","STD",arm);
    }

    cpl_msg_info(cpl_func,"params retrieved from iframelist");
    

    //Calculate CONTINUUM_CONST
    //err = xsh_molecfit_model_spec_data_calcs(data,is_idp,ilist);
    
    //parlist2 is NOW called mergedlist
    cpl_parameterlist* mergedlist = cpl_parameterlist_new();

    cpl_msg_info(cpl_func,"calling mf_wrap_merge_parameterlists");
    err = mf_wrap_merge_parameterlists(ilist, parlist,mergedlist);

    cpl_msg_info(cpl_func,"calling xsh_molecfit_setup_frameset");
    err = xsh_molecfit_calc_setup_frameset(frameset,mergedlist,arm,input_name);


    //err = set_output_filename_suffix(suffix, parlist); // This could be moved to instlist create function mf_inst_config();

    cpl_msg_info(cpl_func,"calling cpl_fits_count_extensions");
    //cpl_size n_ext = cpl_fits_count_extensions(fname);

    mf_wrap_fits * data = NULL;
	data = mf_wrap_fits_load(fname, CPL_FALSE);
    
    //err = mf_wrap_calc_data(data,fname); // parlist2, parameters->use_only_input_pri_ext);
    cpl_msg_info(cpl_func,"AFTER mf_wrap_calc_data() next=%d ; %s",data->n_ext,cpl_error_get_message());

    //mf_wrap_model_parameter parameters = mf_wrap_config_init(frameset, parlist2);
    cpl_msg_info(cpl_func,"calling mf_wrap_config_calc_init");
    //cpl_propertylist_dump(data->v_ext[0].header,stdout);

    molecfit_calctrans_parameter* parameters = mf_wrap_config_calc_init(frameset,mergedlist,data->v_ext[0].header, data->n_ext, combined_suffix);
    //TODO: Do we need to add code from /* Recipe Parameters : Need scientific_header_primary */ in molecfit_calctrans.c????
    cpl_msg_info(cpl_func,"AFTER mf_wrap_config_calc_init() %s",cpl_error_get_message());
        //cpl_frame* input_frame = cpl_frameset_find(frameset,input_name);

    err = mf_wrap_data_convert_to_table( data,
                                                parameters->chip_extensions,
                                                parameters->use_only_input_pri_ext,
                                                parameters->dflux_extension_data,
                                                parameters->mask_extension_data,
                                                parameters->mf_config->parameters->inputs.column_lambda,
                                                parameters->mf_config->parameters->inputs.column_flux,
                                                parameters->mf_config->parameters->inputs.column_dflux,
                                                parameters->mf_config->parameters->inputs.column_mask);
    cpl_msg_info(cpl_func,"AFTER mf_wrap_data_convert_to_table() %s",cpl_error_get_message());



    cpl_msg_info(cpl_func,"calling mf_wrap_calc_molecules");
    //LOAD MODEL_MOLECULES
    cpl_table *molecules = NULL;
    err= mf_wrap_calc_molecules(frameset,&molecules,arm);

    //LOAD CALCTRANS_MAPPING_KERNEL // CALCTRANS_KERNEL_LIBRARY
    mf_wrap_fits *kernel_data    = NULL;
    cpl_table       *mapping_kernel = NULL;

    err= mf_wrap_calc_kernel_library(kernel_data,&mapping_kernel,frameset,parameters,arm);

    //LOAD MAPPING_ATMOSPHERIC
    cpl_table *mapping_atmospheric = NULL;
    err= mf_wrap_calc_mapping_atm(frameset,&mapping_atmospheric,parameters,arm);
    
    //LOAD MAPPING_CONVOLVE
    cpl_table *mapping_convolve = NULL;
    err= mf_wrap_calc_mapping_conv(frameset,&mapping_convolve,parameters,arm);
  /*** Save generic output files */
  if (!err) {

      cpl_msg_info(cpl_func, "Save generic multi-extension output FITS file ('%s','%s') ...", MOLECFIT_TELLURIC_DATA, MOLECFIT_TELLURIC_CORR);

      const char* output_fname = mf_wrap_tag_suffix(MOLECFIT_LBLRTM_RESULTS,arm,CPL_TRUE);
      const char* tag_name = mf_wrap_tag_suffix(MOLECFIT_LBLRTM_RESULTS,arm,CPL_FALSE);
      err     += mf_wrap_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, tag_name, output_fname);

      output_fname = mf_wrap_tag_suffix(MOLECFIT_TELLURIC_DATA,arm,CPL_TRUE);
      tag_name = mf_wrap_tag_suffix(MOLECFIT_TELLURIC_DATA,arm,CPL_FALSE);
      err     += mf_wrap_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, tag_name, output_fname);

      output_fname = mf_wrap_tag_suffix(MOLECFIT_TELLURIC_CORR,arm,CPL_TRUE);
      tag_name = mf_wrap_tag_suffix(MOLECFIT_TELLURIC_CORR,arm,CPL_FALSE);
      err     += mf_wrap_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl, tag_name, output_fname);

      if (kernel_data) {
      output_fname = mf_wrap_tag_suffix(MOLECFIT_CALCTRANS_KERNEL_LIBRARY,arm,CPL_TRUE);
      tag_name = mf_wrap_tag_suffix(MOLECFIT_CALCTRANS_KERNEL_LIBRARY,arm,CPL_FALSE);
          err += mf_wrap_save(frameset, frameset, parlist, RECIPE_NAME, parameters->pl,  tag_name, output_fname);
      }
  }

  cpl_msg_info(cpl_func,"AFTER saving lots of files %s",cpl_error_get_message());
  /*** Execute molecfit (lblrtm) ***/
  if (!err) {

      int null;

      /* Execution extensions */
      cpl_size n_ext;
      if (     parameters->use_only_input_pri_ext) n_ext = 1;
      else if (parameters->chip_extensions    ) n_ext = data->v_ext[0].spectrum_data ? 1 : 2;
      else n_ext = data->n_ext;

      cpl_msg_info(cpl_func,"AFTER Execution extensions %s",cpl_error_get_message());

      for (cpl_size ext = 0; ext < n_ext && !err; ext++) {
          cpl_msg_info(cpl_func,"look over extensions %lld %s",ext,cpl_error_get_message());

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
          cpl_table        *lblrtm_spec_out = NULL;

          cpl_table        *telluric_data   = NULL;
          cpl_vector       *telluric_vec    = NULL;

          cpl_matrix       *kernel_matrix   = NULL;

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

              /* Wrap the data */
              double *telluric_corr_column = cpl_table_get_data_double(telluric_data, MF_COL_OUT_TELLURIC_CORR);

              /* Convert in a spectrum (cpl_vector) */
              cpl_vector *vAux = cpl_vector_wrap(cpl_table_get_nrow(telluric_data), telluric_corr_column);
              telluric_vec = cpl_vector_duplicate(vAux);
              cpl_vector_unwrap(vAux);

              /* Get Kernel results */
              if (parameters->results_convolution[ext]->kernel_resampled_normalized) {
                  kernel_matrix = parameters->results_convolution[ext]->kernel_resampled_normalized;
              }
          }

          /*** SAVE results : LBLRTM_RESULTS, TELLURIC_DATA, TELLURIC_CORR and KERNEL_MOLECFIT results ***/
          if (!err && (parameters->use_only_input_pri_ext || ext > 0)) {

              cpl_msg_info(cpl_func, "Saving %s, %s, %s, %s ... (ext =%lld : lblrtm_results ? %d, convolution_results ? %d)", MOLECFIT_LBLRTM_RESULTS, MOLECFIT_CALCTRANS_KERNEL_LIBRARY, MOLECFIT_TELLURIC_DATA, MOLECFIT_TELLURIC_CORR, ext, lblrtm_spec_out != NULL, telluric_data != NULL);

              const char* output_fname = mf_wrap_tag_suffix(MOLECFIT_LBLRTM_RESULTS,arm,CPL_TRUE);
              err     += mf_wrap_save_mf_results( lblrtm_header,                   output_fname,           CPL_FALSE, NULL,          lblrtm_spec_out, NULL         );
              output_fname = mf_wrap_tag_suffix(MOLECFIT_TELLURIC_DATA,arm,CPL_TRUE);
              err     += mf_wrap_save_mf_results( data->v_ext[ext].header,         output_fname,            CPL_FALSE, NULL,          telluric_data,   NULL         );
              output_fname = mf_wrap_tag_suffix(MOLECFIT_TELLURIC_CORR,arm,CPL_TRUE);
              err     += mf_wrap_save_mf_results( data->v_ext[ext].header,         output_fname,            CPL_FALSE, NULL,          NULL,            telluric_vec );

              if (kernel_data) {
                  output_fname = mf_wrap_tag_suffix(MOLECFIT_CALCTRANS_KERNEL_LIBRARY,arm,CPL_TRUE);
                  err += mf_wrap_save_mf_results( kernel_data->v_ext[ext].header,  output_fname, CPL_FALSE, kernel_matrix, NULL,            NULL         );
              }
          }

          /* Cleanup */
          if (lblrtm_header) cpl_propertylist_delete(lblrtm_header);
          if (telluric_vec ) cpl_vector_delete(telluric_vec);
      }
  }

  /* Cleanup */
  if (parameters         ) molecfit_calctrans_parameter_delete( parameters         );
  if (data               ) mf_wrap_fits_delete(                data               );
  if (molecules          ) cpl_table_delete(                    molecules          );
  if (kernel_data        ) mf_wrap_fits_delete(                kernel_data        );
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

   

 
    //----------------
    //ROUGH OUTLINE OF MOLECFIT_CALCTRANS()
    //-------------------------

    //initial frameset checks 
    //BLOCK: 
  /* Check mandatory TAGS/Parameters */
  //to go into mf_wrap_calc_preconfig_check
      //SCIENCE_CALCTRANS/SCIENCE; 
      //MODEL_MOLECULES
      //ATM_PARAMETERS
      //BEST_FIT_PARAMETERS
      //MAPPING_CONVOLVE
      //TEST KERNEL MAPPING
      //TEST _ATMOSPHERIC: MAPPING 
      //TEST _CONVOLVE: MAPPING

//Separate functions  for each ; TODO: Check if we can just reuse the same function from molecfit_model (if we change parameter == the TAG ==  to these functions)
  /* Load TAG = SCIENCE_CALCTRANS/SCIENCE */

  /* Load TAG = MOLECULES */ //SAME
  /* Load TAG : CALCTRANS_KERNEL_LIBRARY */
  /* Load TAG = MAPPING_ATMOSPHERIC */
  /* Load TAG = MAPPING_CONVOLVE */

//Parameters: molecfit_calctrans_parameter* parameters
//call to mf_wrap_data_convert_to_table

//Saving output files  

//DO THE CALCTRANS STUFF
  /*** Execute molecfit (lblrtm) ***/
  //1. mf_calctrans_lblrtm
  //2. mf_calctrans_convolution

//SAVE RESULTS

    //Load up calctrans *input files*
    //ATM_PARAMETERS_(SCI/STD)_(UVB/VIS/NIR) - molecfit_calctrans_parameters()
    //BEST_FIT_PARAMETERS_(SCI/STD)_(UVB/VIS/NIR) - molecfit_calctrans_parameters()
    //BEST_FIT_MODEL_(SCI/STD)_(UVB/VIS/NIR) - 
    //MODEL_MOLECULES_(UVB/VIS/NIR) - mf_wrap_calc_molecules()
    //KERNEL_LIBRARY_(UVB/VIS/NIR) - mf_wrap_calc_kernel_library()



    
    //Load up calctrans *parameters*
    //DO ON TELLURICCORR SIDE -> mf_wrap_calc.c

    //CALCTRANS_MAPPING_KERNEL
    //MAPPING_ATMOSPHERIC
    //MAPPING_CONVOLVE

    //USE_ONLY_INPUT_PRIMARY_DATA
    //USE_DATA_EXTENSION_AS_DFLUX
    //USE_DATA_EXTENSION_AS_MASK
    //CHIP_EXTENSIONS













    return err;

}


//cpl_error_code xsh_check_and_set_input_tags(frameset){
	//map the XSHOOTER input frame tags to MOLECFIT SCIENCE etc.

	//check the frameset for the following tags =>
	//SCI_SLIT_FLUX_IDP_XXX, SCI_SLIT_FLUX_MERGE1D_XXX,
	//TELL_SLIT_FLUX_MERGE_XXX, TELL_SLIT_MERGE1D_XXX,
	//STD_SLIT_FLUX_IDP_YYY_XXX
	//where XXX = VIS/NIR/ UBB, and YYY= NOD, STARE, OFFSET

	//call a function in Telluriccorr that maps out the input frame names to molecfit defaults();

//   return CPL_ERROR_OK;
//}


/// -----------------------------------INSTRUMENT SIDE END-----------------------------------///
