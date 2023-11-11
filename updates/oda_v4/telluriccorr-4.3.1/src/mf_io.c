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

#define _XOPEN_SOURCE      700         /* For nftw(), stpcpy(), mkdtemp() */
#define _DARWIN_C_SOURCE               /* macOS mkdtemp() is not available if _POSIX_C_SOURCE=200809L (Apple bug report #35851865) */

#include <ctype.h>
#include <netdb.h>
#include <ftw.h>
#include <fcntl.h>

#include "mf_molecules.h"

#include "mf_io.h"

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

/* Create a symbolic link */
static cpl_error_code mf_io_symlink(
    const char               *target,
    const char               *linkpath);

/* Callback for mf_io_rm_rf function */
static int mf_io_nftw_rm_rf(
    const char               *path_name,
    const struct stat        *stat_buf,
    int                      file_type,
    struct FTW               *ftw_buf);

/*  */
static int mf_io_get_socket_connection(
    const char               *host,
    const char               *port);

/*  */
static int mf_io_ftp_reply(
    int                      sockfd,
    char                     **message);

/*  */
static int mf_io_verify_ftp_code(
    char                     *msg,
    int                      length);

/*  */
static int mf_io_send_ftpcmd(
    int                      sockfd,
    const char                *cmd);

/*  */
static int mf_io_send_pasv(
    int                      sockfd,
    const char               *cmd);

/*  */
static char * mf_io_get_ftp_file(
    int                      sockfd,
    int                      *data_length);

/* Trim an input string in-place */
static void mf_io_str_trim(
    char                     *str);           /* 'char *' to trim spaces from start and end */

/*----------------------------------------------------------------------------*/
/**
 *                 Functions
 */
/*----------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/**
 * @defgroup mf_io   Tools for IO and POSIX calls.
 *
 * @brief
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* ---------------------------------------------------------------------------*/
/**
  * @brief Get the absolute current working directory.
  *
  * @return char string
  *
  * @note Must be deleted by the user with cpl_free.
  *
  */
 /* ---------------------------------------------------------------------------*/
 char * mf_io_pwd(void)
 {
     size_t size = MF_LEN_MAX;
     char * buf;
     errno = 0;

     /* if only we could use sane GNU functions instead of this posix crap */
     while (1) {

         buf = cpl_calloc(size, sizeof(char));
         if (getcwd(buf, size) != 0) {

             break;

         } else if (errno == ERANGE) {

             /* increase buffer, repeat */
             errno = 0;
             size *= 2;
             cpl_free(buf);

         } else {

             cpl_free(buf);
             cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                   "Could not determine current working directory: %s",
                                   strerror(errno));
             return NULL;
         }
     }

     return buf;
 }

 /* ---------------------------------------------------------------------------*/
 /**
  * @brief Get value from a environment variable, NULL if not exist.
  *
  * @param env
  *
  * @return Value of the environment variable or NULL if not exist.
  *
  */
 /* ---------------------------------------------------------------------------*/
 const char * mf_io_getenv(
     const char               *env)
 {
     return getenv(env);
 }

 /* Synchronize memory to disk */
 /* ---------------------------------------------------------------------------*/
 /**
  * @brief Synchronize memory to disk.
  *
  * @return result
  *
  */
 /* ---------------------------------------------------------------------------*/
 void mf_io_sync(void)
 {
     sync();
 }

 /* ---------------------------------------------------------------------------*/
 /**
  * @brief Check if the file exist in the disk
  *
  * @param file
  *
  */
 /* ---------------------------------------------------------------------------*/
 cpl_error_code mf_io_access(
     const char               *file)
 {
     return access(file, F_OK);
 }

/* ---------------------------------------------------------------------------*/
/**
 * @brief Wrapper from system_calls
 *
 * @param command         System call for execute
 * @param path            NULL or path for the execution
 * @param runtime         NULL o pointer to double for storage the time spend in the execution
 *
 * @note The return int code is converted in a CPL_ERROR_CODE, be careful doesn't match with the significant
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_io_system(
    const char               *command,
    const char               *path,
    double                   *runtime)
{
  /* Check command */
  if (!command) return CPL_ERROR_NULL_INPUT;

  /* Debug info */
  cpl_msg_debug(cpl_func, "(mf_io        ) %s", command);

  /* Get current directory */
  char *current = (path) ? mf_io_pwd() : NULL;

  /* Change to the work directory */
  cpl_error_code err = CPL_ERROR_NONE;
  if (path && current) err = chdir(path);
  if (err) cpl_msg_error(cpl_func,
               "Cannot change to directory %s", current);

  /* Execute command getting time */
  if (!err) {

      double ts = cpl_test_get_walltime();
      err       = system(command);
      double te = cpl_test_get_walltime();
      if (err) cpl_msg_error(cpl_func,
               "System Call %s failed!", command);

      if (runtime) *runtime = te - ts;
  }

  /* Return to the currrent directory */
  /* (Check for error but do not overwrite any previous error)*/
  cpl_error_code err2 = CPL_ERROR_NONE;
  if (path && current) err2 = chdir(current);
  if (err2) cpl_msg_error(cpl_func,
               "Cannot change to directory %s", current);
  if (!err && err2) err=err2;

  /* Cleanup */
  if (current) cpl_free(current);

  return err;
}
/* ---------------------------------------------------------------------------*/
/**
 * @brief Alternative Wrapper for system_calls that is safer with OpenMP
 *
 * @param command         System call for execute
 * @param path            NULL or path for the execution
 * @param runtime         NULL o pointer to double for storage the time spend in the execution
 *
 * @note The return int code is converted in a CPL_ERROR_CODE, be careful doesn't match with the significant
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_io_systemOpenMP (
    const char               *command,
    const char               *path,
    double                   *runtime)
{
  /* Check command */
  if (!command) return CPL_ERROR_NULL_INPUT;

  /* Debug info */
  cpl_msg_debug(cpl_func, "(mf_io        ) %s", command);

  double ts = cpl_test_get_walltime();
  
  /* Create a new system command that prefixes the old with "cd [path];"*/
  char *newcommand;
  int size;
  size=4;
  size=size+strlen(path);
  size=size+strlen(command);
  newcommand = malloc(sizeof(char)*(size+1));
  newcommand=strcpy(newcommand, "cd ");
  newcommand=strcat(newcommand, path);
  newcommand=strcat(newcommand, ";");
  newcommand=strcat(newcommand, command);
  
  /* Call the new system command*/
  cpl_error_code err = CPL_ERROR_NONE;
  err       = system(newcommand);
  if (err) cpl_msg_error(cpl_func,
           "System Call %s failed!", command);

  double te = cpl_test_get_walltime();

  if (runtime) *runtime = te - ts;

  /* Cleanup */
  free (newcommand);

  return err;
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief Remove a file from the disk
 *
 * @param file
 *
 * @note  It is possible to use also : unlink(file)
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_io_rm(
    const char               *file)
{
    /* Check file */
    if (!file) return CPL_ERROR_NULL_INPUT;

    return remove(file);
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Move one file between paths
 *
 * @param source_path
 * @param dest_path
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_io_mv(
    const char               *source_path,
    const char               *dest_path)
{
  /* Check file */
  if (!source_path || !dest_path) return CPL_ERROR_NULL_INPUT;
  
  int rename_err;
  rename_err=rename(source_path, dest_path);
  
  return rename_err;
  
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Remove a path recursively without follow symbolic links
 *
 * @param path           .
 * @param max_dir_depth  .
 *
 * @note Continue if error or path doesn't exist --> Similar to $rm -rf [path]
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_io_rm_rf(
    char                     *path,
    int                      max_dir_depth)
{
    /* Check path */
    if (!path) return CPL_ERROR_NULL_INPUT;

    return nftw(path, mf_io_nftw_rm_rf, max_dir_depth, FTW_DEPTH | FTW_PHYS);
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Remove a path recursively without follow symbolic links
 *
 * @param config_tmp_path  .
 * @param default_tmp_path .
 *
 * @return cpl_boolean     Indicate if the new folder was created or not
 *
 */
/* ---------------------------------------------------------------------------*/
char * mf_io_mkstemp(
    const char               *config_tmp_path,
    const char               *default_tmp_path)
{
  /* Get generic temporary name for mkstemp */
  char *tmp_path = NULL;
  if (config_tmp_path) {

      /* Defined in the config_tmp_path */
      cpl_msg_info(cpl_func, "(mf_io        ) Using temporary path defined by default/user = %s", config_tmp_path);
      tmp_path = cpl_sprintf("%s", config_tmp_path);

  } else {

      const char *tmpdir_env = mf_io_getenv(MF_IO_TMP_FOLDER_ENV);
      if (tmpdir_env) {

          /* Get the value from the MF_IO_TMP_FOLDER_ENV environment variable */
          cpl_msg_warning(cpl_func, "(mf_io        ) Using temporary path defined by environment variable %s = %s", MF_IO_TMP_FOLDER_ENV, tmpdir_env);
          tmp_path = cpl_sprintf("%s", tmpdir_env);

      } else {

          /* Using the temporary value specified in the call */
          cpl_msg_warning(cpl_func, "(mf_io        ) Using temporary path defined by default in telluriccorr = %s", default_tmp_path);
          tmp_path = cpl_sprintf("%s", default_tmp_path);
      }
  }

  /* Create New temporary directory */
  char *tmp_folder = NULL;
  for (cpl_size attempt = 1; attempt <= MF_IO_TMP_FOLDER_MAX_ATTEMPTS && !tmp_folder; attempt++) {

      /* Name of the temporary file */
      tmp_folder = cpl_sprintf("%s/"MF_IO_TMP_FOLDER_INIT, tmp_path);

      /* Create the temporary file */
      int fd = mkstemp(tmp_folder);
      if (fd < 0) {

          cpl_free(tmp_folder);
          tmp_folder = NULL;
          break;

      } else {

          /* Close the file */
          close(fd);

          /* Delete temporary file */
          mf_io_rm(tmp_folder);

          /* Create and check temporary directory */
          if (mf_io_mkdir(tmp_folder) != CPL_ERROR_NONE) {
              cpl_msg_warning(cpl_func, "Temporary directory creation failed [Name = %s, Attempt = %lld] !", tmp_folder, attempt);
              cpl_free(tmp_folder);
              tmp_folder = NULL;
          }
      }
  }

  /* Cleanup */
  cpl_free(tmp_path);

  return tmp_folder;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Remove a path recursively without follow symbolic links
 *
 * @param new_folder       .
 *
 * @return cpl_boolean     Indicate if the new folder was created or not
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_io_mkdir(
    const char               *new_folder)
{
  /* Check path */
  if (!new_folder) return CPL_ERROR_NULL_INPUT;

  return mkdir(new_folder, S_IRWXU);
}

/*----------------------------------------------------------------------------*/
/**
 * @brief    Retrieve the Earth Orientation Parameters computed by IERS
 *
 * @param    ftp_host     The FTP host to retrieve the data from
 * @param    url_path     The full path to the data file
 * @param    local_dst    The total size of the data buffer returned (returned)
 *
 * @return   cpl_error_code   CPL_ERROR_NONE or the next in case of error:
 *                            - CPL_ERROR_NULL_INPUT if eop_host, data_length or url_path are NULL
 *                            - CPL_ERROR_DATA_NOT_FOUND if the connection to the host cannot be successfully established.
 *                            - CPL_ERROR_DATA_NOT_FOUND if the FTP transaction cannot be fullfilled.
 *
 * @note This function will connect to a given FTP host specified in ftp_host
 *         and the given eop_urlpath and retrieve the ascii file with the EOP data.
 *
 */
/*----------------------------------------------------------------------------*/
cpl_error_code mf_io_curl(
    const char               *ftp_host,
    const char               *url_path,
    const char               *local_dst)
{
    const char ftp_port[] = "21";
    int cmd_socket, data_socket;

    /* Check and dump the input */
    cpl_ensure (ftp_host,    CPL_ERROR_NULL_INPUT, CPL_ERROR_NULL_INPUT);
    cpl_ensure (url_path, CPL_ERROR_NULL_INPUT, CPL_ERROR_NULL_INPUT);
    cpl_msg_debug(cpl_func, "(mf_io        ) Using URL ftp://%s%s", ftp_host, url_path);

    /* Getting the communication socket. */
    cmd_socket = mf_io_get_socket_connection(ftp_host, ftp_port);
    if (cmd_socket == 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Couldn't connect to the host");
        return CPL_ERROR_FILE_IO;
    }

    if (!mf_io_ftp_reply(cmd_socket, NULL)) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "FTP server didn't reply");
        close(cmd_socket);
        return CPL_ERROR_FILE_IO;
    }

    cpl_msg_debug(cpl_func, "(mf_io        ) SEND");
    if (!mf_io_send_ftpcmd(cmd_socket, "USER anonymous\n")) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Failed to send anonymous user");
        close(cmd_socket);
        return CPL_ERROR_FILE_IO;
    }

    if (!mf_io_send_ftpcmd(cmd_socket, "PASS ftp@eso.org\n")) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Failed to send pasword");
        close(cmd_socket);
        return CPL_ERROR_FILE_IO;
    }

    int data_port = mf_io_send_pasv(cmd_socket, "PASV\n");
    if (!data_port) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Failed to get data port");
        close(cmd_socket);
        return CPL_ERROR_FILE_IO;
    }

    /* Getting the data socket in passive mode */
    char data_port_s[256];
    snprintf(data_port_s, 255, "%d", data_port);
    data_socket = mf_io_get_socket_connection(ftp_host, data_port_s);
    if (data_socket == 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Couldn't open ftp data connection");
        return CPL_ERROR_FILE_IO;
    }

    /* Retrieving the file */
    if (!mf_io_send_ftpcmd(cmd_socket, "TYPE I\n")) {
        return CPL_ERROR_FILE_IO;
    }

    char * retr_command = cpl_malloc(strlen(url_path) + 7);
    snprintf(retr_command, strlen(url_path) + 7, "RETR %s\n", url_path);
    if(!mf_io_send_ftpcmd(cmd_socket, retr_command)) {
        close(cmd_socket);
        close(data_socket);
        cpl_free(retr_command);
        return CPL_ERROR_FILE_IO;
    }
    cpl_free(retr_command);

    int  data_length = -1;
    char *data = mf_io_get_ftp_file(data_socket, &data_length);

    /* Close connection and free resources */
    close(cmd_socket);
    close(data_socket);

    cpl_msg_info(cpl_func, "(mf_io    ) Download [ftp://%s%s] in : %s", ftp_host, url_path, local_dst);

    cpl_error_code err = CPL_ERROR_NONE;
    int fd = open(local_dst, O_CREAT | O_WRONLY, S_IRWXU);
    if (fd >= 0) {
        if (write(fd, data, data_length - 1)) {};
        close(fd);
    } else {
        perror("ERROR: ");
        err = CPL_ERROR_FILE_IO;
    }

    /* Cleanup */
    cpl_free(data);

    return err;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Get number of lines in gdas tarball.
 *
 * @param path               File path of tarball.
 *
 * @return long              Number of files in tarball, 0 on error.
 *
 */
/* ---------------------------------------------------------------------------*/
long mf_io_tarball_nfiles(
    const char               *path)
{
    /* Open tarball */
    char *tar_sys = cpl_sprintf("tar -tf \"%s\"  2>/dev/null | wc -l", path);
    cpl_msg_info(cpl_func, "(mf_io        ) Load TAR file: %s (mf_gdas_get_tarball_nlines)", tar_sys);
    FILE *stream = popen(tar_sys, "r");
    cpl_free(tar_sys);
    if (!stream) {
        return 0;
    }

    /* Read number of lines */
    long n_lines = 0;
    char line[MF_LEN_MAX];
    if (fread(line, 1, MF_LEN_MAX, stream) > 0) {
        char *endptr;
        n_lines = strtol(line, &endptr, 10);
    }
    pclose(stream);

    return n_lines;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Read GDAS profile.
 *
 * @param gdas_file_ASCII    Filename of an atmospheric GDAS profile.
 * @param hgt_units          Height unit MF_UNIT_DIST.
 *
 * @return cpl_table         GDAS profile.
 *
 * @note This function reads an atmospheric GDAS profile into a CPL_TABLE.
 *              The output CPL_TABLE must exist before calling this routine.
 *              It will get resized and overwritten.
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_table * mf_io_read_gdas_file_and_create_table(
    const char               *gdas_file_ASCII,
    const char               *hgt_units)
{
    /* open file for reading */
    cpl_msg_info(cpl_func, "(mf_io        ) Load ASCII file: %s", gdas_file_ASCII);
    FILE  *stream = fopen(gdas_file_ASCII, "r");
    if (!stream) {
        cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                              "Could not open gdas_file: %s",
                              gdas_file_ASCII);
        return NULL;

    } else {

        char line[MF_LEN_MAX];

        /* Find header (non-data lines at beginning of file) */
        static char *save;
        long        nhead = 0;
        while (fgets(line, MF_LEN_MAX - 1, stream) != NULL) {

            /* Remove leading and trailing blanks from string */
            mf_io_str_trim(line);

            char   *ptr    = strtok_r(line, " \t", &save);
            double vals[4] = {0., 0., 0., 0.};

            for (cpl_size i = 0; i < 4; i++) {
                if (ptr == NULL) break;
                vals[i] = atof(ptr);
                ptr = strtok_r(NULL, " \t", &save);
            }

            if (vals[0] + vals[1] + vals[2] + vals[3] == 0) {
                nhead++;
            } else {
                // first data line found
                break;
            }
        }
        rewind(stream);

        /* Skip header lines */
        for (cpl_size i = 0; i < nhead; i++) {
            if (fgets(line, MF_LEN_MAX - 1, stream)) {}
        }

        /* Count data lines (excluding comments) */
        cpl_size nrows = 0;
        while (fgets(line, MF_LEN_MAX - 1, stream) != NULL) {
            mf_io_str_trim(line);
            if (line[0] != '#') nrows++;
        }

        /* Initialize GDAS profile tables and tag array */
        cpl_table *gdas_profile = cpl_table_new(1);
        cpl_table_set_size(gdas_profile, nrows);

        cpl_table_new_column(gdas_profile, MF_COL_GDAS_PRESS,  CPL_TYPE_DOUBLE);
        cpl_table_new_column(gdas_profile, MF_COL_GDAS_HEIGHT, CPL_TYPE_DOUBLE);
        cpl_table_new_column(gdas_profile, MF_COL_GDAS_TEMP,   CPL_TYPE_DOUBLE);
        cpl_table_new_column(gdas_profile, MF_COL_GDAS_RELHUM, CPL_TYPE_DOUBLE);

        rewind(stream);

        /* Skip header lines */
        for (cpl_size i = 0; i < nhead; i++) {
            if (fgets(line, MF_LEN_MAX - 1, stream)) {}
        }

        cpl_size i = 0;
        while (fgets(line, MF_LEN_MAX - 1, stream) != NULL) {

            mf_io_str_trim(line);

            /* skip comments */
            if (line[0] == '#') {
                continue;
            }

            char *ptr;

            ptr = strtok_r(line, " \t", &save);
            cpl_table_set_double(    gdas_profile, MF_COL_GDAS_PRESS,  i, atof(ptr)        );

            /* height in km */
            ptr = strtok_r(NULL, " \t", &save);
            if (strcmp(hgt_units, "m") == 0) {
                cpl_table_set_double(gdas_profile, MF_COL_GDAS_HEIGHT, i, atof(ptr) / 1000.);
            } else {
                cpl_table_set_double(gdas_profile, MF_COL_GDAS_HEIGHT, i, atof(ptr)        );
            }

            ptr = strtok_r(NULL, " \t", &save);
            cpl_table_set_double(    gdas_profile, MF_COL_GDAS_TEMP,   i, atof(ptr)        );

            ptr = strtok_r(NULL, " \t", &save);
            cpl_table_set_double(    gdas_profile, MF_COL_GDAS_RELHUM, i, atof(ptr)        );
            i++;
        }

        fclose(stream);

        return gdas_profile;
    }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create execute LNFL configuration.
 *
 * @param data_path         .
 * @param w_dir             Working directory to save the file
 * @param wn_start          Initial wavelength
 * @param wn_end            End     wavelength
 * @param lbl_molecs        Char array when each position is a activation flag for this molecule
 * @param config            Configuration LNFL structure
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_error_code mf_io_write_lnfl_configuration(
    const char               *data_path,
    const char               *w_dir,
    const double             wn_start,
    const double             wn_end,
    const char               *lbl_molecs,
    const mf_io_lnfl_config  *config)
{
  /* Create symbolic link to HITRAN in working directory */
  char *tape1        = cpl_sprintf("%s/%s/%s", data_path, MF_HITRAN_PATH, config->line_db);
  char *w_dir_TAPE1  = cpl_sprintf("%s/%s", w_dir, MF_AER_TAPE1_FILE);
  cpl_error_code err = mf_io_symlink(tape1, w_dir_TAPE1);
  cpl_free(w_dir_TAPE1);
  cpl_free(tape1);
  if (err != CPL_ERROR_NONE) {
      return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                   "could not create symbolic link to %s", MF_AER_TAPE1_FILE);
  }

  /* Create MF_AER_TAPE5 file :
   *  1st line: 72 chars user info (e.g.: "$ f100 format")
   *  2nd line: lower & upper wavenumber (25cm-1 wider than required)
   *            format: F10.3,  F10.3 (e.g.: "   300.      3500.")
   *  3rd line: molecule indicator & hollerith indicator
   *            format: 39I1,3X,     A40
   *            (e.g.: "111111111111111111111111111111111111111   NBLK1 LNOUT")
   *  4th line: fixed: "%%%%%%%%%%%%%%%%%%"
   *  5th line: fixed: "1234567890123456789012345678901234567890"\
   *                   "1234567890123456789012345678901234567890"
   */

  char *w_dir_TAPE5 = cpl_sprintf("%s/%s", w_dir, MF_AER_TAPE5_FILE);
  cpl_msg_info(cpl_func, "(mf_io        ) Write output ASCII file: %s", w_dir_TAPE5);
  FILE *stream = fopen(w_dir_TAPE5, "w");
  cpl_free(w_dir_TAPE5);

  if (!stream) {

      return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                   "Could not open %s for writing",
                                   MF_AER_TAPE5_FILE);
  } else {

      /* NOTE: The number of molecules (MF_MOLEC_NUMBER in mf_molecules.h) need to match with the line 3 --> Now, 47 but can change with new versions of LBLRTM */

      if (   fprintf(stream, "$ created by lnfl\n"                                              ) == EOF    /* 1st line */
          || fprintf(stream, "%10.3f%10.3f\n",                   wn_start, wn_end               ) == EOF    /* 2nd line */
          || fprintf(stream, "%47s    LNOUT F%i\n",              lbl_molecs, config->line_db_fmt) == EOF    /* 3rd line */
          || fputs("%%%%%%%%%%%%%%%%%%\n",                       stream                         ) == EOF    /* 4th line */
          || fputs("1234567890123456789012345678901234567890"
                   "1234567890123456789012345678901234567890\n", stream                         ) == EOF) { /* 5th line */

          fclose(stream);
          return cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                                       "Unexpected file structure");
      }
  }
  fclose(stream);

  return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Create execute LBLRTM configuration.
 *
 * @param w_dir             Working directory to save the file
 * @param tape3
 * @param V1
 * @param V2
 * @param vbar
 * @param angle             Configuration ANGLE for LBLRTM
 * @param emission_spec     Flag : TRUE:Emission / FALSE:Transmission
 * @param lbl_molecs        Char array when each position is a activation flag for this molecule
 * @param config            Configuration LBLRTM structure
 * @param prof              cpl_table profile
 *
 * @return cpl_error_code   .
 */
/* ---------------------------------------------------------------------------*/
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
    const cpl_table          *prof)
{
    /* Create a symbolic link to MF_AER_TAPE3 in the working directory */
    char       *w_dir_TAPE3 = cpl_sprintf("%s/%s", w_dir, MF_AER_TAPE3_FILE);
    cpl_error_code err      = mf_io_symlink(tape3, w_dir_TAPE3);
    cpl_free(w_dir_TAPE3);
    if (err != CPL_ERROR_NONE) {
        err = cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                    "Could not create link from %s to %s/%s",
                                    tape3, w_dir, MF_AER_TAPE3_FILE);
        return err;
    }

    /* Flag recording successful writing of MF_AER_TAPE5 file */
    cpl_boolean success = CPL_TRUE;

    /* Create MF_AER_TAPE5 file */
    char *w_dir_TAPE5  = cpl_sprintf("%s/%s", w_dir, MF_AER_TAPE5_FILE);
    cpl_msg_info(cpl_func, "(mf_io        ) Write output ASCII file: %s", w_dir_TAPE5);
    FILE *stream = fopen(w_dir_TAPE5, "w");
    cpl_free(w_dir_TAPE5);
    if (!stream) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                     "Could not open %s for writing",
                                     MF_AER_TAPE5_FILE);
    } else {

        /*** LBLRTM record 1.1 ***/

        if (fputs("$ created by lblrtm_start1\n", stream) == EOF) {
            success = CPL_FALSE;
        }


        /*** LBLRTM record 1.2 ***/

        /* Format string in the FORTRAN code
         * 4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1, 4X,I1, 4X,I1, 4X,I1, 3X,A2, 4X,I1, 4X,I1,  4X,I1, 1X,I4, 1X,I4,  4X,I1,4X,I1
         * 4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1, 4X,I1, 4X,I1, 4X,I1, 3X,A2, 4X,I1, 4X,I1,  4X,I1, 1X,I4, 1X,I4 */
        if (fprintf(stream, "    %1i    %1i    %1i    %1i    %1i    %1i    %1i"
                    "    %1i    %1i    %1i   %2s    %1i    %1i    %1i %4i %4i    %1i    %1i"
                    "\n", 1, 1, config->icntnm, config->iaersl, 1, 0, 0, 0, 0, 1, " 0", 0, 0,
                    0, config->mpts, config->npts, 0, 0) < 0) {
            success = CPL_FALSE;
        }


        /*** LBLRTM record 1.3 ***/

        /* Format string in the FORTRAN code
         * E10.3,  E10.3,    E10.3,   E10.3,   E10.3,    E10.3,    E10.3,    E10.3,    4X,I1,  5X,E10.3,       3x,I2
         * E10.3,  E10.3,    E10.3,   E10.3,   E10.3,    E10.3,    E10.3,    E10.3,    4X,I1,  5X,E10.3,       3x,I2  */
        if (fprintf(stream, "%10.3e%10.3e%10i%10.3e%10.3e%10.3e%10.3e%10.3e    "
                    "%1i     %10s   %2s\n", V1, V2, config->sample, 0., config->alfal0, config->avmass,
                    config->dptmin, config->dptfac, 0, "", "") < 0) {
            success = CPL_FALSE;
        }


        /*** LBLRTM record 1.4 ***/

        /* Format string in the FORTRAN code
         * E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3    4X,1A
         * E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3    4X,1A */
        if (fprintf(stream, "%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e    %c\n",
                    config->tbound,
                    config->sremis[0], config->sremis[1], config->sremis[2],
                    config->srrefl[0], config->srrefl[1], config->srrefl[2],
                    's') < 0) {
            success = CPL_FALSE;
        }


        /*** LBLRTM record 3.1 ***/

        /* Get string with fit molecules and find the last molecule */
        const char *molec_string = lbl_molecs;
        int last_mol = strrchr(molec_string, '1') - lbl_molecs;

        /* Set Atmospheric MODEL */
        int MODEL = 0;

        /* Format string in the FORTRAN code
         * I5,     I5,    I5,      I5,      I5,    I5,     I5,     I2,   1X, I2, F10.3,  F10.3, F10.3,   10x, F10.3
         * I5,     I5,    I5,      I5,      I5,    I5,     I5,     I2,   1X, I2, F10.3,  F10.3, F10.3,   10x, F10.3 */
        if (fprintf(stream, "%5i%5i%5i%5i%5i%5i%5i%2i %2i%10.3e%10.3e%10.3e"
                    "          %10.3e\n", MODEL, config->itype, 0, config->nozero, config->noprnt,
                    last_mol + 1, config->ipunch, 0, 0, config->re, config->hspace, vbar, config->ref_lat) < 0) {
            success = CPL_FALSE;
        }


        /*** LBLRTM record 3.2 ***/

        /* Format string in the FORTRAN code
         * F10.3, F10.3,   F10.3,   F10.3,  F10.3,    I5, 5X,F10.3
         * F10.3, F10.3,   F10.3,   F10.3,  F10.3,    I5, 5X,F10.3 */
        if (fprintf(stream, "%10.3e%10.3e%10.3e%10.3e%10.3e%5i     %10.3e\n",
                    config->h[0], config->h[1], angle, config->range, config->beta, config->len, config->hobs) < 0) {
            success = CPL_FALSE;
        }


        /*** LBLRTM record 3.3a ***/

        /* Format string in the FORTRAN code
         * F10.3,  F10.3,  F10.3, F10.3, F10.3
         * F10.3,  F10.3,  F10.3, F10.3, F10.3 */
        if (fprintf(stream, "%10.3e%10.3e%10.3e%10.3e%10.3e\n",
                    config->avtrat, config->tdiff[0], config->tdiff[1], config->altd[0], config->altd[1])   < 0) {
            success = CPL_FALSE;
        }


        /*** LBLRTM record 3.4 ***/

        /* Store number of height levels from length of merged profile and number of LBLRTM molecules */
        int immax = cpl_table_get_nrow(prof);

        /* Format string in the FORTRAN code
         * I5,    3A8
         * I5,    3A8 */
        if (fprintf(stream, "%5i%8s%8s%8s\n", immax, "", "", "") < 0) {
            success = CPL_FALSE;
        }


        /*** Get names of molecules in atmospheric profile ***/
        cpl_array *allmolecs  = mf_molecules_create_array();
        cpl_array *atm_molecs = cpl_table_get_column_names(prof);

        /* Loop over all height levels */
        for (int level = 0; level < immax; level++) {

            /*** LBLRTM record 3.5 ***/

            /* Format string in the FORTRAN code
             * E10.3, E10.3, E10.3,   5x,  A1,     A1,  1x, A1,     1x,    39A1
             * E10.3, E10.3, E10.3,   5x,  A1,     A1,  1x, A1,     1x,    39A1 */
            if (fprintf(stream, "%10.3e%10.3e%10.3e     AA   ",
                        cpl_table_get_double(prof, MF_COL_ATM_HGT, level, NULL),
                        cpl_table_get_double(prof, MF_COL_ATM_PRE, level, NULL),
                        cpl_table_get_double(prof, MF_COL_ATM_TEM, level, NULL)) < 0) {
                success = CPL_FALSE;
            }

            /* Loop over required molecules in molec_string / allmolecs */
            for (int allmol = 0; allmol <= last_mol; allmol++) {
                fputc('A', stream);
            }

            if (fputc('\n', stream) == EOF) {
                success = CPL_FALSE;
            }


            /*** LBLRTM record 3.6.1 ... 3.6.n ***/

            /* Format string in the FORTRAN code
             * 8E10.3
             * 8E10.3 */

            /* Counter for molecules in list */
            int mol = 0;

            /* Loop over all LBLRTM molecules */
            for (int current_mol = 0; current_mol <= last_mol; current_mol++) {

                /* String for a single molecule */
                const char *molec = cpl_array_get_string(allmolecs, current_mol);

                /* Either write out value from table or 0. otherwise */
                if (   cpl_table_has_column(prof, molec)
                    && strncmp(molec_string + current_mol, "1", 1) == 0 ){

                    if (fprintf(stream, "%10.3e", cpl_table_get_double(prof, molec, level, NULL)) < 0) {
                        success = CPL_FALSE;
                    }
                    mol++;

                } else {
                    if (fprintf(stream, "%10s", "0.000e+00") < 0) success = CPL_FALSE;
                }

                if (mol > last_mol) {
                    break;
                }

                /* 8 molecules per line */
                if ((current_mol + 1) % 8 == 0) {
                    if (fputc('\n', stream) == EOF) success = CPL_FALSE;
                }
            }

            if ((last_mol + 1) % 8 != 0) {
                if (fputc('\n', stream) == EOF) success = CPL_FALSE;
            }
        }

        cpl_array_delete(atm_molecs);
        cpl_array_delete(allmolecs);


        if (fputs("-1\n", stream) == EOF) {
            success = CPL_FALSE;
        }

        if (fprintf(stream, "$ Transfer to ASCII plotting data\n") < 0 ||
            fprintf(stream, " HI=0 F4=0 CN=0 AE=0 EM=0 SC=0 FI=0 PL=1 TS=0 AM=0 "
                    "MG=0 LA=0 MS=0 XS=0    0    0\n") < 0 ||
            fprintf(stream, "# Plot title not used\n") < 0) {
            success = CPL_FALSE;
        }

        if (emission_spec) {

            /* Prepare input by output TAPE27 (Emission) */
            if (fprintf(stream, "%10.4e%10.4e%10.4e%10.4e%5i%5i%5i%5i%10.3e%2i"
                        "%3i%5i\n", V1, V2, 10.2, config->delv, 1, 0, 12, 0, 1., 0,
                        0, 0) < 0) {
                success = CPL_FALSE;
            }

            if (fprintf(stream, "%10.4g%10.4g%10.3e%10.3e%5i%5i%5i%5i%5i%5i%2i"
                        "   %2i%3i\n", 0., 1.2, 7.02, 0.2, 4, 0, 1, 1, 0, 0,
                        1, 3, 27) < 0) {
                success = CPL_FALSE;
            }

        } else {

            /* Prepare input by output TAPE28 (Transmission) */
            if (fprintf(stream, "%10.4e%10.4e%10.4e%10.4e%5i%5i%5i%5i%10.3e%2i"
                        "%3i%5i\n", V1, V2, 10.2, config->delv, 1, 0, 12, 0, 1., 0,
                        0, 0) < 0) {
                success = CPL_FALSE;
            }

            if (fprintf(stream, "%10.4g%10.4g%10.3e%10.3e%5i%5i%5i%5i%5i%5i%2i"
                        "   %2i%3i\n", 0., 1.2, 7.02, 0.2, 4, 0, 1, 0, 0, 0,
                        1, 3, 28) < 0) {
                success = CPL_FALSE;
            }
        }

        if (fputs("-1\n", stream) == EOF) {
            success = CPL_FALSE;
        }

        if (fputs("% created by lblrtm_start1\n", stream) == EOF) {
            success = CPL_FALSE;
        }
    }
    fclose(stream);

    return (success) ? CPL_ERROR_NONE : cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                                              "Problem writing %s file",
                                                              MF_AER_TAPE5_FILE);
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Finds the wavenumber ranges of a set of LBLRTM output files.
 *
 * @param w_dir_range          directory that contains the LBLRTM output files.
 * @param lblrtm_out_filename  Standard filename (usually MF_AER_TAPE27 or MF_AER_TAPE28 (The name is supplemented by "_" and the file number)
 *
 * @return cpl_error_code    CPL array of wavenumbers
 *
 */
/* ---------------------------------------------------------------------------*/
cpl_array * mf_io_find_klim(
    const char               *w_dir_range,
    const char               *lblrtm_out_filename)
{
    cpl_array *klim_all = cpl_array_new(0, CPL_TYPE_DOUBLE);

    /* Get Number of files */
    cpl_size nfile = 0;
    do {

        char *filename_range_ASCII = cpl_sprintf("%s/%s_%lld", w_dir_range, lblrtm_out_filename, nfile + 1);
        cpl_error_code err = mf_io_access(filename_range_ASCII);
        cpl_free(filename_range_ASCII);

        if (!err) {

            /* File found */
            nfile++;

        } else if (nfile == 0) {

            /* Error: Not found any file in the disk */
            cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                  "File opening failed: %s",
                                  lblrtm_out_filename);
            cpl_array_delete(klim_all);
            return NULL;

        } else {

            /* Not found more files in the disk */
            cpl_array_set_size(klim_all, nfile + 1);
            break;
        }

    } while (CPL_TRUE);

    /* Allocate memory for wavenumber limits */
    double *kmin = cpl_calloc(nfile, sizeof(double));
    double *kmax = cpl_calloc(nfile, sizeof(double));

    /* Read wavenumber limits */
    cpl_boolean linestruct = CPL_TRUE;
    for (cpl_size ifile = 1; ifile <= nfile; ifile++) {

        char *filename_range_ASCII = cpl_sprintf("%s/%s_%lld", w_dir_range, lblrtm_out_filename, ifile);
        cpl_msg_info(cpl_func, "(mf_io        ) Load ASCII file: %s (mf_lblrtm_find_klim)", filename_range_ASCII);
        FILE *stream = fopen(filename_range_ASCII, "r");
        cpl_free(filename_range_ASCII);

        /* Get first line */
        char line[MF_LEN_MAX];
        if (fgets(line, MF_LEN_MAX - 1, stream)) {}

        if (line[0] != '1') {
            linestruct = CPL_FALSE;
        } else {

            cpl_boolean v1line = CPL_FALSE;
            do {
                if (strstr(line, "V1 =") != NULL) {
                    v1line = CPL_TRUE;

                    for (int j = 0; j < 2; j++) {

                        /* Split read line into strings */
                        static char *save;
                        char *str = strtok_r(line, "\n\t =", &save);
                        for (int i = 1; i < 3; i++) {
                            str = strtok_r(NULL, "\n\t =", &save);
                            if (str == NULL) {
                                linestruct = CPL_FALSE;
                            }
                        }

                        /* Convert wavenumber string into double */
                        if (j == 0) {
                            kmin[ifile - 1] = strtod(str, NULL);
                            if (fgets(line, MF_LEN_MAX - 1, stream)) {}
                        } else {
                            kmax[ifile - 1] = strtod(str, NULL);
                        }
                    }
                }

                /* Next lines */
                if (fgets(line, MF_LEN_MAX - 1, stream)) {}

            } while (v1line == CPL_FALSE);

            if (kmin[ifile - 1] == 0 || kmax[ifile - 1] == 0) {
                linestruct = CPL_FALSE;
            }
        }

        fclose(stream);
    }

    /* Create klim_all */
    cpl_array_set(klim_all, 0, kmin[0]);

    /* Compute mean wavenumber limits (overlaps!) and write them into output array "klimall" */
    cpl_array_set(klim_all, 0, kmin[0]);
    for (cpl_size ifile = 1; ifile < nfile; ifile++) {
        cpl_array_set(klim_all, ifile, (kmax[ifile - 1] + kmin[ifile]) / 2);
    }
    cpl_array_set(klim_all, nfile, kmax[nfile - 1]);

    /* Free memory occupied by "kmin" and "kmax" */
    cpl_free(kmin);
    cpl_free(kmax);

    /* Handle file structure errors */
    if (!linestruct) {
        cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                              "Unexpected file structure: %s (wavenumber limits missing or = 0)",
                              lblrtm_out_filename);
        cpl_array_delete(klim_all);
        return NULL;
    }

    /* Test correct order of wavenumbers */
    for (cpl_size ifile = 1; ifile <= nfile; ifile++) {
        if (cpl_array_get(klim_all, ifile - 1, NULL) >= cpl_array_get(klim_all, ifile, NULL)) {
            cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                  "Invalid order of data points: cpl_array *klimall (wavenumber limits of LBLRTM output files");
            cpl_array_delete(klim_all);
            return NULL;
        }
    }

    return klim_all;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Rebins LBLRTM spectra (wrapper output) in wavelength units [mu m] (variable step size possible).
 *
 * @param nrow               .
 * @param lamv               .
 * @param fluxv              .
 * @param spectrum_filename  .
 * @param llim               .
 * @param usampl             .
 * @param jmin               .
 * @param jmax               .
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
    mf_io_lnfl_config        *lnfl_config)
{

    FILE *stream;
    if (!lnfl_config->use_ODA) {
        /* Check file existence */
        cpl_msg_info(cpl_func, "(mf_io        ) Load ASCII file: %s (mf_lblrtm_rebin_spectrum)", spectrum_filename);
        stream = fopen(spectrum_filename, "r");
        if (!stream) {
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                        "File opening failed: %s",
                                        spectrum_filename);
        }

        /* Skip header */
        char line[MF_LEN_MAX];
        cpl_boolean endhead = CPL_FALSE;
        if (fgets(line, MF_LEN_MAX - 1, stream)) {}

        if (line[0] != '1') {
            fclose(stream);
            return cpl_error_set_message(cpl_func, CPL_ERROR_BAD_FILE_FORMAT,
                                        "Unexpected file structure: %s (first character != '1')",
                                        spectrum_filename);
        } else {

            do {
                if (strstr(line, "WAVENUMBER") != NULL) endhead = CPL_TRUE;
                if (fgets(line, MF_LEN_MAX - 1, stream)) {}
            } while (endhead == CPL_FALSE);
        }
    }

    /* Search for limiting pixels and wavelengths of input grid that correspond to the given wavenumber interval
    * pixels : jmin, *jmax;  wavelengths : lmin0, lmax0 */
    double lmin0 =  0.;
    double lmax0 =  0.;
    for (cpl_size j = 0; j < nrow; j++) {

        if (*jmin < 0 && lamv[j] > llim[0]) {

            *jmin = j;

            if (*jmin == 0) lmin0 = 1.5 * lamv[*jmin    ] - 0.5 * lamv[*jmin + 1];
            else           lmin0 =      (lamv[*jmin - 1] +       lamv[*jmin    ]) / 2;
        }

        if (*jmin >= 0 && lamv[j] > llim[1]) {

            *jmax = j - 1;

            if (*jmax == nrow - 1) lmax0 = 1.5 * lamv[*jmax] - 0.5 * lamv[*jmax - 1];
            else                  lmax0 = (lamv[*jmax] + lamv[*jmax + 1]) / 2;

            break;

        } else if (j == nrow - 1 && lamv[j] <= llim[1]) {
            *jmax  = j;
            lmax0 = 1.5 * lamv[*jmax] - 0.5 * lamv[*jmax - 1];
        }
    }

    /* Read wavelengths and fluxes of input file. Average all flux values inside a bin of the output wavelength grid. */
    cpl_boolean empty  =  CPL_FALSE;
    double      lam    =  HUGE_VAL;
    double      lmin   =  0.;
    int         num    =  0;
    double      k      =  0.;
    double      flux   =  0.;
    int         ncol   =  2;

    double* xv;
    double* yv;
    int bivector_size;
    int idx;
    if (lnfl_config->use_ODA) {
        xv=cpl_bivector_get_x_data(bvec);
        yv=cpl_bivector_get_y_data(bvec);
        bivector_size=cpl_bivector_get_size(bvec);
        idx=0;
    }

    for (cpl_size j = *jmax; j >= *jmin; j--) {

        /* Adapt lmin for next bin. (lmin: grid related, variable, lower lambda limit) */
        if (j == *jmin) lmin = lmin0;
        else           lmin = (lamv[j - 1] + lamv[j]) / 2;

        /* Read new data point(s) from file? YES: take saved flux value as first summand */
        if (lam > lmin && ncol == 2) {
            empty = CPL_FALSE;
            if (j == *jmax) {
                fluxv[j] = 0.;
                num      = 0;
            } else {
                fluxv[j] = flux;
                num      = 1;
            }
        } else {
            empty = CPL_TRUE;
            num   = 0;
        }

        if (!lnfl_config->use_ODA) {
            while (empty == CPL_FALSE && (ncol = fscanf(stream, "%le %le", &k, &flux)) == 2) {
                lam = MF_CONV_K_LAM / k;
                if (lam <= lmax0) {
                    if (lam > lmin) {
                        /* Sum up fluxes and count data points */
                        fluxv[j] += flux;
                        num++;
                    } else {
                        break;
                    }
                }
            }
        } else {
            while (empty == CPL_FALSE && idx<bivector_size) {
                //double xval=xv[idx];
                //double yval=yv[idx++];
                k   =xv[idx];
                flux=yv[idx++];
                //cpl_msg_info(cpl_func,"i=%d,%f=%f,%f=%f",idx-1,xval,k,yval,flux);
                lam = MF_CONV_K_LAM / k;
                if (lam <= lmax0) {
                    if (lam > lmin) {
                        /* Sum up fluxes and count data points */
                        fluxv[j] += flux;
                        num++;
                    } else {
                        break;
                    }
                }
            }
        }

        if (num != 0) {

            /* Average fluxes */
            fluxv[j] /= (double)num;

            /* Avoid negative values */
            if (fluxv[j] < 0.) fluxv[j] = 0.;

        } else {

            /* Mark "empty" bins */
            fluxv[j] = -HUGE_VAL;

            if (!(*usampl)) *usampl = CPL_TRUE;
        }

    }

    if (! lnfl_config->use_ODA) {
        fclose(stream);
    }

    return CPL_ERROR_NONE;
}

cpl_vector* mf_io_molecule_abundancies(mf_parameters* params, cpl_array* fitpar) {
    int nmolec = params->config->internal.molecules.n_molec;
    const double *par = cpl_array_get_data_double_const(fitpar);
    cpl_vector* v=cpl_vector_new(nmolec);
    for (cpl_size i = 0; i < nmolec; i++) {
         cpl_vector_set(v,i,par[i]);
    }
    return v;
}


// ===================================================================
// ===================================================================
// -------------------------------
// BEGINNING OF MNB-ODA-INSERTION
// -------------------------------

/*
 * NOTE: MNB Hack For ODA Tables
 * Use a single function to set or get specifc data from ODA Table database which
 * are defined as static variables within this routine.
 * Individual, set get and delete wrapper routines are defined around this

*/

cpl_matrix* mf_io_oda_tableDB(int range, int molecule, double* vec, int nrows, int nmols, int option) {

    static cpl_matrix* TABLE_LIST[MF_PARAMETERS_MAXIMUM_NO_OF_RANGES]={NULL};
    static int         TABLE_NROW[MF_PARAMETERS_MAXIMUM_NO_OF_RANGES]={0};

    const int SET_DIMS  =1;
    const int SET_VECTR =2;
    const int FREE_ALL  =3;

    cpl_matrix* od_matrix=TABLE_LIST[range];

    if (option==SET_DIMS) {

        if (od_matrix!=NULL) cpl_matrix_delete(od_matrix);
        TABLE_LIST[range]=cpl_matrix_new(nrows,nmols);
        TABLE_NROW[range]=nrows;
        return TABLE_LIST[range];


    } else if (option==SET_VECTR) {

        if (od_matrix==NULL) return NULL;

        for (int i=0;i<nrows;i++) {
            if (i>TABLE_NROW[range]) break;
            cpl_matrix_set(od_matrix,i,molecule,vec[i]);
        }
        return od_matrix;

    } else if (option==FREE_ALL) {
        cpl_msg_info(cpl_func,"FREE ALL");
        for (int i=0;i<2;i++) {
            cpl_msg_info(cpl_func,"Test OD table for range %d",i);
            cpl_msg_info(cpl_func,"Test OD TABEL_NROW[%d]=%d ",i,TABLE_NROW[i]);
            if (TABLE_NROW[i]==0) {
                cpl_msg_info(cpl_func,"OD table for range %d NOT listed as allocated so leaving alone",i);
            } else {
                cpl_msg_info(cpl_func,"OD table for range %d listed as allocated, Freeing now",i);
                cpl_matrix_delete(TABLE_LIST[i]);
                TABLE_LIST[i]=NULL;
                TABLE_NROW[i]=0;
            };// endif
        };//end for
        cpl_msg_info(cpl_func,"OPTION FREE ALL LOOP FINISHED");
    } else {
        return TABLE_LIST[range];
    }; // end if FREE_ALL

    cpl_msg_info(cpl_func,"MNB ABOUT TO RETURN NULL");
    return NULL;
}
// -------------------------------------------------------------------

void mf_io_oda_init_tableDB(int range, int nmols, int nrows) {

    const int SET_DIMS  =1;
    //const int SET_VECTR =2;

    //cpl_matrix* ret;
    mf_io_oda_tableDB(range, 0, NULL, nrows, nmols, SET_DIMS);
}
// -------------------------------------------------------------------

void mf_io_oda_set_tableDB(int range, int molecule, double *vec, int nrows) {

    //const int SET_DIMS  =1;
    const int SET_VECTR =2;

    //cpl_matrix* ret;
    mf_io_oda_tableDB(range, molecule, vec, nrows, 0, SET_VECTR);
}
// -------------------------------------------------------------------

cpl_matrix*  mf_io_oda_get_tableDB(int range) {

    cpl_matrix* ret;
    ret=mf_io_oda_tableDB(range, 0, NULL, 0, 0, 0);
    return ret;
}
// -------------------------------------------------------------------


cpl_array* mf_io_klim_from_odatable(int range) {

    cpl_matrix* mat=mf_io_oda_get_tableDB(range);
    if (mat==NULL) return NULL;

    /* The wave idx is the 0th column and we want the min and max value*/
    int wave_idx=0;
    int m=cpl_matrix_get_nrow(mat);
    double v1=cpl_matrix_get(mat,0,  wave_idx);
    double v2=cpl_matrix_get(mat,m-1,wave_idx); /* Remember cpl_matrix indecies start at 0*/

    /* Create and return a 2 element cpl_array with v1 and v2 value*/
    cpl_array* klim=cpl_array_new(2,CPL_TYPE_DOUBLE);
    cpl_array_set(klim,0,v1);
    cpl_array_set(klim,1,v2);
    return klim;
}
// -------------------------------------------------------------------

void mf_io_oda_delete_tableDB(void) {
    const int FREE_ALL  =3;
    cpl_msg_info(cpl_func,"About to free OD Tables");
    mf_io_oda_tableDB(0, 0, NULL, 0, 0, FREE_ALL);
}
// -------------------------------------------------------------------
/* END OF : MNB Hack For ODA Tables*/
// ===================================================================


cpl_array* mf_io_molecstring2Names(char* molec_string) {

    /* =========
     * Purpose:
     * =========
     * From the molec string eg "10000110000010000001...."
     * parse to a list of moecule names to be stored as a
     * string array to be returned
    */

    /* Create a molecule table that lists molecule name strings per idx */
    cpl_array *allMoleculeNames=mf_molecules_create_array();
    int str_size=strlen(molec_string);
    int mol_idxv[str_size];
    int nmols=0;
    for (int i=0;i<str_size;i++) {
        if (molec_string[i]=='0') continue;
        mol_idxv[nmols++]=i;
    }

    cpl_array* moleculeNames=cpl_array_new(nmols,CPL_TYPE_STRING);
    for (int i=0;i<nmols;i++) {
        int idx=mol_idxv[i];
        const char* name=cpl_array_get_string(allMoleculeNames,idx);
        cpl_array_set_string(moleculeNames,i,name);
    }

    /* Cleanup */
    cpl_array_delete(allMoleculeNames);

    return moleculeNames;

}/* end mf_io_molecstring2Names */
// -------------------------------------------------------------------

int oda_mol_paramidx2tableidxV[ODA_MAX_MOL_IDX]; /* Permutation Vector from mol idx in parameter list to the official mol id list*/


cpl_bivector* mf_io_mergeODTables(const int range, cpl_vector* mol_abuns, const char* lblrtm_out_filename) {

    /* ========
     * Purpose:
     * ========
     * Given the list of molecular abundancies (mol_abuns) generate a transmission
     * profile as the exponential of the -ve of the linear combination of the
     * optical depth profiles of the molecules using these abundanciy values and return
     * as a bivector.
     *
     * Hacks as follows:
     * 1) GEN_HYBRID=CPL_TRUE implies create a TAPE28_HYBRID value to compare with any TAPE28
     * file. This is a hard coding debug setting.
    */


    //cpl_boolean GEN_HYBRID=CPL_FALSE;
    cpl_boolean GEN_HYBRID=CPL_FALSE;
    cpl_boolean PERMUTE_MOLECULE_ORDERING=CPL_TRUE;/* MNB This is here for testing purposes*/

    if (mf_io_use_stdlblrtm(CPL_TRUE)) GEN_HYBRID=CPL_TRUE;

    int nmols=cpl_vector_get_size(mol_abuns);
    double* v=cpl_vector_get_data(mol_abuns);

    /* Get the Optical Depth Matrix for this range
     * Note Matrix is defined as follows:
     * No of rows = no of data points
     * Column 0 is the wave number cpl_bivector_wrap_vector
     * Column's 1, to end are the Optical depth values for the molecule associated to tha column
     */
    cpl_matrix* oda_table=mf_io_oda_get_tableDB(range);

    /* In cases where optical depth calculation is not being performed yje oda_table
     * will be NULL so return NULL here
     */
    if (oda_table==NULL) {
        cpl_msg_info(cpl_func,"MNBX-> ODA TABLE is NULL");
        return NULL;
    }

    /* Get mol index permutation vector */
    /* MNB int idx1=oda_mol_paramidx2tableidxV[ODA_MAX_MOL_IDX-1]; */
    /*MNB cpl_msg_info(cpl_func,"idx1=%d",idx1);*/


    /* Allocate the vector componenets of the return bivector*/
    int nrows=cpl_matrix_get_nrow(oda_table);
    cpl_vector* wavnum=cpl_vector_new(nrows);
    cpl_vector* merged=cpl_vector_new(nrows);
    double*  wv=cpl_vector_get_data(wavnum);
    double*  mv=cpl_vector_get_data(merged);

    /* Now combine the optical depths as a linear combination*/
    for (int i=0; i<nrows;i++) {
        wv[i]=cpl_matrix_get(oda_table,i,0);
        mv[i]=0.0;
        for (int j=0; j<nmols;j++) {
            double scale=v[j];
            if (PERMUTE_MOLECULE_ORDERING) {
                int mol_idx=oda_mol_paramidx2tableidxV[j];
                mv[i]=mv[i]+scale*cpl_matrix_get(oda_table,i,mol_idx+1);
            } else {
                mv[i]=mv[i]+scale*cpl_matrix_get(oda_table,i,j+1);
            }
        }
        /* Trnsmission is the exponential of the -ve of the combination of the OD's*/
        mv[i]=exp(-mv[i]);
    }

    /* Define the return bivector asa wrap to these vectors*/
    cpl_bivector * rbv=cpl_bivector_wrap_vectors(wavnum,merged);

    if (GEN_HYBRID) {
        /* Dump a TAPE28_HYBRID file for diagnostic/debug purposes*/
        char *h_filename = cpl_sprintf("%s_HYBD",lblrtm_out_filename);
        FILE *stream = fopen(h_filename, "w");
        for (int j=0; j<nmols;j++) {
            double scale=v[j];
            int i=0;
            fprintf(stream,"#MNBX scale=%g odaTable=%g\n",scale,cpl_matrix_get(oda_table,i,j+1));
        }
        for (int i=0; i<nrows;i++) {
            fprintf(stream,"%f %f\n", wv[i], mv[i]);
        }
        fclose(stream);
        cpl_free(h_filename);
    }

    /* Return the calculated transmission as a bivector*/
    return(rbv);

}/* end mf_io_mergeODTables*/
// -------------------------------------------------------------------

cpl_boolean mf_io_use_odatable(cpl_boolean use_ODA) {

    /* Hack routine to return a flag based on the existance/value of
     * env var ODA_OPTION.
     * If env var does not exist return false
     * If env var exists and has value "STD" return false.
     * Otherwise return true.
     */

    /* Read the ODA_OPTION env var */
    char* oda_option=getenv("ODA_OPTION");

    if (oda_option==NULL) return use_ODA;
    /*if (oda_option==NULL) return CPL_FALSE;*/
    /*if (oda_option==NULL) return CPL_TRUE;*/

    cpl_boolean return_flag=CPL_TRUE;
    if (strcmp(oda_option,"STD"  )==0) return_flag=CPL_TRUE;
    if (strcmp(oda_option,"BOTH" )==0) return_flag=CPL_FALSE;
    if (strcmp(oda_option,"DEBUG")==0) return_flag=CPL_TRUE;

    return return_flag;

}/* end mf_io_use_odatable */
// -------------------------------------------------------------------


cpl_boolean mf_io_use_debug(void) {

    /* Hack routine to return a flag based on the existance/value of
     * env var ODA_OPTION.
     * If env var does not exist return false
     * If env var exists and has value "STD" return false.
     * Otherwise return true.
     */

    /* Read the ODA_OPTION env var */
    char* oda_option=getenv("ODA_OPTION");

    if (oda_option==NULL) return CPL_FALSE;

    cpl_boolean return_flag=CPL_FALSE;
    if (strcmp(oda_option,"DEBUG")==0) return_flag=CPL_TRUE;

    return return_flag;

}/* end mf_io_use_odatable */
// -------------------------------------------------------------------


cpl_boolean mf_io_use_stdlblrtm(cpl_boolean use_ODA) {

    /* Hack routine to return a flag based on the existance/value of
     * env var ODA_OPTION.
     * If env var does not exist return true
     * If env var exists and has value "STD" return false.
     * Otherwise return true.
     */

    /* Read the ODA_OPTION env var */
    char* oda_option=getenv("ODA_OPTION");

    if (oda_option==NULL) return !use_ODA;


    cpl_boolean return_flag=CPL_FALSE;
    if (strcmp(oda_option,"STD"  )==0) return_flag=CPL_FALSE;
    if (strcmp(oda_option,"BOTH" )==0) return_flag=CPL_TRUE;
    if (strcmp(oda_option,"DEBUG")==0) return_flag=CPL_FALSE;

    return return_flag;

}/* end mf_io_use_stdlblrtm */
// -------------------------------------------------------------------


cpl_bivector* mf_io_read_lblrtm_spec(
    const char               *spectrum_filename)
{
    /*
     * Parses LBLRTM output file (TAPE28) and returns the data as a bivector (wavenum,trans)
     * Note:
     * This section of code was extracted from mf_io_read_lblrtm_and_update_spec so that
     * the bivector data could be supplied by other means, eg from the ODa table and so
     * that ODa code could simply extract this data from TAPE28 files.
     */

    /* Check file existence */
    cpl_msg_info(cpl_func, "(mf_io        ) Load ASCII file: %s (mf_lblrtm_rebin_spectrum)", spectrum_filename);
    fpos_t fpos;
    int n_lines=0;
    FILE *stream = fopen(spectrum_filename, "r");
    if (!stream) {
        cpl_msg_info(cpl_func,"UH OH");
        return NULL;
    }

    /* Skip header */
    char line[MF_LEN_MAX];
    cpl_boolean endhead = CPL_FALSE;
    if (fgets(line, MF_LEN_MAX - 1, stream)) {}

    if (line[0] != '1') {
        fclose(stream);
        return NULL;
    } else {

        do {
            if (strstr(line, "WAVENUMBER") != NULL) {
                endhead = CPL_TRUE;
                fgetpos(stream,&fpos);
            }
            if (fgets(line, MF_LEN_MAX - 1, stream)) {}
        } while (endhead == CPL_FALSE);
    }

    // Count the number of lines in this file befoer proceeding
    double dummy_x, dummy_y;
    while ( fscanf(stream, "%le %le", &dummy_x, &dummy_y) == 2) {
        n_lines++;
    }
    fsetpos(stream,&fpos);

    cpl_vector* wv=cpl_vector_new(n_lines);
    cpl_vector* fv=cpl_vector_new(n_lines);


    /* Read wavelengths and fluxes of input file. Average all flux values inside a bin of the output wavelength grid. */
    double      k      =  0.;
    double      flux   =  0.;
    int idx=0;

    while (fscanf(stream, "%le %le", &k, &flux) == 2) {
        cpl_vector_set(wv,idx,k);
        cpl_vector_set(fv,idx,flux);
        idx++;
    }

    fclose(stream);

    cpl_bivector* bvec=cpl_bivector_wrap_vectors(wv,fv);

    return bvec;

}/* end mf_io_read_lblrtm_spec*/
// -------------------------------------------------------------------


cpl_bivector* mf_io_merge_wavefiles(
    const char               *w_dir,
    const char               *lblrtm_out_filename)
{

    /*
     * Finds all the LBLRTM output files associated with a specific range that
     * may have been broken up into wavenumber subsections, eg TAPE28_1,TAPE28_2,
     * ...TAPE28_n where each file has a designated range [nu0,nu1] where some
     * overlap is expected.
     * function mf_io_find_klim determines the number of files and produces an
     * array to deine non overlapping intervals that we want from each file.
     * After calling mf_io_find_klim, this routine extracts the data from
     * each file and into bivectors which it then merges into a single
     * (non overlapping) bivector to be returned.
     */

    /* Get the required wavenumber ranges to extract from each file*/
    cpl_array* klim_all = mf_io_find_klim(w_dir,lblrtm_out_filename);
    cpl_size nfiles=cpl_array_get_size(klim_all)-1;

    cpl_msg_info(cpl_func,"Found %lld %s files in directory %s", nfiles,lblrtm_out_filename,w_dir);
    cpl_msg_info(cpl_func,"KLM range:%f %f", cpl_array_get(klim_all,0,NULL),cpl_array_get(klim_all,nfiles,NULL));
    /*Declare nfile bivectors to contain the data from each file */
    cpl_bivector* bivector_lst[nfiles];

    /*Declare the begin and end indecies for the range that we will use for each bivector */
    int idx0V[nfiles];
    int idx1V[nfiles];

    /* Iterate through each file */
    int npts=0;
    for (cpl_size file_i=1;file_i<=nfiles;file_i++) {

        /*From the klim all array get the [nu0,nu1] range to select from */
        double nu0=cpl_array_get(klim_all,file_i-1,NULL);
        double nu1=cpl_array_get(klim_all,file_i  ,NULL);

        /*Import the data from this file a a bivector (wavenum,trans)*/
        char *filename = cpl_sprintf("%s/%s_%lld",w_dir,lblrtm_out_filename,file_i);
        bivector_lst[file_i-1]=mf_io_read_lblrtm_spec(filename);
        cpl_free(filename);

        /* Get the wavenumber axis */
        cpl_vector* nu_axis=cpl_bivector_get_x(bivector_lst[file_i-1]);

        /* Find the first index of this vector that is within [nu0,nu1] range */
        int n=cpl_vector_get_size(nu_axis);
        cpl_msg_info(cpl_func,"nu range=[%f,%f]",cpl_vector_get(nu_axis,0),cpl_vector_get(nu_axis,n-1));
        idx0V[file_i-1]=0; /*Default value*/
        for (int i=0;i<n;i++) {
            if (cpl_vector_get(nu_axis,i)<nu0) continue;
            idx0V[file_i-1]=i;
            cpl_msg_info(cpl_func,"idx0=%d, val=%f",i,cpl_vector_get(nu_axis,i));
            break;
        } /* end i loop */

        /* Find the last index of this vector that is within [nu0,nu1] range */
        idx1V[file_i-1]=n-1; /*Default value*/
        for (int i=idx0V[file_i-1];i<n;i++) {
            if (cpl_vector_get(nu_axis,i)<nu1) continue;
            idx1V[file_i-1]=i;
            cpl_msg_info(cpl_func,"idx1=%d, val=%f",i,cpl_vector_get(nu_axis,i));
            break;
        } /* end i loop */

        /* Add the number of new points we will need for our merged bivector*/
        npts=npts+idx1V[file_i-1]-idx0V[file_i-1]+1;

    }/* end file_i loop*/

    /* Allocate a return bivector to contain the merged data points*/
    cpl_bivector* bivec_ret=cpl_bivector_new(npts);
    cpl_vector* x_ret = cpl_bivector_get_x(bivec_ret);
    cpl_vector* y_ret = cpl_bivector_get_y(bivec_ret);

    /* Iterate through each file and copy the required sections into the return bivector */
    int ret_idx=0;
    for (cpl_size file_i=1;file_i<=nfiles;file_i++) {
        cpl_vector* nu_axis=cpl_bivector_get_x(bivector_lst[file_i-1]);
        cpl_vector* trans_v=cpl_bivector_get_y(bivector_lst[file_i-1]);
        for (int i=idx0V[file_i-1];i<=idx1V[file_i-1]; i++) {
            double nuval=cpl_vector_get(nu_axis,i);
            double trval=cpl_vector_get(trans_v,i);
            cpl_vector_set(x_ret,ret_idx,nuval);
            cpl_vector_set(y_ret,ret_idx,trval);
            ret_idx++;
        } /* end i loop */
    }/* end file_i loop */

    /* Cleanup*/
    cpl_array_delete(klim_all);
    for (cpl_size file_i=1;file_i<=nfiles;file_i++) {
        cpl_bivector_delete(bivector_lst[file_i-1]);
    }

    return bivec_ret;

} /*end mf_io_merge_wavefiles*/
// -------------------------------------------------------------------


void mf_io_oda_symlink(char* target, char* destination) {

    /* ========
     * Purpose:
     * ========
     * It is late I am tired and I could not get the correct access
     * status for mf_io_symlink to work in mf_lblrtm.c routines
     * so I added an accesible wrapper here. (MNB)
     */

    mf_io_symlink(target,destination);

}/* end mf_io_oda_symlink*/
// -------------------------------------------------------------------

cpl_error_code mf_oda_mol_idx(char** mol, int nmol) {

    /* Initialise the parameter list to table list permutation as th eidentity*/
    for (int i=0; i<ODA_MAX_MOL_IDX; i++) oda_mol_paramidx2tableidxV[i]=i;

    cpl_vector* idx_v=cpl_vector_new(nmol);
    cpl_msg_info(cpl_func,"MOELCUL STRING==========");
    for (int i=0;i<nmol;i++) {
        cpl_msg_info(cpl_func,"STRING=%s",mol[i]);
    }


    cpl_array *allMoleculeNames=mf_molecules_create_array();
    cpl_size n_all_mols=cpl_array_get_size(allMoleculeNames);
    int match_cnt=0;
    for (cpl_size i=0;i<n_all_mols;i++) {
        const char* name=cpl_array_get_string(allMoleculeNames,i);
        cpl_msg_info(cpl_func,"STRING=%s -> %lld",name,i);
        for (int j=0;j<nmol;j++) {
            if (strcmp(name,mol[j])==0) {
                cpl_vector_set(idx_v,match_cnt,j);
                oda_mol_paramidx2tableidxV[match_cnt]=j;
                match_cnt++;
                cpl_msg_info(cpl_func,"Match %s with %s: Map %s to %d", name,mol[j],name,j);
            }
        }
    }

    cpl_array_delete(allMoleculeNames);

    return CPL_ERROR_NONE;
} /*  mf_oda_mol_idx*/
// -------------------------------------------------------------------

// ------------------------
// END OF MNB-ODA-INSERTION
// ------------------------
// ===================================================================
// ===================================================================

/** @cond PRIVATE */

/* ---------------------------------------------------------------------------*/
/**
 * @brief Remove a file from the disk
 *
 * @param file
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code mf_io_symlink(
    const char               *target,
    const char               *linkpath)
{
    /* Check inputs */
     
    if (!target || !linkpath) return CPL_ERROR_NULL_INPUT;

    /* Check if target file exist */
    if (mf_io_access(target) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func,"Cannot find target file %s",target);
	return CPL_ERROR_FILE_NOT_FOUND;
    }

    return symlink(target, linkpath);
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Callback function for nftw(...)
 *
 * @param fpath
 * @param sb
 * @param typeflag
 * @param ftwbuf
 *
 * Remove the file path.
 *
 */
/* ---------------------------------------------------------------------------*/
static int mf_io_nftw_rm_rf(
    const char               *path_name,
    const struct stat        *stat_buf,
    int                      file_type,
    struct FTW               *ftw_buf)
{
    cpl_msg_debug(cpl_func, "(mf_io        ) Remove file = %s ; stat->mode = %i, type = %i, ftw_buf->level = %i",
                  path_name, stat_buf->st_mode, file_type, ftw_buf->level);

    int err = mf_io_rm(path_name);

    if (err) perror(path_name);

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
static int mf_io_get_socket_connection(
    const char               *host,
    const char               *port)
{
    int sockfd;

    /* IP Name resolution */

    /* Set the hints. First we initialize to 0 the structure.
       Only retrieve IPv4 or IPv6 if configured in the system */

    //CPL_DIAG_PRAGMA_PUSH_IGN(-Wmissing-field-initializers);
    //struct addrinfo hints = { 0 };
    //CPL_DIAG_PRAGMA_POP;

    struct addrinfo hints;
    memset(&hints, 0, sizeof(hints));
    hints.ai_flags = AI_ADDRCONFIG;
    hints.ai_socktype = SOCK_STREAM;
    /* Getting the list of IP addresses */
    cpl_msg_debug(cpl_func, "(mf_io        ) Getting IP");
    struct addrinfo * addr_list ;
    if (getaddrinfo(host, port, &hints, &addr_list) != 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Couldn't get address for host");
        return 0;
    }

    /* Connecting to the server for the FTP commands.
       The first address to which we can connect will be used.
       addr_list is a linked list */
    cpl_msg_debug(cpl_func, "(mf_io        ) Connecting to server");
    struct addrinfo *this_addr;
    for(this_addr = addr_list; this_addr != NULL; this_addr = this_addr->ai_next)
    {
        /* Opening the socket */
        if ((sockfd = socket(this_addr->ai_family, this_addr->ai_socktype,
                this_addr->ai_protocol)) == -1) {
            continue;
        }

        if (connect(sockfd, this_addr->ai_addr, this_addr->ai_addrlen) == -1) {
            close(sockfd);
            if(errno == ECONNREFUSED)
                errno = 0; //Reset errno if no remote partner at this address
            continue;
        }
        cpl_msg_debug(cpl_func, "(mf_io        ) Connection established");
        break;
    }

    if (this_addr == NULL)
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Couldn't connect to the host");
        return 0;
    }

    freeaddrinfo(addr_list);

    return sockfd;
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
static int mf_io_ftp_reply(
    int                      sockfd,
    char                     **message)
{
    int  length, n;
    char buffer[MF_IO_BUFFER_SOCKET_LENGTH];
    char * msg = NULL;

    length = 0;
    if(message != NULL)
        *message = NULL;

    while( ( n = recv(sockfd, buffer, MF_IO_BUFFER_SOCKET_LENGTH - 1, 0) ) > 0)
    {
        void * tmp = realloc(msg, length + n + 1);
        if(tmp == NULL)
        {
            free(msg);
            return 0;
        }
        else
        {
            msg =tmp;
        }

        strncpy(msg + length, buffer, n);
        length += n;

        //Check FTP end of message.
        //http://www.tcpipguide.com/free/t_FTPRepliesReplyCodeFormatandImportantReplyCodes-5.htm
        if(msg[length-1]=='\n') //Reply contains full lines
        {
            //Search for the beggining of the last line
            char * eolchar= msg + length - 1;
            while(--eolchar != msg)
                if (*(eolchar - 1) == '\n')
                    break;

            //The FTP code has 3 numbers and a space afterwards if it is the last line
            if(*(eolchar+3) == ' ')
                break;
        }
    }
    if(length == 0)
    {
        free(msg);
        return 0;
    }
    msg[length] = '\0';

    cpl_msg_debug(cpl_func,"(mf_io        ) FTP reply: <<%s>>", msg);

    /* verify */
    int verify = mf_io_verify_ftp_code(msg, length + 1);

    if(message != NULL && verify)
        *message = msg;
    else
        free(msg);
    return verify;
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
static int mf_io_verify_ftp_code(
    char                     *msg,
    int                      length)
{
    char * line = msg;

    //FTP protocol specifies that lines starting with 2xx codes are ok
    //Starting with 3xx are ok but the server expects some extra input
    //Starting with 4xx, 5xx or 6xx it denotes an error.
    //http://www.tcpipguide.com/free/t_FTPRepliesReplyCodeFormatandImportantReplyCodes-2.htm
    while(line[0] == '1' || line[0] == '2' || line[0] == '3')
    {
        line = strchr(line, '\n');
        if(line == NULL || line - msg + 2 == length)
            return 1;
        line = line + 1;
    }

    return 0;
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
static int mf_io_send_ftpcmd(
    int                      sockfd,
    const char                *cmd)
{
    cpl_msg_debug(cpl_func, "(mf_io        ) Sending FTP command <<%s>>", cmd);

    if(write(sockfd, cmd, strlen(cmd))==0)
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Problem during FTP transaction");
        return 0;
    }

    char * msg;
    if(!mf_io_ftp_reply(sockfd, &msg))
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Problem during FTP transaction");
        return 0;
    }
    if (msg != NULL)
        free(msg);

    return 1;
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
static int mf_io_send_pasv(
    int                      sockfd,
    const char               *cmd)
{
    if(write(sockfd, cmd, strlen(cmd))==0)
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Problem during FTP transaction");
        return 0;
    }

    char * msg = NULL;
    char * new_con;
    unsigned int v[6];

    if(!mf_io_ftp_reply(sockfd, &msg))
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
            "Problem during FTP transaction");
        return 0;
    }
    if (msg == NULL)
      return 0;

    new_con = strchr(msg,'(');
    if (new_con == NULL)
        return 0;
    sscanf(new_con+1,"%10u,%10u,%10u,%10u,%10u,%10u",&v[2],&v[3],&v[4],&v[5],&v[0],&v[1]);

    //Get the new port connection of passive mode
    //This is coded in the reply of PASV
    //See http://www.freefire.org/articles/ftpexample.php
    int data_port = v[0] * 256 + v[1];

    free(msg);

    return data_port;
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
static char * mf_io_get_ftp_file(
    int                      sockfd,
    int                      *data_length)
{
    int  length, n;
    char buffer[MF_IO_BUFFER_SOCKET_LENGTH];
    char * msg = NULL;

    length = 0;

    /* Get the data */
    cpl_msg_debug(cpl_func, "(mf_io        ) Get the data");
    while( ( n = recv(sockfd, buffer, MF_IO_BUFFER_SOCKET_LENGTH - 1, 0) ) > 0)
    {
        void * tmp = realloc(msg, length + n + 1);
        if(tmp == NULL)
        {
            free(msg);
            return 0;
        }
        else
        {
            msg =tmp;
        }

        memcpy(msg + length, buffer, n);
        length += n;

        cpl_msg_debug(cpl_func, "(mf_io        ) Received %d bytes so far",length);
    }
    if(length == 0)
        return 0;
    msg[length] = '\0';
    *data_length = length+1;

    return msg;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Trim an input string in-place.
 *
 * @param str                .
 *
 * @return cpl_error_code    CPL_ERROR_NONE is everything is OK.
 *                           If not, these are the errors:
 *                           - .
 *                           - Error in subroutine (see subroutines).
 *
 * @note This function removes all leading and trailing " " characters from str using isspace().
 *       The action is performed in place, the input gets overwritten.
 *
 */
/* ---------------------------------------------------------------------------*/
static void mf_io_str_trim(
    char                     *str)
{
    if (str != NULL && *str != 0) {

        /* length of input string */
        int len = strlen(str);

        /* Set pointers to beginning of the input string: for looping through */
        char *p = str;

        /* Skip over leading spaces */
        char *q = str;
        int i = 0;                /* Counter for start position of string */
        while (isspace(*q)) {
            i++;
            q++;
        }

        /* Skip over trailing spaces */
        q = str + len - 1;
        int j = 0;                /* Counter for end position of string */
        while (isspace(*q)) {
            j++;
            if (*q == *p) break;  /* Avoids valgrind error */
            q--;
        }

        /* j now has end position of string */

        j = len - j - i;          /* Count number of remaining characters */
        q = str + i;              /* Set pointer to beginning of trimmed string */
        for (; j > 0; j--) {
            *p = *q;              /* Copy characters */
            p++;
            q++;
        }

        /* Terminate string */
        *p = '\0';
    }
}

/** @endcond */


/**@}*/
