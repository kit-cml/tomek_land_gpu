/**
 * @file glob_funct.cpp
 * @brief Implementation of global functions used in the simulation framework.
 *
 * This file contains the implementation of various utility functions used throughout the simulation framework,
 * including functions for printing, reading parameters, and handling files and directories.
 */

#include "glob_funct.hpp"

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// To make it more "portable" between OSes.
#if defined _WIN32
  #include <direct.h>
  #define snprintf _snprintf
  #define vsnprintf _vsnprintf
  #define strcasecmp _stricmp
  #define strncasecmp _strnicmp
#else
  #include <dirent.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/**
 * @brief Print formatted output to standard output.
 *
 * @param node_id The ID of the node.
 * @param fmt The format string.
 * @param ... Additional arguments for the format string.
 */
void mpi_printf(unsigned short node_id, const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
}

/**
 * @brief Print formatted output to the specified file stream.
 *
 * @param node_id The ID of the node.
 * @param stream The file stream to print to.
 * @param fmt The format string.
 * @param ... Additional arguments for the format string.
 */
void mpi_fprintf(unsigned short node_id, FILE *stream, const char *fmt, ...)
{
#ifndef _WIN32
  if(mympi::rank == node_id){
    va_list args;
    va_start(args, fmt);
    vfprintf(stream, fmt, args);
    va_end(args);
  }
#else
  va_list args;
  va_start(args, fmt);
  vfprintf(stream, fmt, args);
  va_end(args);
#endif
}

/**
 * @brief Assign parameters from command line arguments and input deck file.
 *
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @param p_param Pointer to the parameter structure to be filled.
 */
void edison_assign_params(int argc, char *argv[], param_t *p_param)
{
  bool is_default = false;
  char buffer[100];
  char key[100];
  char value[100];
  char file_name[150];
  FILE *fp_inputdeck;

  // Parameters from arguments
  for (int idx = 1; idx < argc; idx += 2) {
    if (!strcasecmp(argv[idx], "-input_deck"))
      strncpy(file_name, argv[idx + 1], sizeof(file_name));
    else if (!strcasecmp(argv[idx], "-hill_file"))
      strncpy(p_param->hill_file, argv[idx + 1], sizeof(p_param->hill_file));
    // else if (!strcasecmp(argv[idx], "-cvar_file"))
    //   strncpy(p_param->hill_file, argv[idx + 2], sizeof(p_param->cvar_file));
  }

  fp_inputdeck = fopen(file_name, "r");
  if(fp_inputdeck == NULL){
    fprintf(stderr, "Cannot open input deck file %s!!!\nUse default value as the failsafe.\n", file_name);
    is_default = true;
  }

  // Read input deck line by line and store each line in the buffer
  while (!is_default && fgets(buffer, 100, fp_inputdeck) != NULL) {
    sscanf(buffer, "%s %*s %s", key, value);
    if (strcasecmp(key, "Simulation_Mode") == 0) {
      p_param->simulation_mode = strtod(value, NULL);
    }
    else if (strcasecmp(key, "Celltype") == 0) {
      p_param->celltype = strtod(value, NULL);
    }
    else if (strcasecmp(key, "Is_Dutta") == 0) {
      p_param->is_dutta = strtol(value, NULL, 10);
    }
    else if (strcasecmp(key, "Use_Conductance_Variability") == 0) {
      p_param->is_cvar = strtol(value, NULL, 10);
    }
    else if (strcasecmp(key, "Pace_Find_Steepest") == 0) {
      p_param->find_steepest_start = strtod(value, NULL);
    }
    else if (strcasecmp(key, "GPU_Index") == 0) {
      p_param->gpu_index = strtod(value, NULL);
    }
    else if (strcasecmp(key, "Basic_Cycle_Length") == 0) {
      p_param->bcl = strtod(value, NULL);
    }
    else if (strcasecmp(key, "Number_of_Pacing") == 0) {
      p_param->pace_max = strtod(value, NULL);
    }
    else if (strcasecmp(key, "Time_Step") == 0) {
      p_param->dt = strtod(value, NULL);
    }
  }

  if (!is_default) fclose(fp_inputdeck);
}

/**
 * @brief Create a directory with the specified name.
 *
 * @param dirname The name of the directory to create.
 * @return 0 on success, -1 on failure.
 */
int make_directory(const char* dirname)
{
#if defined _WIN32
  return _mkdir(dirname);
#else
  return mkdir(dirname, 0775);
#endif
}

/**
 * @brief Check if a file exists at the specified path.
 *
 * @param pathname The path of the file to check.
 * @return 0 if the file exists, -1 otherwise.
 */
int is_file_existed(const char* pathname)
{
#if defined _WIN32
  struct _stat buf;
  return _stat(pathname, &buf);
#else
  struct stat st = {0};
  return stat(pathname, &st);
#endif
}
