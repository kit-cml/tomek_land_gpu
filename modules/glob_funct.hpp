#ifndef GLOB_FUNCT_HPP
#define GLOB_FUNCT_HPP

#include <cstdio>

#include "glob_type.hpp"
#include "param.hpp"
#include "../cellmodels/cellmodel.hpp"

/**
 * @brief Print formatted output to standard output.
 *
 * @param node_id The ID of the node.
 * @param fmt The format string.
 * @param ... Additional arguments for the format string.
 */
void mpi_printf(unsigned short node_id, const char *fmt, ...);

/**
 * @brief Print formatted output to the specified file stream.
 *
 * @param node_id The ID of the node.
 * @param stream The file stream to print to.
 * @param fmt The format string.
 * @param ... Additional arguments for the format string.
 */
void mpi_fprintf(unsigned short node_id, FILE *stream, const char *fmt, ...);

/**
 * @brief Assign parameters from command line arguments and input deck file.
 *
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @param p_param Pointer to the parameter structure to be filled.
 */
void edison_assign_params(int argc, char *argv[], param_t *p_param);

/**
 * @brief Create a directory with the specified name.
 *
 * @param dirname The name of the directory to create.
 * @return 0 on success, -1 on failure.
 */
int make_directory(const char* dirname);

/**
 * @brief Check if a file exists at the specified path.
 *
 * @param pathname The path of the file to check.
 * @return 0 if the file exists, -1 otherwise.
 */
int is_file_existed(const char* pathname);

#endif
