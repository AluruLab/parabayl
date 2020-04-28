/***
 *  $Id: read_exp.hpp 255 2009-06-15 14:38:15Z zola $
 **
 *  File: read_exp.hpp
 *  Created: Dec 16, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 */

#ifndef READ_EXP_HPP
#define READ_EXP_HPP

#include <jaz/string_add.hpp>
#include <jaz/sys_tools.hpp>
#include <mpix/MPI_env.hpp>

#include <fstream>
#include <set>
#include <sstream>

#include "AppConfig.hpp"
#include "InputData.hpp"


bool read_exp(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
	      InputData& data);

#endif // READ_EXP_HPP
