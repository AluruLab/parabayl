/***
 *  $Id: parabayl.hpp 557 2010-03-20 22:41:27Z zola $
 **
 *  File: parabayl.hpp
 *  Created: May 22, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *          Olga Nikolova <olga.nikolova@gmail.com>
 */

#ifndef PARA_BAY_L_HPP
#define PARA_BAY_L_HPP

#include <jaz/bit_util.hpp>
#include <mpix/MPI_env.hpp>

#include "AppConfig.hpp"
#include "BitCounter.hpp"
#include "InputData.hpp"
#include "MDL_score.hpp"
#include "BDe_score.hpp"
#include "bit_config.hpp"
#include "iomanip.hpp"

#include <iostream>
#include <fstream>
#include <limits>


bool parabayl(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
	      InputData& data);

#endif // PARA_BAY_L_HPP
