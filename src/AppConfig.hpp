/***
 *  $Id: AppConfig.hpp 566 2010-04-01 21:32:04Z olia $
 **
 *  File: AppConfig.hpp
 *  Created: May 22, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *          Olga Nikolova <olga.nikolova@gmail.com>
 */

#ifndef APP_CONFIG_HPP
#define APP_CONFIG_HPP

#include <string>

#include "config.h"
#include "iomanip.hpp"


struct AppConfig {

  AppConfig() : input(""), output(""), deg_bound(64), sync_step(0), mem_report(true) { }
  
  // name of the input file
  std::string input;
  
  // name of the output file
  std::string output;
  
  // in-degree bound
  int deg_bound;
  
  // sync step
  long int sync_step;
  
  // mem report
  bool mem_report;
  
}; // struct AppConfig

#endif // APP_CONFIG_HPP
