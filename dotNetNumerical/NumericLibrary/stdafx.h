// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>



// TODO: reference additional headers your program requires here
#include <stdio.h>
#include <tchar.h>
#include<math.h>
#include<float.h>
#include <stdexcept>

// TODO: reference additional headers your program requires here
#include<map>
#include<vector>
#include<iostream>
#include<memory>

#include <boost/preprocessor/stringize.hpp>

#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>
using namespace boost::timer;

#include <Eigen/Dense>
using Eigen::Matrix;
using Eigen::MatrixXd;