/**
 * @file common.h
 * @author Caleb Johnston
 * @version 0.1
 *
 * @section LICENSE
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef __FDL_COMMON_H
#define __FDL_COMMON_H
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>

/**
 * fdl is used for all members of the fluid dynamics library.
 */
namespace fdl {

typedef boost::numeric::ublas::vector<float> Vector;

typedef enum medium_T { FLUID, SOLID, SMOKE, AIR } medium_T;
typedef enum interp_T { LINEAR, RK2, CATMULLROM } interp_T;
// typedef enum fileFormat_T { POV_RAY, BLENDER, YAFARAY, PPM, PBRT, PNG } fileFormat_T;

const int DIMENSIONS = 3;
const float EPSILON = 1e-5; // use FLT_EPSILON in <cfloat> instead? 
const double PI = M_PI;
const double E = M_E;

}

#endif	// __FDL_COMMON_H