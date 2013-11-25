/**
 * @file pngexporter.h
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

#ifndef __FDL_PNG_EXPORTER_H
#define __FDL_PNG_EXPORTER_H

#include <string>

#include "core/grid.hpp"
#include "logger/logger.h"
#include "io/exporterbase.h"

namespace fdl {
	
class PngExporter : public ExporterBase {
public:
	PngExporter(std::string prefix="density_export_");
	~PngExporter();
	// virtual void write(const fdl::Grid& grid);
	
private:	
	virtual void write();
	
	void exportDensity(int counter, std::string prefix, float* field, int xRes, int yRes, int zRes);
	int writePNG(const char* filename, unsigned char** rowsp, int w, int h);

};

}	// namespace fdl

#endif	// __FDL_PNG_EXPORTER_H