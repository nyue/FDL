/**
 * @file exporter.h
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

#ifndef __FDL_EXPORTER_BASE_H
#define __FDL_EXPORTER_BASE_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/thread/thread.hpp>

#include "core/grid.hpp"
#include "core/common.h"

namespace fdl {
	
class ExporterBase {
public:
	ExporterBase(std::string prefix="density_export_");
	~ExporterBase();
	bool isWritingFile() const { return m_writingFile; }
	virtual long int start(fdl::Grid& grid);
	virtual void stop();
	
	
protected:
	// boost::thread m_writeThread;
	boost::mutex m_writeMutex; 
	std::ofstream* m_filestream;
	std::string m_filenamePrefix;
	int m_filenameCounter;
	bool m_writingFile;
	bool m_isCancelled;
	
	Grid* m_grid;
	// fdl::Grid m_grid;
	virtual void write() = 0;
	
};	// class Exporter
};	// namespace fdl

#endif	// __FDL_EXPORTER_BASE_H