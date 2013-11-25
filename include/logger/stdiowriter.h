/**
 * @file stdiowriter.h
 * @author M. Hank Kiedrowski
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

#ifndef __FDL_StdOutWriter_H
#define __FDL_StdOutWriter_H

#include "logger/logwriter.h"

namespace boost {
	class mutex;
}

namespace fdl {
	
class StdOutWriter : public LogWriter {
	
public:
	StdOutWriter();
	~StdOutWriter();

	void write(Logger::LEVEL level, const std::string& identity, const std::string& message);
	void setFormat(const char* newFormat);
	
private:
	static boost::mutex m_coutMutex;
	const char* m_format;
	
};

}	// namespace fdl

#endif	//__FDL_StdOutWriter_H