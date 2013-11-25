/**
 * @file logger.h
 * @author M. Hank Kiedrowski, Caleb Johnston
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

#ifndef __FDL_Logger_H
#define __FDL_Logger_H

#define LOG(level) \
fdl::Logger::canLog(level) && fdl::Logger().log(level)

#define DEV() \
fdl::Logger::canLog(fdl::Logger::DEV) && fdl::Logger().log(fdl::Logger::DEV)

#define DEBUG() \
fdl::Logger::canLog(fdl::Logger::DEBUG) && fdl::Logger().log(fdl::Logger::DEBUG)

#define INFO() \
fdl::Logger::canLog(fdl::Logger::INFO) && fdl::Logger().log(fdl::Logger::INFO)

#define ERROR() \
fdl::Logger::canLog(fdl::Logger::ERROR) && fdl::Logger().log(fdl::Logger::ERROR)

#define LOGIdent(level, identity) \
fdl::Logger::canLog(level) && ecs::Logger().log(level, identity)

#include <sstream>
#include <string>
#include <list>
#include <vector>

namespace fdl {
	
class LogWriter;

class Logger  {
	
public:	
	Logger();
	~Logger();
	
	enum LEVEL {
		
		// Error messages; something has gone wrong.
		ERROR = 1,
		
		// Possible error, possible misconfiguration
		WARN = 2,
		
		// Messages that may be of interest; verbose action logging.
		INFO = 4,
		
		// Debugging messages
		DEBUG = 8,
		
		// Development messages, should not be used in production
		DEV = 16
		
	};
		
	// if level is not set default will be used; FAIL through INFO
	static const int LEVEL_FLOOR = Logger::ERROR;
	static const int ALL = 255;
	static const int DEFAULT_REPORTING_LEVEL = 31;
	
	static bool canLog(Logger::LEVEL level);
	static const int setLevel(const int level); // sets logging level
	static const int getLevel(); // returns the current level
	static void pushLevel(Logger::LEVEL level); // pushes the passed log level onto the level stack
	static void popLevel(); // pops the top level off of the log level stack
	
	std::ostringstream& log(Logger::LEVEL level); // logs a message using registered LogWriters(s)
	std::ostringstream& log(Logger::LEVEL level, const std::string& identity);
	
	// Identity is used to label the use of the log, typically this should be the application name
	static void setIdentity(const std::string& identity);
	static const std::string getIdentity();
	
	static void registerWriter(LogWriter* writer); // registers the passed LogWriter
	static void unRegisterWriter(LogWriter* writer); // un registers the passed LogWriter
	
	static std::string currentTime(const char* format="%d-%b-%Y %H:%M:%S");
	static std::string loggerLevelAsString(Logger::LEVEL level);
	static Logger::LEVEL stringAsLoggerLevel(std::string level_string);
	
	// returns the passed level OR'd with all lower levels
	static int levelAndBelow(Logger::LEVEL level);
	
protected:
	static int _level;
	static std::vector<int> _levelStack;
	static const int _maxWriters = 10;
	static std::list<LogWriter *> _registeredWriters;
	static std::string _identity;
	std::string* m_local_identity;
	
	std::ostringstream m_message;
	Logger::LEVEL m_messageLevel;
	
	void writeMessage(); // sends the current message to the RegisteredWriters
	
};

}	// namespace fdl

#endif	//__FDL_Logger_H