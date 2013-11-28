/*
 *  SyslogWriter.cpp
 *  XMPPBase
 *
 *  Created by M. Hank Kiedrowski on 8/24/10.
 *  Copyright 2010 Clockwork Active Media Systems. All rights reserved.
 *
 */

#include "logger/syslogwriter.h"

#include <boost/thread/mutex.hpp>
#ifndef WIN32
#include <syslog.h>
#endif // WIN32
#include <stdarg.h>


namespace fdl {
	
boost::mutex SyslogWriter::m_syslogOutMutex;

/**
 * constructor, initializes superclass LogWriter with "Syslog Writer"
 *
 */
SyslogWriter::SyslogWriter( ) :
LogWriter( "Syslog Writer" )
{
}

/**
 * destructor
 *
 */
SyslogWriter::~SyslogWriter( ) {
}

/**
 * Writes a message to the log, given a logging level and logger identity.
 *
 * @param level the logging level
 * @param identity string value to be attached to the logged output
 * @param message the string value to be logged
 *
 */
void SyslogWriter::write( Logger::LEVEL level, const std::string& identity, const std::string& message ) {
	
	boost::mutex::scoped_lock lock( m_syslogOutMutex );
	
	std::ostringstream m;
	m << "[" << Logger::loggerLevelAsString( level ) << "] " << message;
#ifndef WIN32
	openlog( identity.c_str(), LOG_CONS | LOG_PID | LOG_NDELAY, LOG_LOCAL5 );
	setlogmask( LOG_UPTO( LOG_NOTICE ) );
	syslog(mapLevel( level ), "%s", m.str().c_str() );
#endif // WIN32
}

/**
 * Maps an internal logging level to a syslog level and returns it
 *
 * @param level to be mapped to syslog
 * 
 * @return an integer syslog level
 *
 */
int SyslogWriter::mapLevel( const Logger::LEVEL& level ) {
#ifdef WIN32
	return 0;
#else // WIN32
	switch ( level ) {
		case Logger::ERROR:
			return LOG_ERR;
		case Logger::WARN:
			return LOG_WARNING;
		case Logger::INFO:
		case Logger::DEV:
		case Logger::DEBUG:
			return LOG_NOTICE;
	}
	return LOG_DEBUG;
#endif // WIN32
}

}	// namespace fdl