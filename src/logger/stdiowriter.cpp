/*
 *  StdOutWriter.cpp
 *  ECSLogger
 *
 *  Created by Jerry Belich on 8/11/10.
 *  Copyright 2010 Clockwork. All rights reserved.
 *
 */

#include "logger/stdiowriter.h"

#include <iostream>

#include <boost/thread/mutex.hpp>

namespace fdl {
	
boost::mutex StdOutWriter::m_coutMutex;

/**
 * Constructor. Initializes the format property.
 *
 */
StdOutWriter::StdOutWriter() : LogWriter( "Standard Output" )
{
	// set default format
	m_format = "%d-%b-%Y %H:%M:%S";
}

/**
 * Destructor.
 *
 */
StdOutWriter::~StdOutWriter()
{
}


/**
 * Accessor method for the logging format property.
 *
 * @param newFormat the new value for the logging format
 *
 */
void StdOutWriter::setFormat(const char* newFormat)
{
	m_format = newFormat;
}

/**
 * Writes a message to stdout given a logging level and logger identity.
 *
 * @param level the logging level
 * @param identity string value to be attached to the logged output
 * @param message the string value to be logged
 *
 */
void StdOutWriter::write( Logger::LEVEL level, const std::string& identity, const std::string& message )
{
	boost::mutex::scoped_lock lock(m_coutMutex);
	std::cout << "[" << Logger::currentTime(m_format) << " " << Logger::loggerLevelAsString( level ) << "] " <<  message << std::endl;
}

}	// namespace fdl