/*
 *  Logger.cpp
 *  ECSLogger
 *
 *  Created by M. Hank Kiedrowski on 8/4/10.
 *  Copyright 2010 Clockwork Active Media Systems. All rights reserved.
 *
 */


#include "logger/logger.h"
#include "logger/logwriter.h"

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <iostream>
#include <locale>

namespace fdl {
	
int Logger::_level = Logger::DEFAULT_REPORTING_LEVEL;
std::vector<int> Logger::_levelStack( 1, Logger::DEFAULT_REPORTING_LEVEL );
std::list<LogWriter *> Logger::_registeredWriters;
std::string Logger::_identity = "fdlLogger";

/**
 * Constructor.
 *
 */
Logger::Logger() :
m_local_identity( NULL ) {
}

/**
 * Destructor.
 *
 */
Logger::~Logger() {
	writeMessage();
	if( m_local_identity != NULL ) {
		delete m_local_identity;
	}
}

/**
 * Sets the new identity of the logger
 *
 * @param identity new string value for the logger
 *
 */
void Logger::setIdentity( const std::string& identity ) {
	Logger::_identity = identity;
}

/**
 * Returns the identity value of the logger
 * 
 * @return the string name of the logger
 *
 */
const std::string Logger::getIdentity() {
	return Logger::_identity;
}

/**
 * Sets the level of the logger
 *
 * @param level the new logger level
 * 
 * @return the old level
 *
 */
const int Logger::setLevel( const int level ) {
	
	int old_level = Logger::_level;
	int l;
	
	if ( level >= Logger::LEVEL_FLOOR ) {
		l = level;
	}
	else {
		l = Logger::DEFAULT_REPORTING_LEVEL;
	}
	
	Logger::_level = l;
	
	// clear the level stack and set the current value
	Logger::_levelStack.clear();
	Logger::_levelStack.push_back( Logger::_level );
	return old_level;
}

/**
 * Returns the logger level
 * 
 * @return level
 *
 */
const int Logger::getLevel() {
	return Logger::_level;
}


/**
 * Registers a logwriter, thus adding it to the logging process
 *
 * @param writer to be registered
 *
 */
void Logger::registerWriter( LogWriter * writer ) {
	std::list<LogWriter *>::iterator itr;
	
	for (itr = _registeredWriters.begin(); itr != _registeredWriters.end(); itr++) {
		if (writer->type() == (*itr)->type()) {
			return;
		}
	}
	Logger::_registeredWriters.push_back( writer );
}


/**
 * Unregisters a logwriter, thus removing it from the logging process
 *
 * @param writer to be unregistered
 *
 */
void Logger::unRegisterWriter( LogWriter * writer ) {
	
	Logger::_registeredWriters.remove( writer );
}


/**
 * Adds a level to the logger level stack.
 *
 * @param level to be added
 *
 */
void Logger::pushLevel( Logger::LEVEL level ) {

	Logger::_levelStack.push_back( Logger::_level );
	
	Logger::_level = Logger::_level | level;
	
}


/**
 * Removes the last logger level on the stack.
 *
 */
void Logger::popLevel() {
	
	if( Logger::_levelStack.empty() == false ) {
		Logger::_level = Logger::_levelStack.back();
		Logger::_levelStack.pop_back();
	}
}


/**
 * Determines whether or not the level is currently being logged.
 *
 * @param compare logger level to check
 * 
 * @return a boolean flag
 *
 */
bool Logger::canLog( Logger::LEVEL compare ) {
	bool isWithin = 0;
	
	if (!(_level >= Logger::LEVEL_FLOOR) || !(compare >= Logger::LEVEL_FLOOR)) {
		return false;
	}
	
	isWithin = ( ( _level & compare ) != 0 );
	
	return isWithin;
}


/**
 * Stores the input message to be logged at the next write step.
 *
 * @param level a logger level mask
 *
 */
std::ostringstream& Logger::log( Logger::LEVEL level ) {
	m_messageLevel = level;

	return m_message;
}


/**
 * Assigns the local identity of a given message and stores it to be logged at the next write step.
 *
 * @param level a logger level
 * @param identity for the message
 *
 */
std::ostringstream& Logger::log( Logger::LEVEL level, const std::string& identity ) {
	
	if( m_local_identity == NULL ) {
		m_local_identity = new std::string( identity );
	}
	return log( level );
}


/**
 * Iterates through the logwriter stack and calls write() on each of them.
 *
 */
void Logger::writeMessage() {
	
	for ( std::list<LogWriter *>::iterator it = Logger::_registeredWriters.begin(); it != Logger::_registeredWriters.end(); it++ ) {
		if( m_local_identity != NULL ) {
			( *it )->write( m_messageLevel, *m_local_identity, m_message.str() );
		}
		else {
			( *it )->write( m_messageLevel, getIdentity(), m_message.str() );
		}
	}
	
}


/**
 * Returns the current time using the Boost dateTime library given an input format.
 *
 * @param format the dateTime format for the time
 * 
 * @return a string representation of the time using the input format
 *
 */
std::string Logger::currentTime(const char* format) {
	
	std::ostringstream dateTime;
	const boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	boost::posix_time::time_facet* const f = new boost::posix_time::time_facet(format);
	dateTime.imbue(std::locale(dateTime.getloc(),f));
	dateTime << now;
	return dateTime.str();
}


/**
 * Returns a string representation of the logger level
 *
 * @param level the logging level
 * 
 * @return a string for the level
 *
 */
std::string Logger::loggerLevelAsString( Logger::LEVEL level ) {
	
	std::ostringstream debugLevel;
	debugLevel.flags( std::ios::right );

	switch ( level ) {
		case Logger::ERROR:
			debugLevel << "ERROR";
			break;
		case Logger::WARN:
			debugLevel << "WARN";
			break;
		case Logger::INFO:
			debugLevel << "INFO";
			break;
		case Logger::DEBUG:
			debugLevel << "DEBUG";
			break;
		case Logger::DEV:
			debugLevel << "DEV";
			break;
		case Logger::ALL:
			debugLevel << "ALL";
			break;
		default:
			break;
	}
	return debugLevel.str();
}


/**
 * Returns a logger level given a string representation of an existing level
 *
 * @param a string for the level
 *
 * @return the logging level
 *
 */
Logger::LEVEL Logger::stringAsLoggerLevel( std::string level_string ) {
	
	Logger::LEVEL l;
		
	if ( "ERROR" == level_string ) {
		l = Logger::ERROR;
	}
		
	if ( "WARN" == level_string ) {
		l = Logger::WARN;
	}
		
	if ( "INFO" == level_string ) {
		l = Logger::INFO;
	}
		
	if ( "DEBUG" == level_string ) {
		l = Logger::DEBUG;
	}
		
	if ( "DEV" == level_string ) {
		l = Logger::DEV;
	}
			
	return l;
	
}

/**
 * Computes a logger level bit mask that will cover the parameter level and everything more important
 *
 * @param level a logger level mask
 *
 * @return the logging level
 *
 */
int Logger::levelAndBelow( Logger::LEVEL level ) {
	
	int l = level;
	
	if ( Logger::ERROR <= level ) {
		l = l | Logger::ERROR;
	}
	
	if ( Logger::WARN <= level ) {
		l = l | Logger::WARN;
	}
	
	if ( Logger::INFO <= level ) {
		l = l | Logger::INFO;
	}
	
	if ( Logger::DEBUG <= level ) {
		l = l | Logger::DEBUG;
	}
	
	if ( Logger::DEV <= level ) {
		l = l | Logger::DEV;
	}
	
	return l;
	
}

}	// namespace fdl