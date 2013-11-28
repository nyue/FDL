#include <boost/bind.hpp>

#include "logger/logger.h"
#include "io/exporterbase.h"

namespace fdl {

/**
 * Constructor
 *
 * @param prefix string filename prefix to apply to the files.
 * @param format file format to write to
 *
 */
ExporterBase::ExporterBase(std::string prefix)
{
	// m_filestream = new std::ofstream();
	m_isCancelled = m_writingFile = false;
	m_filenamePrefix = prefix;
	m_filenameCounter = 0;
}

/**
 * Destructor.
 *
 */
ExporterBase::~ExporterBase()
{
	// m_filestream->close();
	delete m_filestream;	// not used now..
}

/**
 * Dispatches thread to write the file to disk.
 *
 */
long int ExporterBase::start(fdl::Grid& grid)
{
	m_writeMutex.lock();
	
	m_grid = grid.copy();
		
	boost::thread t1(boost::bind(&ExporterBase::write,this));
	t1.join();
	
	m_writeMutex.unlock();
	return 0;
}

/**
 * NOT YET IMPLEMENTED
 *
 */
void ExporterBase::stop()
{
	//m_isCancelled = true;
}

}	// namespace fdl