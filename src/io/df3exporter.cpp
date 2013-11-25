#include <stdio.h>
#include <algorithm>

#include "io/df3exporter.h"
#include "logger/logger.h"

namespace fdl {

/**
 * Constructor.
 *
 */	
Df3Exporter::Df3Exporter(std::string prefix) : ExporterBase(prefix)
{
}

/**
 * Destructor.
 *
 */
Df3Exporter::~Df3Exporter()
{
}

/**
 * Writes a file for the input grid.
 *
 * @param grid a reference to the grid to pull numbers from
 *
 */
void Df3Exporter::write()
{
	// m_writeMutex.lock();
	
	m_writingFile = true;
	
	int numCells = m_grid->getNumberOfGridCells();
	int sizeX = m_grid->getGridSizeX();
	int sizeY = m_grid->getGridSizeY();
	int sizeZ = m_grid->getGridSizeZ();
	
	float* density = m_grid->getDensityArray();
	exportDensity(m_filenameCounter, m_filenamePrefix, density, sizeX, sizeY, sizeZ);

	m_filenameCounter++;
	
	m_writingFile = false;
	
	// m_writeMutex.unlock();
}


/**
 * <strong><em>UNTESTED!!</em></strong><br/>
 * Exports density field as a DF3 file. DF3s can be rendered using <a href="http://www.povray.org/">POV-Ray.</a>
 * 
 *
 * @param counter number extension for the output file
 * @param prefix string filename to write to
 * @param field array of float densities
 * @param xRes resolution of the grid in the x dimension
 * @param yRes resolution of the grid in the y dimension
 * @param zRes resolution of the grid in the z dimension
 *
 */
void Df3Exporter::exportDensity(int counter, std::string prefix, float* field, int xRes, int yRes, int zRes)
{
	char buffer[256];
	sprintf(buffer,"%04i", counter);
	std::string number = std::string(buffer);
	
	int i,j,k, index;
	float value = 0;
	float range = 0;
	float min = 1e32;
	float max = -1e32;
	FILE* fptr;

	// Calculate the bounds
	for (i=0;i<xRes;i++) {
		for (j=0;j<yRes;j++) {
			for (k=0;k<zRes;k++) {
				index = (k * yRes * xRes) + (j * xRes) + i;
				max = std::max(max,field[index]);
				min = std::min(min,field[index]);
			}
		}
	}
	if (min >= max) {
		max = min + 1;
		min -= 1;
	}

	// Write it to a file
	std::string filename = prefix + number + std::string(".df3");
	DEV() << "Writing DF3 file " << filename;
	
	if((fptr = fopen(filename.c_str(), "w")) == NULL){
		ERROR() << " Couldn't write file " << filename << "!";
		return;
	}
	
	// Start the file
	fputc(xRes >> 8,fptr);
	fputc(xRes & 0xff,fptr);
	fputc(yRes >> 8,fptr);
	fputc(yRes & 0xff,fptr);
	fputc(zRes >> 8,fptr);
	fputc(zRes & 0xff,fptr);
	
	range = max - min;
	for (k=0;k<zRes;k++) {
		for (j=0;j<yRes;j++) {
			for (i=0;i<xRes;i++) {
				index = (k * yRes * xRes) + (j * xRes) + i;
				value = 255 * (field[index] - min)/range;
				fputc((int)value, fptr);
				
				if(m_isCancelled){
					fclose(fptr);
					return;
				}
			}
		}
	}
	
	fclose(fptr);
}

}	// namespace fdl