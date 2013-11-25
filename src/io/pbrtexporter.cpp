#include <png.h>
#include "io/pbrtexporter.h"
#include "logger/logger.h"

namespace fdl {

/**
 * Constructor.
 *
 */	
PbrtExporter::PbrtExporter(std::string prefix) : ExporterBase(prefix)
{
	
}

/**
 * Destructor.
 *
 */
PbrtExporter::~PbrtExporter()
{
	
}

/**
 * Writes a file for the input grid.
 *
 * @param grid a reference to the grid to pull numbers from
 *
 */
// void PbrtExporter::write(const fdl::Grid& grid)
void PbrtExporter::write()
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
 * Exports a 3D density field as a PBRT file. A <a href="http://www.cse.ohio-state.edu/~parent/classes/782/labs/PBRT/pbrtscene.html">PBRT format</a>
 * is used by the <a href="http://www.pbrt.org">PBRT renderer.</a>
 *
 * @param counter number extension for the output file
 * @param prefix string filename to write to
 * @param field array of float densities
 * @param xRes resolution of the grid in the x dimension
 * @param yRes resolution of the grid in the y dimension
 * @param zRes resolution of the grid in the z dimension
 *
 */
void PbrtExporter::exportDensity(int counter, std::string prefix, float* field, int xRes, int yRes, int zRes)
{
	char buffer[256];
	sprintf(buffer,"%04i", counter);
	std::string number = std::string(buffer);

	std::string filenamePbrt = prefix + number + std::string(".pbrt.gz");
	DEV() << "Writing PBRT " << filenamePbrt;

	//float* field = new float[xRes * yRes * zRes];
	float maxDensVal = fabs(field[0]);
	float targetNorm = 0.5;
	for (int i = 0; i < xRes * yRes * zRes; i++) {
		if(fabs(field[i])>maxDensVal) maxDensVal = fabs(field[i]);
		//field[i] = 0.0;
	}
	if(maxDensVal>0.0) {
		for (int i = 0; i < xRes * yRes * zRes; i++) {
			field[i] = fabs(field[i]) / maxDensVal * targetNorm;
		}
	}

	std::fstream fout;
	fout.open(filenamePbrt.c_str(), std::ios::out);

	int maxRes = (xRes > yRes) ? xRes : yRes;
	maxRes = (maxRes > zRes) ? maxRes : zRes;

	const float xSize = 1.0 / (float)maxRes * (float)xRes;
	const float ySize = 1.0 / (float)maxRes * (float)yRes;
	const float zSize = 1.0 / (float)maxRes * (float)zRes;

	gzFile file;
	file = gzopen(filenamePbrt.c_str(), "wb1"); 
	if (file == NULL) {
		ERROR() << " Couldn't write file " << filenamePbrt << "!";
		return;
	}

	// write file
	gzprintf(file, "Volume \"volumegrid\" \n");
	gzprintf(file, " \"integer nx\" %i\n", xRes);
	gzprintf(file, " \"integer ny\" %i\n", yRes);
	gzprintf(file, " \"integer nz\" %i\n", zRes);
	gzprintf(file, " \"point p0\" [ 0.0 0.0 0.0 ] \"point p1\" [%f %f %f ] \n", xSize, ySize, zSize);
	gzprintf(file, " \"float density\" [ \n");
	for (int i = 0; i < xRes * yRes * zRes; i++){
		gzprintf(file, "%f ", field[i]);
		
		if(m_isCancelled){
			gzclose(file);
			delete[] field;
			return;
		}
	}
	gzprintf(file, "] \n \n");

	gzclose(file);
	delete[] field;
}

}	// namespace fdl