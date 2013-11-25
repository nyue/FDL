#include <png.h>

#include "io/pngexporter.h"
#include "logger/logger.h"

namespace fdl {

/**
 * Constructor.
 *
 */	
PngExporter::PngExporter(std::string prefix) : ExporterBase(prefix)
{
	
}

/**
 * Destructor.
 *
 */
PngExporter::~PngExporter()
{
	
}

/**
 * Writes a file for the input grid.
 *
 * @param grid a reference to the grid to pull numbers from
 *
 */
// void PngExporter::write(const fdl::Grid& grid)
void PngExporter::write()
{
	m_writingFile = true;
	
	int numCells = m_grid->getNumberOfGridCells();
	int sizeX = m_grid->getGridSizeX();
	int sizeY = m_grid->getGridSizeY();
	int sizeZ = m_grid->getGridSizeZ();
	float* density = m_grid->getDensityArray();
	exportDensity(m_filenameCounter, m_filenamePrefix, density, sizeX, sizeY, sizeZ);

	m_filenameCounter++;
	
	m_writingFile = false;
}

/**
 * Exports a 3D density field as a 2D PNG file by accumulating the density values along the z dimension.
 *
 * @param counter number extension for the output file
 * @param prefix string filename to write to
 * @param field array of float densities
 * @param xRes resolution of the grid in the x dimension
 * @param yRes resolution of the grid in the y dimension
 * @param zRes resolution of the grid in the z dimension
 *
 */
void PngExporter::exportDensity(int counter, std::string prefix, float* field, int xRes, int yRes, int zRes)
{
	char buffer[256];
	sprintf(buffer,"%04i", counter);
	std::string number = std::string(buffer);
	int totalSize = xRes * yRes * zRes;

	unsigned char pngbuf[xRes*yRes*4];
	unsigned char* rows[yRes];
	for (int j=0; j<yRes; ++j) {
		for (int i=0; i<xRes; ++i) {
			float val = 0;
			for (int k=0; k<zRes; ++k) {
				int index = i + yRes*(j + zRes*k);
				val += field[index];
			}
			if(val>100.0) val = 100.0;
			if(val<0.0) val = 0.0;
			val /= 100.0;

			pngbuf[(j*xRes + i)*4+0] = (unsigned char)(val*255.0);	// R
			pngbuf[(j*xRes + i)*4+1] = (unsigned char)(val*255.0);	// G
			pngbuf[(j*xRes + i)*4+2] = (unsigned char)(val*255.0);	// B
			pngbuf[(j*xRes + i)*4+3] = 255;							// A
			rows[j] = &pngbuf[(yRes-j-1) * xRes*4];
		}
	}
	std::string filenamePNG = prefix + number + std::string(".png");

	DEV() << "Writing " << filenamePNG;
	writePNG(filenamePNG.c_str(), rows, xRes, yRes);
}

/**
 * Writes a PNG of chars
 *
 * @param filename for output png
 * @param rowsp array of chars to write
 * @param w width of the output png image
 * @param h height of the output png image
 *
 */
int PngExporter::writePNG(const char* filename, unsigned char** rowsp, int w, int h)
{
	// defaults 
	const int colortype = PNG_COLOR_TYPE_RGBA;
	const int bitdepth = 8;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep* rows = rowsp;

	FILE* file = NULL;
	file = fopen(filename, "wb");
	if (file == NULL) {
		ERROR() << "PngExporter:\tcould not open for writing!";
		return -1;
	}

	if(!png_ptr) {
		png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		if (png_ptr == NULL) {
			fclose(file);
			ERROR() << "PngExporter:\tcould not create png write struct!";
			return -1;
		}
	}
	if(!info_ptr) {
		info_ptr = png_create_info_struct(png_ptr);
		if (info_ptr == NULL) {
			fclose(file);
			ERROR() << "PngExporter:\tcould not create png info struct!";
			return -1;
		}
	}

	if (setjmp(png_jmpbuf(png_ptr))){
		ERROR() << "PngExporter:\tFAILED.";
		fclose(file);
		png_destroy_write_struct(&png_ptr, &info_ptr);
		return -1;
	}
	//init IO
	png_init_io(png_ptr, file);
	//write header
	png_set_IHDR(png_ptr, info_ptr, w, h, bitdepth, colortype, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	// write info
	png_write_info(png_ptr, info_ptr);
	// write image
	png_write_image(png_ptr, rows);
	// write end
	png_write_end(png_ptr, NULL);
	// write destroy structs
	png_destroy_write_struct(&png_ptr, &info_ptr);

	fclose(file);
	
	return 0;
}

}