#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include "omp.h"
#include "tiffio.h"

#include "base/Global.h"
#include "base/SOA3D.h"
#include "base/Vector3D.h"
#include "base/Transformation3D.h"
#include "volumes/PlacedVolume.h"
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedBox.h"
#include "management/GeoManager.h"
#include "base/Stopwatch.h"

#undef VERBOSE // silence navigator output

#include "navigation/SimpleNavigator.h"

using namespace vecgeom;

// external functions from libCMS geometry library

extern LogicalVolume * lvol0;
extern Transformation3D * idtrans;
extern void GenerateTransformations_part0();
extern void GenerateTransformations_part1();
extern void GenerateTransformations_part2();
extern void GenerateTransformations_part3();
extern void CreateLogicalVolumes();
extern void GeneratePlacedVolumes_part0();
extern void GeneratePlacedVolumes_part1();
extern void GeneratePlacedVolumes_part2();
extern void GeneratePlacedVolumes_part3();
extern void GeneratePlacedVolumes_part4();

VPlacedVolume const * generateDetector() {
	GenerateTransformations_part0();
	GenerateTransformations_part1();
	GenerateTransformations_part2();
	GenerateTransformations_part3();
	CreateLogicalVolumes();
	GeneratePlacedVolumes_part0();
	GeneratePlacedVolumes_part1();
	GeneratePlacedVolumes_part2();
	GeneratePlacedVolumes_part3();
	GeneratePlacedVolumes_part4();
	return lvol0->Place( idtrans );
}

int write_amira(const char* fname, unsigned char *data, uint32_t nx, uint32_t ny, uint16_t nz)
{
	FILE* output;

	if (!(output = fopen(fname, "w"))) {
		fprintf(stderr, "%s\n", strerror(errno));
		return 1;
	}

	fprintf(output,"# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n\n");
	fprintf(output,"define Lattice %d %d %d\n\n", nx, ny, nz);
	fprintf(output,"Parameters {\n   Content \"%dx%dx%d byte, uniform coordinates\",\n", nx, ny, nz);
	fprintf(output,"   BoundingBox 0.000 %.3f 0.000 %.3f 0.000 %.3f\n   CoordType \"uniform\",\n}\n\n",
		(float) nx, (float) ny, (float) nz);
	fprintf(output,"Lattice { byte Data } @1\n\n@1\n");
	fwrite(data, sizeof(unsigned char), nx*ny*nz, output);
	fclose(output);
	return 0;
}

int write_tiff(const char* fname, unsigned char *data, uint32_t nx, uint32_t ny, uint16_t nz)
{
    float xres, yres;
    uint16_t spp, bpp, photo, res_unit;
    TIFF *output = TIFFOpen(fname, "w");

    if (!output) {
		fprintf (stderr, "Can't open tiff file for writing\n");
		return 1;
    }

    spp = 1; /* Samples per pixel */
    bpp = 8; /* Bits per sample */
    photo = PHOTOMETRIC_MINISBLACK;

    for (uint16_t page = 0; page < nz; page++)
    {
        TIFFSetField(output, TIFFTAG_IMAGEWIDTH, nx / spp);
        TIFFSetField(output, TIFFTAG_IMAGELENGTH, ny);
        TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, bpp);
        TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, spp);
        TIFFSetField(output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(output, TIFFTAG_PHOTOMETRIC, photo);
        TIFFSetField(output, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        /* It is good to set resolutions too (but it is not nesessary) */
        xres = yres = 100;
        res_unit = RESUNIT_INCH;
        TIFFSetField(output, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(output, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(output, TIFFTAG_RESOLUTIONUNIT, res_unit);

        /* We are writing single page of the multipage file */
        TIFFSetField(output, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        /* Set the page number */
        TIFFSetField(output, TIFFTAG_PAGENUMBER, page, nz);

        for (uint32_t j = 0; j < ny; j++)
            TIFFWriteScanline(output, &data[(page*ny + j) * nx], j, 0);

        TIFFWriteDirectory(output);
    }

    TIFFClose(output);
    return 0;
}

int main(int argc, char *argv[])
{
	if(argc < 2) {
		fprintf(stderr, "Usage: %s VolumeName MaxLayers", argv[0]);
		return 0;
	}

	GeoManager &geoManager = GeoManager::Instance();

	geoManager.SetWorld(generateDetector());

	VPlacedVolume *foundvolume  = geoManager.FindPlacedVolume(argv[1]);
	LogicalVolume *foundlvolume = geoManager.FindLogicalVolume(argv[1]);

	if (!foundlvolume) {
		fprintf(stderr, "Volume %s not found!\n", argv[1]);
		return 1;
	}

	long int layers = strtol(argv[2], NULL, 10);

	if (layers < 0 || layers > 1 << 20) {
		fprintf(stderr, "Number of layers out of range\n");
		return 1;
	}

	Vector3D<Precision> min, max;

	foundvolume->Extent(min, max);

	Vector3D<Precision> size = max - min;

	printf("Extent of bounding box: [%.3f, %.3f, %.3f]\n",
		size[0], size[1], size[2]);

	UnplacedBox worldbox = UnplacedBox(size);

	LogicalVolume lworld = LogicalVolume("world", &worldbox);

	lworld.PlaceDaughter(foundvolume);

	VPlacedVolume *world = lworld.Place();

	geoManager.SetWorld(world);

	Vector3D<Precision> origin = 0.5 * (max + min);

	double max_size = std::max(size[0], std::max(size[1], size[2]));

	int nx = layers * size[0] / max_size;
	int ny = layers * size[1] / max_size;
	int nz = layers * size[2] / max_size;

	double dx = size[0] / nx;
	double dy = size[1] / ny;
	double dz = size[2] / nz;

	Stopwatch timer;

	geoManager.CloseGeometry();

	timer.Start();

	// create volumetric image here
	int n_threads = omp_get_max_threads();
	int max_depth = geoManager.getMaxDepth();

	NavigationState * nstates[n_threads];

	for(size_t i = 0; i < n_threads; ++i) {
		nstates[i] = NavigationState::MakeInstance(max_depth);
	}

	unsigned char *data = new unsigned char[nx*ny*nz];

	// shift by half a pixel width to fall in the center of each voxel
	min = min + 0.5 * Vector3D<Precision>(dx, dy, dz);

#pragma omp parallel for shared(data) collapse(3) schedule(dynamic)
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				size_t thread_id = omp_get_thread_num();
				SimpleNavigator nav;
				NavigationState *nav_state = nstates[thread_id];
				Vector3D<Precision> p = min + Vector3D<Precision>(i*dx, j*dy, k*dz);

				nav_state->Clear();

				const VPlacedVolume *volume = nav.LocatePoint(geoManager.GetWorld(), p, *nav_state, true);

				if (volume != NULL && volume->GetLogicalVolume()->GetDaughtersp()->size() == 0) {
					data[(k*ny + j)*nx + i] = 255 * nav_state->GetCurrentLevel() / max_depth;
				}
			}
		}
	}

	for(size_t i = 0; i < n_threads; ++i) {
		NavigationState::ReleaseInstance(nstates[i]);
	}

	timer.Stop();

	char filename[64];

	strncpy(filename, argv[1], 58);
	strncat(filename, ".tiff", 6);
	write_tiff(filename, data, nx, ny, nz);

	printf("Elapsed time %.3f\n", timer.Elapsed());

	return 0;
}
