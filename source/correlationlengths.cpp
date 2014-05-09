#include <iostream>
#include <cstdlib>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/porsol/common/Rock.hpp>
#include <dune/grid/CpGrid.hpp>

using namespace Opm;
using namespace Dune;
using namespace std;

typedef GridInterfaceEuler<Dune::CpGrid> GridInterface;

int position2cell(GridInterface gridinterface, GridInterface::Vector pos)
{
    typedef GridInterface::CellIterator CI;
    typedef CI::FaceIterator FI;
    int c = -1;
    bool succeed = false;
    for (CI ci = gridinterface.cellbegin(); ci != gridinterface.cellend(); ++ci) {
	c = ci->index();
	for (FI fi = ci->facebegin(); fi != ci->faceend(); ++fi) {
	    GridInterface::Vector v = pos - fi->centroid();
	    if (fi->normal()*v > 0) {
		break;
	    }
	    succeed = true;
	}
	if (succeed) break;
    }
    return c; 
}

int main(int argc, char** argv)
{
    // Process input parameters
    if (argc != 2) {
        cout << "Usage: correlationlengths filename.grdecl" << endl;
        exit(1);
    }
    const char* ECLIPSEFILENAME = argv[1];
    
    // Parse grdecl file
    ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
    if (eclipsefile.fail()) {
        cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
        exit(1);
    }
    eclipsefile.close();
    EclipseGridParser * eclParser_p;
    try {
        eclParser_p = new EclipseGridParser(ECLIPSEFILENAME);
    }
    catch (...) {
        cout << "Error: Filename " << ECLIPSEFILENAME << " does not look like an eclipse grid file." << endl;
        exit(1);
    }
    EclipseGridParser& eclParser = *eclParser_p;
    EclipseGridInspector eclInspector(eclParser);

    // Create grid
    CpGrid grid;
    grid.processEclipseFormat(eclParser, 0, false);
    GridInterface gridinterf(grid);

    // Load rock data
    Rock<3> rock;
    rock.init(eclParser, grid.globalCell());

    double porevol = 0.0;
    double volume  = 0.0;
    for (CI ci = gridinterface.cellbegin(); ci != gridinterface.cellend(); ++ci) {
	c = ci->index();
	volume  += cell_iter->volume();
	porevol += cell_iter->volume()*rock.porosity(c)
    }
    meanporo = porevol/volume;

    array<double,6> gridlimits = eclInspector.getGridLimits();
    array<double,3> size;
    size[0] = gridlimits[1]-gridlimits[0];
    size[1] = gridlimits[3]-gridlimits[2];
    size[2] = gridlimits[5]-gridlimits[4];

    const int N = 100;
    const double lmin = 10;
    vector<double> lvec;
    vector<double> corrvec;
    for (double l = size[0]/2; l > 10; l -= lmin) {
	double sumpor = 0.0;
	double varpor = 0.0;
	int nsumpor = 0;
	for (int i=0; i<N; ++i) {
	    GridInterface::Vector pos;
	    pos[0] = ((double) rand() / (RAND_MAX))*(length[0]-l) + gridlimits[0];
	    pos[1] = ((double) rand() / (RAND_MAX))*length[1] + gridlimits[2];
	    pos[2] = ((double) rand() / (RAND_MAX))*length[2] + gridlimits[4];
	    int cell = position2cell(gridinterf, pos);
	    if (cell == -1) {
		cerr << "Could not find cell index at position [" 
		     << pos[0] << ", " << pos[1] << ", " << pos[2] << "]" << endl;
	    }
	    
	    GridInterface::Vector next_pos = pos;
	    next_pos[0] = pos[0] + l;
	    int next_cell = position2cell(gridinterf, next_pos);
	    if (next_cell == -1) {
		cerr << "Could not find cell index at position [" 
		     << next_pos[0] << ", " << next_pos[1] << ", " << next_pos[2] << "]" << endl;
	    }

	    poro1 = rock.porosity(cell);
	    poro2 = rock-porosity(next_cell);

	    sumpor += (poro1 - meanpor) * (poro2 - meanpor);
	    varpor += (poro1 - meanpor) * (poro1 - meanpor);
	    nsumpor++;
	       
	}
	varpor = varpor/nsumpor;
	sumpor = sumpor/(nsumpor*varpor);
	lvec.push_back(l);
	corrvec.push_back(sumpor);
    }


    return 0;
}


