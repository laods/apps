#include <iostream>
#include <cstdlib>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/Rock.hpp>
#include <dune/grid/CpGrid.hpp>

using namespace Opm;
using namespace Dune;
using namespace std;

typedef GridInterfaceEuler<CpGrid> GridInterface;
typedef GridInterface::CellIterator CI;
typedef CI::FaceIterator FI;

int position2cell(CI ibegin, CI iend, GridInterface::Vector pos)
{
    int c = -1;
    bool succeed = false;
    for (CI ci = ibegin; ci != iend; ++ci) {
	c = ci->index();
	for (FI fi = ci->facebegin(); fi != ci->faceend(); ++fi) {
	    GridInterface::Vector v = pos - fi->centroid();
	    if (fi->normal()*v > 0) {
		succeed = false;
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
    for (CI ci = gridinterf.cellbegin(); ci != gridinterf.cellend(); ++ci) {
	int c = ci->index();
	volume  += ci->volume();
	porevol += ci->volume()*rock.porosity(c);
    }
    double meanporo = porevol/volume;

    std::array<double,6> gridlimits = eclInspector.getGridLimits();
    std::array<double,3> size;
    size[0] = gridlimits[1]-gridlimits[0];
    size[1] = gridlimits[3]-gridlimits[2];
    size[2] = gridlimits[5]-gridlimits[4];

    const int N = 1000;
    const double lmin = 10;
    vector<double> lvec;
    vector<double> corrvec;
    for (double l = size[0]/2; l > 10; l -= lmin) {
	double sumpor = 0.0;
	double varpor = 0.0;
	int nsumpor = 0;
	for (int i=0; i<N; ++i) {
	    GridInterface::Vector pos;
	    pos[0] = ((double) rand() / (RAND_MAX))*(size[0]-l) + gridlimits[0];
	    pos[1] = ((double) rand() / (RAND_MAX))*size[1] + gridlimits[2];
	    pos[2] = ((double) rand() / (RAND_MAX))*size[2] + gridlimits[4];
	    int cell = position2cell(gridinterf.cellbegin(), gridinterf.cellend(), pos);
	    if (cell == -1) {
		cerr << "Could not find cell index at position [" 
		     << pos[0] << ", " << pos[1] << ", " << pos[2] << "]" << endl;
	    }
	    
	    GridInterface::Vector next_pos = pos;
	    next_pos[0] = pos[0] + l;
	    int next_cell = position2cell(gridinterf.cellbegin(), gridinterf.cellend(), next_pos);
	    if (next_cell == -1) {
		cerr << "Could not find cell index at position [" 
		     << next_pos[0] << ", " << next_pos[1] << ", " << next_pos[2] << "]" << endl;
	    }

	    double poro1 = rock.porosity(cell);
	    double poro2 = rock.porosity(next_cell);

	    sumpor += (poro1 - meanporo) * (poro2 - meanporo);
	    varpor += (poro1 - meanporo) * (poro1 - meanporo);
	    nsumpor++;
	   
	    /*
	    cout << "Position [" << pos[0] << ", " << pos[1] << ", " << pos[2] << "], "
		 << "(i " << cell << ") "
		 << "Next position [" << next_pos[0] << ", " << next_pos[1] << ", " << next_pos[2] << "], "
		 << "(i " << next_cell << ") "
		 << "Length " << l << ", poro1 " << poro1 << ", poro2 " << poro2 << endl;
	
	    */
	}
	varpor = varpor/nsumpor;
	sumpor = sumpor/(nsumpor*varpor);
	lvec.push_back(l);
	corrvec.push_back(sumpor);
    }

    stringstream outputtmp;
    
    outputtmp << "#########################################################################" << endl
	      << "# Correlation lengths for model " << ECLIPSEFILENAME << endl
	      << "# Length  Corr" << endl;
    for (int i=0; i<lvec.size(); ++i) {
	outputtmp << lvec[i] << "\t" << corrvec[i] << endl;
    }
    
    cout << outputtmp.str();
    
    ofstream outfile;
    outfile.open("output_corrlengths.txt", ios::out | ios::trunc);
    outfile << outputtmp.str();
    outfile.close();

    return 0;
}


