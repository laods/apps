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
	for (FI fi = ci->facebegin(); fi != ci->faceend(); ++fi) {
	    GridInterface::Vector v = pos - fi->centroid();
	    if (fi->normal()*v > 0) {
		succeed = false;
		break;
	    }
	    succeed = true;
	}
	if (succeed) {
	    c = ci->index();
	    break;	    
	}
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
    const double lmin [3] = {20, 20 ,1};
    vector< vector<double> > lvec;
    vector< vector<double> > corrvec;
    for (int dir=0; dir<3; ++dir) {
	cout << "Correlation direction " << dir << ":" << endl;
	vector<double> tmp;
	lvec.push_back(tmp);
	corrvec.push_back(tmp);
	for (double l = size[dir]/2; l > lmin[dir]; l -= lmin[dir]) {
	    cout << "Calculating correlation for l=" << l << "... ";
	    double sumpor = 0.0;
	    double varpor = 0.0;
	    int nsumpor = 0;
	    for (int i=0; i<N; ++i) {
		GridInterface::Vector pos;
		for (int d=0; d<3; ++d) {
		    if (d == dir) {
			pos[d] = ((double) rand() / (RAND_MAX))*(size[d]-l) + gridlimits[2*d];
		    }
		    else {
		    	pos[d] = ((double) rand() / (RAND_MAX))*size[d] + gridlimits[2*d];  
		    }
		}
		int cell = position2cell(gridinterf.cellbegin(), gridinterf.cellend(), pos);
		if (cell == -1) {
		    cerr << "Could not find cell index at position [" 
			<< pos[0] << ", " << pos[1] << ", " << pos[2] << "]" << endl;
		}
		
		GridInterface::Vector next_pos = pos;
		next_pos[dir] = pos[dir] + l;
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
	    
	    cout << "Correlation: " << sumpor << endl;
	    
	    lvec[dir].push_back(l);
	    corrvec[dir].push_back(sumpor);
	  
	}
	cout << endl;
    }
    
    stringstream outputtmp, outputtmpx, outputtmpy, outputtmpz;
    
    outputtmp << "#########################################################################" << endl
	      << "# Correlation lengths for model " << ECLIPSEFILENAME << endl
	      << "#########################################################################" << endl;
	      
    outputtmpx << "# Lengthx  Corrx" << endl;
    outputtmpy << "# Lengthy  Corry" << endl;
    outputtmpz << "# Lengthz  Corrz" << endl;

    for (int i=0; i<lvec[0].size(); ++i) {
	outputtmpx << lvec[0][i] << "\t" << corrvec[0][i] << endl;
    }
    for (int i=0; i<lvec[1].size(); ++i) {
	outputtmpy << lvec[1][i] << "\t" << corrvec[1][i] << endl;
    }
    for (int i=0; i<lvec[2].size(); ++i) {
	outputtmpz << lvec[2][i] << "\t" << corrvec[2][i] << endl;
    }
    
    cout << outputtmp.str() << outputtmpx.str() << outputtmpy.str() << outputtmpz.str();
    
    ofstream outfilex, outfiley, outfilez;
    
    outfilex.open("output_corrlengths_x.txt", ios::out | ios::trunc);
    outfilex << outputtmp.str() << outputtmpx.str();
    outfilex.close();
    
    outfiley.open("output_corrlengths_y.txt", ios::out | ios::trunc);
    outfiley << outputtmp.str() << outputtmpy.str();
    outfiley.close();    
    
    outfilez.open("output_corrlengths_z.txt", ios::out | ios::trunc);
    outfilez << outputtmp.str() << outputtmpz.str();
    outfilez.close();    
    
    cout << "Output written to files!" << endl;

    return 0;
}


