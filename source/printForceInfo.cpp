/* Arguments:
 *   1) eclipse file
 *   2) surface tension (sigma) in N/m
 *   3) Pressure drop in Pa
 */

#include "config.h"

#include <iostream>
#include <fstream>

#include <dune/grid/CpGrid.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/porsol/common/GridInterfaceEuler.hpp>
#include <opm/porsol/common/Rock.hpp>
#include <opm/core/utility/MonotCubicInterpolator.hpp>

using namespace std;


class RockInterface
{
public:
    RockInterface() {};
    void read(const char* file_list);
    void calcFracFlow(vector<double> visc, int points);
    void calcSatAtVL(double fracflow);
    const double cpAtVL(int satnum);
    
private:
    vector<Opm::MonotCubicInterpolator> cp_;
    vector<Opm::MonotCubicInterpolator> krw_;
    vector<Opm::MonotCubicInterpolator> kro_;
    vector<Opm::MonotCubicInterpolator> fw_;
    vector<Opm::MonotCubicInterpolator> fw_inv_;
    vector<double> sat_VL_;
    int nSatnum_;
};

void RockInterface::read(const char* file_list)
{
    // Open file_list
    ifstream file_list_stream(file_list, ios::in);
    if (file_list_stream.fail()) {
        cerr << "Error: Filename " << file_list << " not found or not readable." << endl;
        exit(1);
    }
    
    // Read rockfile one by one
    string rockfile;
    nSatnum_ = 0;
    getline(file_list_stream, rockfile);
    while ( getline(file_list_stream, rockfile) ) {
        Opm::MonotCubicInterpolator rock_cp(rockfile.c_str(), 1, 4);
        cp_.push_back(rock_cp);
        Opm::MonotCubicInterpolator rock_krw(rockfile.c_str(), 1, 2);
        krw_.push_back(rock_krw);
        Opm::MonotCubicInterpolator rock_kro(rockfile.c_str(), 1, 3);
        kro_.push_back(rock_kro);
        ++nSatnum_;
    }
    file_list_stream.close();
}

void RockInterface::calcFracFlow(vector<double> visc, int points = 100)
{
    for (int r=0; r<nSatnum_; ++r) {
        double swir = krw_[r].getMinimumX().first;
        double swor = krw_[r].getMaximumX().first;
        vector<double> fracflowVec(points, 0.0);
        vector<double> s(points, 0.0);
        for (int i=0; i<points; ++i) {
            s[i] = swir + i*(swor-swir)/(points-1);
            double mobW = krw_[r].evaluate(s[i]) / visc[0];
            double mobO = kro_[r].evaluate(s[i]) / visc[1];
            fracflowVec[i] = mobW / (mobW + mobO);
        }
        Opm::MonotCubicInterpolator fracflow(s, fracflowVec);
        fw_.push_back(fracflow);
        if (fracflow.isStrictlyMonotone()) {
            fw_inv_.push_back(Opm::MonotCubicInterpolator(fracflow.get_fVector(), fracflow.get_xVector()));
        }
        else { 
            cerr << "Frac flow function for rock " << r << " is not invertible!" << endl; 
        }
    }
}

void RockInterface::calcSatAtVL(double fracflow = 0.5)
{
    for (int r=0; r<nSatnum_; ++r) {
        sat_VL_.push_back(fw_inv_[r].evaluate(fracflow));
    }
}

const double RockInterface::cpAtVL(int satnum)
{
    assert(satnum <= nSatnum_);
    return cp_[satnum-1].evaluate(sat_VL_[satnum-1]);
}


int main(int varnum, char** vararg) {
    
    // Read input
    const char* ECLIPSEFILENAME(vararg[1]);
    const int dir = atoi(vararg[2]);
    const double sigma = atof(vararg[3]);
    const double deltap = atof(vararg[4]);
    const string section_name = string(vararg[5]);
    const char* ROCKLIST(vararg[6]);
    const double viscW = atof(vararg[7]);
    const double viscO = atof(vararg[8]);
    const double satur_for_visc = atof(vararg[9]);
     
    // Cap pressure data rocks and map poro and satnum
    // OBS!!! HARD-CODED!
    bool isParallel = false;
    double Pc_avg[20] = { };
    double Pb[20] = { };
    if (section_name == "SN74" || section_name == "sn74") {
	Pc_avg = {23771,86713,1987,50694,28470,66404,26280,58528,4547,79796,61050,10623};
	Pb = {3038,12617,679,3712,5526,8420,4949,9626,1440,14447,8060,3386};
    }
    else if (section_name == "SN82" || section_name == "sn82") {
	Pc_avg = {23771,86713,50694,20412,2882,31040,4547,79796,20369,10623};
	Pb = {3038,12617,3712,1630,1153,10224,1440,14447,1573,3386};
    }
    else if (section_name == "SN91" || section_name == "sn91") {
	Pc_avg = {23771,86713,1987,50694,78336,3516,9222,66404,6282,5675,20412,26280,4547,79796,7158};
	Pb = {3038,12617,679,3712,17123,1327,1869,8420,1188,1401,1630,4949,1440,14447,1775};
    }
    else if (section_name == "Parallel" || section_name == "parallel") {
	Pc_avg = {1.2286e4,2.7473e5};
	Pb = {4554.9,1.0185e5};
	isParallel = true;
    }
    else {
	cerr << "Model type (argument 5: '" << section_name << "') unknown. Either SN74, SN82, SN91 or Parallel.";
	exit(1);
    }

    // Test if filename exists and is readable
    ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
    if (eclipsefile.fail()) {
        cerr << "Error: Filename " << ECLIPSEFILENAME << " not found or not readable." << endl;
        exit(1);
    }
    eclipsefile.close();
    
    // Parse grdecl file
    //std::cout << "Parsing grid file '" << ECLIPSEFILENAME << "' ..." << std::endl;
    Opm::ParserPtr parser(new Opm::Parser);
    Opm::DeckConstPtr deck(parser->parseFile(ECLIPSEFILENAME));

    Dune::CpGrid grid;
    grid.processEclipseFormat(deck, 0, false);
    
    typedef Opm::GridInterfaceEuler<Dune::CpGrid> GridInterface;
    GridInterface gridinterf(grid);

    std::vector<int> satnum_from_deck = deck->getKeyword("SATNUM")->getIntData();
    std::vector<int> satnum(grid.globalCell().size());
    for (int i=0; i<satnum.size(); ++i) {
      satnum[i] = satnum_from_deck[grid.globalCell()[i]];
    }
    
    Opm::Rock<3> rock;
    rock.init(deck, grid.globalCell());
   
    // Read and process rock curves
    RockInterface ri;
    ri.read(ROCKLIST);
    vector<double> visc;
    visc.push_back(viscW);
    visc.push_back(viscO);
    ri.calcFracFlow(visc);
    ri.calcSatAtVL(satur_for_visc);
    
    double accumulated_porevol = 0.0;
    double nabla_Pc_J = 0.0;
    double nabla_Pc_avg = 0.0;
    double nabla_Pb = 0.0;
    double nabla_Pc_VL = 0.0;
    double min_coord = 1e16;
    double max_coord = -1e16;
    
    typedef GridInterface::CellIterator CI;
    typedef CI::FaceIterator FI;
    int i = 0;
    for (CI cell_iter = gridinterf.cellbegin(); cell_iter != gridinterf.cellend(); ++cell_iter) {
        int cell_index = cell_iter->index();
        const double cell_perm = rock.permeability(cell_index)(dir,dir);
        double cell_poro = rock.porosity(cell_index);
        const double cell_porevol = cell_iter->volume()*cell_poro;
        const double cell_lambda = sqrt(cell_perm/cell_poro);
	const int cell_satnum = satnum[cell_index];
        const double cell_Pc_avg = Pc_avg[cell_satnum];
        const double cell_Pb = Pb[cell_satnum];
        const double cell_Pc_VL = ri.cpAtVL(cell_satnum);
        for (FI face_iter = cell_iter->facebegin(); face_iter != cell_iter->faceend(); ++face_iter) {
            if ( ! face_iter->boundary()) {
                int neighbour_index = face_iter->neighbourCellIndex();
                const double neighbour_perm = rock.permeability(neighbour_index)(dir,dir);
                double neighbour_poro = rock.porosity(neighbour_index);
                const double neighbour_porevol = face_iter->neighbourCellVolume()*neighbour_poro;
                const double neighbour_lambda = sqrt(neighbour_perm/neighbour_poro);
		const int neighbour_satnum = satnum[neighbour_index];
                const double neighbour_Pc_avg = Pc_avg[neighbour_satnum];
                const double neighbour_Pb = Pb[neighbour_satnum];
                const double neighbour_Pc_VL = ri.cpAtVL(neighbour_satnum);
                const double length = (cell_iter->centroid() - face_iter->neighbourCell().centroid()).two_norm();
                const double porevol = cell_porevol + neighbour_porevol;
                accumulated_porevol += porevol;
                nabla_Pc_J += abs(1/cell_lambda - 1/neighbour_lambda)*porevol/length;
                nabla_Pc_avg += abs(cell_Pc_avg - neighbour_Pc_avg)*porevol/length;
                nabla_Pb += abs(cell_Pb - neighbour_Pb)*porevol/length;
                nabla_Pc_VL += abs(cell_Pc_VL - neighbour_Pc_VL)*porevol/length;
            }
            else {
                min_coord = min(min_coord, face_iter->centroid()[dir]);
                max_coord = max(max_coord, face_iter->centroid()[dir]);
            }
        }
	++i;
    }
    nabla_Pc_J = sigma*nabla_Pc_J/accumulated_porevol;
    nabla_Pc_avg = nabla_Pc_avg/accumulated_porevol;
    nabla_Pb = nabla_Pb/accumulated_porevol;
    nabla_Pc_VL = nabla_Pc_VL/accumulated_porevol;
    
    double model_length = max_coord - min_coord;
    double nabla_visc_pressure = deltap/model_length;
    
    if (false) {
      cout << scientific;
      cout << endl;
      cout << "Pressure drop:        " << deltap << endl;
      cout << "||nabla p||           " << nabla_visc_pressure << "\t(or " << 1/model_length << "*DeltaP)" << endl;
      cout << endl;
      cout << "With |J|=1:" << endl;
      cout << "||nabla Pc||          " << nabla_Pc_J << endl;
      cout << "Ca:                   " << nabla_visc_pressure/nabla_Pc_J << "\t(or " << 1/(model_length*nabla_Pc_J) << "*DeltaP)" << endl;
      cout << "Ca = 1 when Delta p = " << nabla_Pc_J * model_length << endl;
      cout << endl;
      cout << "With Pc_avg:" << endl;
      cout << "||nabla Pc||          " << nabla_Pc_avg << endl;
      cout << "Ca:                   " << nabla_visc_pressure/nabla_Pc_avg << "\t(or " << 1/(model_length*nabla_Pc_avg) << "*DeltaP)" << endl;
      cout << "Ca = 1 when Delta p = " << nabla_Pc_avg * model_length << endl;
      cout << endl;
      cout << "With Pb:" << endl;
      cout << "||nabla Pc||          " << nabla_Pb << endl;
      cout << "Ca:                   " << nabla_visc_pressure/nabla_Pb << "\t(or " << 1/(model_length*nabla_Pb) << "*DeltaP)" << endl;
      cout << "Ca = 1 when Delta p = " << nabla_Pb * model_length << endl;
      cout << endl;
      cout << "With Pc at VL:" << endl;
      cout << "||nabla Pc||          " << nabla_Pc_VL << endl;
      cout << "Ca:                   " << nabla_visc_pressure/nabla_Pc_VL << "\t(or " << 1/(model_length*nabla_Pc_VL) << "*DeltaP)" << endl;
      cout << "Ca = 1 when Delta p = " << nabla_Pc_VL * model_length << endl;
      cout << endl;
    }
    else {
      cout << nabla_visc_pressure/nabla_Pc_VL << endl;
   }
}
