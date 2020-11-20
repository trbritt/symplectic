#include <iostream>
#include <fstream>
#include "particle.hpp"
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include "matrixMulCUBLAS.hpp"
#include <omp.h>
#include "mapping.hpp"
void write_current_to_file(std::string dir, Set& s, int tstep);
void write_current_to_file(std::string dir, Set& s, int tstep)
{
    static const auto bar_length = 70;
    std::string filename = dir;
    filename.append("/data/beam_" + boost::lexical_cast<std::string>(tstep) + ".txt");
    std::ofstream fh(filename);
    fh << "###############################" << std::endl;
    fh << "# ID, X, Y, Z, GBx, GBy, GBz" << std::endl;
    for (int i = 0; i < (int)s.size(); ++i) {
        if (i % 1000 == 0) std::cout << i << "/" << s.size() << "\r" << std::flush;
        fh << i + 1 << "," << s[i].x << "," << s[i].y << "," << s[i].z << "," << s[i].GB_x << "," << s[i].GB_y << "," << s[i].GB_z << std::endl;
    }
    fh << "###############################" << std::endl;
    fh.close();
}