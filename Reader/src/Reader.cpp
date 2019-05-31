#include "Reader.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>


Domain Reader::Read(const char* filename) const
{
    std::vector<std::vector<bool> > points;

    std::ifstream f(filename);
    if (!f) 
        throw std::invalid_argument("Invalid file given");

    int m, n;
    f >> m >> n;

    points.resize(m, std::vector<bool>(n, false));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double tmp;
            f >> tmp;
            points[i][j] = bool(tmp);
        }
    }

    Domain dom(points);
    return dom;
}