#include "Writer.h"
#include <fstream>


void Writer::Write(const char* filename, const Domain& dom, const std::vector<double>& result) const
{
    int nx = dom.GetPoints().size();
    int ny = dom.GetPoints()[0].size();

    std::vector<std::vector<double> > output(nx, std::vector<double>(ny, 0.0));

    for (int i = 0; i < result.size(); i++)
    {
        std::pair<int, int> coordinates = dom.GetInterior()[i];
        
        output[coordinates.first][coordinates.second] = result[i];
    }

    std::ofstream ofs(filename, std::ios::out);
    if (!ofs)
        throw std::invalid_argument("Invalid filename in Writer::Write");

    for (int i = 0; i < output.size(); i++)
    {
        for (int j = 0; j < output[0].size(); j++)
        {
            ofs << output[i][j] << ", ";
        }
        ofs << std::endl;
    }

    ofs.close();
}