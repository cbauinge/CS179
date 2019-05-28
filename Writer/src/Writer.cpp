#include "Writer.h"
#include <fstream>


void Writer::Write(const char* filename, const Domain& dom, const std::vector<double>& result, const std::vector<double>& bc) const
{
    int nx = dom.GetPoints().size();
    int ny = dom.GetPoints()[0].size();

    std::vector<std::vector<double> > output(nx, std::vector<double>(ny, 0.0));
    if (result.size() != dom.GetInteriorWOBoundary().size())
        throw std::invalid_argument("Result vector does not have the same size as the interior points");

    for (int i = 0; i < result.size(); i++)
    {
        std::pair<int, int> coordinates = dom.GetInteriorWOBoundary()[i];
        
        output[coordinates.first][coordinates.second] = result[i];
    }

    for (int i = 0; i < bc.size(); i++)
    {
        int coordinatex = dom.GetBoundary().GetOrderedPoint(i).i;
        int coordinatey = dom.GetBoundary().GetOrderedPoint(i).j;

        output[coordinates.first][coordinates.second] = bc[i];
    }

    std::ofstream ofs(filename, std::ios::out);
    if (!ofs)
        throw std::invalid_argument("Invalid filename in Writer::Write");

    for (int i = 0; i < output.size(); i++)
    {
        for (int j = 0; j < output[0].size(); j++)
        {
            ofs << output[i][j];
            if (j < output[0].size()-1)
                ofs << ", ";
        }
        if (i < output.size()-1)
            ofs << std::endl;
    }

    ofs.close();
}