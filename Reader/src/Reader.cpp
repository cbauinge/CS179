#include "Reader.h"

#include <vector>
#include <fstream>
#include <cmath>


Domain Reader::Read(const char* filename) const
{
    std::vector<std::vector<bool> > points;
    BMPFileHeader file_header;

    std::ifstream f(filename, std::ios::binary);

    if (!f) 
        throw "Invalid file given";

    f.read((char*)&file_header, sizeof(BMPFileHeader));

    if (file_header.bitsPerPixel != 1) 
    {
        f.close();
        throw "Invalid bitmap loaded, not monochrome";
    }

    int height = file_header.height;
    int width = file_header.width;

    // Lines are aligned on a 4-byte boundary
    int lineSize = (width / 8 + (width / 8) % 4) * 2;
    int fileSize = lineSize * height;

    std::vector<unsigned char> rawFile(fileSize);
    points.resize(file_header.height, std::vector<bool>(width, 0));

    // Skip to where the actual image data is
    f.seekg(file_header.offset);

    // Read in all of the file
    f.read((char*)&rawFile[0], fileSize);

    // Decode the actual boolean values of the pixesl
    int row;
    int reverseRow; // Because bitmaps are stored bottom to top for some reason
    int columnByte;
    int columnBit;

    for (row = 0, reverseRow = height - 1; row < height; ++row, --reverseRow) 
    {
        columnBit = 0;
        for (columnByte = 0; columnByte < std::ceil((width / 8.0)); ++columnByte) 
        {
            int rawPos = (row * lineSize) + columnByte;

            for (int k = 7; k >= 0 && columnBit < width; --k, ++columnBit) {
                points[reverseRow][columnBit] = (rawFile[rawPos] >> k) & 1;
            }
        }
    }

    f.close();

    Domain dom(points);
    return dom;
}