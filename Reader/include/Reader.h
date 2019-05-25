#ifndef READER_H
#define READER_H


#include "Domain.h"

////For the reader, I followed the post here
////https://stackoverflow.com/questions/49215933/reading-a-monochrome-bitmap-in-c-requires-reading-every-other-line

class Reader
{
    friend class ReaderTest;

public:
    Domain Read(const char* filename) const;


private:
    uint32_t make_stride_aligned(uint32_t align_stride);


    uint32_t row_stride{ 0 };

    #pragma pack(1)
    struct BMPFileHeader {
        char magic[2];          // 0-1
        uint32_t fileSize;      // 2-5
        uint32_t reserved;      // 6-9
        uint32_t offset;        // 10-13
        uint32_t headerSize;    // 14-17
        uint32_t width;         // 18-21
        uint32_t height;        // 22-25
        uint16_t bitsPerPixel;  // 26-27
        uint16_t bitDepth;      // 28-29
    };
    #pragma pack()

};


#endif /* REAER_H */