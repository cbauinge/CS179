#ifndef READER_H
#define READER_H


#include "Domain.h"

///@brief Reader class which opens and reads the bitmap and returns a domain.
class Reader
{
    friend class ReaderTest;

public:
    Domain Read(const char* filename) const;

};


#endif /* READER_H */