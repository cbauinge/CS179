#ifndef WRITER_H
#define WRITER_H

#include "Domain.h"

///class that writes the result to a file than can be read by Matlab.
///Very simple first version without much options.
class Writer
{
public:
    ///function that takes a filename, the domain and the result and writes the result in csv 
    ///format for matlab to be able to read it as a matrix for visualization.
    void Write(const char* filename, const Domain& dom, const std::vector<double>& result) const;
};


#endif /* WRITER_H */