#ifndef FUNCTIONS_INCLUDED
#define FUNCTIONS_INCLUDED

#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "vcf.h"

void printUsage();

std::vector<std::string> parseString(std::string line);

void getRefAltCounts(BamTools::BamReader& reader, const BamTools::BamRegion &r, int& nReadsRef, int& nReadsAlt, vcfRecord record);

std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems);

#endif // FUNCTIONS_INCLUDED