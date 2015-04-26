#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <sstream>
#include "vcf.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "functions.h"

using namespace std;
using namespace BamTools;

void printUsage()
{
	cout << "Usage for Ast:" << endl;
	cout << "\t--vcf, -i\tVCF file containing genotypes for one individual" << endl
		<< "\t--bam, -b\tBAM file" << endl
		<< "\t--out-vcf\tOutput VCF file" << endl;
}

vector<string> parseString(string line)
{
	stringstream ss(line);
	vector<string> outVector;
	string outString;
	while ( ss >> outString )
		outVector.push_back(outString);
	return outVector;
}

void getRefAltCounts(BamReader& reader, const BamRegion &r, int& nReadsRef, int& nReadsAlt, vcfRecord record)
{
	// Return the number of reads aligning reference and alternate allele
	if (!reader.IsOpen())
	{
		cout << "Can't access BAM file" << endl;
		return;
	}
	if (!reader.HasIndex())
	{
		cout << "Can't access BAM file index" << endl;
		return;
	}
	
	
	BamAlignment al;
	assert(reader.SetRegion(r));
		
	while ( reader.GetNextAlignment(al) )
	{
		// int baseAlignment = al.Position;
		
		bool found = false;
		int index = 0;
		int indexPosition = al.Position;
		for ( size_t i = 0; i < al.CigarData.size(); i++ )
		{
			switch ( al.CigarData.at(i).Type )
			{
				case 'H' : break; // Hard clip
				case 'S' : break; // Soft clip
				case 'M' : // Match or mismatch
				case 'D' :
				case 'N' :
				case 'P' : 
					if ( indexPosition + al.CigarData.at(i).Length > (record.pos-1) )
					{
						index += (record.pos-1) - indexPosition;
						found = true;
					}
					else
					{
						index += al.CigarData.at(i).Length;
						indexPosition += al.CigarData.at(i).Length;
					}
					break;
				case 'I' :
					index += al.CigarData.at(i).Length;
					break;
			}
			if ( found )
				break;
		}

		if ( !found || index >= al.AlignedBases.size())
		{
			cerr << "Problem with finding variant in read:\n";
			for ( unsigned int i = 0; i < al.AlignedBases.size(); i++ )
				cout << al.AlignedBases.at(i) << " ";
			cout << endl;
			for ( unsigned int i = 0; i < al.CigarData.size(); i++ )
				cout << al.CigarData.at(i).Type << al.CigarData.at(i).Length << " ";
			cout << endl;
			cout << "Read position: " << al.RefID+1 << " " << al.Position+1 << endl;
			cout << "Variant position: " << record.chr << " " << record.pos << endl;
			exit(1);
		}
		
		if ( al.AlignedBases.at(index) == record.ref.at(0) )
			nReadsRef++;
		else if ( al.AlignedBases.at(index) == record.alt.at(0) )
			nReadsAlt++;
	}
}

std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}