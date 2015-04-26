#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <sstream>
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "functions.h"
#include "vcf.h"

using namespace std;
using namespace BamTools;

void printUsage();
vector<string> parseString(string line);

int main(int argc, char* argv[])
{
	
	// Read in arguments and variables
	if ( argc < 6 )
	{
		cerr << "You must give the variant, bam, and output files" << endl;
		printUsage();
		exit(1);
	}
	
	string  bamFileName, vcfFileName, vcfOutName;
	int k = 1;
	while ( k < argc )
	{
		if ( strcmp(argv[k], "--vcf") == 0 || strcmp(argv[k], "-i") == 0 )
			vcfFileName = argv[k+1];
		else if ( strcmp(argv[k], "--bam") == 0 || strcmp(argv[k], "-b") == 0 )
			bamFileName = argv[k+1];
		else if ( strcmp(argv[k], "--out-vcf") == 0 )
			vcfOutName = argv[k+1];
		else
		{
			cout << "Unkown argument " << argv[k] << endl;
			printUsage();
			exit(1);
		}
		k += 2;
	}

	
	// Open the VCF or variant input file
	cout << "Setting up input and output VCF files... ";
	
	vcfFileInStream vcfIn(vcfFileName);
	vcfFileOutStream vcfOut(vcfOutName);
	
	vcfHeader header = vcfIn.readHeader();
	vcfOut.writeHeader(header);
	cout << "Done" << endl;
	
	
	
	// Open Bam file
	cout << "Opening BAM file... ";
	BamReader reader;
	assert(reader.Open(bamFileName));
	assert(reader.OpenIndex(bamFileName + ".bai"));
	cout << "Done" << endl;
	
	

	// Loop through every variant
	cout << "Looping through variants... ";
	while (true)
	{
		// Store variant information in VcfRecord format		
		if (vcfIn.atEnd())
			break;
		
		vcfRecord record = vcfIn.readRecord();		
		
		
		// Get genotype information
		int a1, a2; // Hold genotype of each allele.
		size_t split;
		string geno_info = record.samples.at(0).value.at(0);
		if ( geno_info.find("|") != string::npos )
			split = geno_info.find_first_of('|');
		else
			split = geno_info.find_first_of('/');

		a1 = round(atof(geno_info.substr(0,split).c_str()));
		a2 = round(atof(geno_info.substr(split+1,  geno_info.find_first_of(':') - split - 1).c_str()));

		record.format.value.push_back("AC");
		
		// If this person is homozygous
		if ( a1 ==  a2 )
		{	
			record.samples.at(0).value.push_back("NA,NA");
			vcfOut.writeRecord(record);
			continue; // Move to next variant
		}
				
		int nReadsRef = 0;
		int nReadsAlt = 0;
		
		// 0 based, half open interval
		BamRegion region(record.chr-1, record.pos-1, record.chr-1, record.pos);
				
		getRefAltCounts(reader, region, nReadsRef, nReadsAlt, record);
				
		stringstream sRef;
		sRef << nReadsRef;
		stringstream sAlt;
		sAlt << nReadsAlt;
		
		string app;
		if ( geno_info.find('|') != string::npos && a1 > a2)
		{
			app = sAlt.str() + "," + sRef.str();
		}
		else
		{
			app = sRef.str() + "," + sAlt.str();
		}
		
		record.samples.at(0).value.push_back(app);
		
		vcfOut.writeRecord(record);

	}
	
	cout << "Finished" << endl;
		
	vcfIn.close();
	vcfOut.close();
	reader.Close();

	return 0;
}


