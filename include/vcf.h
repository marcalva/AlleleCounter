#ifndef VCF_INCLUDED
#define VCF_INCLUDED

#include <string>
#include <vector>

struct infoField{
	// Each component is separated by a semicolon in the VCF format
	std::vector<std::string> key;
	std::vector<std::string> value;
};

struct formatField{
	// Each component is separated by a colon in the VCF format
	std::vector<std::string> value;
};

struct sampleField{
	// Each component is separated by a colon in the VCF format
	// Must have the same length as formatField.
	std::vector<std::string> value;
};

struct vcfHeader
{
	// Meta lines that start with ## 
	std::vector<std::string> meta;
	
	// Header lne beginning with #CHROM
	std::string head;
};

struct vcfRecord{
	int chr;
	int pos;
	std::string id;
	std::string ref;
	std::string alt;
	int qual;
	std::string filter;
	infoField info;
	formatField format;
	std::vector<sampleField> samples;
};

class vcfFileInStream{
public:
	std::string file_name;
	
	vcfFileInStream(std::string fn);
	vcfRecord readRecord();
	vcfHeader readHeader();
	bool atEnd();
	void close();
	bool good();

private:
	std::ifstream ifs;
};

class vcfFileOutStream{
public:
	std::string file_name;
	
	vcfFileOutStream(std::string fn);
	void writeRecord(vcfRecord v);
	void writeHeader(vcfHeader h);
	void close();
	bool good();
	
private:
	std::ofstream ofs;
};

#endif // VCF_INCLUDED