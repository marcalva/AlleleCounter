#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <sstream>
#include "functions.h"
#include "vcf.h"

using namespace std;

/*
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
*/

vcfFileInStream::vcfFileInStream(std::string fn)
{
	ifs.open(fn.c_str()); 
}

vcfRecord vcfFileInStream::readRecord()
{
	vcfRecord r;
	
	if (!ifs.good())
		return r;
	
	char line_array[100000];
	ifs.getline(line_array, 100000, '\n');
	
	std::string line = line_array;
	
	std::vector<std::string> line_vector;
	
	// MAKE SURE THAT VCF RECORD IS TAB-DELIMITED
	split(line, '\t', line_vector);
	
	if (line_vector.size() < 10)
	{
		std::cout << "There are less than 10 columns in this VCF record. Please correct file format." << std::endl;
		std::cout << "Number of columns: " << line_vector.size() << endl;
		for (int i = 0 ; i < line_vector.size(); i++)
			cout << line_vector.at(i) << "\t";
		cout << endl;
		std::exit(1);
	}
	
	r.chr = std::stoi(line_vector.at(0));
	r.pos = std::stoi(line_vector.at(1));
	r.id = line_vector.at(2);
	r.ref = line_vector.at(3);
	r.alt = line_vector.at(4);
	r.qual = std::stoi(line_vector.at(5));
	r.filter = line_vector.at(6);
	
	std::vector<std::string> iField;
	split(line_vector.at(7), ';', iField);
	
	if (! ( line_vector.at(7) == "" || line_vector.at(7) == "." ) ) // If the info field is not empty
	{
		for (unsigned long i = 0; i < iField.size(); i++)
		{
			std::vector<std::string> v;
			split(iField.at(i), '=', v);
			
			// If this info field is a flag
			if ( v.size() == 1 )
			{
				r.info.key.push_back(v.at(0));
				r.info.value.push_back("");
			}
			else // If the info field is a key value pair
			{
				r.info.key.push_back(v.at(0));
				r.info.value.push_back(v.at(1));
			}
		}
	}
	
	split(line_vector.at(8), ':', r.format.value);
	
	for (int i = 9; i < line_vector.size(); i++)
	{
		sampleField s;
		split(line_vector.at(i), ':', s.value);
		r.samples.push_back(s);
	}

	return r;
	
}

vcfHeader vcfFileInStream::readHeader()
{
	vcfHeader h;
	
	if (!ifs.good())
		return h;
	
	// Go to beginning of file
	ifs.seekg(0, ifs.beg);
	
	char line_array[100000];
	std::string line;
	do
	{
		ifs.getline(line_array, 100000, '\n');
		
		line = line_array;
		if (line.substr(0, 2) == "##")
			h.meta.push_back(line);
		
	} while (line.substr(0, 2) == "##");
	
	assert(line.substr(0, 6) == "#CHROM");
	
	h.head = line;
	
	return h;
	
}

bool vcfFileInStream::atEnd()
{
	ifs.peek();
	return ifs.eof();
}

void vcfFileInStream::close()
{
	ifs.close();
}

bool vcfFileInStream::good()
{
	return ifs.good();
};

vcfFileOutStream::vcfFileOutStream(std::string fn)
{
	ofs.open(fn.c_str());
}

void vcfFileOutStream::writeRecord(vcfRecord v)
{
	if (!ofs.good())
	{
		std::cout << "Attempting to write to non-good output file." << std::endl;
		std::exit(1);
	}
	
	std::string output_line = std::to_string(v.chr) + "\t" + \
		std::to_string(v.pos) + "\t" + \
			v.id + "\t" + v.ref + "\t" + v.alt + "\t" + \
				std::to_string(v.qual) + "\t" + v.filter + "\t";

	if (v.info.key.size() > 0)
	{
		for (int i = 0; i < v.info.key.size(); i++)
		{
			if ( v.info.value.at(i) == "" ) // If key is a flag
				output_line += v.info.key.at(i);
			else // If key-value pair
				output_line += v.info.key.at(i) + "=" + output_line += v.info.value.at(i);
			
			// Don't place ; delimiter if at the last value.
			if (i < v.info.key.size() - 1)
				output_line += ";";
		}
	}
	else
		output_line += ".";
		
	output_line += "\t";
	
	for (int i = 0; i < v.format.value.size(); i++)
	{
		output_line += v.format.value.at(i);
		
		if (i < v.format.value.size() - 1)
			output_line += ":";
	}
	
	for (int i = 0; i < v.samples.size(); i++)
	{
		output_line += "\t";
		
		for (int j = 0; j < v.samples.at(i).value.size(); j++)
		{
			output_line += v.samples.at(i).value.at(j);
			
			if (j < v.samples.at(i).value.size() - 1)
				output_line += ":";
		}
	}
		
	output_line += "\n";
	
	ofs << output_line;
	
}

void vcfFileOutStream::writeHeader(vcfHeader h)
{
	if (!ofs.good())
	{
		std::cout << "Attempting to write to non-good output file." << std::endl;
		std::exit(1);
	}
	
	ofs.seekp(0, ofs.beg);
	
	for (int i = 0; i < h.meta.size(); i++)
	{
		ofs << h.meta.at(i) << "\n";
	}
	
	ofs << h.head << "\n";
}

void vcfFileOutStream::close()
{
	ofs.close();
}

bool vcfFileOutStream::good()
{
	return ofs.good();
}