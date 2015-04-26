CC = g++
FLAGS = -v -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -std=c++11 -stdlib=libc++
LD_LIBS = -L ThirdParty/lib -I include -I ThirdParty/include \
		-lbamtools -lz

AlleleCount : src/main.cpp
	$(CC) $(FLAGS) src/main.cpp src/functions.cpp src/vcf.cpp -o $@ $(LD_LIBS)
