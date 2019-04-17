all: leath.cc
	g++ -std=c++14 leath.cc -o leath

clean:
	%(RM) leath
