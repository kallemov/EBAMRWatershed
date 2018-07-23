#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <zlib.h>

using namespace std;

int
main (int arc, char**argv)
{
    if (arc<2) 
    {
	cout<<"Please enter file name of dem input file to make gzipped arc"<<endl;
	exit(1);
    }
 
    char* filename1= argv[1];
    ifstream DEM_file;
    DEM_file.open(filename1); //the input file
    if (!DEM_file.good())
    {
	cout<<"Bad DEM_file in readDEM"<<endl;
    }
    string line;
    string zipname(filename1);
    zipname += ".gz";
    gzFile outfile = gzopen(zipname.c_str(), "wb");

    char astring[1024];
    int ncols, nrows;
    double hx,hy;

    DEM_file>>astring;
    DEM_file>>ncols;
    gzwrite(outfile, (char*) &ncols, sizeof(int));
    DEM_file>>nrows;
    gzwrite(outfile, (char*) &nrows, sizeof(int));

    DEM_file>>astring;
    DEM_file>>hx;
    gzwrite(outfile, (char*) &hx, sizeof(double));
    DEM_file>>hy;
    gzwrite(outfile, (char*) &hy, sizeof(double));
    double **M = new double *[ncols];
    for (int i = 0; i < ncols;i++)
    {
	M[i] = new double[nrows];
    }
    for (int j=0;j<nrows;j++)
    {
	for (int i=0;i<ncols;i++)
	{
	    
	    DEM_file>>M[i][j];    
	}
    }

    for (int i=0;i<ncols;i++)
    {
    	gzwrite(outfile, (char*) M[i], sizeof(double)*nrows);
    }
    cout<<"Number of elements="<<sizeof(double)*nrows*ncols<<endl;
    gzclose(outfile);
    
    return 0;
}


