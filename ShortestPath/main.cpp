/**
*   usage : pass a file name of a list of P-Center instance file names by the command line arguments.
*           or use the standard input to a list of instance files.
*           each line should contain only one file name in the list.
*           check log.csv for result.
*/


#include <iostream>
#include <fstream>
#include "Graph.h"
#include "Timer.h"
#include "RangeRand.h"
#include "PCenter.h"
#include "RandSelect.h"

using namespace std;

void readInstanceList( istream &is, vector<string> &filename )
{
    string buf;
    do {    // read all instance file names
        is >> buf;
        filename.push_back( buf );
    } while (!is.eof());
}

Graph::ArcList readPmedInstance( const string &fname, int &nodeNum, int &arcNum, int &pNum )
{
    ifstream ifs( fname );
    UndirectedGraph::ArcList arcList;
    UndirectedGraph::Arc arc;
    ifs >> nodeNum >> arcNum >> pNum;

    do {
        ifs >> arc.startVertex >> arc.endVertex >> arc.dist;
        arcList.push_back( arc );
    } while (!ifs.eof());

    ifs.close();
    return arcList;
}

int main( int argc, char **argv )
{
    const int runTime = 10;

    ofstream ofs( "log.csv", ios::app );
    //PCenter::initResultSheet( ofs ); // call if log.csv is not exist

    vector<string> filename;

    if (argc == 2) {    // use command line argument as a file contains a list of instance file names
        ifstream ifs( argv[1] );
        readInstanceList( ifs, filename );
        ifs.close();
    } else if (argc == 1) {    // use standard input to get each file name
        readInstanceList( cin, filename );
    } else {
        cout << "invalid arguments" << endl;
        return -1;
    }

    for (vector<string>::iterator iter = filename.begin(); iter != filename.end(); iter++) {
        // for each instance, run some times for judging average performance
        int nodeNum, arcNum, pNum;
        UndirectedGraph::ArcList arcList( readPmedInstance( *iter, nodeNum, arcNum, pNum ) );
        UndirectedGraph g( arcList, nodeNum, 1 );

        for (int i = 1; i <= runTime; i++) {
            PCenter pc( g, pNum, i * 2500 );

            pc.solve( 50*i, 5*i );
            //pc.basicSolve();
            pc.printResult( cout );
            if (!pc.check()) {
                ofs << "[LogicError] " << endl;
            }

            pc.appendResultToSheet( *iter, ofs );
        }
    }

    ofs.close();
    return 0;
}