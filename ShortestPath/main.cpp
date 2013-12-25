/**
*   usage : pass a file name of a list of P-Center instance file names by the command line arguments.
*           or use the standard input to a list of instance files.
*           each line should contain only one file name in the list.
*           check log.csv for result.
*/


#include <iostream>
#include <fstream>
#include <cfloat>
#include "Graph.h"
#include "Timer.h"
#include "RangeRand.h"
#include "PCenter.h"
#include "RandSelect.h"
#include "Double.h"

using namespace std;

void readInstanceList( istream &is, vector<string> &filename )
{
    string buf;
    do {    // read all instance file names
        is >> buf;
        filename.push_back( buf );
    } while (!is.eof());
}

TopologicalGraph<>::ArcList readPmedInstance( const string &fname, unsigned &nodeNum, unsigned &arcNum, unsigned &pNum )
{
    ifstream ifs( fname );
    UndirectedGraph<>::ArcList arcList;
    UndirectedGraph<>::Arc arc;
    ifs >> nodeNum >> arcNum >> pNum;

    do {
        ifs >> arc.startVertex >> arc.endVertex >> arc.dist;
        arcList.push_back( arc );
    } while (!ifs.eof());

    ifs.close();
    return arcList;
}

int solve_pmed( int argc, char **argv, ofstream &csvFile )
{
    const int runTime = 8;

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
        unsigned nodeNum, arcNum, pNum;
        UndirectedGraph<unsigned>::ArcList arcList( readPmedInstance( *iter, nodeNum, arcNum, pNum ) );
        UndirectedGraph<unsigned> g( arcList, nodeNum, 1 );

        for (int i = 1; i <= runTime; i++) {
            PCenter<unsigned> pc( g, pNum, i * 20 * nodeNum );

            pc.solve( i*nodeNum / 2, nodeNum / 2 );

            pc.printResult( cout );
            if (!pc.check()) {
                csvFile << "[LogicError] ";
            }

            pc.appendResultToSheet( *iter, csvFile );
        }
    }
}


GeometricalGraph::PointList read_rl1323Instance( const string &fname )
{
    ifstream ifs( fname );

    GeometricalGraph::Coord x, y;
    GeometricalGraph::PointList pl;

    do {
        ifs >> x >> y;
        pl.push_back( GeometricalGraph::Point( x, y ) );
    } while (!ifs.eof());

    ifs.close();
    return pl;
}

int solve_rl1323( ofstream &csvFile )
{
    const string fname( "rl1323.tsp" );
    const int nodeNum = 1323;

    const int runTime = 4;

    GeometricalGraph gg( read_rl1323Instance( fname ) );
    UndirectedGraph<Double> dug( gg );

    for (int pnum = 10; pnum <= 100; pnum += 10) {
        // for each instance, run some times for judging average performance

        for (int i = 1; i <= runTime; i++) {
            PCenter<Double> pc( dug, pnum, i * 500000 );

            pc.solve( i*nodeNum / 2, nodeNum / 2 );

            pc.printResult( cout );
            if (!pc.check()) {
                csvFile << "[LogicError] ";
            }

            ostringstream fn( fname );
            fn << '(' << pnum << ')';
            pc.appendResultToSheet( fn.str(), csvFile );
        }
    }

    return -1;
}

int main( int argc, char **argv )
{

    ofstream ofs( "log.csv", ios::app );
    //PCenter::initResultSheet( ofs ); // call if log.csv is not exist

    //solve_pmed( argc, argv, ofs );
    solve_rl1323( ofs );

    ofs.close();
    return 0;
}