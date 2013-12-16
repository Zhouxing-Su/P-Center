#include <iostream>
#include <fstream>
#include "Graph.h"
#include "Timer.h"
#include "RangeRand.h"
#include "PCenter.h"
#include "RandSelect.h"

using namespace std;

Graph::ArcList readPmedInstance( const char *fname, int &nodeNum, int &arcNum, int &pNum )
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
    ofstream ofs( "log.csv" );
    PCenter::initResultSheet( ofs );

    for (int i = 1; i <= argc; i++) {
        int nodeNum, arcNum, pNum;
        UndirectedGraph::ArcList arcList( readPmedInstance( argv[i], nodeNum, arcNum, pNum ) );

        UndirectedGraph g( arcList, nodeNum, 1 );
        PCenter pc( g, pNum );

        pc.solve( 10000000 );
        if (!pc.check()) {

            return -1;
        }

        pc.appendResultToSheet( argv[i], ofs );
        //pc.printResult( cout );
    }

    ofs.close();
    return 0;
}