/**
*   usage :
*           (single center is not considered here)

(example code)
====================================================================
#include <iostream>
#include "Graph.h"

using namespace std;

int main( int argc, char **argv )
{
return 0;
}
====================================================================
*/
#ifndef PCENTER_H
#define PCENTER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <ctime>
#include "Graph.h"
#include "RangeRand.h"
#include "RandSelect.h"
#include "Timer.h"

class PCenter
{
public:
    struct Solution
    {
        Graph::Distance serveRadius;
        Graph::VertexSet center;
        int iterCount;
        double duration;
    };

    const int pnum;

    PCenter( UndirectedGraph &ug, int pnum, int maxIterCount );
    ~PCenter();

    void solve( int tabuTenureBase, int tabuTenureAmplitude );
    void BasicSolve();
    void greedyBasicSolve();
    bool check() const; // check the result by shortestDist

    void printResult( std::ostream &os ) const;
    static void initResultSheet( std::ofstream &csvFile );
    void appendResultToSheet( const std::string &instanceFileName, std::ofstream &csvFile ) const;

private:
    class CenterSwap
    {
    public:
        CenterSwap() {}
        CenterSwap( int oc, int nc ) : oldCenter( oc ), newCenter( nc ) {}

        int oldCenter;
        int newCenter;
    };

    class ClosestCenterQueue    // each vertex will have one
    {
    public:
        ClosestCenterQueue() {}
        ClosestCenterQueue( int c1, Graph::Distance d1, int c2, Graph::Distance d2 )
        {
            center[0] = c1;
            center[1] = c2;
            dist[0] = d1;
            dist[1] = d2;
        }

        std::array<int, 2> center;   // index of the center
        std::array<Graph::Distance, 2> dist;    // distance to the center
    };

    typedef std::vector<ClosestCenterQueue> ClosestCenterTable;
    typedef std::vector< std::vector<int> > TabuTable;

    class TabuTenureCalculator
    {
    public:
        TabuTenureCalculator( int ttb, int tta )
            : tabuTenureBase( ttb ), tabuTenureAmplitude( tta ), rr( -tta, tta )
        {
        }

        int operator()( int iterCount )
        {
            return iterCount + rr() + tabuTenureBase;
        }

        const int tabuTenureBase;
        const int tabuTenureAmplitude;
        RangeRand rr;
    };

    void genInitSolution();
    void initClosestCenter( int firstCenter, int secondCenter );

    // find one of the farthest vertices from the center set randomly.
    int findFarthestVertex( ClosestCenterTable &closestCenter ) const;         // available after initClosestCenter() is called
    // find one of the longest serve arcs randomly
    Graph::Arc findLongestServeArc( ClosestCenterTable &closestCenter ) const; // available after initClosestCenter() is called
    // find the set of the farthest vertices from the center set.
    Graph::VertexSet findFarthestVertices( ClosestCenterTable &closestCenter ) const;        // available after initClosestCenter() is called
    // find the set of the longest serve arcs
    Graph::ArcSet findLongestServeArcs( ClosestCenterTable &closestCenter ) const;  // available after initClosestCenter() is called

    // update the closest center queue on each vertex (will not update the center set)
    void addCenter( int newCenter, ClosestCenterTable &closestCenter );    // available after initClosestCenter() is called
    // update the closest center queue on each vertex (will not update the center set)
    void removeCenter( int oldCenter ); // available after initClosestCenter() is called

    UndirectedGraph graph;
    Graph::VertexSet center;
    ClosestCenterTable closestCenter;

    TabuTable tabu;

    int maxIterCount;
    Solution bestSolution;
    Timer timer;
    std::string solvingAlgorithm;
};



#endif