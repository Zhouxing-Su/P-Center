/**
*   usage :
*           (single center is not considered here)

*   Problems:
* 1. choose update or not randomly when the distance to the new center is the same as 
*    the second nearest center, or the distance to the new center is not less than 
*    the nearest center and the distance to the first and second nearest centers are the same?

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
        TopologicalGraph::Distance serveRadius;
        TopologicalGraph::VertexSet center;
        int iterCount;
        double duration;
    };

    const unsigned pnum;

    PCenter( UndirectedGraph &ug, unsigned pnum, int maxIterCount );
    ~PCenter();

    void solve( int tabuTenureBase, int tabuTenureAmplitude );
    void tabuSolve( int tabuTenureBase, int tabuTenureAmplitude );
    void basicSolve();
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
        ClosestCenterQueue( int c1, TopologicalGraph::Distance d1, int c2, TopologicalGraph::Distance d2 )
        {
            center[0] = c1;
            center[1] = c2;
            dist[0] = d1;
            dist[1] = d2;
        }

        std::array<int, 2> center;   // index of the center
        std::array<TopologicalGraph::Distance, 2> dist;    // distance to the center
    };

    typedef std::vector<ClosestCenterQueue> ClosestCenterTable;
    typedef std::vector< std::vector<int> > TabuTable;

    class TabuTenureCalculator
    {
    public:
        TabuTenureCalculator( int ttb, int tta )
            : tabuTenureBase( ttb ), tabuTenureAmplitude( tta ), threshold( ttb * 8 ), punishment( 1 ), rr( -tta, tta )//, ofs( "TabuTenure.csv" )
        {
        }

        //~TabuTenureCalculator()
        //{
        //    ofs.close();
        //}

        const int tabuTenureBase;
        const int tabuTenureAmplitude;
        const int threshold;
        static const int attenuation = 16;

        int operator()( int iterCount, bool isBest )
        {
            if (isBest) {
                punishment /= attenuation;
            } else {
                ++punishment;   // (punishment <= 1 || punishment >= threshold)
                if (punishment < threshold) {
                    punishment *= 9;
                    punishment /= 8;
                }
            }
            int tenure = iterCount + rr() + tabuTenureBase + punishment;
            //ofs << tenure << '\n';
            return tenure;
        }

        int operator()( int iterCount )
        {
            return iterCount + rr() + tabuTenureBase;
        }

        void reset()
        {
            punishment = 1;
        }

        std::ofstream ofs;
        int punishment;
        RangeRand rr;
    };

    void genInitSolution();
    void initClosestCenter( int firstCenter, int secondCenter );

    // find one of the farthest vertices from the center set randomly.
    int findFarthestVertex( ClosestCenterTable &closestCenter ) const;         // available after initClosestCenter() is called
    // find one of the longest serve arcs randomly
    TopologicalGraph::Arc findLongestServeArc( ClosestCenterTable &closestCenter ) const; // available after initClosestCenter() is called
    // find the set of the farthest vertices from the center set.
    TopologicalGraph::VertexSet findFarthestVertices( ClosestCenterTable &closestCenter ) const;        // available after initClosestCenter() is called
    // find the set of the longest serve arcs
    TopologicalGraph::ArcSet findLongestServeArcs( ClosestCenterTable &closestCenter ) const;  // available after initClosestCenter() is called

    // update the closest center queue on each vertex (will not update the center set)
    void addCenter( int newCenter, ClosestCenterTable &closestCenter );    // available after initClosestCenter() is called
    // update the closest center queue on each vertex (will not update the center set)
    void removeCenter( int oldCenter ); // available after initClosestCenter() is called
    // select a pair of (oldCenter,newCenter)
    CenterSwap getRandSwap() const;
    // Random Remove, Greedy Add perturbation, return new minRadius
    TopologicalGraph::Distance perturbRRGA( int perturbStrength );

    UndirectedGraph graph;
    TopologicalGraph::VertexSet center;
    ClosestCenterTable closestCenter;

    TabuTable tabu;

    int maxIterCount;
    Solution bestSolution;
    Timer timer;
    std::string solvingAlgorithm;
};



#endif