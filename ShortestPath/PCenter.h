/**
*   usage :
*           (single center is not considered here)

*   Problems:
* 1. choose update or not randomly when the distance to the new center is the same as
*    the second nearest center, or the distance to the new center is not less than
*    the nearest center and the distance to the first and second nearest centers are the same?
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

template <typename T_DIST>
class PCenter
{
public:
    struct Solution
    {
        typename TopologicalGraph<T_DIST>::Distance serveRadius;
        typename TopologicalGraph<T_DIST>::VertexSet center;
        int iterCount;
        double duration;
    };

    const unsigned pnum;

    PCenter( UndirectedGraph<T_DIST> &ug, unsigned pnum, int maxIterCount );
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
        ClosestCenterQueue( int c1, typename TopologicalGraph<T_DIST>::Distance d1, int c2, typename TopologicalGraph<T_DIST>::Distance d2 )
        {
            center[0] = c1;
            center[1] = c2;
            dist[0] = d1;
            dist[1] = d2;
        }

        std::array<int, 2> center;   // index of the center
        std::array<typename TopologicalGraph<T_DIST>::Distance, 2> dist;    // distance to the center
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
    typename TopologicalGraph<T_DIST>::Arc findLongestServeArc( ClosestCenterTable &closestCenter ) const; // available after initClosestCenter() is called
    // find the set of the farthest vertices from the center set.
    typename TopologicalGraph<T_DIST>::VertexSet findFarthestVertices( ClosestCenterTable &closestCenter ) const;        // available after initClosestCenter() is called
    // find the set of the longest serve arcs
    typename TopologicalGraph<T_DIST>::ArcSet findLongestServeArcs( ClosestCenterTable &closestCenter ) const;  // available after initClosestCenter() is called

    // update the closest center queue on each vertex (will not update the center set)
    void addCenter( int newCenter, ClosestCenterTable &closestCenter );    // available after initClosestCenter() is called
    // update the closest center queue on each vertex (will not update the center set)
    void removeCenter( int oldCenter ); // available after initClosestCenter() is called
    // select a pair of (oldCenter,newCenter)
    CenterSwap getRandSwap() const;
    // Random Remove, Greedy Add perturbation, return new minRadius
    typename TopologicalGraph<T_DIST>::Distance perturbRRGA( int perturbStrength );

    UndirectedGraph<T_DIST> graph;
    typename TopologicalGraph<T_DIST>::VertexSet center;
    ClosestCenterTable closestCenter;

    TabuTable tabu;

    int maxIterCount;
    Solution bestSolution;
    Timer timer;
    std::string solvingAlgorithm;
};












using namespace std;

template <typename T_DIST>
PCenter<T_DIST>::PCenter( UndirectedGraph<T_DIST> &ug, unsigned pn, int mic )
: pnum( pn ), graph( ug ), closestCenter( ug.vertexAllocNum, ClosestCenterQueue() ),
tabu( ug.vertexAllocNum, vector<int>( ug.vertexAllocNum, 0 ) ), maxIterCount( mic )
{
    graph.getDistSeqTable();
}


template <typename T_DIST>
PCenter<T_DIST>::~PCenter()
{
}

template <typename T_DIST>
void PCenter<T_DIST>::solve( int tabuTenureBase, int tabuTenureAmplitude )
{
    ostringstream ss;
    ss << "perturb(RRGA)+tabu search(B=" << tabuTenureBase << "&A=" << tabuTenureAmplitude << ')';
    solvingAlgorithm = ss.str();

    TabuTenureCalculator getTabuTenure( tabuTenureBase, tabuTenureAmplitude );
    unsigned noImproveCount = 0;

    genInitSolution();

    RandSelect rs( 2 );
    for (int iterCount = 0; iterCount < maxIterCount; iterCount++) {
        bool isSwaped = false;
        CenterSwap centerSwap;
        TopologicalGraph<T_DIST>::Distance minRadius = TopologicalGraph<T_DIST>::MAX_DISTANCE;

        TopologicalGraph<T_DIST>::Arc longestServeArc = findLongestServeArc( closestCenter );

        int longestEnd = longestServeArc.endVertex;
        TopologicalGraph<T_DIST>::Distance longestDist = longestServeArc.dist;
        // try each vertex whose distance to longestEnd is shorter than longestDist
        for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
            int newCenter = graph.nthClosestVertex( longestEnd, i );
            if (graph.distance( longestEnd, newCenter ) < longestDist) {
                // find the best swap between center i and non-center vertices
                ClosestCenterTable tmpCCT( closestCenter );
                addCenter( newCenter, tmpCCT );
                TopologicalGraph<T_DIST>::Distance radiusAfterAdd = tmpCCT[findFarthestVertex( tmpCCT )].dist[0];
                // calculate new radius for removing each center except the newly added one
                for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = center.begin(); iter != center.end(); iter++) {
                    // when *iter is removed
                    int removedCenter = *iter;
                    TopologicalGraph<T_DIST>::Distance radiusAfterRemove = radiusAfterAdd;
                    for (int k = graph.minVertexIndex; k <= graph.maxVertexIndex; k++) {
                        if (tmpCCT[k].center[0] == removedCenter) {
                            TopologicalGraph<T_DIST>::Distance newDist = tmpCCT[k].dist[1];
                            if (radiusAfterRemove < newDist) {
                                radiusAfterRemove = newDist;
                            }
                        }
                    }
                    // check if the swap between the candidate and the old is better
                    if (radiusAfterRemove < minRadius) {
                        if (radiusAfterRemove < bestSolution.serveRadius
                            || iterCount > tabu[removedCenter][newCenter]) {
                            centerSwap = CenterSwap( removedCenter, newCenter );
                            minRadius = radiusAfterRemove;
                            rs.reset( 2 );
                            isSwaped = true;
                        }
                    } else if (radiusAfterRemove == minRadius) {
                        if (rs.isSelected()) {
                            if (radiusAfterRemove < bestSolution.serveRadius
                                || iterCount > tabu[removedCenter][newCenter]) {
                                centerSwap = CenterSwap( removedCenter, newCenter );
                                isSwaped = true;
                            }
                        }
                    }
                }
            } else {
                break;
            }
        }

        if (!isSwaped || noImproveCount > (4 * graph.vertexNum)) {    // "Random Remove, Greedy Add" perturbation
            minRadius = perturbRRGA( pnum / 4 );
            getTabuTenure.reset();
            noImproveCount = 0;
        } else {    // commit the swap
            center.erase( centerSwap.oldCenter );
            center.insert( centerSwap.newCenter );
            addCenter( centerSwap.newCenter, closestCenter );
            removeCenter( centerSwap.oldCenter );
        }
        // record if it is the best solution
        if (minRadius < bestSolution.serveRadius) {
            timer.record();
            bestSolution.duration = timer.getDuration();
            bestSolution.iterCount = iterCount;
            bestSolution.serveRadius = minRadius;
            bestSolution.center = center;
            noImproveCount = 0;
        } else {
            ++noImproveCount;
        }
        // update tabu
        if (isSwaped) {
            tabu[centerSwap.oldCenter][centerSwap.newCenter] = getTabuTenure( iterCount, (noImproveCount == 0) );
        }
    }
}

template <typename T_DIST>
void PCenter<T_DIST>::tabuSolve( int tabuTenureBase, int tabuTenureAmplitude )
{
    ostringstream ss;
    ss << "tabu search(B=" << tabuTenureBase << "&A=" << tabuTenureAmplitude << ')';
    solvingAlgorithm = ss.str();

    TabuTenureCalculator getTabuTenure( tabuTenureBase, tabuTenureAmplitude );
    genInitSolution();

    RandSelect rs( 2 );
    for (int iterCount = 0; iterCount < maxIterCount; iterCount++) {
        bool isSwaped = false;
        CenterSwap centerSwap;
        TopologicalGraph<T_DIST>::Distance minRadius = TopologicalGraph<T_DIST>::MAX_DISTANCE;

        TopologicalGraph<T_DIST>::Arc longestServeArc = findLongestServeArc( closestCenter );

        int longestEnd = longestServeArc.endVertex;
        TopologicalGraph<T_DIST>::Distance longestDist = longestServeArc.dist;
        // try each vertex whose distance to longestEnd is shorter than longestDist
        for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
            int newCenter = graph.nthClosestVertex( longestEnd, i );
            if (graph.distance( longestEnd, newCenter ) < longestDist) {
                // find the best swap between center i and non-center vertices
                ClosestCenterTable tmpCCT( closestCenter );
                addCenter( newCenter, tmpCCT );
                TopologicalGraph<T_DIST>::Distance radiusAfterAdd = tmpCCT[findFarthestVertex( tmpCCT )].dist[0];
                // calculate new radius for removing each center except the newly added one
                for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = center.begin(); iter != center.end(); iter++) {
                    // when *iter is removed
                    int removedCenter = *iter;
                    TopologicalGraph<T_DIST>::Distance radiusAfterRemove = radiusAfterAdd;
                    for (int k = graph.minVertexIndex; k <= graph.maxVertexIndex; k++) {
                        if (tmpCCT[k].center[0] == removedCenter) {
                            TopologicalGraph<T_DIST>::Distance newDist = tmpCCT[k].dist[1];
                            if (radiusAfterRemove < newDist) {
                                radiusAfterRemove = newDist;
                            }
                        }
                    }
                    // check if the swap between the candidate and the old is better
                    if (radiusAfterRemove < minRadius) {
                        if (radiusAfterRemove < bestSolution.serveRadius
                            || iterCount > tabu[removedCenter][newCenter]) {
                            centerSwap = CenterSwap( removedCenter, newCenter );
                            minRadius = radiusAfterRemove;
                            rs.reset( 2 );
                            isSwaped = true;
                        }
                    } else if (radiusAfterRemove == minRadius) {
                        if (rs.isSelected()) {
                            if (radiusAfterRemove < bestSolution.serveRadius
                                || iterCount > tabu[removedCenter][newCenter]) {
                                centerSwap = CenterSwap( removedCenter, newCenter );
                                isSwaped = true;
                            }
                        }
                    }
                }
            } else {
                break;
            }
        }

        if (!isSwaped) {    // do strong perturbation
            centerSwap = getRandSwap();
            //do strong perturbation !!!
        }
        // commit the swap
        center.erase( centerSwap.oldCenter );
        center.insert( centerSwap.newCenter );
        addCenter( centerSwap.newCenter, closestCenter );
        removeCenter( centerSwap.oldCenter );
        if (!isSwaped) {    // update minRadius
            minRadius = closestCenter[findFarthestVertex( closestCenter )].dist[0];
        }
        // record if it is the best solution
        if (minRadius < bestSolution.serveRadius) {
            timer.record();
            bestSolution.duration = timer.getDuration();
            bestSolution.iterCount = iterCount;
            bestSolution.serveRadius = minRadius;
            bestSolution.center = center;
            tabu[centerSwap.oldCenter][centerSwap.newCenter] = getTabuTenure( iterCount, true );
        } else {
            tabu[centerSwap.oldCenter][centerSwap.newCenter] = getTabuTenure( iterCount, false );
        }
    }
}

template <typename T_DIST>
void PCenter<T_DIST>::basicSolve()
{
    solvingAlgorithm = "basic local search";
    genInitSolution();

    RandSelect rs( 2 );
    for (int iterCount = 0; iterCount < maxIterCount; iterCount++) {
        CenterSwap centerSwap;
        TopologicalGraph<T_DIST>::Distance minRadius = TopologicalGraph<T_DIST>::MAX_DISTANCE;

        TopologicalGraph<T_DIST>::Arc longestServeArc = findLongestServeArc( closestCenter );

        int longestEnd = longestServeArc.endVertex;
        TopologicalGraph<T_DIST>::Distance longestDist = longestServeArc.dist;
        // try each vertex whose distance to longestEnd is shorter than longestDist
        for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
            int newCenter = graph.nthClosestVertex( longestEnd, i );
            if (graph.distance( longestEnd, newCenter ) < longestDist) {
                // find the best swap between center i and non-center vertices
                ClosestCenterTable tmpCCT( closestCenter );
                addCenter( newCenter, tmpCCT );
                TopologicalGraph<T_DIST>::Distance radiusAfterAdd = tmpCCT[findFarthestVertex( tmpCCT )].dist[0];
                // calculate new radius for removing each center except the newly added one
                for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = center.begin(); iter != center.end(); iter++) {
                    // when *iter is removed
                    int removedCenter = *iter;
                    TopologicalGraph<T_DIST>::Distance radiusAfterRemove = radiusAfterAdd;
                    for (int k = graph.minVertexIndex; k <= graph.maxVertexIndex; k++) {
                        if (tmpCCT[k].center[0] == removedCenter) {
                            TopologicalGraph<T_DIST>::Distance newDist = tmpCCT[k].dist[1];
                            if (radiusAfterRemove < newDist) {
                                radiusAfterRemove = newDist;
                            }
                        }
                    }
                    // check if the swap between the candidate and the old is better
                    if (radiusAfterRemove < minRadius) {
                        centerSwap = CenterSwap( removedCenter, newCenter );
                        minRadius = radiusAfterRemove;
                        rs.reset( 2 );
                    } else if (radiusAfterRemove == minRadius) {
                        if (rs.isSelected()) {
                            centerSwap = CenterSwap( removedCenter, newCenter );
                        }
                    }
                }
            } else {
                break;
            }
        }
        // commit the swap
        center.erase( centerSwap.oldCenter );
        center.insert( centerSwap.newCenter );
        addCenter( centerSwap.newCenter, closestCenter );
        removeCenter( centerSwap.oldCenter );
        // record if it is the best solution
        if (minRadius < bestSolution.serveRadius) {
            timer.record();
            bestSolution.duration = timer.getDuration();
            bestSolution.iterCount = iterCount;
            bestSolution.serveRadius = minRadius;
            bestSolution.center = center;
        }
    }
}

template <typename T_DIST>
void PCenter<T_DIST>::greedyBasicSolve()
{
    solvingAlgorithm = "greedy local search";
    genInitSolution();

    RandSelect rs( 2 );
    for (int iterCount = 0; iterCount < maxIterCount; iterCount++) {
        CenterSwap centerSwap;
        TopologicalGraph<T_DIST>::Distance minRadius = TopologicalGraph<T_DIST>::MAX_DISTANCE;

        TopologicalGraph<T_DIST>::ArcSet lsa = findLongestServeArcs( closestCenter );

        for (TopologicalGraph<T_DIST>::ArcSet::iterator iter = lsa.begin(); iter != lsa.end(); iter++) {
            int longestEnd = iter->endVertex;
            TopologicalGraph<T_DIST>::Distance longestDist = iter->dist;
            // try each vertex whose distance to longestEnd is shorter than longestDist
            for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
                int newCenter = graph.nthClosestVertex( longestEnd, i );
                if (graph.distance( longestEnd, newCenter ) < longestDist) {
                    // find the best swap between center i and non-center vertices
                    ClosestCenterTable tmpCCT( closestCenter );
                    addCenter( newCenter, tmpCCT );
                    TopologicalGraph<T_DIST>::Distance radiusAfterAdd = tmpCCT[findFarthestVertex( tmpCCT )].dist[0];
                    // calculate new radius for removing each center
                    for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = center.begin(); iter != center.end(); iter++) {
                        // when *iter is removed
                        int removedCenter = *iter;
                        TopologicalGraph<T_DIST>::Distance radiusAfterRemove = radiusAfterAdd;
                        for (int k = graph.minVertexIndex; k <= graph.maxVertexIndex; k++) {
                            if (tmpCCT[k].center[0] == removedCenter) {
                                TopologicalGraph<T_DIST>::Distance newDist = tmpCCT[k].dist[1];
                                if (radiusAfterRemove < newDist) {
                                    radiusAfterRemove = newDist;
                                }
                            }
                        }
                        // check if the swap between the candidate and the old is better
                        if (radiusAfterRemove < minRadius) {
                            centerSwap = CenterSwap( removedCenter, newCenter );
                            minRadius = radiusAfterRemove;
                            rs.reset( 2 );
                        } else if (radiusAfterRemove == minRadius) {
                            if (rs.isSelected()) {
                                centerSwap = CenterSwap( removedCenter, newCenter );
                            }
                        }
                    }
                } else {
                    break;
                }
            }
        }

        // commit the swap
        center.erase( centerSwap.oldCenter );
        center.insert( centerSwap.newCenter );
        addCenter( centerSwap.newCenter, closestCenter );   // add first may pop the being removed center out of the
        removeCenter( centerSwap.oldCenter );   // closestCenter queue which will left less work to removeCenter()
        // record if it is the best solution
        if (minRadius < bestSolution.serveRadius) {
            timer.record();
            bestSolution.duration = timer.getDuration();
            bestSolution.iterCount = iterCount;
            bestSolution.serveRadius = minRadius;
            bestSolution.center = center;
        }
    }
}

template <typename T_DIST>
bool PCenter<T_DIST>::check() const
{
    if (bestSolution.center.size() != pnum) {
        return false;
    }

    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        TopologicalGraph<T_DIST>::Distance minRadius = TopologicalGraph<T_DIST>::MAX_DISTANCE;
        for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = bestSolution.center.begin(); iter != bestSolution.center.end(); iter++) {
            if (minRadius > graph.distance( i, *iter )) {
                minRadius = graph.distance( i, *iter );
            }
        }

        if (minRadius > bestSolution.serveRadius) {
            return false;
        }
    }

    return true;
}

template <typename T_DIST>
void PCenter<T_DIST>::printResult( ostream &os ) const
{
    os << "The max serving radius is : " << bestSolution.serveRadius << endl;
    os << "The indexes of the vertices which are chosed as centers are :\n";
    for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = bestSolution.center.begin(); iter != bestSolution.center.end(); iter++) {
        os << *iter << "|";
    }
    os << endl;
}

template <typename T_DIST>
void PCenter<T_DIST>::initResultSheet( std::ofstream &csvFile )
{
    csvFile << "Date, Instance, Algorithm, TotalIter, Duration, IterCount, ServingRadius, Centers" << endl;
}

template <typename T_DIST>
void PCenter<T_DIST>::appendResultToSheet( const string &instanceFileName, ofstream &csvFile ) const
{
    csvFile << Timer::getLocalTime() << ", " << solvingAlgorithm << ", " << instanceFileName << ", " << maxIterCount << ", "
        << bestSolution.duration << ", " << bestSolution.iterCount << ", " << bestSolution.serveRadius << ", ";
    for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = bestSolution.center.begin(); iter != bestSolution.center.end(); iter++) {
        csvFile << *iter << "|";
    }
    csvFile << endl;
}




template <typename T_DIST>
void PCenter<T_DIST>::genInitSolution()
{
    // select a vertex as center randomly
    RangeRand viRand( graph.minVertexIndex, graph.maxVertexIndex );
    int firstCenter = viRand();
    center.insert( firstCenter );

    // select one of the longest arc with length longestDist randomly
    RandSelect rs( 2 );
    int longestEnd = graph.nthClosestVertex( firstCenter, graph.maxVertexIndex );
    TopologicalGraph<T_DIST>::Distance longestDist = graph.distance( firstCenter, longestEnd );
    for (int i = (graph.maxVertexIndex - 1); i >= graph.minVertexIndex; i--) {
        if (graph.distance( firstCenter, graph.nthClosestVertex( firstCenter, i ) ) == longestDist) {
            if (rs.isSelected()) {
                longestEnd = graph.nthClosestVertex( firstCenter, i );
            }
        } else {
            break;
        }
    }

    // select one of the vertex whose distance to longestEnd is shorter than longestDist randomly
    int secondCenter = graph.findVertexWithinRadius( longestEnd, longestDist );
    center.insert( secondCenter );

    // select other (pnum-2) centers
    initClosestCenter( firstCenter, secondCenter );
    for (int i = pnum - 2; i > 0; i--) {
        int fv = findFarthestVertex( closestCenter );
        int newCenter = graph.findVertexWithinRadius( fv, closestCenter[fv].dist[0] );
        center.insert( newCenter );
        addCenter( newCenter, closestCenter );
    }

    // init the min max min serve radius
    timer.record();
    bestSolution.serveRadius = closestCenter[findFarthestVertex( closestCenter )].dist[0];
    bestSolution.center = center;
    bestSolution.duration = timer.getDuration();
    bestSolution.iterCount = 0;
}

template <typename T_DIST>
void PCenter<T_DIST>::initClosestCenter( int firstCenter, int secondCenter )
{
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        TopologicalGraph<T_DIST>::Distance d1 = graph.distance( firstCenter, i );
        TopologicalGraph<T_DIST>::Distance d2 = graph.distance( secondCenter, i );
        if (d1 < d2) {
            closestCenter[i] = ClosestCenterQueue( firstCenter, d1, secondCenter, d2 );
        } else {
            closestCenter[i] = ClosestCenterQueue( secondCenter, d2, firstCenter, d1 );
        }
    }
}

template <typename T_DIST>
int PCenter<T_DIST>::findFarthestVertex( ClosestCenterTable &cct ) const
{
    typename TopologicalGraph<T_DIST>::Distance maxDist = TopologicalGraph<T_DIST>::MIN_DISTANCE;
    int farthestVertex;
    RandSelect rs( 2 );
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        if (maxDist < cct[i].dist[0]) {
            maxDist = cct[i].dist[0];
            farthestVertex = i;
            rs.reset( 2 );
        } else if (maxDist == cct[i].dist[0]) {
            if (rs.isSelected()) {
                farthestVertex = i;
            }
        }
    }

    return farthestVertex;
}

template <typename T_DIST>
typename TopologicalGraph<T_DIST>::Arc PCenter<T_DIST>::findLongestServeArc( ClosestCenterTable &cct ) const
{
    int fv = findFarthestVertex( cct );

    return TopologicalGraph<T_DIST>::Arc( cct[fv].center[0], fv, cct[fv].dist[0] );
}

template <typename T_DIST>
typename TopologicalGraph<T_DIST>::VertexSet PCenter<T_DIST>::findFarthestVertices( ClosestCenterTable &cct ) const
{
    typename TopologicalGraph<T_DIST>::Distance maxDist = TopologicalGraph<T_DIST>::MIN_DISTANCE;
    typename TopologicalGraph<T_DIST>::VertexSet farthestVertices;
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        if (maxDist < cct[i].dist[0]) {
            maxDist = cct[i].dist[0];
            farthestVertices.clear();
            farthestVertices.insert( i );
        } else if (maxDist == cct[i].dist[0]) {
            farthestVertices.insert( i );
        }
    }

    return farthestVertices;
}

template <typename T_DIST>
typename TopologicalGraph<T_DIST>::ArcSet PCenter<T_DIST>::findLongestServeArcs( ClosestCenterTable &cct ) const
{
    typename TopologicalGraph<T_DIST>::VertexSet fvs = findFarthestVertices( cct );
    typename TopologicalGraph<T_DIST>::ArcSet longestServeArcs;

    for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = fvs.begin(); iter != fvs.end(); iter++) {
        longestServeArcs.insert( TopologicalGraph<T_DIST>::Arc( cct[*iter].center[0], *iter, cct[*iter].dist[0] ) );
    }

    return longestServeArcs;
}

template <typename T_DIST>
void PCenter<T_DIST>::addCenter( int newCenter, ClosestCenterTable &cct )
{
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        TopologicalGraph<T_DIST>::Distance newDist = graph.distance( newCenter, i );
        if (newDist < cct[i].dist[0]) {
            cct[i].center[1] = cct[i].center[0];
            cct[i].dist[1] = cct[i].dist[0];
            cct[i].center[0] = newCenter;
            cct[i].dist[0] = newDist;
        } else if (newDist < cct[i].dist[1]) {
            cct[i].center[1] = newCenter;
            cct[i].dist[1] = newDist;
        }
    }
}

template <typename T_DIST>
void PCenter<T_DIST>::removeCenter( int oldCenter )
{
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        bool changed = false;
        if (closestCenter[i].center[0] == oldCenter) {
            closestCenter[i].center[0] = closestCenter[i].center[1];
            closestCenter[i].dist[0] = closestCenter[i].dist[1];
            changed = true;
        }

        if (changed || closestCenter[i].center[1] == oldCenter) {
            // locate the closest center and add it to the queue
            for (int k = graph.minVertexIndex; k <= graph.maxVertexIndex; k++) {
                TopologicalGraph<T_DIST>::VertexSet::iterator iter = center.find( graph.nthClosestVertex( i, k ) );
                if ((iter != center.end()) && (*iter != closestCenter[i].center[0]) && (*iter != oldCenter)) {
                    closestCenter[i].center[1] = *iter;
                    closestCenter[i].dist[1] = graph.distance( *iter, i );
                    break;
                }
            }
        }
    }
}

template <typename T_DIST>
typename PCenter<T_DIST>::CenterSwap PCenter<T_DIST>::getRandSwap() const
{
    RandSelect selectOld( 1 );
    RandSelect selectNew( 1 );

    CenterSwap cs;
    vector<bool> isCenter( graph.vertexAllocNum, false );

    for (TopologicalGraph<T_DIST>::VertexSet::iterator iter = center.begin(); iter != center.end(); iter++) {
        isCenter[*iter] = true;
    }

    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        if (isCenter[i]) {
            if (selectOld.isSelected()) {
                cs.oldCenter = i;
            }
        } else {
            if (selectNew.isSelected()) {
                cs.newCenter = i;
            }
        }
    }

    return cs;
}

template <typename T_DIST>
typename TopologicalGraph<T_DIST>::Distance PCenter<T_DIST>::perturbRRGA( int perturbStrength )
{
    // remove some centers
    for (int i = perturbStrength; i > 0; i--) {
        RangeRand rr( 0, center.size() - 1 );
        TopologicalGraph<T_DIST>::VertexSet::iterator iter = center.begin();
        advance( iter, rr() );
        removeCenter( *iter );
        center.erase( iter );
    }
    // select some new centers
    for (int i = perturbStrength; i > 0; i--) {
        int fv = findFarthestVertex( closestCenter );
        int newCenter = graph.findVertexWithinRadius( fv, closestCenter[fv].dist[0] );
        center.insert( newCenter );
        addCenter( newCenter, closestCenter );
    }
    // recalculate the min radius
    return closestCenter[findFarthestVertex( closestCenter )].dist[0];
}




#endif