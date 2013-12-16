#include "PCenter.h"

using namespace std;

PCenter::PCenter( UndirectedGraph &ug, int pn )
: graph( ug ), pnum( pn ), closestCenter( graph.vertexAllocNum, ClosestCenterQueue() )
{
    graph.getDistSeqTable();
}


PCenter::~PCenter()
{
}


void PCenter::solve( int maxIterCount )
{
    genInitSolution();

    RandSelect rs( 2 );
    for (int iterCount = 0; iterCount < maxIterCount; iterCount++) {
        CenterSwap centerSwap;
        Graph::Distance minRadius = Graph::MAX_DISTANCE;

        Graph::Arc longestServeArc = findLongestServeArc( closestCenter );

        int longestEnd = longestServeArc.endVertex;
        Graph::Distance longestDist = longestServeArc.dist;
        // try each vertex whose distance to longestEnd is shorter than longestDist
        for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
            int newCenter = graph.nthClosestVertex( longestEnd, i );
            if (graph.distance( longestEnd, newCenter ) < longestDist) {
                // find the best swap between center i and non-center vertices
                ClosestCenterTable tmpCCT( closestCenter );
                addCenter( newCenter, tmpCCT );
                Graph::Distance radiusAfterAdd = closestCenter[findFarthestVertex( tmpCCT )].dist[0];
                // calculate new radius for removing each center
                for (Graph::VertexSet::iterator iter = center.begin(); iter != center.end(); iter++) {
                    // when *iter is removed
                    int removedCenter = *iter;
                    Graph::Distance radiusAfterRemove = radiusAfterAdd;
                    for (int k = graph.minVertexIndex; k <= graph.maxVertexIndex; k++) {
                        if (closestCenter[k].center[0] == removedCenter) {
                            Graph::Distance newDist = closestCenter[k].dist[1];
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

void PCenter::greedySolve( int maxIterCount )
{
    genInitSolution();

    RandSelect rs( 2 );
    for (int iterCount = 0; iterCount < maxIterCount; iterCount++) {
        CenterSwap centerSwap;
        Graph::Distance minRadius = Graph::MAX_DISTANCE;

        Graph::ArcSet lsa = findLongestServeArcs( closestCenter );

        for (Graph::ArcSet::iterator iter = lsa.begin(); iter != lsa.end(); iter++) {
            int longestEnd = iter->endVertex;
            Graph::Distance longestDist = iter->dist;

            for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
                int newCenter = graph.nthClosestVertex( longestEnd, i );
                if (graph.distance( longestEnd, newCenter ) < longestDist) {
                    // find the best swap between center i and non-center vertices
                    ClosestCenterTable tmpCCT( closestCenter );
                    addCenter( newCenter, tmpCCT );
                    Graph::Distance radiusAfterAdd = closestCenter[findFarthestVertex( tmpCCT )].dist[0];
                    // calculate new radius for removing each center
                    for (Graph::VertexSet::iterator iter = center.begin(); iter != center.end(); iter++) {
                        // when *iter is removed
                        int removedCenter = *iter;
                        Graph::Distance radiusAfterRemove = radiusAfterAdd;
                        for (int k = graph.minVertexIndex; k <= graph.maxVertexIndex; k++) {
                            if (closestCenter[k].center[0] == removedCenter) {
                                Graph::Distance newDist = closestCenter[k].dist[1];
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

bool PCenter::check() const
{
    for (int i = graph.minVertexIndex; i < graph.maxVertexIndex; i++) {
        Graph::Distance minRadius = Graph::MAX_DISTANCE;
        for (Graph::VertexSet::iterator iter = bestSolution.center.begin(); iter != bestSolution.center.end(); iter++) {
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

void PCenter::printResult( ostream &os ) const
{
    os << "The max serving radius is : " << bestSolution.serveRadius << endl;
    os << "The indexes of the vertices which are chosed as centers are :\n";
    for (Graph::VertexSet::iterator iter = bestSolution.center.begin(); iter != bestSolution.center.end(); iter++) {
        os << *iter << "|";
    }
}

void PCenter::initResultSheet( std::ofstream &csvFile )
{
    csvFile << "Instance, " << "Duration, " << "IterCount, " << "ServingRadius, " << "Centers" << endl;
}

void PCenter::appendResultToSheet( const char *instanceFileName, ofstream &csvFile ) const
{
    csvFile << instanceFileName << ", " << bestSolution.duration << ", "
        << bestSolution.iterCount << ", " << bestSolution.serveRadius << ", ";
    for (Graph::VertexSet::iterator iter = bestSolution.center.begin( ); iter != bestSolution.center.end( ); iter++) {
        csvFile << *iter << "|";
    }
    csvFile << endl;
}




void PCenter::genInitSolution()
{
    // select a vertex as center randomly
    RangeRand viRand( graph.minVertexIndex, graph.maxVertexIndex );
    int firstCenter = viRand();
    center.insert( firstCenter );

    // select one of the longest arc with length longestDist randomly
    RandSelect rs( 2 );
    int longestEnd = graph.nthClosestVertex( firstCenter, graph.maxVertexIndex );
    Graph::Distance longestDist = graph.distance( firstCenter, longestEnd );
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
        addCenter( newCenter, closestCenter );
        center.insert( newCenter );
    }

    // init the min max min serve radius
    bestSolution.serveRadius = closestCenter[findFarthestVertex( closestCenter )].dist[0];
}

void PCenter::initClosestCenter( int firstCenter, int secondCenter )
{
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        Graph::Distance d1 = graph.distance( firstCenter, i );
        Graph::Distance d2 = graph.distance( secondCenter, i );
        if (d1 < d2) {
            closestCenter[i] = ClosestCenterQueue( firstCenter, d1, secondCenter, d2 );
        } else {
            closestCenter[i] = ClosestCenterQueue( secondCenter, d2, firstCenter, d1 );
        }
    }
}

int PCenter::findFarthestVertex( ClosestCenterTable &cct ) const
{
    Graph::Distance maxDist = Graph::MIN_DISTANCE;
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

Graph::Arc PCenter::findLongestServeArc( ClosestCenterTable &cct ) const
{
    int fv = findFarthestVertex( cct );

    return Graph::Arc( cct[fv].center[0], fv, cct[fv].dist[0] );
}

Graph::VertexSet PCenter::findFarthestVertices( ClosestCenterTable &cct ) const
{
    Graph::Distance maxDist = Graph::MIN_DISTANCE;
    Graph::VertexSet farthestVertices;
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

Graph::ArcSet PCenter::findLongestServeArcs( ClosestCenterTable &cct ) const
{
    Graph::VertexSet fvs = findFarthestVertices( cct );
    Graph::ArcSet longestServeArcs;

    for (Graph::VertexSet::iterator iter = fvs.begin(); iter != fvs.end(); iter++) {
        longestServeArcs.insert( Graph::Arc( cct[*iter].center[0], *iter, cct[*iter].dist[0] ) );
    }

    return longestServeArcs;
}

void PCenter::addCenter( int newCenter, ClosestCenterTable &cct )
{
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        Graph::Distance newDist = graph.distance( newCenter, i );
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

void PCenter::removeCenter( int oldCenter )
{
    for (int i = graph.minVertexIndex; i <= graph.maxVertexIndex; i++) {
        if (closestCenter[i].center[0] == oldCenter) {
            closestCenter[i].center[0] = closestCenter[i].center[1];
            closestCenter[i].dist[0] = closestCenter[i].dist[1];
        }

        if (closestCenter[i].center[1] == oldCenter) {
            // locate the sequence number of old center 
            // (there is no available center to be added into the queue before it)
            int k = graph.minVertexIndex;
            for (; k <= graph.maxVertexIndex; k++) {
                if (graph.nthClosestVertex( i, k ) == oldCenter) {
                    break;
                }
            }
            // locate the closest center and add it to the queue
            for (; k <= graph.maxVertexIndex; k++) {
                Graph::VertexSet::iterator iter = center.find( graph.nthClosestVertex( i, k ) );
                if (iter != center.end()) {
                    closestCenter[i].center[1] = *iter;
                    closestCenter[i].dist[1] = graph.distance( *iter, i );
                    break;
                }
            }
        }
    }
}