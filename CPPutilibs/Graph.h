/**
*   usage : template for Geometrical Graph and Topological Graph.
*           some basic related operation on graph is implemented.
*/

#ifndef GRAPH_H
#define GRAPH_H



#include <vector>
#include <set>
#include <array>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <typeinfo>
#include "RandSelect.h"

class Graph
{
public:
    const unsigned vertexNum;

protected:
    Graph( unsigned vn ) : vertexNum( vn ) {}
    ~Graph() {}
};


class GeometricalGraph : public Graph
{
public:
    typedef double Distance;
    typedef double Coord;
    class Point
    {
    public:
        Point() {}
        Point( Coord cx, Coord cy ) : x( cx ), y( cy ) {}

        Coord x;
        Coord y;
    };
    typedef std::vector<Point> PointList;

    GeometricalGraph( const PointList &pl ) : Graph( pl.size() ), pointList( pl ) {}
    ~GeometricalGraph() {}

    const Point& point( int i ) const
    {
        return pointList[i];
    }

private:
    PointList pointList;
};


template <typename T_DIST = unsigned>
class TopologicalGraph : public Graph
{
public:
    // the distance or weight on edges
    typedef std::set<int> VertexSet;
    typedef T_DIST Distance;
    struct Arc
    {
        class CmpClass
        {
        public:
            bool operator() ( const Arc& lhs, const Arc& rhs ) const
            {
                return (lhs.startVertex < rhs.startVertex);
            }
        };

        Arc() {}
        Arc( int sv, int ev, Distance d )
            : startVertex( sv ), endVertex( ev ), dist( d )
        {
        }
        ~Arc() {}

        int startVertex;
        int endVertex;
        Distance dist;
    };
    typedef std::vector<Arc> ArcList;
    typedef std::set<Arc, typename Arc::CmpClass> ArcSet;
    // the No. N element on line M is the distance from vertex N to vertex M
    typedef std::vector< std::vector<Distance> > DistanceMatrix;
    // for each line N, it records the index of vertices sorted by their distance to vertex N with ascending order
    typedef std::vector< std::vector<int> > DistSeqTable;


    static const Distance MAX_DISTANCE;
    static const Distance MIN_DISTANCE;
    static const int DEFAULT_MIN_VERTEX_INDEX = 1;
    static const int MULTIPLICATION = 10000;


    const int minVertexIndex;
    const int maxVertexIndex;
    const int vertexAllocNum;


    // get the shortest distance between start and end
    Distance distance( int start, int end ) const
    {
        return shortestDist[start][end];
    }

    // get the Nth closest vertex from start
    // shiftedN = (N-1) + minVertexIndex;
    int nthClosestVertex( int start, int shiftedN ) const
    {
        return distSeq[start][shiftedN];
    }

    // find a vertex whose distance to start is shorter the radius randomly
    // return -1 if can not find any one
    int findVertexWithinRadius( int start, Distance radius ) const;
    // find the number of vertices whose distance to start is shorter the radius
    // this number can be used to get all vertices within radius by distSeq
    int findVertexNumWithinRadius( int start, Distance radius ) const;

    const DistanceMatrix& getShortestPath();
    void printShortestDist( std::ostream &os ) const;
    const DistSeqTable& getDistSeqTable();
    void printDistSeqTable( std::ostream &os ) const;
    bool multiplied() const
    {
        return isMultiplied;
    }

protected:
    TopologicalGraph( unsigned vertexNumber, int minVertexIndex );
    virtual ~TopologicalGraph();


    // essential elements ( must be initialized in constructor )
    DistanceMatrix adjMat;          // use elements within [minVertexIndex,maxVertexIndex]

    // optional elements ( can be generated later by calling member methods )
    DistanceMatrix shortestDist;    // use elements within [minVertexIndex,maxVertexIndex]
    DistSeqTable distSeq;           // use elements within [minVertexIndex,maxVertexIndex]
    bool shortestDistSolved;
    bool distSeqSolved;
    bool isMultiplied;

private:
    const DistanceMatrix& getShortestPathByDijkstra();
    const DistanceMatrix& getShortestPathByFloyd();
    const DistSeqTable& getDistSeqTableBySTLsort();
    const DistSeqTable& getDistSeqTableByInsertSort();
};

template <typename T_DIST = unsigned>
class UndirectedGraph : public TopologicalGraph<T_DIST>
{
public:
    // get a symmetrical adjMat
    UndirectedGraph( const ArcList &arcList, unsigned vertexNumber, int minVertexIndex = DEFAULT_MIN_VERTEX_INDEX );
    // making T_DIST double is recommanded, or the distances will got precision loss
    // ( the shortestDist will be generated automatically )
    UndirectedGraph( const GeometricalGraph& gg );
    ~UndirectedGraph();
};

template <typename T_DIST = unsigned>
class DirectedGraph : public TopologicalGraph<T_DIST>
{
public:
    // get a asymmetrical adjMat
    DirectedGraph( const ArcList &arcList, unsigned vertexNumber, int minVertexIndex = DEFAULT_MIN_VERTEX_INDEX );
    ~DirectedGraph();
};










// GeometricalGraph =======================


// TopologicalGraph =======================

template <typename T_DIST>
const typename TopologicalGraph<T_DIST>::Distance TopologicalGraph<T_DIST>::MAX_DISTANCE = (std::numeric_limits<T_DIST>::max() / static_cast<T_DIST>(2));

template <typename T_DIST>
const typename TopologicalGraph<T_DIST>::Distance TopologicalGraph<T_DIST>::MIN_DISTANCE = 0;



template <typename T_DIST>
TopologicalGraph<T_DIST>::TopologicalGraph( unsigned vn, int mvi )
: Graph( vn ), minVertexIndex( mvi ), maxVertexIndex( vn + mvi - 1 ), vertexAllocNum( vn + mvi ), isMultiplied(false),
adjMat( vn + mvi, std::vector<Distance>( vn + mvi, MAX_DISTANCE ) ), shortestDistSolved( false ), distSeqSolved( false )
{
    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        adjMat[i][i] = 0;
    }
}

template <typename T_DIST>
TopologicalGraph<T_DIST>::~TopologicalGraph()
{
}

template <typename T_DIST>
int TopologicalGraph<T_DIST>::findVertexWithinRadius( int start, typename TopologicalGraph::Distance radius ) const
{
    RandSelect rs( 1 );
    int vertex = -1;
    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        if (distance( start, nthClosestVertex( start, i ) ) < radius) {
            if (rs.isSelected()) {
                vertex = nthClosestVertex( start, i );
            }
        } else {
            break;
        }
    }

    return vertex;
}

template <typename T_DIST>
int TopologicalGraph<T_DIST>::findVertexNumWithinRadius( int start, typename TopologicalGraph::Distance radius ) const
{
    int i = minVertexIndex;
    for (; i <= maxVertexIndex; i++) {
        if (distance( start, nthClosestVertex( start, i ) ) >= radius) {
            break;
        }
    }

    return i;
}

template <typename T_DIST>
const typename TopologicalGraph<T_DIST>::DistanceMatrix& TopologicalGraph<T_DIST>::getShortestPath()
{
    if (!shortestDistSolved) {
        getShortestPathByFloyd();   /// put your own specified algorithm here

        shortestDistSolved = true;
        distSeqSolved = false;
    }

    return shortestDist;
}

template <typename T_DIST>
const typename TopologicalGraph<T_DIST>::DistSeqTable& TopologicalGraph<T_DIST>::getDistSeqTable()
{
    if (!distSeqSolved) {
        getShortestPath();  // prerequisite  
        /// put your own specified algorithm below
        getDistSeqTableBySTLsort();

        distSeqSolved = true;
    }

    return distSeq;
}

template <typename T_DIST>
void TopologicalGraph<T_DIST>::printDistSeqTable( std::ostream &os ) const
{
    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        for (int j = minVertexIndex; j < maxVertexIndex; j++) {
            os << shortestDist[i][j] << ',';
        }
        os << shortestDist[i][maxVertexIndex] << '\n';
    }
}

template <typename T_DIST>
void TopologicalGraph<T_DIST>::printShortestDist( std::ostream &os ) const
{
    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        for (int j = minVertexIndex; j < maxVertexIndex; j++) {
            if (shortestDist[i][j] == MAX_DISTANCE) {
                os << "x,";
            } else {
                os << shortestDist[i][j] << ',';
            }
        }

        if (shortestDist[i][maxVertexIndex] == MAX_DISTANCE) {
            os << "x\n";
        } else {
            os << shortestDist[i][maxVertexIndex] << '\n';
        }
    }
}


template <typename T_DIST>
const typename TopologicalGraph<T_DIST>::DistanceMatrix& TopologicalGraph<T_DIST>::getShortestPathByDijkstra()
{
    enum VertexState { IN_SET, OUT_OF_SET };

    shortestDist = adjMat;

    // find the shortest path start from each vertex
    for (int startVertex = minVertexIndex; startVertex <= maxVertexIndex; startVertex++) {
        // create an array of flag to record whether a vertex is in the final set
        std::vector<VertexState> vertexState( vertexAllocNum, VertexState::OUT_OF_SET );
        vertexState[startVertex] = VertexState::IN_SET;
        // loop (vertexNum-1) times to add every vertex into the final set
        for (int k = minVertexIndex; k < maxVertexIndex; k++) {
            int closestVertex = 0;    // vertex with index of 0 is not used, so the distance will always be MAX_DISTANCE
            // find the closest vertex outside the final set
            for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
                if (vertexState[i] == VertexState::OUT_OF_SET &&
                    shortestDist[startVertex][i] < shortestDist[startVertex][closestVertex]) {
                    closestVertex = i;
                }
            }
            // add it into the final set
            vertexState[closestVertex] = VertexState::IN_SET;
            // update current shortest path from start vertex to vertexs outside of the final set
            for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
                Distance newDist = shortestDist[startVertex][closestVertex] + shortestDist[closestVertex][i];
                if (vertexState[i] == VertexState::OUT_OF_SET &&
                    newDist < shortestDist[startVertex][i]) {
                    shortestDist[startVertex][i] = newDist;
                }
            }
        }
    }


    return shortestDist;
}

template <typename T_DIST>
const typename TopologicalGraph<T_DIST>::DistanceMatrix& TopologicalGraph<T_DIST>::getShortestPathByFloyd()
{
    shortestDist = adjMat;

    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        for (int j = minVertexIndex; j <= maxVertexIndex; j++) {
            for (int k = minVertexIndex; k <= maxVertexIndex; k++) {
                Distance newDist = shortestDist[j][i] + shortestDist[i][k];
                if (newDist < shortestDist[j][k]) {
                    shortestDist[j][k] = newDist;
                }
            }
        }
    }

    return shortestDist;
}

template <typename T_DIST>
const typename TopologicalGraph<T_DIST>::DistSeqTable& TopologicalGraph<T_DIST>::getDistSeqTableBySTLsort()
{
    class DistCmp
    {
    public:
        DistCmp( const TopologicalGraph<T_DIST> &g, int sv ) : graph( g ), startVertex( sv ) {}
        ~DistCmp() {}

        bool operator() ( const int &a, const int &b )
        {
            return graph.distance( startVertex, a ) < graph.distance( startVertex, b );
        }

    private:
        const TopologicalGraph<T_DIST> &graph;
        int startVertex;
    };

    // init distSeq
    std::vector<int> seq( vertexAllocNum );
    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        seq[i] = i;
    }
    distSeq = DistSeqTable( vertexAllocNum, seq );

    // do sort
    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        DistCmp cmp( *this, i );
        sort( distSeq[i].begin() + minVertexIndex, distSeq[i].end(), cmp );
    }

    return distSeq;
}

template <typename T_DIST>
const typename TopologicalGraph<T_DIST>::DistSeqTable& TopologicalGraph<T_DIST>::getDistSeqTableByInsertSort()
{
    // init distSeq
    std::vector<int> seq( vertexAllocNum );
    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        seq[i] = i;
    }
    distSeq = DistSeqTable( vertexAllocNum, seq );

    // do sort
    for (int i = minVertexIndex; i <= maxVertexIndex; i++) {
        //  do insertion sort on each vertex
        for (int j = minVertexIndex; j <= maxVertexIndex; j++) {
            int k = j;
            for (Distance d = distance( i, j ); k > 0; k--) {
                if (d < distance( i, distSeq[i][k - 1] )) {
                    distSeq[i][k] = distSeq[i][k - 1];
                } else {
                    distSeq[i][k] = j;
                    break;
                }
            }
            if (k == 0) {
                distSeq[i][0] = j;
            }
        }
    }

    return distSeq;
}



// UndirectedGraph =======================
template <typename T_DIST>
UndirectedGraph<T_DIST>::UndirectedGraph( const ArcList &arcList, unsigned vn, int mvi )
: TopologicalGraph( vn, mvi )
{
    for (ArcList::const_iterator iter = arcList.begin(); iter != arcList.end(); iter++) {
        adjMat[iter->startVertex][iter->endVertex] = iter->dist;
        adjMat[iter->endVertex][iter->startVertex] = iter->dist;
    }
}

template <typename T_DIST>
UndirectedGraph<T_DIST>::UndirectedGraph( const GeometricalGraph &gg )
: TopologicalGraph( gg.vertexNum, 0 )
{
    if (typeid(T_DIST) == typeid(unsigned)) {
        isMultiplied = true;
        // calculate distance between each pair of points
        for (unsigned i = 0; i < gg.vertexNum; i++) {
            for (unsigned j = 0; j < i; j++) {
                GeometricalGraph::Coord dx = gg.point( i ).x - gg.point( j ).x;
                GeometricalGraph::Coord dy = gg.point( i ).y - gg.point( j ).y;
                adjMat[i][j] = adjMat[j][i] = MULTIPLICATION * sqrt( dx*dx + dy*dy );
            }
            adjMat[i][i] = 0;
        }
    } else {
        // calculate distance between each pair of points
        for (unsigned i = 0; i < gg.vertexNum; i++) {
            for (unsigned j = 0; j < i; j++) {
                GeometricalGraph::Coord dx = gg.point( i ).x - gg.point( j ).x;
                GeometricalGraph::Coord dy = gg.point( i ).y - gg.point( j ).y;
                adjMat[i][j] = adjMat[j][i] = sqrt( dx*dx + dy*dy );
            }
            adjMat[i][i] = 0;
        }
    }
    shortestDist = adjMat;
    shortestDistSolved = true;
}

template <typename T_DIST>
UndirectedGraph<T_DIST>::~UndirectedGraph()
{
}


// DirectedGraph =======================
template <typename T_DIST>
DirectedGraph<T_DIST>::DirectedGraph( const ArcList &arcList, unsigned vn, int mvi )
: TopologicalGraph( vn, mvi )
{
    for (ArcList::const_iterator iter = arcList.begin(); iter != arcList.end(); iter++) {
        adjMat[iter->startVertex][iter->endVertex] = iter->dist;
    }
}

template <typename T_DIST>
DirectedGraph<T_DIST>::~DirectedGraph()
{
}



#endif