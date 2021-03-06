#include "Tester.h"

using namespace std;


int main( int argc, char *argv[] )
{
    testDouble();
    //testTimer();
    //testRangeRand();
    //testRandSelect();
    //testGraph();

    return 0;
}



void testDouble()
{
    Double d( 10.000001 );
    Double d1( 10 );
    double d2( 10 );
    double dd( 9.999999 );

    cout << d1 + d2 << endl;
    cout << d2 + d1 << endl;
    cout << d + d1 << endl;
    cout << (d1 == d) << ',' << (d == d1) << endl;
    cout << (d1 == d2) << ',' << (d2 == d1) << endl;
    cout << (d == d2) << ',' << (d2 == d) << endl;
    cout << (dd == d2) << ',' << (d2 == dd) << endl;
    cout << &d << endl;
}

void testTimer()
{
    cout << Timer::getLocalTime() << endl;

    Timer timer;

    timer.record( Timer::INFO, string( "1000 loop start" ) );
    for (int i = 0; i < 1000; i++) {}    // some time-consuming procedure writen by you.
    timer.record( Timer::INFO, string( "1000 loop end" ) );
    timer.printAll( cout );

    timer.reset();
    timer.recordAndPrint( Timer::INFO, string( "5000 loop start" ), cout );
    for (int i = 0; i < 5000; i++) {}    // some time-consuming procedure writen by you.
    timer.recordAndPrint( Timer::INFO, string( "5000 loop end" ), cout );
}

void testRangeRand()
{
    RangeRand r( 1, 10 );

    for (int i = 0; i < 20; i++) {
        cout << r() << endl;
    }
}

void testRandSelect()
{
    int result = 0;

    RandSelect rs;
    for (int i = 0; i < 10; i++) {
        if (rs.isSelected()) {
            result = i;
        }
    }
    cout << result << endl;

    rs.reset( 2 );    // select 1 out of 11 elements (0-9 and the last selected number)
    for (int i = 0; i < 10; i++) {
        if (rs.isSelected()) {
            result = i;
        }
        cout << result << endl;
    }
}

void testGraph()
{
    // ====== usage for an GeometricalGraph ======
    const int pointNum = 20;
    const int coordRange = 10;
    RangeRand plrr( -coordRange, coordRange );
    RangeRand vrr( 0, pointNum - 1 );

    // generate a pointList
    GeometricalGraph::PointList pointList;
    for (int i = 0; i < pointNum; i++) {
        pointList.push_back( GeometricalGraph::Point( plrr(), plrr() ) );
    }

    GeometricalGraph gg( pointList );



    // ====== usage for an UndDirectedGraph ======
    UndirectedGraph<double> ug( gg );
    ug.getShortestPath();    // can be leave out
    ug.getDistSeqTable();
    ug.printShortestDist( cout );
    ug.printDistSeqTable( cout );

    int vertex = ug.nthClosestVertex( vrr(), vrr() );
    vertex = ug.findVertexWithinRadius( vrr(), coordRange );
    int vnum = ug.findVertexNumWithinRadius( vrr(), coordRange );
    TopologicalGraph<double>::Distance dd = ug.distance( vrr(), vrr() );



    // ====== usage for an DirectedGraph ======
    const unsigned nodeNum = 10;
    const unsigned arcNum = 20;
    const int minVertexIndex = 1;
    RangeRand alrr( minVertexIndex, nodeNum + minVertexIndex - 1 );
    RangeRand drr( 1, 255 );

    // generate an arcList for an directed graph with unsigned distance randomly
    TopologicalGraph<>::ArcList arcList;
    for (int i = 0; i < arcNum; i++) {
        int startVertex = alrr();
        int endVertex;
        do {
            endVertex = alrr();
        } while (endVertex != startVertex);
        arcList.push_back( TopologicalGraph<>::Arc( startVertex, endVertex, drr() ) );
    }

    DirectedGraph<> dg( arcList, nodeNum, minVertexIndex );
    dg.getShortestPath();    // can be leave out
    dg.getDistSeqTable();
    dg.printShortestDist( cout );
    dg.printDistSeqTable( cout );

    vertex = dg.nthClosestVertex( alrr(), alrr() );
    vertex = dg.findVertexWithinRadius( alrr(), drr() );
    vnum = dg.findVertexNumWithinRadius( alrr(), drr() );
    TopologicalGraph<>::Distance ud = dg.distance( alrr(), alrr() );
}