#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <assert.h>

//#include <gmp.h>
#include <gmpxx.h>

#include "DGtal/base/BasicTypes.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

using namespace DGtal;

#define SIZE 0.1
#define BORDER 1
#define EPSILON 1e-6

// a/b
typedef mpq_class Rational;
// x,y
typedef std::pair<Rational,Rational> RationalPoint;
//ax + by + c = 0
typedef std::tuple<Rational,Rational,Rational> Line;
//(a,b,c) => a*a+b*bb=c*c, sin = a/c, cos = b/c
typedef std::tuple<int,int,int> PythagoreanTriple;

std::ostream& operator<<(std::ostream& os, const RationalPoint p);

Line getLineFromPoints(RationalPoint p1, RationalPoint p2);
double getRealValue(Rational rp);
Rational min(Rational r1, Rational r2);
Rational max(Rational r1, Rational r2);

bool isEqual(RationalPoint p1, RationalPoint p2);
bool isSameVector(const std::vector<int>& v1, const std::vector<int>& v2);
bool isSameVector(const std::vector<Rational>& v1, const std::vector<Rational>& v2);
bool isSameVector(const std::vector<RationalPoint>& v1, const std::vector<RationalPoint>& v2);
std::vector<int> permuteCircularVector(const std::vector<int>& v, int k);
bool isSameCircularVector(const std::vector<int>& v1, const std::vector<int>& v2);
std::vector<RationalPoint> permuteCircularVector(const std::vector<RationalPoint>& v, int k);
bool isSameCircularVector(const std::vector<RationalPoint>& v1, const std::vector<RationalPoint>& v2);

int getUnion(std::vector<Z2i::Point> v1, std::vector<Z2i::Point> v2);
int getIntersection(std::vector<Z2i::Point> v1, std::vector<Z2i::Point> v2);
int getDifference(std::vector<Z2i::Point> v1, std::vector<Z2i::Point> v2);

int getSymetricDifference(std::vector<Z2i::Point> v1, std::vector<Z2i::Point> v2);
int getSymetricDifference(Z2i::DigitalSet s1, Z2i::DigitalSet s2);

double distance(RationalPoint p1, RationalPoint p2);
double distance(Z2i::Point p1, Z2i::Point p2);
double distance(Z2i::RealPoint p1, Z2i::RealPoint p2);

//Return area(face) > if face belongs to polygon, and -area(face) otherwise
double getDistanceArea(const std::vector<RationalPoint>& polygon, const std::vector<RationalPoint>& face);
double getIntersectArea(const std::vector<RationalPoint>& polygon1, const std::vector<RationalPoint>& polygon2, const std::vector<RationalPoint>& face);
double getIntersectArea(const std::vector<std::vector<RationalPoint> >& polygon1, const std::vector<std::vector<RationalPoint> >& polygon2, const std::vector<std::vector<RationalPoint> >& faces);

RationalPoint getBaryCenter(std::vector<RationalPoint> vP);
Z2i::RealPoint getBaryCenter(std::vector<Z2i::Point> vP);

Line getLine(std::pair<RationalPoint, RationalPoint> p);

bool isConvex(std::vector<RationalPoint> points);

bool isCollinear(RationalPoint p1, RationalPoint p2, RationalPoint p3);

std::vector<RationalPoint> sortCollinearPoints(std::vector<RationalPoint>);

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(RationalPoint p, RationalPoint q, RationalPoint r);

//returns true if line segment p1q1 and p2q2 intersect
bool isIntersect(RationalPoint , RationalPoint , RationalPoint , RationalPoint );
RationalPoint intersection(RationalPoint , RationalPoint , RationalPoint , RationalPoint );
std::vector<RationalPoint> intersection(const std::vector<RationalPoint>& face, RationalPoint sp1, RationalPoint sp2);
std::vector<RationalPoint> intersection(const std::vector<RationalPoint>& face, std::pair<RationalPoint,RationalPoint> line);
std::vector<RationalPoint> intersection(const std::vector<RationalPoint>& face, std::vector<std::pair<RationalPoint,RationalPoint> > lines);

// Returns true if the point p lies inside the polygon
bool isInsidePolygon(const std::vector<RationalPoint>& polygon, RationalPoint p);
// Returns true if the point p is completly inside the polygon
bool isInterieurPolygon(const std::vector<RationalPoint>& polygon, RationalPoint p);
// Returns true if the point p lies inside the pixel
bool isInsidePixel(DGtal::Z2i::Point px, RationalPoint p);
//Returns area of a polygon given order vertices
double areaPolygon(const std::vector<RationalPoint>& polygon);


double fullAngle(std::pair<RationalPoint,RationalPoint> v1, std::pair<RationalPoint, RationalPoint> v2);
double accuteAngle(std::pair<RationalPoint,RationalPoint> v1, std::pair<RationalPoint, RationalPoint> v2);

bool isAllTrue(std::vector<bool> vec);
bool isAllFalse(std::vector<bool> vec);
int getFirstTrue(std::vector<bool> vec);
int getFirstFalse(std::vector<bool> vec);

int findVertex(const std::vector<RationalPoint>& vector, RationalPoint v);
int findVertex(const std::vector<Z2i::Point>& vector, Z2i::Point v);
bool containElement(const std::vector<int>& vector, int element);
bool containElement(const std::vector<RationalPoint>& vector, RationalPoint element);
bool containElement(const std::vector<Z2i::Point>& vector, Z2i::Point element);

PythagoreanTriple convertAngle2Pythagore(double B, double e);

RationalPoint getFaceCenter(const std::vector<RationalPoint>& face);
std::vector<RationalPoint> getFaceCenter(const std::vector<std::vector<RationalPoint> >& face);

Z2i::RealPoint getFaceCenter(const std::vector<Z2i::RealPoint>& face);
std::vector<Z2i::RealPoint> getFaceCenter(const std::vector<std::vector<Z2i::RealPoint> >& face);

std::pair<RationalPoint,RationalPoint> getBoundingBox(const std::vector<RationalPoint>& points);
std::pair<RationalPoint,RationalPoint> getBoundingBox(const std::vector<std::pair<RationalPoint,RationalPoint> >& points);
std::pair<Z2i::Point,Z2i::Point> getBoundingBox(const std::vector<Z2i::Point>& points);
std::pair<Z2i::RealPoint,Z2i::RealPoint> getBoundingBox(const std::vector<Z2i::RealPoint>& points);

//Topology handling
std::vector<Z2i::Point> getDigitalSetPoint(Z2i::DigitalSet shape_set);
Z2i::DigitalSet getDigitalSet(std::vector<Z2i::Point> point_set);

bool isSimplePoint(Z2i::DigitalSet shape_set, Z2i::Point point, int adjacency = 4); //Adjacency : objet connectivity
bool isSimplePoint(std::vector<Z2i::Point> point_set, Z2i::Point point, int adjacency = 4);
Z2i::DigitalSet getSimplifiedPoint(Z2i::DigitalSet shape_set, int adjacency = 4);
std::vector<Z2i::Point> getSimplifiedPoint(std::vector<Z2i::Point> point_set, int adjacency = 4);
//(Equiv Equation 17): remove at most the simple points
Z2i::DigitalSet getSimplifiedPoint(Z2i::DigitalSet shape_set, std::vector<Z2i::Point> fixed_point, int adjacency = 4);
//Equation 17: min (|X| - |Y|)
Z2i::DigitalSet getMinimizePoint(Z2i::DigitalSet shape_set, std::vector<Z2i::Point> fixed_point, int adjacency = 4);
//Equation 15: min ( |Xt \ Y|) + |Y \ Xt| )
Z2i::DigitalSet getMinimizeDifference(Z2i::DigitalSet shape_set, std::vector<Z2i::Point> fixed_point, int adjacency = 4);
//Equation 18: min ( distance (barycenter(Xt), barycenter(Y) )
Z2i::DigitalSet getMinimizeDistance(Z2i::DigitalSet shape_set, std::vector<Z2i::Point> fixed_point, int adjacency = 4);
//TODO: Equation 16: min ( area(|Xt \ Y|) + area(|Y \ Xt|) )
Z2i::DigitalSet getMinimizeArea(Z2i::DigitalSet shape_set, std::vector<Z2i::Point> fixed_point, int adjacency = 4);

//method = 0 => getSimplifiedPoint (=getMinimizePoint)
//method = 1 => getMinimizePoint
//method = 2 => getMinimizeDifference
//method = 3 => getMinimizeDistance
std::vector<Z2i::Point> getMinimizeDigitization(std::vector<Z2i::Point> point_set, std::vector<Z2i::Point> fixed_point, int adjacency = 4, int method=0);

std::vector<Z2i::DigitalSet> getConnectedComponent(Z2i::DigitalSet shape_set, int adjacency = 4);
std::vector<std::vector<Z2i::Point> > getConnectedComponent(std::vector<Z2i::Point> point_set, int adjacency = 4);
std::vector<int> getIdConnectedComponent(std::vector<Z2i::Point> point_set, int adjacency = 4);

int getBetii0(Z2i::DigitalSet shape_set, int adjacency = 4); //Connected component of objet
int getBetii1(Z2i::DigitalSet shape_set, int adjacency = 8); //Hole = connected component of complement - 1
int getEulerCharacteristic(Z2i::DigitalSet shape_set, int adjacency = 4);

int getBetii0(std::vector<Z2i::Point> point_set, int adjacency = 4);
int getBetii1(std::vector<Z2i::Point> point_set, int adjacency = 8);
int getEulerCharacteristic(std::vector<Z2i::Point> point_set, int adjacency = 4);

bool isHomotopopyEquivalence(Z2i::DigitalSet shape1, Z2i::DigitalSet shape2, int adjacency = 4);
bool isHomotopopyEquivalence(std::vector<Z2i::Point> shape1, std::vector<Z2i::Point> shape2, int adjacency = 4);

void writePoints(std::vector<Z2i::Point> pts, std::string filename);

std::vector<std::vector<RationalPoint> > sortAreaFaces(std::vector<std::vector<RationalPoint> > faces);
std::vector<size_t> sortAreaFaces(std::vector<double> vecArea);

//#include "Functions.ih" // Includes inline functions.
#endif // FUNCTIONS_H
