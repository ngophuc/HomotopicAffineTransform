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

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/images/ImageSelector.h"

#include "DGtal/io/colormaps/GradientColorMap.h"

using namespace std;
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

RationalPoint getIntegerPoint(RationalPoint p);

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

bool containEdge(const std::vector<RationalPoint>& face, RationalPoint e1, RationalPoint e2);
bool containEdge(const std::vector<int>& face, std::pair<int,int> e);
int findEdege(std::vector<std::pair<int, int> > edges, std::pair<int,int> e);

std::vector<int> findContainEdge(const std::vector<std::vector<RationalPoint> >& face, RationalPoint e1, RationalPoint e2);
std::vector<int> findContainEdge(const std::vector<std::vector<int> >& face, std::pair<int,int> e);

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

//#include "Functions.ih" // Includes inline functions.
#endif // FUNCTIONS_H
