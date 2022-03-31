#include "Functions.h"

RationalPoint getIntegerPoint(RationalPoint p)
{
  Z2i::Point pp = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second))); //pixel center
  RationalPoint pr = RationalPoint(Rational(pp[0]),Rational(pp[1]));
  return pr;
}

//Equation of line passing through two points : ax+by+c=0
//(y1−y2)x+(x2−x1)y+x1y2−x2y1=0
Line getLineFromPoints(RationalPoint p1, RationalPoint p2)
{
  Rational a = p1.second - p2.second;
  Rational b = p2.first - p1.first;
  Rational c = p1.first*p2.second - p2.first*p1.second;
  Line l = std::make_tuple(a,b,c);
  return l;
}

double getRealValue(Rational rp)
{
  //return static_cast<double>(rp.get_num())/rp.get_den();
  return rp.get_d();
}

Rational min(Rational r1, Rational r2)
{
  if(r1<r2)
    return r1;
  return r2;
}

Rational max(Rational r1, Rational r2)
{
  if(r1<r2)
    return r2;
  return r1;
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
  
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
// See https://www.geeksforgeeks.org/orientation-3-ordered-points/
int orientation(RationalPoint p, RationalPoint q, RationalPoint r)
{
  Rational v1 = (q.second - p.second);
  Rational v2 = (r.first - q.first);
  Rational v3 = (q.first - p.first);
  Rational v4 = (r.second - q.second);
  Rational val = v1 * v2 - v3 * v4;
  //Rational val = (q.second - p.second) * (r.first - q.first) - (q.first - p.first) * (r.second - q.second);
  
  if (fabs(getRealValue(val))<EPSILON)//if (val == Rational(0))
    return 0;  // colinear
  
  return (val > 0)? 1: 2; // clock or counterclock wise
}

int orientationTmp(RationalPoint p, RationalPoint q, RationalPoint r)
{
  Rational v1 = (q.second - p.second);
  Rational v2 = (r.first - q.first);
  Rational v3 = (q.first - p.first);
  Rational v4 = (r.second - q.second);
  Rational val = v1 * v2 - v3 * v4;
  //Rational val = (q.second - p.second) * (r.first - q.first) - (q.first - p.first) * (r.second - q.second);
  std::cout<<"VAL ="<<getRealValue(val);
  if (val == Rational(0))
    return 0;  // colinear
  
  return (val > 0)? 1: 2; // clock or counterclock wise
}

bool isCollinear(Rational x1, Rational y1, Rational x2, Rational y2, Rational x3, Rational y3)
{
  // Calculation the area of triangle
  Rational a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
  if (a == Rational(0))
    return true;
  return false;
}

bool isCollinear(RationalPoint p1, RationalPoint p2, RationalPoint p3)
{
  return orientation(p1,p2,p3)==0;
}

std::vector<RationalPoint> sortCollinearPoints(std::vector<RationalPoint> vecP)
{
  std::vector<RationalPoint> vecP_sorted;
  assert(vecP.size() > 1);
  RationalPoint p1 = vecP.at(0);
  RationalPoint p2 = vecP.at(1);
  Rational dx = p2.first-p1.first;
  Rational dy = p2.second-p1.second;
  RationalPoint d = RationalPoint(dx,dy);
  //Line : p = p1 + alpha * d
  Rational alpha_min = 0; //p1
  //find the most-left point
  for(int it=2; it<vecP.size(); it++) {
    RationalPoint p = vecP.at(it);
    assert(isCollinear(p1, p2, p)==true);
    assert(d.first!=0 || d.second!=0);
    Rational alpha;
    if(d.first!=Rational(0))
      alpha = (p.first - p1.first) / d.first;
    else if(d.second!=Rational(0))
      alpha = (p.second - p1.second) / d.second;
    if(alpha < alpha_min)
      p1 = vecP.at(it);
  }
  //Most-left point is p1, thus the line is p = p1 + alpha * d
  std::vector<Rational> vecAlpha;
  for(int it=0; it<vecP.size(); it++) {
    Rational alpha;
    assert(d.first!=0 || d.second!=0);
    if(d.first!=Rational(0))
      alpha = (vecP.at(it).first - p1.first) / d.first;
    else if(d.second!=Rational(0))
      alpha = (vecP.at(it).second - p1.second) / d.second;
    vecAlpha.push_back(alpha);
  }
  std::vector<size_t> alpha_sorted = sort_indexes(vecAlpha);
  for(int it=0; it<alpha_sorted.size(); it++)
  vecP_sorted.push_back(vecP.at(alpha_sorted.at(it)));
  return vecP_sorted;
}

// Given three points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(RationalPoint p, RationalPoint q, RationalPoint r)
{
  if(!isCollinear(p, q, r)) //verify collinearity
    return false;
  if(fabs(getRealValue(q.first)-getRealValue(p.first))<EPSILON && fabs(getRealValue(q.second)-getRealValue(p.second))<EPSILON) //q = p
    return true;
  if(fabs(getRealValue(q.first)-getRealValue(r.first))<EPSILON && fabs(getRealValue(q.second)-getRealValue(r.second))<EPSILON) //q = r
    return true;
  if(fabs(getRealValue(q.first)-getRealValue(r.first))<EPSILON &&
     fabs(getRealValue(p.first)-getRealValue(r.first))<EPSILON &&
     q.second <= std::max(p.second, r.second) && q.second >= min(p.second, r.second) ) //vertical line
    return true;
  if(fabs(getRealValue(q.second)-getRealValue(r.second))<EPSILON &&
     fabs(getRealValue(p.second)-getRealValue(r.second))<EPSILON &&
     q.first <= std::max(p.first, r.first) && q.first >= std::min(p.first, r.first) ) //horizontal line
    return true;
  if (q.first <= std::max(p.first, r.first) && q.first >= std::min(p.first, r.first) &&
      q.second <= std::max(p.second, r.second) && q.second >= min(p.second, r.second))
    return true;
  return false;
}

// The main function that returns true if line segment 'p1q1' and 'p2q2' intersect.
bool isIntersect(RationalPoint p1, RationalPoint q1, RationalPoint p2, RationalPoint q2)
{
  // Find the four orientations needed for general and
  // special cases
  int o1 = orientation(p1, q1, p2);
  int o2 = orientation(p1, q1, q2);
  int o3 = orientation(p2, q2, p1);
  int o4 = orientation(p2, q2, q1);
  
  // General case
  if (o1 != o2 && o3 != o4)
    return true;
  
  // Special Cases
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1
  if (o1 == 0 && onSegment(p1, p2, q1)) return true;
  
  // p1, q1 and q2 are colinear and q2 lies on segment p1q1
  if (o2 == 0 && onSegment(p1, q2, q1)) return true;
  
  // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  if (o3 == 0 && onSegment(p2, p1, q2)) return true;
  
  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  if (o4 == 0 && onSegment(p2, q1, q2)) return true;
  
  return false; // Doesn't fall in any of the above cases
}

RationalPoint intersection(RationalPoint A, RationalPoint B, RationalPoint C, RationalPoint D)
{
  // If the lines are parallel then return a pair of INT_MAX
  RationalPoint p = std::make_pair(Rational(INT_MAX), Rational(INT_MAX));
  if(isIntersect(A, B, C, D)) {
    // Line AB represented as a1x + b1y = c1
    Rational a1 = B.second - A.second;
    Rational b1 = A.first - B.first;
    Rational c1 = a1*(A.first) + b1*(A.second);
    
    // Line CD represented as a2x + b2y = c2
    Rational a2 = D.second - C.second;
    Rational b2 = C.first - D.first;
    Rational c2 = a2*(C.first)+ b2*(C.second);
    
    Rational determinant = a1*b2 - a2*b1;
    // The lines are parallel
    if (determinant == Rational(0))
      return p;
    Rational x = (b2*c1 - b1*c2)/determinant;
    Rational y = (a1*c2 - a2*c1)/determinant;
    p = std::make_pair(x, y);
  }
  return p;
}

std::vector<RationalPoint> intersection(const std::vector<RationalPoint>& face, RationalPoint sp1, RationalPoint sp2)
{
  std::vector<RationalPoint> vecIntersection;
  for(int it=0; it<face.size()-1; it++) {
    RationalPoint p1 = face.at(it);
    RationalPoint p2 = face.at(it+1);
    if(isIntersect(p1, p2, sp1, sp2)) {
      RationalPoint p = intersection(p1,p2,sp1,sp2);
      vecIntersection.push_back(p);
    }
  }
  return vecIntersection;
}

std::vector<RationalPoint> intersection(const std::vector<RationalPoint>& face, std::pair<RationalPoint,RationalPoint> line)
{
  RationalPoint p1 = line.first;
  RationalPoint p2 = line.second;
  return intersection(face,p1,p2);
}

std::vector<RationalPoint> intersection(const std::vector<RationalPoint>& face, std::vector<std::pair<RationalPoint,RationalPoint> > lines)
{
  std::vector<RationalPoint> vecIntersection;
  for(int it=0; it<lines.size(); it++) {
    std::vector<RationalPoint> p = intersection(face, lines.at(it));
    for(int it_bis=0; it_bis<p.size(); it++)
    vecIntersection.push_back(p.at(it_bis));
  }
  return vecIntersection;
}

//see https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
bool isInsidePolygon(const std::vector<RationalPoint>& polygon, RationalPoint p)
{
  int n = polygon.size();
  // There must be at least 3 vertices in polygon[]
  if (n < 3)  return false;
  
  // Create a point for line segment from p to infinite (horizontal direction)
  RationalPoint extreme1 = RationalPoint(Rational(INT_MAX), Rational(p.second));
  // Degenerate case when the intersection point is a vertex of polygon, thus lauch two rays
  // Create a point for line segment from p to infinite (diagonal direction)
  //RationalPoint extreme2 = RationalPoint(Rational(INT_MAX), Rational(INT_MAX));
  // Create a point for line segment from p to infinite (vertical direction)
  RationalPoint extreme2 = RationalPoint(Rational(p.first), Rational(INT_MAX));
  
  RationalPoint extreme = extreme1;
  bool belong = true;
  //Find the ray that doest pass by any vertices of polygon
  do {
    belong = false;
    for(size_t it=0; it<polygon.size(); it++) {
      // Check p is a vertex of polygon, then it is inside
      if(polygon.at(it)==p)
        return true;
      else {
        // Check if 'polygon[i]' belongs to the line segment from 'p' to 'extreme',
        // if it is then launch another ray using dichotomy search
        if(isCollinear(p, polygon.at(it), extreme)) {
          extreme.first = (extreme.first + extreme2.first)/2;
          extreme.second = (extreme.second + extreme2.second)/2;
          belong = true;
        }
      }
    }
  } while (belong);
  
  // Count intersections of the above line with sides of polygon
  int count = 0, i = 0;
  do {
    int next = (i+1)%n;
    
    // Check if the line segment from 'p' to 'extreme' intersects
    // with the line segment from 'polygon[i]' to 'polygon[next]'
    bool intersect = isIntersect(polygon.at(i), polygon.at(next), p, extreme);
    if (intersect)
    {
      // If the point 'p' is colinear with line segment 'i-next',
      // then check if it lies on segment. If it lies, return true,
      // otherwise false
      if (orientation(polygon[i], p, polygon[next]) == 0)
        return onSegment(polygon[i], p, polygon[next]);
      count++;
    }
    i = next;
  } while (i != 0);
  // Return true if count is odd, false otherwise
  return (count&1);  // Same as (count%2 == 1)
}

bool isInsidePixel(DGtal::Z2i::Point px, RationalPoint p)
{
  DGtal::Z2i::RealPoint pr(getRealValue(p.first), getRealValue(p.second));
  double dx = px[0] - pr[0];
  double dy = px[1] - pr[1];
  if(fabs(dx)<=0.5 && fabs(dy)<=0.5)
    return true;
  return false;
}

double areaPolygon(const std::vector<RationalPoint>& polygon)
{
  Rational area(0);
  int n = polygon.size();
  int j= n - 1;
  for(int i=0; i<n; i++) {
    Rational xj = polygon.at(j).first;
    Rational xi = polygon.at(i).first;
    Rational yj = polygon.at(j).second;
    Rational yi = polygon.at(i).second;
    area += (xj+xi)*(yj-yi);
    j = i;
  }
  area = area/2;
  return abs(getRealValue(area));
}

double fullAngle(std::pair<RationalPoint,RationalPoint> v1, std::pair<RationalPoint, RationalPoint> v2)
{
  Z2i::RealPoint line1Start = Z2i::RealPoint(getRealValue(v1.first.first),getRealValue(v1.first.second));
  Z2i::RealPoint line1End = Z2i::RealPoint(getRealValue(v1.second.first),getRealValue(v1.second.second));
  Z2i::RealPoint line2Start = Z2i::RealPoint(getRealValue(v2.first.first),getRealValue(v2.first.second));
  Z2i::RealPoint line2End = Z2i::RealPoint(getRealValue(v2.first.first),getRealValue(v2.first.second));
  double angle1 = atan2(line1Start[1]-line1End[1], line1Start[0]-line1End[0]);
  double angle2 = atan2(line2Start[1]-line2End[1], line2Start[0]-line2End[0]);
  double result = (angle2-angle1);
  if (result<0)
    result+=2*M_PI;
  return result;
}

double accuteAngle(std::pair<RationalPoint,RationalPoint> v1, std::pair<RationalPoint, RationalPoint> v2)
{
  double result = fullAngle(v1, v2);
  if (result>M_PI)
    result-=M_PI;
  return result;
}

bool isAllTrue(std::vector<bool> vec)
{
  for(int it=0; it<vec.size(); it++)
  if(vec.at(it)==false)
    return false;
  return true;
}

bool isAllFalse(std::vector<bool> vec)
{
  for(int it=0; it<vec.size(); it++)
  if(vec.at(it)==true)
    return false;
  return true;
}

int getFirstTrue(std::vector<bool> vec)
{
  for(int it=0; it<vec.size(); it++)
  if(vec.at(it)==true)
    return it;
  return -1;
}

int getFirstFalse(std::vector<bool> vec)
{
  for(int it=0; it<vec.size(); it++)
  if(vec.at(it)==false)
    return it;
  return -1;
}

//Find and return the index of a vertex in the list
int findVertex(const std::vector<RationalPoint>& vector, RationalPoint v)
{
  for(int it=0; it<vector.size(); it++) {
    //if(vector.at(it)==v)
    if(fabs(getRealValue(vector.at(it).first)-getRealValue(v.first))<EPSILON && fabs(getRealValue(vector.at(it).second)-getRealValue(v.second))<EPSILON)
      return it;
  }
  return -1;
}

int findVertex(const std::vector<Z2i::Point>& vector, Z2i::Point v)
{
  for(int it=0; it<vector.size(); it++)
  if(vector.at(it)==v)
    return it;
  return -1;
}

bool containEdge(const std::vector<RationalPoint>& face, RationalPoint e1, RationalPoint e2)
{
  if(findVertex(face, e1) && findVertex(face, e2))
    return true;
  return false;
}

bool containEdge(const std::vector<int>& face, std::pair<int,int> e)
{
  if(containElement(face, e.first) && containElement(face, e.second))
    return true;
  return false;
}

int findEdege(std::vector<std::pair<int, int> > edges, std::pair<int,int> e)
{
  for(size_t it=0; it<edges.size(); it++) {
    if(edges.at(it).first==e.first && edges.at(it).second==e.second)
      return  it;
    if(edges.at(it).first==e.second && edges.at(it).second==e.first)
      return  it;
  }
  return -1;
}

std::vector<int> findContainEdge(const std::vector<std::vector<RationalPoint> >& face, RationalPoint e1, RationalPoint e2)
{
  std::vector<int> f;
  for(size_t it=0; it<face.size(); it++) {
    if(containEdge(face.at(it), e1, e2))
      f.push_back(it);
  }
  return f;
}

std::vector<int> findContainEdge(const std::vector<std::vector<int> >& face, std::pair<int,int> e)
{
  std::vector<int> f;
  for(size_t it=0; it<face.size(); it++) {
    if(containEdge(face.at(it), e))
      f.push_back(it);
  }
  return f;
}

bool containElement(const std::vector<int>& vector, int element)
{
  auto it = std::find(vector.begin(), vector.end(), element);
  if (it != vector.end())
    return true;
  return false;
}

bool containElement(const std::vector<RationalPoint>& vector, RationalPoint element)
{
  auto it = std::find(vector.begin(), vector.end(), element);
  if (it != vector.end())
    return true;
  return false;
}

bool containElement(const std::vector<Z2i::Point>& vector, Z2i::Point element)
{
  auto it = std::find(vector.begin(), vector.end(), element);
  if (it != vector.end())
    return true;
  return false;
}

// Function that returns true if it is possible
// to form a polygon with the given sides
bool isConvex(std::vector<double> a)
{
  int n = a.size();
  // Sum stores the sum of all the sides
  // and maxS stores the length of
  // the largest side
  double sum = 0, maxS = 0;
  for (int i = 0; i < n; i++) {
    sum += a[i];
    maxS = std::max(a[i], maxS);
  }
  
  // If the length of the largest side
  // is less than the sum of the
  // other remaining sides
  if ((sum - maxS) > maxS)
    return true;
  return false;
}

bool isEqual(RationalPoint p1, RationalPoint p2)
{
  return (fabs(getRealValue(p1.first) - getRealValue(p2.first)) < EPSILON) && (fabs(getRealValue(p1.second) - getRealValue(p2.second)) < EPSILON);
  //return p1.first==p2.first && p1.second==p2.second;
}

bool isSameVector(const std::vector<int>& v1, const std::vector<int>& v2)
{
  if(v1.size() != v2.size())
    return false;
  for(int it=0; it<v1.size(); it++)
  if(v1.at(it) != v2.at(it))
    return false;
  return true;
}

bool isSameVector(const std::vector<Rational>& v1, const std::vector<Rational>& v2)
{
  if(v1.size() != v2.size())
    return false;
  for(int it=0; it<v1.size(); it++)
  if(v1.at(it) != v2.at(it))
    return false;
  return true;
}

bool isSameVector(const std::vector<RationalPoint>& v1, const std::vector<RationalPoint>& v2)
{
  if(v1.size() != v2.size())
    return false;
  for(int it=0; it<v1.size(); it++)
  if(!isEqual(v1.at(it), v2.at(it)))
    return false;
  return true;
}

std::vector<int> permuteCircularVector(const std::vector<int>& v, int k)
{
  int start = k % v.size();
  if(start==0)
    return v;
  std::vector<int> res;
  for(int it=start; it<v.size(); it++)
  res.push_back(v.at(it));
  for(int it=0; it<start; it++)
  res.push_back(v.at(it));
  return res;
}

bool isSameCircularVector(const std::vector<int>& v1, const std::vector<int>& v2)
{
  if(v1.size() != v2.size())
    return false;
  for(int it=0; it<v2.size(); it++) {
    std::vector<int> tmp = permuteCircularVector(v2, it);
    if(isSameVector(v1, tmp))
      return true;
  }
  return false;
}

std::vector<RationalPoint> permuteCircularVector(const std::vector<RationalPoint>& v, int k)
{
  int start = k % v.size();
  if(start==0)
    return v;
  std::vector<RationalPoint> res;
  for(int it=start; it<v.size(); it++)
  res.push_back(v.at(it));
  for(int it=0; it<start; it++)
  res.push_back(v.at(it));
  return res;
}
bool isSameCircularVector(const std::vector<RationalPoint>& v1, const std::vector<RationalPoint>& v2)
{
  if(v1.size() != v2.size())
    return false;
  for(int it=0; it<v2.size(); it++) {
    std::vector<RationalPoint> tmp = permuteCircularVector(v2, it);
    if(isSameVector(v1, tmp))
      return true;
  }
  return false;
}

int getUnion(std::vector<Z2i::Point> v1, std::vector<Z2i::Point> v2)
{
  std::vector<Z2i::Point> v;
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  std::set_union(v1.begin(),v1.end(), v2.begin(),v2.end(), back_inserter(v));
  return v.size();
}

int getIntersection(std::vector<Z2i::Point> v1, std::vector<Z2i::Point> v2)
{
  std::vector<Z2i::Point> v;
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  std::set_intersection(v1.begin(),v1.end(), v2.begin(),v2.end(), back_inserter(v));
  return v.size();
}

int getDifference(std::vector<Z2i::Point> v1, std::vector<Z2i::Point> v2)
{
  std::vector<Z2i::Point> v;
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(v, v.begin()));
  return v.size();
}

int getSymetricDifference(std::vector<Z2i::Point> v1, std::vector<Z2i::Point> v2)
{
  std::vector<Z2i::Point> v;
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(v, v.begin()));
  return v.size();
}

int getSymetricDifference(Z2i::DigitalSet s1, Z2i::DigitalSet s2)
{
  std::vector<Z2i::Point> v1;
  for(auto it=s1.begin(); it!=s1.end(); it++)
  v1.push_back(*it);
  std::vector<Z2i::Point> v2;
  for(auto it=s2.begin(); it!=s2.end(); it++)
  v2.push_back(*it);
  return getSymetricDifference(v1,v2);
}

double distance(RationalPoint p1, RationalPoint p2)
{
  Rational dx = p1.first - p2.first;
  Rational dy = p1.second - p2.second;
  Rational dd = dx*dx + dy*dy;
  double d = sqrt(getRealValue(dd));
  return d;
}

double distance(Z2i::Point p1, Z2i::Point p2)
{
  int dx = p1[0] - p2[0];
  int dy = p1[1] - p2[1];
  int dd = dx*dx + dy*dy;
  return sqrt(dd);
}

double distance(Z2i::RealPoint p1, Z2i::RealPoint p2)
{
  double dx = p1[0] - p2[0];
  double dy = p1[1] - p2[1];
  double dd = dx*dx + dy*dy;
  return sqrt(dd);
}

double getDistanceArea(const std::vector<RationalPoint>& polygon, const std::vector<RationalPoint>& face)
{
  RationalPoint c = getBaryCenter(face);
  double area = areaPolygon(face);
  if(!isInsidePolygon(polygon, c))
    return -area;
  return area;
}

double getIntersectArea(const std::vector<RationalPoint>& polygon1, const std::vector<RationalPoint>& polygon2, const std::vector<RationalPoint>& face)
{
  double d1 = getDistanceArea(polygon1, face);
  double d2 = getDistanceArea(polygon2, face);
  if(d1>0 && d2>0) //face belongs to poly1 and poly2
    return  d1;
  return 0;
}

double getIntersectArea(const std::vector<std::vector<RationalPoint> >& polygon1, const std::vector<std::vector<RationalPoint> >& polygon2, const std::vector<std::vector<RationalPoint> >& faces)
{
  double d=0,d1=0,d2=0;
  int count=0;
  bool check=false;
  for(int it3=0; it3<faces.size(); it3++) {
    check=false;
    for(int it2=0; it2<polygon2.size() && !check; it2++) {
      if(distance(getBaryCenter(polygon2.at(it2)), getBaryCenter(faces.at(it3)))<1) {//the two face are not far to each other
        d2 = getDistanceArea(polygon2.at(it2), faces.at(it3));
        if(d2>0) {
          for(int it1=0; it1<polygon1.size() && !check; it1++) {
            if(distance(getBaryCenter(polygon1.at(it1)), getBaryCenter(faces.at(it3)))<1) {//the two face are not far to each other
              d1 = getDistanceArea(polygon1.at(it1), faces.at(it3));
              if(d1>0) {
                check=true;
                d+=d1;
                count++;
              }
            }
          }
        }
      }
      
    }
  }
  assert(count<=faces.size());
  return d;
}

RationalPoint getBaryCenter(std::vector<RationalPoint> vP)
{
  Rational x=0, y=0;
  for(int it=0; it<vP.size(); it++) {
    RationalPoint p = vP.at(it);
    x += p.first;
    y += p.second;
  }
  RationalPoint c = RationalPoint(x/vP.size(),y/vP.size());
  return c;
}

Z2i::RealPoint getBaryCenter(std::vector<Z2i::Point> vP)
{
  double x=0, y=0;
  for(int it=0; it<vP.size(); it++) {
    x += vP.at(it)[0];
    y += vP.at(it)[1];
  }
  Z2i::RealPoint c = Z2i::RealPoint(x/vP.size(),y/vP.size());
  return c;
}

Line getLine(std::pair<RationalPoint, RationalPoint> p)
{
  RationalPoint p1 = p.first;
  RationalPoint p2 = p.second;
  return getLineFromPoints(p1, p2);
}

// Return the cross product AB x BC.
// The cross product is a vector perpendicular to AB
// and BC having length |AB| * |BC| * Sin(theta) and
// with direction given by the right-hand rule.
// For two vectors in the X-Y plane, the result is a
// vector with X and Y components 0 so the Z component
// gives the vector's length and direction.
double CrossProductLength(double Ax, double Ay, double Bx, double By, double Cx, double Cy)
{
  // Get the vectors' coordinates.
  double BAx = Ax - Bx;
  double BAy = Ay - By;
  double BCx = Cx - Bx;
  double BCy = Cy - By;
  // Calculate the Z coordinate of the cross product.
  return (BAx * BCy - BAy * BCx);
}

bool isConvex(std::vector<RationalPoint> points)
{
  // For each set of three adjacent points A, B, C,
  // find the cross product AB · BC. If the sign of
  // all the cross products is the same, the angles
  // are all positive or negative (depending on the
  // order in which we visit them) so the polygon
  // is convex.
  bool got_negative = false;
  bool got_positive = false;
  int num_points = points.size();
  int B, C;
  for (int A = 0; A < num_points; A++)
  {
    B = (A + 1) % num_points;
    C = (B + 1) % num_points;
    double Ax = getRealValue(points.at(A).first);
    double Ay = getRealValue(points.at(A).second);
    double Bx = getRealValue(points.at(B).first);
    double By = getRealValue(points.at(B).second);
    double Cx = getRealValue(points.at(C).first);
    double Cy = getRealValue(points.at(C).second);
    double cross_product = CrossProductLength(Ax,Ay,Bx,By,Cx,Cy);
    if (cross_product < 0)
      got_negative = true;
    else if (cross_product > 0)
      got_positive = true;
    if (got_negative && got_positive)
      return false;
  }
  // If we got this far, the polygon is convex.
  return true;
}

//Find the smallest positive integer V such that [XV] + 1 < YV
//Return V or -1 if V is not exist
int findSmallestValue(double X, double Y)
{
  int start=0;
  int end=ceil(fabs(1.0/(Y-X)));//1e6;
  for(int i=start; i<=end; i++)
  if((int(X*i)+1)<Y*i)
    return i;
  return -1;
}

//Approximage angle by Pythagore triangle (a,b,c) with a<b<c and a*a+b*b=c*c
//sin(B) = a/c and cos (B) = b/c
PythagoreanTriple convertAngle2Pythagore(double B, double e)
{
  PythagoreanTriple pt = std::make_tuple(0, 1, 1); //triple corresponding to 0 degree;
  if(B==0)
    return pt;
  int a, b, c;
  int quart=0;
  double angle=B;
  while(angle>M_PI/2){angle-=M_PI/2;quart++;};//Get <= PI/2 angle
  double aa,bb;
  if(fabs(angle-M_PI/2)<e) { //FIXME: Angle=kPI/2 !!!!!
    //std::make_tuple(1, 0, 1); //triple corresponding to 90 degree
    aa=1;
    bb=0;
    c=1;
  }
  else {
    double X=tan(angle-e)+1.0/cos(angle-e);
    double Y=tan(angle+e)+1.0/cos(angle+e);
    double Xp=1.0/tan(angle+e)+1.0/sin(angle+e);
    double Yp=1.0/tan(angle-e)+1.0/sin(angle-e);
    int V=findSmallestValue(X,Y);
    int Vp=findSmallestValue(Xp,Yp);
    if(V==-1 || Vp==-1)
      return pt; //error here, no triple is found !
    int U=int(X*V)+1;
    int Up=int(Xp*Vp)+1;
    if((U*U+V*V)<(Up*Up+Vp*Vp)) {//quart=0
      aa=2*U*V;
      bb=U*U-V*V;
      c=U*U+V*V;
      if(fabs(asin(aa/c)-angle)>e) {
        int tmp=aa;
        aa=bb;
        bb=tmp;
      }
    }
    else {
      aa=2*Up*Vp;
      bb=Up*Up-Vp*Vp;
      c=Up*Up+Vp*Vp;
    }
  }
  
  if(quart%4==0) {
    a=aa;//sinx
    b=bb;//cosx
  }
  else if(quart%4==1) {
    a=bb;//sin(x+pi/2)=cosx
    b=-aa;//cos(x+pi/2)=-sinx
    
  }
  else if(quart%4==2) {
    a=-aa;//sin(x+pi)=-sinx
    b=-bb;//cos(x+pi/2)=-cosx
    
  }
  else { //if(quart%4==3)
    a=-bb;//sin(x+pi/2)=-cosx
    b=aa;//cos(x+pi/2)=sinx
  }
  return std::make_tuple(a, b, c);
}

Z2i::RealPoint getFaceCenter(const std::vector<Z2i::RealPoint>& face)
{
  double x=0, y=0;
  int closed = 0;
  if(face.front()==face.back()) //Ingore last vertex as the face is closed
    closed = 1;
  for(int it=0; it<face.size()-closed; it++) {
    Z2i::RealPoint p = face.at(it);
    x += p[0];
    y += p[1];
  }
  Z2i::RealPoint c = Z2i::RealPoint(x/(face.size()-closed),y/(face.size()-closed));
  return c;
}

std::vector<Z2i::RealPoint> getFaceCenter(const std::vector<std::vector<Z2i::RealPoint> >& face)
{
  std::vector<Z2i::RealPoint> centers;
  for(int it=0; it<face.size(); it++) {
    Z2i::RealPoint c = getFaceCenter(face.at(it));
    centers.push_back(c);
  }
  return centers;
}

RationalPoint getFaceCenter(const std::vector<RationalPoint>& face)
{
  Rational x=0, y=0;
  int closed = 0;
  if(isEqual(face.front(), face.back())) //Ingore last vertex as the face is closed
    closed = 1;
  for(int it=0; it<face.size()-closed; it++) {
    RationalPoint p = face.at(it);
    x += p.first;
    y += p.second;
  }
  RationalPoint c = RationalPoint(x/(face.size()-closed),y/(face.size()-closed));
  return c;
}

std::vector<RationalPoint> getFaceCenter(const std::vector<std::vector<RationalPoint> >& face)
{
  std::vector<RationalPoint> centers;
  for(int it=0; it<face.size(); it++) {
    RationalPoint c = getFaceCenter(face.at(it));
    centers.push_back(c);
  }
  return centers;
}

std::pair<RationalPoint,RationalPoint> getBoundingBox(const std::vector<RationalPoint>& points)
{
  if(points.size()==0)
    return std::make_pair(RationalPoint(0,0), RationalPoint(0,0));
  Rational minx=points.front().first, maxx=points.front().second;
  Rational miny=points.front().first, maxy=points.front().second;
  for(int it=1; it<points.size(); it++) {
    //std::cout<<"Point="<<points.at(it)<<std::endl;
    if(minx>points.at(it).first)
      minx=points.at(it).first;
    if(maxx<points.at(it).first)
      maxx=points.at(it).first;
    if(miny>points.at(it).second)
      miny=points.at(it).second;
    if(maxy<points.at(it).second)
      maxy=points.at(it).second;
    //std::cout<<"minx="<<minx<<", miny="<<miny<<" and maxx="<<maxx<<", maxy="<<maxy<<std::endl;
  }
  RationalPoint p1 = RationalPoint(minx,miny);
  RationalPoint p2 = RationalPoint(maxx,maxy);
  return std::make_pair(p1, p2);
}

std::pair<RationalPoint,RationalPoint> getBoundingBox(const std::vector<std::pair<RationalPoint,RationalPoint> >& points)
{
  std::vector<RationalPoint> vecPoints;
  for(int it=0; it<points.size(); it++) {
    vecPoints.push_back(points.at(it).first);
    vecPoints.push_back(points.at(it).second);
  }
  return getBoundingBox(vecPoints);
}

std::pair<Z2i::Point,Z2i::Point> getBoundingBox(const std::vector<Z2i::Point>& points)
{
  int minx=points.front()[0], maxx=points.front()[0];
  int miny=points.front()[1], maxy=points.front()[1];
  for(int it=1; it<points.size(); it++) {
    //std::cout<<"Point="<<points.at(it)<<std::endl;
    if(minx>points.at(it)[0])
      minx=points.at(it)[0];
    if(maxx<points.at(it)[0])
      maxx=points.at(it)[0];
    if(miny>points.at(it)[1])
      miny=points.at(it)[1];
    if(maxy<points.at(it)[1])
      maxy=points.at(it)[1];
    //std::cout<<"minx="<<minx<<", miny="<<miny<<" and maxx="<<maxx<<", maxy="<<maxy<<std::endl;
  }
  Z2i::Point p1 = Z2i::Point(minx,miny);
  Z2i::Point p2 = Z2i::Point(maxx,maxy);
  return std::make_pair(p1, p2);
}

std::pair<Z2i::RealPoint,Z2i::RealPoint> getBoundingBox(const std::vector<Z2i::RealPoint>& points)
{
  double minx=points.front()[0], maxx=points.front()[0];
  double miny=points.front()[1], maxy=points.front()[1];
  for(int it=1; it<points.size(); it++) {
    if(minx>points.at(it)[0])
      minx=points.at(it)[0];
    if(maxx<points.at(it)[0])
      maxx=points.at(it)[0];
    if(miny>points.at(it)[1])
      miny=points.at(it)[1];
    if(maxy<points.at(it)[1])
      maxy=points.at(it)[1];
  }
  Z2i::RealPoint p1 = Z2i::Point(minx,miny);
  Z2i::Point p2 = Z2i::Point(maxx,maxy);
  return std::make_pair(p1, p2);
}

std::vector<Z2i::Point> getDigitalSetPoint(Z2i::DigitalSet shape_set)
{
  std::vector<Z2i::Point> point_set;
  for(auto it=shape_set.begin(); it!=shape_set.end(); it++)
  point_set.push_back(*it);
  return point_set;
}

Z2i::DigitalSet getDigitalSet(std::vector<Z2i::Point> point_set)
{
  std::pair<Z2i::Point, Z2i::Point> bb = getBoundingBox(point_set);
  Z2i::Point p1 = bb.first - Z2i::Point(1,1);
  Z2i::Point p2 = bb.second + Z2i::Point(1,1);
  Z2i::Domain domain( p1, p2 );
  Z2i::DigitalSet shape_set( domain );
  for(int it=0; it<point_set.size(); it++)
  shape_set.insertNew(point_set.at(it));
  return shape_set;
}

std::ostream& operator<<(std::ostream& os, const RationalPoint p)
{
  os<<"("<<p.first<<","<<p.second<<")";
  return os;
}
