#include "AffineTransform.h"

AffineTransform::AffineTransform (Rational v11, Rational v12, Rational v21, Rational v22, Rational x, Rational y)
{
  this->a11 = v11;
  this->a12 = v12;
  this->a21 = v21;
  this->a22 = v22;
  this->tx = x;
  this->ty = y;
}

AffineTransform::AffineTransform (double v11, double v12, double v21, double v22)
{
  this->a11 = Rational(v11);
  this->a12 = Rational(v12);
  this->a21 = Rational(v21);
  this->a22 = Rational(v22);
  this->tx = Rational(0);
  this->ty = Rational(0);
}

AffineTransform::AffineTransform (double v11, double v12, double v21, double v22, double x, double y)
{
  this->a11 = Rational(v11);
  this->a12 = Rational(v12);
  this->a21 = Rational(v21);
  this->a22 = Rational(v22);
  this->tx = Rational(x);
  this->ty = Rational(y);
}

AffineTransform::AffineTransform (double v11, double v12, double v21, double v22, int t1x, int t1y, int t2x, int t2y)
{
  assert(t1y!=0 && t2y!=0);
  this->a11 = Rational(v11);
  this->a12 = Rational(v12);
  this->a21 = Rational(v21);
  this->a22 = Rational(v22);
  this->tx = Rational(t1x,t1y);
  this->ty = Rational(t2x,t2y);
}

AffineTransform AffineTransform::setCoeff(Rational v11, Rational v12, Rational v21, Rational v22)
{
  this->a11 = v11;
  this->a12 = v12;
  this->a21 = v21;
  this->a22 = v22;
  return *this;
}

AffineTransform AffineTransform::setCoeff(int v11, int v12, int v21, int v22)
{
  this->a11 = Rational(v11);
  this->a12 = Rational(v12);
  this->a21 = Rational(v21);
  this->a22 = Rational(v22);
  return *this;
}

AffineTransform AffineTransform::setCoeff(double v11, double v12, double v21, double v22)
{
  this->a11 = Rational(v11);
  this->a12 = Rational(v12);
  this->a21 = Rational(v21);
  this->a22 = Rational(v22);
  return *this;
}

AffineTransform AffineTransform::setTranslation(int x, int y)
{
  this->tx = Rational(x);
  this->ty = Rational(y);
  return *this;
}

AffineTransform AffineTransform::setTranslation(double x, double y)
{
  this->tx = Rational(x);
  this->ty = Rational(y);
  return *this;
}

AffineTransform AffineTransform::setTranslation(Rational x, Rational y)
{
  this->tx = x;
  this->ty = y;
  return *this;
}

AffineTransform AffineTransform::setTranslation(int t1x, int t1y, int t2x, int t2y)
{
  this->tx = Rational(t1x,t1y);
  this->ty = Rational(t2x,t2y);
  return *this;
}

RationalPoint AffineTransform::transform(RationalPoint p)
{
  RationalPoint tp;
  Rational x = p.first;
  Rational y = p.second;
  tp.first = this->a11*x + this->a12*y + this->tx;
  tp.second = this->a21*x + this->a22*y + this->ty;
  //std::cout<<p<<" and "<<tp<<std::endl;
  return tp;
}

RationalPoint AffineTransform::transform(Z2i::Point p)
{
  RationalPoint rp = RationalPoint(Rational(p[0]),Rational(p[1]));
  return transform(rp);
}

RationalPoint transformPoint(RationalPoint p, AffineTransform t)
{
  return t.transform(p);
}

//Ax + By + C = 0
//x * Cos(Theta) + y * Sin(Theta) - p = 0
//p = -C/Sqrt(A^2 + B^2)
std::pair<RationalPoint,RationalPoint> AffineTransform::transform(RationalPoint p1, RationalPoint p2)
{
  RationalPoint tp1 = this->transform(p1);
  RationalPoint tp2 = this->transform(p2);
  std::pair<RationalPoint,RationalPoint> tp = std::make_pair(tp1, tp2);
  return tp;
  //Line l = getLineFromPoints(tp1, tp2);
  //return l;
}

std::pair<RationalPoint,RationalPoint> AffineTransform::transform(std::pair<RationalPoint,RationalPoint> p)
{
  RationalPoint p1 = p.first;
  RationalPoint p2 = p.second;
  return this->transform(p1,p2);
}

std::pair<RationalPoint,RationalPoint> transformSegment(RationalPoint p1, RationalPoint p2, AffineTransform t)
{
  return t.transform(p1,p2);
}

std::pair<RationalPoint,RationalPoint> transformSegment(std::pair<RationalPoint,RationalPoint> p, AffineTransform t)
{
  return t.transform(p);
}

Complex transform(const Complex& c, AffineTransform t)
{
  Complex tc(c);
  //Domain
  Z2i::Point bl = tc.domain.first; //bottom-left
  Z2i::Point tr = tc.domain.second;//top-right
  Z2i::Point br = Z2i::Point(tr[0],bl[1]); //bottom-right
  Z2i::Point tl = Z2i::Point(bl[0],tr[1]);//top-left
  RationalPoint ttl = t.transform(tl);
  RationalPoint tbr = t.transform(br);
  RationalPoint ttr = t.transform(tr);
  RationalPoint tbl = t.transform(bl);
  Rational minX = std::min(ttl.first,std::min(tbr.first,std::min(ttr.first,tbl.first)));
  Rational maxX = std::max(ttl.first,std::max(tbr.first,std::max(ttr.first,tbl.first)));
  Rational minY = std::min(ttl.second,std::min(tbr.second,std::min(ttr.second,tbl.second)));
  Rational maxY = std::max(ttl.second,std::max(tbr.second,std::max(ttr.second,tbl.second)));
  Z2i::Point trp1 = Z2i::Point(floor(getRealValue(minX)-0.5),floor(getRealValue(minY)-0.5));
  Z2i::Point trp2 = Z2i::Point(ceil(getRealValue(maxX)+0.5),ceil(getRealValue(maxY)+0.5));
  tc.domain.first = trp1;
  tc.domain.second = trp2;
  //Transform grid
  for(int it=0; it<c.vertical_lines.size(); it++) {
    std::pair<RationalPoint,RationalPoint> tl = transformSegment(c.vertical_lines.at(it), t);
    tc.vertical_lines.at(it) = tl;
  }
  for(int it=0; it<c.horizontal_lines.size(); it++) {
    std::pair<RationalPoint,RationalPoint> tl = transformSegment(c.horizontal_lines.at(it), t);
    tc.horizontal_lines.at(it) = tl;
  }
  //Transform list_vertex
  for(int it=0; it<c.list_vertex.size(); it++) {
    RationalPoint p = c.list_vertex.at(it);
    RationalPoint tp = transformPoint(p,t);
    tc.list_vertex.at(it) = tp;
  }
  //Cell and star operators
  tc.cell_vertex = c.cell_vertex;
  tc.cell_edge = c.cell_edge;
  tc.cell_face_vertex = c.cell_face_vertex;
  tc.cell_face_edge = c.cell_face_edge;
  tc.star_vertex_edge = c.star_vertex_edge;
  tc.star_vertex_face = c.star_vertex_face;
  tc.star_edge = c.star_edge;
  tc.star_face = c.star_face;
  
  return tc;
}

std::ostream& operator<<(std::ostream& os, const AffineTransform& t)
{
  os << "Print AffineTransform "<<std::endl;
  os << "Affine matrix: a11=" <<t.a11 <<", a12=" <<t.a12<<", a21=" <<t.a21 <<", a22=" <<t.a22<<std::endl;
  os << "Translation: tx=" <<t.tx <<", ty=" <<t.ty<<std::endl;
  return os;
}
