#include "Complex.h"
#include "Graph.hpp"

//Get edges of a face being (non) boundary wrt object's cells
bool isSameEdge(std::pair<int,int> edge1, std::pair<int,int> edge2)
{
  if(edge1.first==edge2.first && edge1.second==edge2.second)
    return true;
  if(edge1.first==edge2.second && edge1.second==edge2.first)
    return true;
  return false;
}

bool isBelongFace(std::vector<int> face, std::pair<int,int> edge)
{
  std::pair<int,int> e;
  for(int it=0; it<face.size()-1; it++) {
    e = std::make_pair(face.at(it), face.at(it+1));
    if(isSameEdge(e, edge))
      return true;
  }
  return false;
}

std::vector<int> getBelongFaces(std::vector<std::vector<int> > faces, std::pair<int,int> edge)
{
  std::vector<int> f;
  for(int it=0; it<faces.size(); it++)
  if(isBelongFace(faces.at(it), edge))
    f.push_back(it);
  return f;
}

std::vector<int> getBelongFaces(std::vector<std::vector<int> > faces, int id_face, std::pair<int,int> edge)
{
  std::vector<int> f;
  for(int it=0; it<faces.size(); it++) {
    if(id_face!=it && isBelongFace(faces.at(it), edge))
      f.push_back(it);
  }
  return f;
}

Complex::Complex(const Complex & c)
{
  this->domain = c.domain;
  //Lines
  for(int it=0; it<c.horizontal_lines.size(); it++)
  this->horizontal_lines.push_back(c.horizontal_lines.at(it));
  for(int it=0; it<c.vertical_lines.size(); it++)
  this->vertical_lines.push_back(c.vertical_lines.at(it));
  //Complex
  for(int it=0; it<c.list_vertex.size(); it++)
  this->list_vertex.push_back(c.list_vertex.at(it));
  for(int it=0; it<c.vertex.size(); it++)
  this->vertex.push_back(c.vertex.at(it));
  for(int it=0; it<c.edge.size(); it++)
  this->edge.push_back(c.edge.at(it));
  for(int it=0; it<c.face.size(); it++)
  this->face.push_back(c.face.at(it));
  //Cell and star operator
  this->cell_vertex = c.cell_vertex;
  this->cell_edge = c.cell_edge;
  this->cell_face_vertex = c.cell_face_vertex;
  this->cell_face_edge = c.cell_face_edge;
  this->star_vertex_edge = c.star_vertex_edge;
  this->star_vertex_face = c.star_vertex_face;
  this->star_edge = c.star_edge;
  this->star_face = c.star_face;
}

Complex::Complex(const std::vector<std::vector<RationalPoint> >& faces)
{
  Rational minx=faces.at(0).front().first, maxx=faces.at(0).front().first;
  Rational miny=faces.at(0).front().second, maxy=faces.at(0).front().second;
  for(int it=0; it<faces.size(); it++) {
    this->addFace(faces.at(it),true,true);
    for(int it_bis=0; it_bis<faces.at(it).size(); it_bis++) {
      //std::cout<<"Vertex="<<faces.at(it).at(it_bis)<<std::endl;
      if(minx>faces.at(it).at(it_bis).first)
        minx=faces.at(it).at(it_bis).first;
      if(maxx<faces.at(it).at(it_bis).first)
        maxx=faces.at(it).at(it_bis).first;
      if(miny>faces.at(it).at(it_bis).second)
        miny=faces.at(it).at(it_bis).second;
      if(maxy<faces.at(it).at(it_bis).second)
        maxy=faces.at(it).at(it_bis).second;
    }
  }
  Z2i::Point p1 = Z2i::Point(floor(getRealValue(minx)-0.5)-BORDER,floor(getRealValue(miny)-0.5)-BORDER);
  Z2i::Point p2 = Z2i::Point(ceil(getRealValue(maxx)+0.5)+BORDER,ceil(getRealValue(maxy)+0.5)+BORDER);
  std::pair<Z2i::Point,Z2i::Point> d = std::make_pair(p1,p2);
  this->domain = d;
}

Complex::Complex(const std::vector<Z2i::Point>& points)
{
  if(points.size()!=0) {
    Z2i::Point p1, p2;
    int minx=points.front()[0], maxx=points.front()[0], miny=points.front()[1], maxy=points.front()[1];
    for(int it=0; it<points.size(); it++) {
      this->addCubicalFace(points.at(it));
      if(minx>points.at(it)[0])
        minx=points.at(it)[0];
      if(maxx<points.at(it)[0])
        maxx=points.at(it)[0];
      if(miny>points.at(it)[1])
        miny=points.at(it)[1];
      if(maxy<points.at(it)[1])
        maxy=points.at(it)[1];
    }
    this->init(Z2i::Point(minx-BORDER,miny-BORDER), Z2i::Point(maxx+BORDER,maxy+BORDER));
  }
}

Complex::Complex(const Z2i::Point point)
{
  this->addCubicalFace(point);
  this->init(Z2i::Point(point[0]-BORDER,point[1]-BORDER), Z2i::Point(point[0]+BORDER,point[1]+BORDER));
}

Complex::Complex(const std::vector<RationalPoint>& vertex, const std::vector<std::pair<RationalPoint,RationalPoint> >& edge, const std::vector<std::vector<RationalPoint> >& face)
{
  for(int it=0; it<vertex.size(); it++)
  this->addVertex(vertex.at(it));
  for(int it=0; it<edge.size(); it++)
  this->addEdge(edge.at(it));
  for(int it=0; it<face.size(); it++)
  this->addFace(face.at(it));
  
  Rational minx=this->list_vertex.front().first, maxx=this->list_vertex.front().first;
  Rational miny=this->list_vertex.front().second, maxy=this->list_vertex.front().second;
  for(int it=1; it<this->list_vertex.size(); it++) {
    if(minx>this->list_vertex.at(it).first)
      minx=this->list_vertex.at(it).first;
    if(maxx<this->list_vertex.at(it).first)
      maxx=this->list_vertex.at(it).first;
    if(miny>this->list_vertex.at(it).second)
      miny=this->list_vertex.at(it).second;
    if(maxy<this->list_vertex.at(it).second)
      maxy=this->list_vertex.at(it).second;
    
  }
  Z2i::Point p1 = Z2i::Point(floor(getRealValue(minx)-0.5)-BORDER,floor(getRealValue(miny)-0.5)-BORDER);
  Z2i::Point p2 = Z2i::Point(ceil(getRealValue(maxx)+0.5)+BORDER,ceil(getRealValue(maxy)+0.5)+BORDER);
  std::pair<Z2i::Point,Z2i::Point> d = std::make_pair(p1,p2);
  this->domain = d;
}

//init with limits of the grid space
void Complex::init(Z2i::Point p1, Z2i::Point p2)
{
  int w = p2[0]-p1[0];
  int h = p2[1]-p1[1];
  assert(w>0 && h>0);
  this->domain = std::make_pair(p1, p2);
  Rational half(1,2);
  RationalPoint pp1, pp2;
  //Line l;
  std::pair<RationalPoint,RationalPoint> l;
  //Horizontal lines : y = k + 1/2 =>  0x + y - (k+1/2) = 0
  for(int y=p1[1]-1; y<=p2[1]; y++) {
    //l = std::make_tuple(Rational(1),Rational(0),-Rational(y)-half);
    pp1 = RationalPoint(Rational(p1[0])-half,Rational(y)+half);
    pp2 = RationalPoint(Rational(p2[0])+half,Rational(y)+half);
    l = std::make_pair(pp1, pp2);
    horizontal_lines.push_back(l);
  }
  //Vertical lines : x = k + 1/2 =>  x + 0y - (k+1/2) = 0
  ///std::cout<<"Vertical lines"<<std::endl;
  for(int x=p1[0]-1; x<=p2[0]; x++) {
    //l = std::make_tuple(Rational(0),Rational(1),-Rational(x)-half);
    pp1 = RationalPoint(Rational(x)+half,Rational(p1[1])-half);
    pp2 = RationalPoint(Rational(x)+half,Rational(p2[1])+half);
    l = std::make_pair(pp1, pp2);
    vertical_lines.push_back(l);
  }
}

void Complex::initGrid(Z2i::Point p1, Z2i::Point p2)
{
  this->init(p1,p2);
  for(int x=p1[0]; x<=p2[0]; x++)
  for(int y=p1[1]; y<=p2[1]; y++) {
    Z2i::Point p = Z2i::Point(x,y);
    this->addCubicalFace(p);
    //this->idCC.push_back(0);//background CC
  }
}

void Complex::initGridFaces(Z2i::Point p1, Z2i::Point p2)
{
  this->domain = std::make_pair(p1, p2);
  for(int x=p1[0]; x<=p2[0]; x++)
  for(int y=p1[1]; y<=p2[1]; y++) {
    Z2i::Point p = Z2i::Point(x,y);
    this->addCubicalFace(p);
    //this->idCC.push_back(0);//background CC
  }
}

void Complex::initGridLines(std::vector<std::pair<RationalPoint,RationalPoint> > vlines, std::vector<std::pair<RationalPoint,RationalPoint> > hlines)
{
  std::pair<RationalPoint,RationalPoint> bb1=getBoundingBox(hlines);
  std::pair<RationalPoint,RationalPoint> bb2=getBoundingBox(vlines);
  RationalPoint rp1 = min(bb1.first, bb2.first);
  RationalPoint rp2 = min(bb1.second, bb2.second);
  Rational minX = rp1.first;
  Rational maxX = rp2.first;
  Rational minY = rp1.second;
  Rational maxY = rp2.second;
  Z2i::Point p1 = Z2i::Point(floor(getRealValue(minX)-0.5),floor(getRealValue(minY)-0.5));
  Z2i::Point p2 = Z2i::Point(ceil(getRealValue(maxX)+0.5),ceil(getRealValue(maxY)+0.5));
  this->domain = std::make_pair(p1, p2);
  this->vertical_lines = vlines;
  this->horizontal_lines = hlines;
}

void Complex::addGridLines(std::vector<std::pair<RationalPoint,RationalPoint> > vlines, std::vector<std::pair<RationalPoint,RationalPoint> > hlines)
{
  std::pair<RationalPoint,RationalPoint> bb1 = getBoundingBox(hlines);
  std::pair<RationalPoint,RationalPoint> bb2 = getBoundingBox(vlines);
  RationalPoint rp1 = min(bb1.first, bb2.first);
  RationalPoint rp2 = min(bb1.second, bb2.second);
  Rational minX = rp1.first;
  Rational maxX = rp2.first;
  Rational minY = rp1.second;
  Rational maxY = rp2.second;
  std::pair<Z2i::Point,Z2i::Point> bb = domain;
  int v1 = floor(getRealValue(minX)-0.5);
  int v2 = floor(getRealValue(minY)-0.5);
  Z2i::Point p1 = Z2i::Point(std::min(v1,domain.first[0]),std::min(v1,domain.first[1]));
  v1 = ceil(getRealValue(maxX)+0.5);
  v2 = ceil(getRealValue(maxY)+0.5);
  Z2i::Point p2 = Z2i::Point(std::max(v1,domain.second[0]),std::max(v1,domain.second[1]));
  this->domain = std::make_pair(p1, p2);
  for(size_t it=0; it<vlines.size(); it++) {
    this->vertical_lines.push_back(vlines.at(it));
  }
  for(size_t it=0; it<hlines.size(); it++) {
    this->vertical_lines.push_back(hlines.at(it));
  }
}

std::vector<Line> Complex::getVerticalLines() const
{
  std::vector<Line> lines;
  for(int it=0; it<vertical_lines.size(); it++) {
    Line l = getLine(vertical_lines.at(it));
    lines.push_back(l);
  }
  return lines;
}

std::vector<Line> Complex::getHorizontalLines() const
{
  std::vector<Line> lines;
  for(int it=0; it<horizontal_lines.size(); it++) {
    Line l = getLine(horizontal_lines.at(it));
    lines.push_back(l);
  }
  return lines;
}

std::vector<Line> Complex::getGridLines() const
{
  std::vector<Line> vlines = getVerticalLines();
  std::vector<Line> hlines = getHorizontalLines();
  std::vector<Line> lines;
  for(int it=0; it<horizontal_lines.size(); it++)
  lines.push_back(hlines.at(it));
  for(int it=0; it<vertical_lines.size(); it++)
  lines.push_back(vlines.at(it));
  return lines;
}

std::pair<int,int> Complex::getDomainSize() const
{
  Z2i::Point p1 = domain.first;
  Z2i::Point p2 = domain.second;
  double w = p2[0]-p1[0] + 1; // +1 border included
  double h = p2[1]-p1[1] + 1; // +1 border included
  return std::make_pair(w,h);
}

int Complex::getBorderSize() const
{
  Z2i::Point p1 = domain.first;
  Z2i::Point p2 = domain.second;
  double w = p2[0]-p1[0] + 1; // +1 border included
  double h = p2[1]-p1[1] + 1; // +1 border included
  return std::max(w,h);
}

int Complex::computeCellVertex (int idV)
{
  if(idV < vertex.size() && containElement(vertex, idV))
    return idV;
  return -1;
}

std::pair<int,int> Complex::computeCellEdge (int idE)
{
  std::pair<int,int> vertices (-1, -1);
  if(idE < edge.size()) {
    std::pair<int,int> e = edge.at(idE);
    int e1 = computeCellVertex(e.first);
    int e2 = computeCellVertex(e.second);
    vertices.first = e1;
    vertices.second = e2;
    if(e1==-1 && e2!=-1)
      vertices.first = e2;
  }
  return vertices;
}

std::vector<int> Complex::computeCellFaceVertex (int idF)
{
  std::vector<int> vertices;
  if(idF < face.size()) {
    std::vector<int> aFace = face.at(idF);
    for(int it=0; it<aFace.size()-1; it++) {
      int id = computeCellVertex(aFace.at(it));
      if(id!=-1)
        vertices.push_back(id);
    }
  }
  return vertices;
}

std::vector<int> Complex::computeCellFaceEdge (int idF)
{
  std::vector<int> edges;
  if(idF < face.size()) {
    std::vector<int> aFace = face.at(idF);
    for(int it=0; it<aFace.size()-1; it++) {
      std::pair<int, int> e = std::make_pair(aFace.at(it), aFace.at(it+1));
      int id = findEdege(this->edge, e);
      if(id!=-1)
        edges.push_back(id);
    }
  }
  return edges;
}

std::vector<int> Complex::computeStarVertexEdge (int idV)
{
  std::vector<int> edges;
  for(size_t it=0; it<edge.size(); it++) {
    if(edge.at(it).first==idV || edge.at(it).second==idV)
      edges.push_back(it);
  }
  return edges;
}

std::vector<int> Complex::computeStarVertexFace (int idV)
{
  std::vector<int> faces;
  for(size_t it=0; it<face.size(); it++) {
    if(containElement(face.at(it), idV))
      faces.push_back(it);
  }
  return faces;
}

std::pair<int, int> Complex::computeStarEdge (int idE)
{
  std::vector<int> f = findContainEdge(face, edge.at(idE));
  assert(f.size()<3);
  std::pair<int, int> faces (-1, -1);
  if(f.size()==2) {
    faces.first = f.front();
    faces.second = f.back();
  }
  if(f.size()==1) {
    faces.first = f.front();
    faces.second = -1;
  }
  return faces;
}

int Complex::computeStarFace (int idF)
{
  return idF;
}

void Complex::updateAllCell()
{
  for(size_t it=0; it<vertex.size(); it++)
    cell_vertex.at(it) = computeCellVertex(it);
  for(size_t it=0; it<edge.size(); it++)
    cell_edge.at(it) = computeCellEdge(it);
  for(size_t it=0; it<face.size(); it++) {
    cell_face_vertex.at(it) = computeCellFaceVertex(it);
    cell_face_edge.at(it) = computeCellFaceEdge(it);
  }
}
void Complex::updateAllStar()
{
  for(size_t it=0; it<vertex.size(); it++) {
    star_vertex_edge.at(it) = computeStarVertexEdge(it);
    star_vertex_face.at(it) = computeStarVertexFace(it);
  }
  for(size_t it=0; it<edge.size(); it++)
    star_edge.at(it) = computeStarEdge(it);
  for(size_t it=0; it<face.size(); it++)
    star_face.at(it) = computeStarFace(it);
}


void Complex::addVertex(RationalPoint v, bool updateCellStar)
{
  int id = findVertex(list_vertex, v);
  if(id==-1) { //the vertex is new
    list_vertex.push_back(v);
    id = list_vertex.size()-1;
    vertex.push_back(id);
    cell_vertex.push_back(id); //cell of vertex is itself
    star_vertex_edge.push_back(std::vector<int>());//empty vector
    star_vertex_face.push_back(std::vector<int>());//empty vector
  }
  if(!containElement(vertex,id)) { //the vertex is not in the list yet
    vertex.push_back(id);
    cell_vertex.push_back(id); //cell of vertex is itself
    star_vertex_edge.push_back(std::vector<int>());//empty vector
    star_vertex_face.push_back(std::vector<int>());//empty vector
  }
  //Compute the cell and star operator
  if(updateCellStar) {
    for(size_t id=0; id<vertex.size(); id++) {
      cell_vertex.at(id) = computeCellVertex(id);
      star_vertex_edge.at(id) = computeStarVertexEdge(id);
      star_vertex_face.at(id) = computeStarVertexFace(id);
    }
  }
  assert(vertex.size()==cell_vertex.size());
  assert(vertex.size()==star_vertex_edge.size());
  assert(vertex.size()==star_vertex_face.size());
}

int isExistEdge(const std::vector<std::pair<int,int> >& vector, std::pair<int,int> e)
{
  std::pair<int,int> e_inv = std::make_pair(e.second, e.first);
  auto it = std::find(vector.begin(), vector.end(), e);
  auto it_inv = std::find(vector.begin(), vector.end(), e_inv);
  if (it != vector.end())
    return it-vector.begin();
 if(it_inv != vector.end())
   return it_inv-vector.begin();
  return -1;
}


void Complex::addEdge(RationalPoint e1, RationalPoint e2, bool addVertex, bool updateCellStar)
{
  int id1 = findVertex(list_vertex, e1);
  int id2 = findVertex(list_vertex, e2);
  if (id1==-1) {
    list_vertex.push_back(e1);
    id1 = list_vertex.size()-1;
  }
  if (id2==-1) {
    list_vertex.push_back(e2);
    id2 = list_vertex.size()-1;
  }
  assert(id1!=id2);
  std::pair<int, int> e = std::make_pair(id1, id2);
  int id = isExistEdge(edge,e);
  if(id==-1) { //the edge is new
    edge.push_back(e);
    id = edge.size() - 1;
    cell_edge.push_back(std::make_pair(-1,-1)); //empty
    star_edge.push_back(std::make_pair(-1,-1)); //empty
    if(addVertex) {
      this->addVertex(e1);//vertex.push_back(id1);
      this->addVertex(e2);//vertex.push_back(id2);
    }
  }
  //Compute the cell and star operator
  if(updateCellStar) {
    for(size_t id=0; id<edge.size(); id++) {
      cell_edge.at(id) = computeCellEdge(id);
      star_edge.at(id) = computeStarEdge(id);
    }
  }
}

void Complex::addEdge(std::pair<RationalPoint, RationalPoint> e, bool addVertex, bool updateCellStar)
{
  RationalPoint e1 = e.first;
  RationalPoint e2 = e.second;
  this->addEdge(e1, e2, addVertex, updateCellStar);
}

bool isExistFace(const std::vector<std::vector<int> >& vector, std::vector<int> f)
{
  for(int it=0; it<vector.size(); it++)
  if(isSameVector(vector.at(it), f))
    return true;
  return false;
}

void Complex::addFace(const std::vector<RationalPoint>& f, bool addEdge, bool addVertex, bool updateCellStar)
{
  std::vector<int> aFace;
  for(int it=0; it<f.size(); it++) {
    RationalPoint p = f.at(it);
    int id = findVertex(list_vertex, p);
    if (id==-1) {
      list_vertex.push_back(p);
      id = list_vertex.size()-1;
    }
    aFace.push_back(id);
  }
  
  //Verify if the face is a closed polygon, it is not then close it
  if(!isEqual(f.front(),f.back()))
    aFace.push_back(aFace.front());
  
  //Verify if the face is not in the list yet
  int id = -1;
  if(!isExistFace(face,aFace)) { //this is a new face
    face.push_back(aFace);
    id = face.size() - 1;
    if(addVertex) {
      for(int it=0; it<aFace.size(); it++)
      this->addVertex(list_vertex.at(aFace.at(it)));
    }
    if(addEdge) {
      for(int it=0; it<aFace.size()-1; it++) {
        RationalPoint e1 = list_vertex.at(aFace.at(it));
        RationalPoint e2 = list_vertex.at(aFace.at(it+1));
        this->addEdge( e1, e2);
      }
    }
    cell_face_vertex.push_back(std::vector<int>());//empty
    cell_face_edge.push_back(std::vector<int>());//empty
    star_face.push_back(-1);//empty
  }

  //Compute the cell and star operator
  if(updateCellStar) {
    for(size_t id=0; id<face.size(); id++) {
      cell_face_vertex.at(id) = computeCellFaceVertex(id);
      cell_face_edge.at(id) = computeCellFaceEdge(id);
      star_face.at(id) = computeStarFace(id);
    }
  }
  assert(face.size()==cell_face_vertex.size());
  assert(face.size()==cell_face_edge.size());
  assert(face.size()==star_face.size());
}

//Add a face by its center
void Complex::addCubicalFace(Z2i::Point p)
{
  Rational half(1,2);
  Rational x(p[0]);
  Rational y(p[1]);
  Rational v1x = x - half;
  Rational v1y = y + half;
  RationalPoint v1(std::make_pair(v1x,v1y));
  Rational v2x = x + half;
  Rational v2y = y + half;
  RationalPoint v2(std::make_pair(v2x,v2y));
  Rational v3x = x + half;
  Rational v3y = y - half;
  RationalPoint v3(std::make_pair(v3x,v3y));
  Rational v4x = x - half;
  Rational v4y = y - half;
  RationalPoint v4(std::make_pair(v4x,v4y));
  
  addVertex(v1);
  addVertex(v2);
  addVertex(v3);
  addVertex(v4);
  addEdge(v1, v2);
  addEdge(v2, v3);
  addEdge(v3, v4);
  addEdge(v4, v1);
  std::vector<RationalPoint> f;
  f.push_back(v1);
  f.push_back(v2);
  f.push_back(v3);
  f.push_back(v4);
  f.push_back(v1);
  addFace(f);
}

RationalPoint Complex::getVertex(int index) const
{
  return this->list_vertex.at(index);
}
std::vector<RationalPoint> Complex::getVertex(const std::vector<int>& index) const
{
  std::vector<RationalPoint> points;
  for(int it=0; it<index.size(); it++)
  points.push_back(this->getVertex(it));
  return points;
}
std::pair<RationalPoint, RationalPoint> Complex::getEdgeVertices(int index) const
{
  std::pair<int, int> e = this->edge.at(index);
  RationalPoint e1 = this->list_vertex.at(e.first);
  RationalPoint e2 = this->list_vertex.at(e.second);
  return std::make_pair(e1,e2);
}
std::vector<std::pair<RationalPoint, RationalPoint> > Complex::getEdgeVertices(const std::vector<int>& index) const
{
  std::vector<std::pair<RationalPoint, RationalPoint> > edges;
  for(int it=0; it<index.size(); it++)
  edges.push_back(this->getEdgeVertices(it));
  return edges;
}

std::vector<RationalPoint> Complex::getFaceVertices(int index) const
{
  std::vector<int> f = this->face.at(index);
  std::vector<RationalPoint> face;
  for(int it=0; it<f.size(); it++)
    face.push_back(this->list_vertex.at(f.at(it)));
  return face;
}

double Complex::getFaceArea(int index) const
{
  return areaPolygon(this->getFaceVertices(index));
}

std::vector<std::vector<RationalPoint> > Complex::getFaceVertices(const std::vector<int>& index) const
{
  std::vector<std::vector<RationalPoint> > faces;
  for(int it=0; it<index.size(); it++)
  faces.push_back(this->getFaceVertices(index.at(it)));
  return faces;
}

std::vector<std::vector<RationalPoint> > Complex::getFaceVertices() const
{
  std::vector<std::vector<RationalPoint> > faces;
  for(int it=0; it<this->face.size(); it++)
  faces.push_back(this->getFaceVertices(it));
  return faces;
}

int Complex::getIdFace(RationalPoint p)
{
  for(size_t it=0; it<this->face.size(); it++)
    if(isInsidePolygon(this->getFaceVertices(it), p))
      return it;
  return -1;
}

RationalPoint Complex::getFaceCenter(const std::vector<RationalPoint>& face) const
{
  double x=0, y=0;
  int closed = 0;
  if(isEqual(face.front(), face.back())) //Ingore last vertex as the face is closed
    closed = 1;
  for(int it=0; it<face.size() - closed; it++) {
    RationalPoint p = face.at(it);
    x += getRealValue(p.first);
    y += getRealValue(p.second);
  }
  Z2i::RealPoint c = Z2i::RealPoint(x/(face.size()-closed),y/(face.size()-closed));
  return RationalPoint(Rational(c[0]), Rational(c[1]));
}

RationalPoint Complex::getFaceCenter(int index) const
{
  std::vector<int> f = face.at(index);
  std::vector<RationalPoint> fp;
  for(int it=0; it<f.size(); it++) {
    RationalPoint p = list_vertex.at(f.at(it));
    fp.push_back(p);
  }
  return this->getFaceCenter(fp);
}

std::vector<int> Complex::getNeighbourFace(int index) const
{
  std::vector<std::pair<int,int> > vE = getNonBoundaryEdges(index);
  std::vector<int> vF;
  for(size_t it=0; it<vE.size(); it++) {
    std::pair<int,int> e = vE.at(it);
    //find all face contain e
    for(size_t it_bis = 0; it_bis<this->face.size(); it_bis++) {
      if(it_bis!= index && isBelongFace(this->face.at(it_bis), e))
        vF.push_back(it_bis);
    }
  }
  return vF;
}

std::vector<std::vector<int> > Complex::getAllNeighbourFace() const
{
  std::vector<std::vector<int> > vN;
  for(size_t it = 0; it<this->face.size(); it++) {
    std::vector<int> vF = getNeighbourFace(it);
    vN.push_back(vF);
  }
  return vN;
}

std::vector<RationalPoint> Complex::getFaceCenter() const
{
  std::vector<RationalPoint> center;
  for(int it=0; it<face.size(); it++) {
    RationalPoint p = this->getFaceCenter(it);
    center.push_back(p);
  }
  return center;
}

std::vector<RationalPoint> intersectionVertex(const std::vector<RationalPoint>& face, const Complex& grid)
{
  //Compute the intersections of the face wrt vertical and horizontal lines
  std::vector<std::pair<RationalPoint,RationalPoint> > segment;
  std::vector<RationalPoint> vecVertex;
  
  std::vector<RationalPoint> vertical;
  for(int it_bis=0; it_bis<grid.vertical_lines.size(); it_bis++) { //For each vertical line
    vertical = intersection(face,grid.vertical_lines.at(it_bis));
    if(vertical.size()!=0) {
      segment.push_back(grid.vertical_lines.at(it_bis));
      for(int it_tris=0; it_tris<vertical.size(); it_tris++)
      if(findVertex(vecVertex, vertical.at(it_tris))==-1) //verify if the vertex is not in the list yet
        vecVertex.push_back(vertical.at(it_tris));
    }
  }
  std::vector<RationalPoint> horizontal;
  for(int it_bis=0; it_bis<grid.horizontal_lines.size(); it_bis++) { //For each horizontal line
    horizontal = intersection(face,grid.horizontal_lines.at(it_bis));
    if(horizontal.size()!=0) {
      segment.push_back(grid.horizontal_lines.at(it_bis));
      for(int it_tris=0; it_tris<horizontal.size(); it_tris++) {
        if(findVertex(vecVertex, horizontal.at(it_tris))==-1) //verify if the vertex is not in the list yet
          vecVertex.push_back(horizontal.at(it_tris));
      }
    }
  }
  //Intersection between horizontal and vertical segments inside face
  if(vecVertex.size()>0 /*&& vecVertex.size()<8*/) {
    std::vector<RationalPoint> hv;
    for(int i=0; i<segment.size()-1;i++)
    for(int j=i+1; j<segment.size();j++) {
      RationalPoint p11=segment.at(i).first;
      RationalPoint p12=segment.at(i).second;
      RationalPoint p21=segment.at(j).first;
      RationalPoint p22=segment.at(j).second;
      if(isIntersect(p11,p12,p21,p22))
        hv.push_back(intersection(p11,p12,p21,p22));
    }
    for(int it_bis=0; it_bis<hv.size(); it_bis++) {
      //std::cout<<it_bis<<":"<<hv.at(it_bis)<<"=>"<<isInsidePolygon(face, hv.at(it_bis))<<" and "<<findVertex(vecVertex, hv.at(it_bis))<<std::endl;
      if(isInsidePolygon(face, hv.at(it_bis)) && findVertex(vecVertex, hv.at(it_bis))==-1) //verify if the vertex is not in the list yet
        vecVertex.push_back(hv.at(it_bis));
    }
  }
  //Add vertex's face
  for(int it_bis=0; it_bis<face.size(); it_bis++) {
    if(findVertex(vecVertex, face.at(it_bis))==-1) //verify if the vertex is not in the list yet
      vecVertex.push_back(face.at(it_bis));
  }
  return vecVertex;
}

std::vector<std::pair<RationalPoint,RationalPoint> > intersectionEdge(const std::vector<RationalPoint>& face, const Complex& grid)
{
  //Compute the intersections of each face wrt vertical and horizontal lines
  std::vector<std::pair<RationalPoint,RationalPoint> > segment;
  std::vector<RationalPoint> vertex;
  std::vector<RationalPoint> vertical;
  for(int it_bis=0; it_bis<grid.vertical_lines.size(); it_bis++) { //For each vertical line
    vertical = intersection(face,grid.vertical_lines.at(it_bis));
    if(vertical.size()!=0) {
      segment.push_back(grid.vertical_lines.at(it_bis));
      for(int it_tris=0; it_tris<vertical.size(); it_tris++)
      if(findVertex(vertex, vertical.at(it_tris))==-1) //verify if the vertex is not in the list yet
        vertex.push_back(vertical.at(it_tris));
    }
  }
  std::vector<RationalPoint> horizontal;
  for(int it_bis=0; it_bis<grid.horizontal_lines.size(); it_bis++) { //For each horizontal line
    horizontal = intersection(face,grid.horizontal_lines.at(it_bis));
    if(horizontal.size()!=0) {
      segment.push_back(grid.horizontal_lines.at(it_bis));
      for(int it_tris=0; it_tris<horizontal.size(); it_tris++)
      if(findVertex(vertex, horizontal.at(it_tris))==-1) //verify if the vertex is not in the list yet
        vertex.push_back(horizontal.at(it_tris));
    }
  }
  //Intersection between horizontal and vertical segments inside face
  if(vertex.size()>0 /* && vertex.size()<8 */) {
    std::vector<RationalPoint> hv;
    for(int i=0; i<segment.size()-1;i++)
    for(int j=i+1; j<segment.size();j++) {
      RationalPoint p11=segment.at(i).first;
      RationalPoint p12=segment.at(i).second;
      RationalPoint p21=segment.at(j).first;
      RationalPoint p22=segment.at(j).second;
      if(isIntersect(p11,p12,p21,p22))
        hv.push_back(intersection(p11,p12,p21,p22));
    }
    for(int it_bis=0; it_bis<hv.size(); it_bis++) {
      if(isInsidePolygon(face, hv.at(it_bis)) && findVertex(vertex, hv.at(it_bis))==-1) //verify if the vertex is not in the list yet
        vertex.push_back(hv.at(it_bis));
    }
  }
  //Add vertex's face
  for(int it_bis=0; it_bis<face.size(); it_bis++)
  if(findVertex(vertex, face.at(it_bis))==-1) //verify if the vertex is not in the list yet
    vertex.push_back(face.at(it_bis));
  //Add edge's face
  for(int it_bis=0; it_bis<face.size()-1; it_bis++) {
    RationalPoint p1 = face.at(it_bis);
    RationalPoint p2 = face.at(it_bis+1);
    segment.push_back(std::make_pair(p1, p2));
  }
  
  //Sort points of each segment for edges
  std::vector<std::pair<RationalPoint,RationalPoint> > vecEdge;
  for(int it_bis=0; it_bis<segment.size(); it_bis++) { //Each segment
    //Find all intersections (vertices) belong to the segment
    std::vector<RationalPoint> vecP;
    RationalPoint p1 = segment.at(it_bis).first;
    RationalPoint p2 = segment.at(it_bis).second;
    for(int it_tris=0; it_tris<vertex.size(); it_tris++) { //For all vertices belonging to face
      RationalPoint p = vertex.at(it_tris);
      if(isCollinear(p1,p,p2) && onSegment(p1,p,p2))
        vecP.push_back(p);
    }
    
    if(vecP.size()>1) { //consider only segment with at least two points
      std::vector<RationalPoint> vecP_sorted = sortCollinearPoints(vecP);
      for(int it_tris=0; it_tris<vecP_sorted.size()-1; it_tris++) {
        RationalPoint p1 = vecP_sorted.at(it_tris);
        RationalPoint p2 = vecP_sorted.at(it_tris+1);
        std::pair<RationalPoint,RationalPoint> e = std::make_pair(p1, p2);
        vecEdge.push_back(e);
      }
    }
  }
  return vecEdge;
}

bool containCollinearPoint(const std::vector<RationalPoint>& vertex, const std::vector<int>& cycle)
{
  //non of three points are collinear
  for(int it=0; it<cycle.size()-2; it++) {
    RationalPoint p1 = vertex.at(cycle.at(it));
    RationalPoint p2 = vertex.at(cycle.at(it+1));
    RationalPoint p3 = vertex.at(cycle.at(it+2));
    if(isCollinear(p1,p2,p3))
      return true;
  }
  return false;
}
bool containInsidePoint(const std::vector<RationalPoint>& vertex, const std::vector<int>& cycle)
{
  std::vector<bool> belong_cycle(vertex.size(),false);
  std::vector<RationalPoint> polygon;
  for(int it=0; it<cycle.size(); it++) {
    belong_cycle.at(cycle.at(it))=true;
    polygon.push_back(vertex.at(cycle.at(it)));
  }
  //The polygon does not contain any other point
  for(int it=0; it<belong_cycle.size(); it++) {
    if(belong_cycle.at(it)==false) {
      RationalPoint p = vertex.at(it);
      if(isInsidePolygon(polygon,p)==true)//if(isInterieurPolygon(polygon,p)==true)
        return true;
    }
  }
  return false;
}
bool isValidCycle(const std::vector<RationalPoint>& vertex, const std::vector<int>& cycle)
{
  std::vector<RationalPoint> polygon;
  for(int it=0; it<cycle.size(); it++)
  polygon.push_back(vertex.at(cycle.at(it)));
  assert(polygon.front()==polygon.back());//closed polygon !
  if(isConvex(polygon)==false) //The polygon must be convex
    return false;
  if(containCollinearPoint(vertex,cycle)==true) //doest not contain collinear points
    return false;
  //FIXME
  //if(containInsidePoint(vertex,cycle)==true) //doest not contain any other point
  //    return false;
  return true;
}

//return index of including cycle and -1 if not exist
int isIncludedCycle(std::vector<std::vector<int> > vecCycles, std::vector<int> aCycle) {
  int count, id;
  for(size_t it=0; it<vecCycles.size(); it++) {
    if(vecCycles.at(it).size()>aCycle.size()) { //must have a bigger size to include the cycle
      count=0;
      for(size_t it_bis=0; it_bis<aCycle.size()-1; it_bis++) //count the number of containing vertices
      if(!containElement(vecCycles.at(it), aCycle.at(it_bis)))
        count++;
      if(count==aCycle.size()-1)
        return it;
    }
  }
  return -1;
}

//verfiy if c1 is sous sequence of c2
bool isIncludedCycle(std::vector<int> c1, std::vector<int> c2) {
  if(c1.size()>c2.size()) //must have a bigger size to include the cycle
    return false;
  int count=0;
  for(size_t it=0; it<c1.size(); it++) //count the number of containing vertices
  if(containElement(c2, c1.at(it)))
    count++;
  if(count==c1.size())
    return true;
  return false;
}

std::vector<std::vector<RationalPoint> > intersectionFace(const std::vector<RationalPoint>& face, const Complex& grid, bool saveFile)
{
  //Get all vertices and edges of a face
  std::vector<RationalPoint> vertex = intersectionVertex(face, grid);
  std::vector<std::pair<RationalPoint,RationalPoint> > edge = intersectionEdge(face, grid);
  if(saveFile) {
    Board2D aBoard;
    grid.drawGrid(aBoard);
    //Draw face
    std::vector<std::vector<RationalPoint> > cells;
    cells.push_back(face);
    Complex complex(cells);
    complex.drawComplex(aBoard, Color::Red);
    //Draw intersection vertices
    aBoard.setPenColor(Color::Blue);
    for(size_t it=0; it<vertex.size(); it++)
    aBoard.fillCircle(getRealValue(vertex.at(it).first), getRealValue(vertex.at(it).second), 0.05);
    //Draw intersection edges
    aBoard.setPenColor(Color::Green);
    aBoard.setLineWidth(1);
    for(size_t it=0; it<edge.size(); it++) {
      RationalPoint p1 = edge.at(it).first;
      RationalPoint p2 = edge.at(it).second;
      aBoard.drawLine(getRealValue(p1.first), getRealValue(p1.second),getRealValue(p2.first), getRealValue(p2.second));
    }
    aBoard.saveSVG("../Illustration/intersetionGrid.svg");
  }
  
  //Build the graph from the edges
  int aObjects[vertex.size()];
  for(int it=0; it<vertex.size(); it++)
  aObjects[it] = it;
  size_t aEdges[2*edge.size()];
  for(int it=0; it<edge.size(); it++) {
    RationalPoint e1 = edge.at(it).first;
    RationalPoint e2 = edge.at(it).second;
    int id1 = findVertex(vertex, e1);
    int id2 = findVertex(vertex, e2);
    assert(id1!=-1 && id2!=-1);
    aEdges[2*it]=size_t(id1);
    aEdges[2*it+1]=size_t(id2);
  }
  
  graph::Graph<int> myGraph(aObjects, vertex.size(), aEdges, edge.size());
  //myGraph.computeFundamentalCycles();
  myGraph.computeAllCycles();
  //myGraph.dump(std::cout);
  
  typedef std::list<const int*> NodePath;
  typedef graph::HalfAdjacencyMatrix CycleMatrix;
  typedef std::vector<CycleMatrix> CycleArray;
  CycleArray allCycles = myGraph.getAllCycles();
  std::vector<std::vector<int> > vecCycles;
  for (const CycleMatrix& cycleMatrix : allCycles) {
    std::vector<int> cycle;
    NodePath path = myGraph.cycleMatrix2nodePath(cycleMatrix);
    for (const int* obj : path)
      cycle.push_back(*obj);
    vecCycles.push_back(cycle);
  }
  //check for inclusion of cycles
  std::vector<bool> includingCycle (vecCycles.size(),false);
  for(int it=0; it<vecCycles.size()-1; it++) {
    if(isValidCycle(vertex, vecCycles.at(it)) && !includingCycle.at(it)) {
      for(int it_bis=it+1; it_bis<vecCycles.size(); it_bis++) {
        if(isIncludedCycle(vecCycles.at(it), vecCycles.at(it_bis)))
          includingCycle.at(it_bis)=true;
      }
    }
  }
  //associate cycles to faces
  std::vector<std::vector<RationalPoint> > vecFace;
  //std::cout<<"Print cycles:"<<std::endl;
  for(int it=0; it<vecCycles.size(); it++) {
    std::vector<RationalPoint> face;
    if(isValidCycle(vertex, vecCycles.at(it)) && !includingCycle.at(it)) {//TODO: symplify isValidCycle with includingCycle = true
      for(int it_bis=0; it_bis<vecCycles.at(it).size(); it_bis++) {
        //std::cout<<vecCycles.at(it).at(it_bis)<<" ";
        int id=vecCycles.at(it).at(it_bis);
        face.push_back(vertex.at(id));
      }
      //std::cout<<std::endl;
    }
    if(face.size()!=0)
      vecFace.push_back(face);
  }
  return vecFace;
}

std::vector<std::vector<RationalPoint> > intersection(const Complex& c, std::pair<RationalPoint,RationalPoint> line)
{
  std::vector<std::vector<RationalPoint> > vecFace;
  for(int it=0; it<c.face.size(); it++) {
    std::vector<RationalPoint> face;
    for(int it_bis=0; it_bis<c.face.at(it).size(); it_bis++)
    face.push_back(c.list_vertex.at(c.face.at(it).at(it_bis)));
    vecFace.push_back(face);
  }
  //Compute intersections of each face wrt to given lines
  std::vector<std::vector<RationalPoint> > vecIntersection;
  for(int it=0; it<vecFace.size(); it++) {
    std::vector<RationalPoint> p = intersection(vecFace.at(it), line);
    vecIntersection.push_back(p);
  }
  return vecIntersection;
}

std::vector<std::vector<RationalPoint> > intersectionVertex(const Complex& c, const Complex& grid)
{
  std::vector<std::vector<RationalPoint> > vecFace;
  for(int it=0; it<c.face.size(); it++) {
    std::vector<RationalPoint> face;
    for(int it_bis=0; it_bis<c.face.at(it).size(); it_bis++)
    face.push_back(c.list_vertex.at(c.face.at(it).at(it_bis)));
    vecFace.push_back(face);
  }
  //Compute the intersections of each face wrt vertical and horizontal lines
  std::vector<std::vector<RationalPoint> > vecVertex;
  for(int it=0; it<vecFace.size(); it++) { //For each face
    std::vector<RationalPoint> vertex = intersectionVertex(vecFace.at(it), grid);
    if(vertex.size()!=0)
      vecVertex.push_back(vertex);
  }
  return vecVertex;
}

std::vector<std::vector<std::pair<RationalPoint,RationalPoint> > > intersectionEdge(const Complex& c, const Complex& grid)
{
  std::vector<std::vector<RationalPoint> > vecFace;
  for(int it=0; it<c.face.size(); it++) {
    std::vector<RationalPoint> face;
    for(int it_bis=0; it_bis<c.face.at(it).size(); it_bis++)
    face.push_back(c.list_vertex.at(c.face.at(it).at(it_bis)));
    vecFace.push_back(face);
  }
  //Compute edges for each face
  std::vector<std::vector<std::pair<RationalPoint,RationalPoint> > > vecEdge;
  for(int it=0; it<vecFace.size(); it++) { //For each face
    std::vector<std::pair<RationalPoint,RationalPoint> > edge = intersectionEdge(vecFace.at(it), grid);
    vecEdge.push_back(edge);
  }
  return vecEdge;
}

//return true if p belongs to an edge face
bool onEdgeFace(const std::vector<RationalPoint>& face, RationalPoint p)
{
  //Compute all edge face
  std::vector<std::pair<RationalPoint, RationalPoint> > edgeFace;
  for(int it=0; it<face.size(); it++)
  edgeFace.push_back(std::make_pair(face.at(it), face.at(it+1)));
  for(int it=0; it<edgeFace.size(); it++) { //for each edge face
    RationalPoint p1 = edgeFace.at(it).first;
    RationalPoint p2 = edgeFace.at(it).second;
    if (onSegment(p1, p2, p))
      return true;
  }
  return false;
}

//return true if p1 and p2 belongs to the same edge face
bool onEdgeFace(const std::vector<RationalPoint>& face, RationalPoint p1, RationalPoint p2)
{
  //Compute all edge face
  std::vector<std::pair<RationalPoint, RationalPoint> > edgeFace;
  for(int it=0; it<face.size(); it++)
  edgeFace.push_back(std::make_pair(face.at(it), face.at(it+1)));
  for(int it=0; it<edgeFace.size(); it++) { //for each edge face
    RationalPoint e1 = edgeFace.at(it).first;
    RationalPoint e2 = edgeFace.at(it).second;
    if (onSegment(e1, e2, p1) && onSegment(e1, e2, p2))
      return true;
  }
  return false;
}

std::vector<std::vector<std::vector<RationalPoint> > > intersectionFace(const Complex& c, const Complex& grid)
{
  //Get faces of the complex c
  std::vector<std::vector<RationalPoint> > vecFace;
  for(int it=0; it<c.face.size(); it++) {
    std::vector<RationalPoint> face;
    for(int it_bis=0; it_bis<c.face.at(it).size(); it_bis++)
      face.push_back(c.list_vertex.at(c.face.at(it).at(it_bis)));
    vecFace.push_back(face);
  }
  //Compute edges for each face
  std::vector<std::vector<std::vector<RationalPoint> > > vecCell;
  for(int it=0; it<vecFace.size(); it++) { //For each face
    std::vector<std::vector<RationalPoint> > cells = intersectionFace(vecFace.at(it), grid);
    vecCell.push_back(cells);
  }
  return vecCell;
}

void embedComplexInGrid(const Complex& c, const Complex& grid, std::vector<bool>& belongFace)
{
  for(int it=0; it<grid.face.size(); it++)
    belongFace.push_back(false);
  
  //Get faces of the complex c
  std::vector<std::vector<RationalPoint> > vecFace;
  for(int it=0; it<c.face.size(); it++) {
    std::vector<RationalPoint> face;
    for(int it_bis=0; it_bis<c.face.at(it).size(); it_bis++)
      face.push_back(c.list_vertex.at(c.face.at(it).at(it_bis)));
    vecFace.push_back(face);
  }
  //Verify if any face of the grid belong to the complex C
  for(int it=0; it<grid.face.size(); it++) {
    RationalPoint p = grid.getFaceCenter(it);
    int res = c.isInsideComplex(p);
    if(res != -1)
      belongFace.at(it) = true;
  }
}

std::vector<std::pair<int,int> > Complex::getBoundaryEdges(int id_face) const
{
  std::vector<std::pair<int,int> > e;
  std::vector<int> f;
  std::vector<int> face = this->face.at(id_face);
  std::pair<int,int> edge;
  for(int it=0; it<face.size()-1; it++) {
    edge = std::make_pair(face.at(it), face.at(it+1));
    f = getBelongFaces(this->face, id_face, edge);
    if(f.size()==0)
      e.push_back(edge);
  }
  return  e;
}

std::vector<std::pair<int,int> > Complex::getNonBoundaryEdges(int id_face) const
{
  std::vector<std::pair<int,int> > e;
  std::vector<int> f;
  std::vector<int> face = this->face.at(id_face);
  std::pair<int,int> edge;
  for(int it=0; it<face.size()-1; it++)
  {
    edge = std::make_pair(face.at(it), face.at(it+1));
    f = getBelongFaces(this->face, id_face, edge);
    if(f.size()!=0)
      e.push_back(edge);
  }
  return  e;
}

//Get vertices of a face being (non) boundary wrt object's cells
bool isBelongFace(std::vector<int> face, int vertex)
{
  for(int it=0; it<face.size()-1; it++)//-1 : closed face
  if(vertex==face.at(it))
    return true;
  return false;
}

std::vector<int> getBelongFaces(std::vector<std::vector<int> > faces, int vertex)
{
  std::vector<int> f;
  for(int it=0; it<faces.size(); it++)
  if(isBelongFace(faces.at(it), vertex))
    f.push_back(it);
  return f;
}

std::vector<int> getBelongFaces(std::vector<std::vector<int> > faces, int id_face, int vertex)
{
  std::vector<int> f;
  for(int it=0; it<faces.size(); it++)
  if(id_face!=it && isBelongFace(faces.at(it), vertex))
    f.push_back(it);
  return f;
}

std::vector<int> Complex::getBoundaryVertices(int id_face) const
{
  std::vector<int> v;
  std::vector<int> f;
  std::vector<int> face = this->face.at(id_face);
  int vertex;
  for(int it=0; it<face.size()-1; it++) {//-1 : closed face
    vertex = face.at(it);
    f = getBelongFaces(this->face, id_face, vertex);
    if(f.size()==0) //vertex does not belong to any other face
      v.push_back(vertex);
  }
  return  v;
}

std::vector<int> Complex::getNonBoundaryVertices(int id_face) const
{
  std::vector<int> v;
  std::vector<int> f;
  std::vector<int> face = this->face.at(id_face);
  int vertex;
  for(int it=0; it<face.size(); it++)
  {
    vertex = face.at(it);
    f = getBelongFaces(this->face, id_face, vertex);
    if(f.size()!=0) //vertex belongs to a face different than id_face
      v.push_back(vertex);
  }
  return  v;
}

bool Complex::isSimpleFace(int id_face) const
{
  int v = getBoundaryVertices(id_face).size();
  int e = getBoundaryEdges(id_face).size();
  //std::cout<<"v="<<v<<" and e="<<e<<std::endl;
  int x = e-v;
  return x == 1;
}

bool Complex::isSimpleFace(const std::vector<RationalPoint>& face) const
{
  Complex tmp = *this;
  tmp.addFace(face, true, true);
  if(tmp.face.size()==this->face.size()) //new face is added
    return isSimpleFace(tmp.face.size()-1);
  else
    assert(false); //TODO return the simplicity of the face
}

bool Complex::isBorderVertex(int id_vertex) const
{
  std::vector<int> f = getBelongFaces(this->face, id_vertex);
  return  f.size()==1; //vertex belong to 1 face
}

bool Complex::isBorderEdge(int id_edge) const
{
  std::pair<int,int> edge = this->edge.at(id_edge);
  std::vector<int> f = getBelongFaces(this->face, edge);
  return  f.size()==1; //edge belong to 1 face
}

bool Complex::isBorderFace(int id_face) const
{
  if(this->getBoundaryEdges(id_face).size()!=0 || this->getBoundaryVertices(id_face).size()!=0)
    return true;
  return false;
}

int Complex::isInsideComplex(RationalPoint p) const
{
  for(size_t it=0; it<this->face.size(); it++) {
    if(isInsidePolygon(this->getFaceVertices(it), p))
      return it;
  }
  return -1;
}

Board2D drawVertex(const RationalPoint vertex, Board2D& board, double size, Color c)
{
  double x = getRealValue(vertex.first);
  double y = getRealValue(vertex.second);
  
  board.setPenColor(c);
  board.fillCircle(x, y, size);
  return board;
}

Board2D drawEdge(const std::pair<RationalPoint,RationalPoint>& edge, Board2D& board, double size, Color c)
{
  RationalPoint rp1 = edge.first;
  RationalPoint rp2 = edge.second;
  
  Rational rp1x = rp1.first;
  Rational rp1y = rp1.second;
  double x1 = getRealValue(rp1x);
  double y1 = getRealValue(rp1y);
  
  Rational rp2x = rp2.first;
  Rational rp2y = rp2.second;
  double x2 = getRealValue(rp2x);
  double y2 = getRealValue(rp2y);
  
  board.setPenColor(c);
  board.setLineWidth( size );
  board.setLineStyle(LibBoard::Arc::SolidStyle);
  double dx = x2-x1;
  double dy = y2-y1;
  double xx1 = x1 + 2*size*dx/10;
  double yy1 = y1 + 2*size*dy/10;
  double xx2 = x1 + (1 - 2*size/10) * dx;
  double yy2 = y1 + (1 - 2*size/10) * dy;
  board.drawLine(xx1, yy1, xx2, yy2);
  
  return board;
}

Board2D drawFace(const std::vector<RationalPoint>& face, Board2D& board, double scale, Color c)
{
  std::vector<LibBoard::Point> points;
  //std::cout<<"Face"<<std::endl;
  for(int it=0; it<face.size(); it++) {
    RationalPoint rp = face.at(it);
    Rational rx = rp.first;
    Rational ry = rp.second;
    double x = getRealValue(rx);
    double y = getRealValue(ry);
    //std::cout<<rp;
    points.push_back(LibBoard::Point(x,y));
  }
  std::cout<<std::endl;
  //Set transparence for face
  c.alpha(100);
  board.setPenColor(c);
  
  LibBoard::Polyline polyline(points,true,c,c,scale);
  polyline.scale(scale);
  std::vector<LibBoard::Point> points_scale;
  for (int it=0; it<points.size(); it++)
  points_scale.push_back(polyline[it]);
  board.fillPolyline(points_scale);
  
  return board;
}

Board2D drawPoint(const Z2i::Point vertex, Board2D& board, double size, Color c)
{
  board.setPenColor(c);
  board.fillCircle(vertex[0], vertex[1], size);
  return board;
}

Board2D Complex::drawVertex(Board2D& board, int index, Color c) const
{
  if(index>=0 && index<list_vertex.size()) {
    RationalPoint rp = list_vertex.at(index);
    Rational rx = rp.first;
    Rational ry = rp.second;
    double x = getRealValue(rx);
    double y = getRealValue(ry);
    
    board.setPenColor(c);
    board.fillCircle(x, y, 4*SIZE/this->getBorderSize());
  }
  return board;
}

Board2D Complex::drawEdge(Board2D& board, int index, Color c) const
{
  if(index>=0 && index<edge.size()) {
    std::pair<int,int> e = edge.at(index);
    int e1 = e.first;
    int e2 = e.second;
    
    RationalPoint rp1 = list_vertex.at(e1);
    RationalPoint rp2 = list_vertex.at(e2);
    
    Rational rp1x = rp1.first;
    Rational rp1y = rp1.second;
    double x1 = getRealValue(rp1x);
    double y1 = getRealValue(rp1y);
    
    Rational rp2x = rp2.first;
    Rational rp2y = rp2.second;
    double x2 = getRealValue(rp2x);
    double y2 = getRealValue(rp2y);
    
    board.setPenColor(c);
    board.setLineWidth( 100*SIZE/this->getBorderSize() );
    board.setLineStyle(LibBoard::Arc::SolidStyle);
    
    double dx = x2-x1;
    double dy = y2-y1;
    double xx1 = x1 + 10*SIZE*dx/this->getBorderSize();
    double yy1 = y1 + 10*SIZE*dy/this->getBorderSize();
    double xx2 = x1 + (1 - 10*SIZE/this->getBorderSize()) * dx;
    double yy2 = y1 + (1 - 10*SIZE/this->getBorderSize()) * dy;
    board.drawLine(xx1, yy1, xx2, yy2);
  }
  return board;
}

Board2D Complex::drawFace(Board2D& board, int index, Color c) const
{
  if(index>=0 && index<face.size()) {
    std::vector<int> aFace = face.at(index);
    std::vector<LibBoard::Point> points;
    for(int it=0; it<aFace.size(); it++) {
      RationalPoint rp = list_vertex.at(aFace.at(it));
      Rational rx = rp.first;
      Rational ry = rp.second;
      double x = getRealValue(rx);
      double y = getRealValue(ry);
      points.push_back(LibBoard::Point(x,y));
    }
    
    //Set transparence for face
    c.alpha(100);
    board.setPenColor(c);
    
    LibBoard::Polyline polyline(points,true,c,c,SIZE);
    polyline.scale(1- 20*SIZE/this->getBorderSize()); //10*SIZE/this->getBorderSize()
    std::vector<LibBoard::Point> points_scale;
    for (int it=0; it<points.size(); it++)
    points_scale.push_back(polyline[it]);
    board.fillPolyline(points_scale);
  }
  return board;
}

Board2D Complex::drawVertex(Board2D& board, Color c) const
{
  for(int it=0; it<vertex.size(); it++)
    this->drawVertex(board, vertex.at(it), c);
  return board;
}

Board2D Complex::drawEdge(Board2D& board, Color c) const
{
  for(int it=0; it<edge.size(); it++)
    this->drawEdge(board, it, c);
  return board;
}

Board2D Complex::drawFace(Board2D& board, Color c) const
{
  for(int it=0; it<face.size(); it++)
    this->drawFace(board, it, c);
  return board;
}

Board2D Complex::drawGrid(Board2D& board, Color c) const
{
  board.setPenColor(c);
  board.setLineWidth( 20*SIZE/this->getBorderSize() );
  board.setLineStyle(LibBoard::Arc::DashStyle);
  for(int it=0; it<vertical_lines.size();it++) {
    RationalPoint p1 = vertical_lines.at(it).first;
    double x1 = getRealValue(p1.first);
    double y1 = getRealValue(p1.second);
    RationalPoint p2 = vertical_lines.at(it).second;
    double x2 = getRealValue(p2.first);
    double y2 = getRealValue(p2.second);
    board.drawLine(x1, y1, x2, y2);
  }
  for(int it=0; it<horizontal_lines.size();it++) {
    RationalPoint p1 = horizontal_lines.at(it).first;
    double x1 = getRealValue(p1.first);
    double y1 = getRealValue(p1.second);
    RationalPoint p2 = horizontal_lines.at(it).second;
    double x2 = getRealValue(p2.first);
    double y2 = getRealValue(p2.second);
    board.drawLine(x1, y1, x2, y2);
  }
  return board;
}

Board2D Complex::drawBorder(Board2D& board, Color c) const
{
  RationalPoint minBB = RationalPoint(INT_MAX, INT_MAX);
  RationalPoint maxBB = RationalPoint(INT_MIN, INT_MIN);
  if(this->vertical_lines.size()!=0) {
    std::pair<RationalPoint,RationalPoint> bb = getBoundingBox(this->vertical_lines);
    minBB = min(bb.first, minBB);
    maxBB = max(bb.second, maxBB);
  }
  if(this->horizontal_lines.size()!=0) {
    std::pair<RationalPoint,RationalPoint> bb = getBoundingBox(this->horizontal_lines);
    minBB = min(bb.first, minBB);
    maxBB = max(bb.second, maxBB);
  }
  if(this->list_vertex.size()!=0) {
    std::pair<RationalPoint,RationalPoint> bb = getBoundingBox(this->list_vertex);
    minBB = min(bb.first, minBB);
    maxBB = max(bb.second, maxBB);
  }
  
  board.setPenColor(c);
  board.setLineWidth( 40*SIZE/this->getBorderSize() );
  board.setLineStyle(LibBoard::Arc::SolidStyle);
  
  double x = getRealValue(minBB.first);
  double y = getRealValue(maxBB.second);
  double w = getRealValue(maxBB.first) - getRealValue(minBB.first);
  double h = getRealValue(maxBB.second) - getRealValue(minBB.second);
  board.drawRectangle(x, y, w, h);
  return board;
}

Board2D Complex::drawCenter(Board2D& board, Color c) const
{
  board.setPenColor(c);
  std::vector<RationalPoint> vecP = this->getFaceCenter();
  for(int it=0; it<vecP.size(); it++) {
    double x = getRealValue(vecP.at(it).first);
    double y = getRealValue(vecP.at(it).second);
    board.fillCircle(x, y, 5*SIZE/this->getBorderSize());
  }
  return board;
}

Board2D Complex::drawComplex(Board2D& board, Color c, bool drawBorder, bool drawDomain) const
{
  Z2i::Point p1 = domain.first;
  Z2i::Point p2 = domain.second;
  Z2i::Domain d( p1, p2 );
  if(drawDomain)
    //board << d; //Draw domain by DGtal
    this->drawGrid(board,c);
  if(drawBorder)
    this->drawBorder(board,c);
  this->drawFace(board, c);
  this->drawEdge(board, c);
  this->drawVertex(board, c);
  return board;
}

Board2D Complex::drawPixel(Board2D& board, Color c, bool drawBorder, bool drawDomain) const
{
  Z2i::Point p1 = domain.first;
  Z2i::Point p2 = domain.second;
  Z2i::Domain d( p1, p2 );
  if(drawDomain)
    //board << d; //Draw domain by DGtal
    this->drawGrid(board,c);
  if(drawBorder)
    this->drawBorder(board,c);
  board.setPenColor(c);
  for(int it=0; it<this->face.size(); it++) {
    std::vector<RationalPoint> aFace = this->getFaceVertices(it);
    assert(aFace.size()-1==4);
    
    std::vector<LibBoard::Point> points;
    for(int it_bis=0; it_bis<aFace.size(); it_bis++) {
      RationalPoint rp = aFace.at(it_bis);
      Rational rx = rp.first;
      Rational ry = rp.second;
      double x = getRealValue(rx);
      double y = getRealValue(ry);
      points.push_back(LibBoard::Point(x,y));
    }
    board.fillPolyline(points);
  }
  return board;
}

std::ostream& operator<<(std::ostream& os, const Complex& c)
{
  os << "Print complex "<<std::endl<<"vertices("<<c.list_vertex.size()<<")"<<std::endl;
  for (int it=0; it<c.list_vertex.size(); it++)
  os << c.list_vertex.at(it).first<<","<<c.list_vertex.at(it).second<<std::endl;
  os <<"Vertex ("<<c.vertex.size()<<")"<<std::endl;
  for (int it=0; it<c.vertex.size(); it++)
  os << c.vertex.at(it)<<" ";
  os << std::endl;
  os <<"Edge("<<c.edge.size()<<")"<<std::endl;
  for (int it=0; it<c.edge.size(); it++)
  os <<"("<< c.edge.at(it).first<<","<<c.edge.at(it).second<<") ";
  os << std::endl;
  os <<"Face("<<c.face.size()<<")"<<std::endl;
  for (int it=0; it<c.face.size(); it++) {
    for (int it_bis=0; it_bis<c.face.at(it).size(); it_bis++)
    os << c.face.at(it).at(it_bis)<<" ";
    os << std::endl;
  }
  return os;
}
