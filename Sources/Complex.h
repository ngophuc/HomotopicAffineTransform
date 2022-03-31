#ifndef COMPLEX_H
#define COMPLEX_H

#include "Functions.h"
#include "DGtal/io/boards/Board2D.h"

class Complex {
public:
  std::vector<std::pair<RationalPoint,RationalPoint> > horizontal_lines;//Vector Horizontal lines
  std::vector<std::pair<RationalPoint,RationalPoint> > vertical_lines; //Vector Vertical lines
  std::vector<RationalPoint> list_vertex; //list of all points
  std::vector<int> vertex; //list of vertices of complex
  std::vector<std::pair<int,int> > edge; //list of edges of complex
  std::vector<std::vector<int> > face; //list of (closed) faces of complex
  std::pair<Z2i::Point,Z2i::Point> domain;
  
  std::vector<int> cell_vertex; //cell of vertices (= vertex)
  std::vector<std::pair<int,int> > cell_edge; //cell of edges (= edge id)
  std::vector<std::vector<int> > cell_face_vertex; //cell of faces (=set of vertices by id)
  std::vector<std::vector<int> > cell_face_edge; //cell of faces (=set of edges by id)
  
  std::vector<std::vector<int> > star_vertex_edge; //star of vertices (= set of edges containing the vertex)
  std::vector<std::vector<int> > star_vertex_face; //star of vertices (= set of faces containing the vertex)
  std::vector<std::pair<int, int> > star_edge; //star of edges (= 2 faces sharing the edge)
  std::vector<int> star_face; //star of faces (= face)
  
  Complex(){}
  Complex(const Complex &);
  Complex(const std::vector<std::vector<RationalPoint> >& faces);
  Complex(const std::vector<Z2i::Point>& points);
  Complex(const Z2i::Point point);
  Complex(const std::vector<RationalPoint>& vertex, const std::vector<std::pair<RationalPoint,RationalPoint> >& edge, const std::vector<std::vector<RationalPoint> >& face);
  
  //Getters
  std::pair<Z2i::Point,Z2i::Point> getDomain() const {return this->domain;};
  std::vector<Line> getVerticalLines() const;
  std::vector<Line> getHorizontalLines() const;
  std::vector<Line> getGridLines() const;
  std::pair<int,int> getDomainSize() const;
  int getBorderSize() const;
  
  //Other getters
  RationalPoint getVertex(int index) const;
  std::vector<RationalPoint> getVertex(const std::vector<int>& index) const;
  std::pair<RationalPoint, RationalPoint> getEdgeVertices(int index) const;
  std::vector<std::pair<RationalPoint, RationalPoint> > getEdgeVertices(const std::vector<int>& index) const;
  std::vector<RationalPoint> getFaceVertices(int index) const;
  double getFaceArea(int index) const;
  std::vector<std::vector<RationalPoint> > getFaceVertices(const std::vector<int>& index) const;
  std::vector<std::vector<RationalPoint> > getFaceVertices() const;
  
  RationalPoint getFaceCenter(const std::vector<RationalPoint>& face) const;
  RationalPoint getFaceCenter(int index) const;
  std::vector<RationalPoint> getFaceCenter() const;
  
  std::vector<int> getNeighbourFace(int index) const;
  std::vector<std::vector<int> > getAllNeighbourFace() const;
  
  //get id of face that p belongs to
  int getIdFace(RationalPoint p);
  
  //Cell and star operators
  int computeCellVertex (int idV);
  std::pair<int,int> computeCellEdge (int idE);
  std::vector<int> computeCellFaceVertex (int idF);
  std::vector<int> computeCellFaceEdge (int idF);
  
  std::vector<int> computeStarVertexEdge (int idV);
  std::vector<int> computeStarVertexFace (int idV);
  std::pair<int, int> computeStarEdge (int idE);
  int computeStarFace (int idF);
  
  void updateAllCell();
  void updateAllStar();
  
  //init with limits of the grid space
  void init(Z2i::Point p1, Z2i::Point p2);
  void initGrid(Z2i::Point p1, Z2i::Point p2);
  void initGridFaces(Z2i::Point p1, Z2i::Point p2);
  void initGridLines(std::vector<std::pair<RationalPoint,RationalPoint> > vlines, std::vector<std::pair<RationalPoint,RationalPoint> > hlines);
  void addGridLines(std::vector<std::pair<RationalPoint,RationalPoint> > vlines, std::vector<std::pair<RationalPoint,RationalPoint> > hlines);
  void addVertex(RationalPoint v, bool updateCellStar = false);
  void addEdge(RationalPoint e1, RationalPoint e2, bool addVertex = false, bool updateCellStar = false);
  void addEdge(std::pair<RationalPoint, RationalPoint> e, bool addVertex = false, bool updateCellStar = false);
  void addFace(const std::vector<RationalPoint>& f, bool addEdge = false, bool addVertex = false, bool updateCellStar = false);
  //Add a face by its center
  void addCubicalFace(Z2i::Point p);
  
  //verify if a face is simple wrt object's cells: simple if Euleur=NonBoundary(V-E)=1
  bool isSimpleFace(int id_face) const;
  //verify if a new face (added) is simple wrt object's cells: simple if Euleur=NonBoundary(V-E)=1
  bool isSimpleFace(const std::vector<RationalPoint>& face) const;
  //verify if a face is border wrt object's cells: it contains boundary vertices / edges
  bool isBorderVertex(int id_vertex) const;
  bool isBorderEdge(int id_edge) const;
  bool isBorderFace(int id_face) const;
  //verify if a face (its center) belongs to object's cells (complex), return index of face, -1 otherwise
  int isInsideComplex(RationalPoint p) const;
  
  //Get edges/vertices of a face being (non) boundary wrt object's cells
  std::vector<std::pair<int,int> > getBoundaryEdges(int id_face) const;
  std::vector<int> getBoundaryVertices(int id_face) const;
  std::vector<std::pair<int,int> > getNonBoundaryEdges(int id_face) const;
  std::vector<int> getNonBoundaryVertices(int id_face) const;
  
  //Drawing functions
  Board2D drawVertex(Board2D& board, int index, Color c = Color::Red) const;
  Board2D drawEdge(Board2D& board, int index, Color c = Color::Green) const;
  Board2D drawFace(Board2D& board, int index, Color c = Color::Blue) const;
  
  Board2D drawVertex(Board2D& board, Color c = Color::Red) const;
  Board2D drawEdge(Board2D& board, Color c = Color::Green) const;
  Board2D drawFace(Board2D& board, Color c = Color::Blue) const;
  Board2D drawGrid(Board2D& board, Color c = Color::Black) const;
  Board2D drawBorder(Board2D& board, Color c = Color::Black) const;
  Board2D drawCenter(Board2D& board, Color c = Color::Black) const;
  Board2D drawComplex(Board2D& board, Color cv, Color ce, Color cf, bool drawDomain=false) const;
  Board2D drawComplex(Board2D& board, Color c, bool drawBorder=false, bool drawDomain=false) const;
  Board2D drawPixel(Board2D& board, Color c = Color::Black, bool drawBorder=false, bool drawDomain=false) const;
  
  friend std::ostream& operator<<(std::ostream& os, const Complex& c);
};

std::vector<RationalPoint> intersectionVertex(const std::vector<RationalPoint>& face, const Complex& grid);
std::vector<std::pair<RationalPoint,RationalPoint> > intersectionEdge(const std::vector<RationalPoint>& face, const Complex& grid);
std::vector<std::vector<RationalPoint> > intersectionFace(const std::vector<RationalPoint>& face, const Complex& grid, bool saveFile = false);

std::vector<std::vector<RationalPoint> > intersection(const Complex& c, std::pair<RationalPoint,RationalPoint> line);
std::vector<std::vector<RationalPoint> > intersectionVertex(const Complex& c, Complex grid);
std::vector<std::vector<std::pair<RationalPoint,RationalPoint> > > intersectionEdge(const Complex& c, const Complex& grid);
std::vector<std::vector<std::vector<RationalPoint> > > intersectionFace(const Complex& c, const Complex& grid);

void embedComplexInGrid(const Complex& c, const Complex& grid, std::vector<bool>& belongFace);
//Function phi (h2 -> f2) and gamma (h2->g2)

Board2D drawVertex(const RationalPoint vertex, Board2D& board, double size, Color c = Color::Red);
Board2D drawEdge(const std::pair<RationalPoint,RationalPoint>& edge, Board2D& board, double size, Color c = Color::Green);
Board2D drawFace(const std::vector<RationalPoint>& face, Board2D& board, double scale, Color c = Color::Blue);
Board2D drawPoint(const Z2i::Point vertex, Board2D& board, double size, Color c = Color::Black);

#endif // COMPLEX_H
