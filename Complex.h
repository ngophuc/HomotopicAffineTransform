#ifndef COMPLEX_H
#define COMPLEX_H

#include "Functions.h"
#include "DGtal/io/boards/Board2D.h"

class Complex {
public:
  //std::vector<Line> horizontal_lines;//Vector Horizontal lines
  //std::vector<Line> vertical_lines; //Vector Vertical lines
  std::vector<std::pair<RationalPoint,RationalPoint> > horizontal_lines;//Vector Horizontal lines
  std::vector<std::pair<RationalPoint,RationalPoint> > vertical_lines; //Vector Vertical lines
  std::vector<RationalPoint> list_vertex; //list of all points
  std::vector<int> vertex; //list of vertices of complex
  std::vector<std::pair<int,int> > edge; //list of edges of complex
  std::vector<std::vector<int> > face; //list of (closed) faces of complex
  std::vector<int> idConnectedComponent; //list of id of connected component of face
  std::pair<Z2i::Point,Z2i::Point> domain;
  
  Complex(){}
  Complex(const Complex &);
  Complex(const std::vector<std::vector<RationalPoint> >& faces);
  Complex(const std::vector<Z2i::Point>& points, int adjacency = 4);
  Complex(const Z2i::Point point);
  Complex(const std::vector<RationalPoint>& vertex, const std::vector<std::pair<RationalPoint,RationalPoint> >& edge, const std::vector<std::vector<RationalPoint> >& face);
  //Constructor of complex given a domain (Limits of the grid space)
  //Complex(Z2i::Point p1, Z2i::Point p2); 
  //Constructor of complex given a face
  //Complex(Z2i::Point p);
  
  std::pair<Z2i::Point,Z2i::Point> getDomain() const {return this->domain;};
  std::vector<Line> getVerticalLines() const;
  std::vector<Line> getHorizontalLines() const;
  std::vector<Line> getGridLines() const;
  std::pair<int,int> getDomainSize() const;
  int getBorderSize() const;
  
  //init with limits of the grid space
  void init(Z2i::Point p1, Z2i::Point p2);
  void initGrid(Z2i::Point p1, Z2i::Point p2);
  void initGridFaces(Z2i::Point p1, Z2i::Point p2);
  void initGridLines(std::vector<std::pair<RationalPoint,RationalPoint> > vlines, std::vector<std::pair<RationalPoint,RationalPoint> > hlines);
  void addGridLines(std::vector<std::pair<RationalPoint,RationalPoint> > vlines, std::vector<std::pair<RationalPoint,RationalPoint> > hlines);
  void addVertex(RationalPoint v);
  void addEdge(RationalPoint e1, RationalPoint e2, bool addVertex = false);
  void addEdge(std::pair<RationalPoint, RationalPoint> e, bool addVertex = false);
  void addFace(const std::vector<RationalPoint>& f, bool addEdge = false, bool addVertex = false);
  //Add a face by its center
  void addCubicalFace(Z2i::Point p);
  //remove complex elements
  void removeVertex(int idV);
  void removeEdge(int idE);
  void removeFace(int idF);
  
  RationalPoint getVertex(int index) const;
  std::vector<RationalPoint> getVertex(const std::vector<int>& index) const;
  std::pair<RationalPoint, RationalPoint> getEdgeVertices(int index) const;
  std::vector<std::pair<RationalPoint, RationalPoint> > getEdgeVertices(const std::vector<int>& index) const;
  std::vector<RationalPoint> getFaceVertices(int index) const;
  double getFaceArea(int index) const;
  std::vector<std::vector<RationalPoint> > getFaceVertices(const std::vector<int>& index) const;
  std::vector<std::vector<RationalPoint> > getFaceVertices() const;
  int getFaceConnectedComponent(int index) const;
  std::vector<int> getFaceConnectedComponent() const;
  
  RationalPoint getFaceCenter(const std::vector<RationalPoint>& face) const;
  RationalPoint getFaceCenter(int index) const;
  std::vector<RationalPoint> getFaceCenter() const;
  
  std::vector<int> getNeighbourFace(int index) const;
  std::vector<std::vector<int> > getAllNeighbourFace() const;
  
  //get id of face that p belongs to
  int getIdFace(RationalPoint p);
  //get id of connected component of face
  int getIndexConnectedComponent(const std::vector<RationalPoint>& face) const;
  
  //Function phi and gamma (h2 (id_face) -> g2 (grid))
  std::vector<RationalPoint> mappingGrid(int id_face, const Complex& grid); //return the face of grid corresponding to the id-face
  int mappingGridIndex(int id_face, const Complex& grid); //return index of the face (of grid)
  std::vector<std::vector<RationalPoint> > mappingGrid(const Complex& grid); // find all mapping
  std::vector<std::vector<RationalPoint> > mappingGridArea(const Complex& grid, double area =  0.1); // find all mapping with area(face) > area
  std::vector<int> mappingGridIndex(const Complex& grid); //return index of the face (of grid)
  
  //Function Phi and Gamma (inverse of phi and gamma respectively) (g2 (face) -> {h2}(*this in return))
  std::vector<std::vector<RationalPoint> > inverseMapping(const std::vector<RationalPoint>& face); //return the set of this->face belonging to face
  std::vector<int> inverseMappingIndex(const std::vector<RationalPoint>& face); //return the indices of this->faces
  //Function Phi and Gamma (g2 (id_face of grid) -> {h2}(*this in return))
  std::vector<std::vector<RationalPoint> > inverseMapping(const Complex& grid, int id_face); //return this->face
  std::vector<int> inverseMappingIndex(const Complex& grid, int id_face); //return the index of this->face
  
  std::vector<std::vector<std::vector<RationalPoint> > > inverseMapping(const std::vector<std::vector<RationalPoint> >& face); // find all inverse mapping wrt to face
  std::vector<std::vector<int> > inverseMappingIndex(const std::vector<std::vector<RationalPoint> >& face); // find all inverse mapping wrt to face
  
  //bool allFace = true => save all faces of GridF
  //bool allFace = false (default) => save only faces that contain at least one inverse mapping
  std::vector<std::vector<std::vector<RationalPoint> > > inverseMapping(const Complex& grid, bool allFace=false); // find all inverse mapping wrt to all faces of grid
  std::vector<std::vector<int> > inverseMappingIndex(const Complex& grid, bool allFace=false); // find all inverse mapping wrt to all faces of grid
  
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
  
  //Gaussian digitization
  std::vector<Z2i::Point> gaussianDigitization(const Complex& grid);
  //Majority vote digitization
  std::vector<Z2i::Point> majorityVoteDigitization(const Complex& grid);
  
  //Get face of grid containing cell, return face id of grid
  Z2i::Point digitization(int id_face, const Complex& grid);
  
  //find the face of grid that contain face centers of the current complex and return those faces (of grid)
  std::vector<std::vector<RationalPoint> > digitizationFace(const Complex& grid);
  //find the points of grid that is the digitization face centers of the current complex and return those points (of grid)
  std::vector<Z2i::Point> digitization(const Complex& grid);
  
  
  Complex getSortAreaFaces();
  Complex getSortAreaFaces(Complex& complexRef);
  
  Board2D drawVertex(Board2D& board, int index, Color c = Color::Red) const;
  Board2D drawEdge(Board2D& board, int index, Color c = Color::Green) const;
  Board2D drawFace(Board2D& board, int index, Color c = Color::Blue) const;
  
  Board2D drawVertex(Board2D& board, Color c = Color::Red) const;
  Board2D drawEdge(Board2D& board, Color c = Color::Green) const;
  Board2D drawFace(Board2D& board, Color c = Color::Blue) const;
  Board2D drawGrid(Board2D& board, Color c = Color::Black) const;
  Board2D drawBorder(Board2D& board, Color c = Color::Black) const;
  Board2D drawCenter(Board2D& board, Color c = Color::Black) const;
  //Board2D drawComplex(Board2D& board, Color cv = Color::Red, Color ce = Color::Green, Color cf = Color::Blue, bool drawDomain=false);
  Board2D drawComplex(Board2D& board, Color cv, Color ce, Color cf, bool drawDomain=false) const;
  Board2D drawComplex(Board2D& board, Color c, bool drawBorder=false, bool drawDomain=false) const;
  Board2D drawPixel(Board2D& board, Color c = Color::Black, bool drawBorder=false, bool drawDomain=false) const;
  Board2D drawConectedComponent(Board2D& board, Color c = Color::Black, bool drawBorder=false, bool drawDomain=false) const;
  
  friend std::ostream& operator<<(std::ostream& os, const Complex& c);
};

std::vector<RationalPoint> intersectionVertex(const std::vector<RationalPoint>& face, const Complex& grid);
std::vector<std::pair<RationalPoint,RationalPoint> > intersectionEdge(const std::vector<RationalPoint>& face, const Complex& grid);
//std::vector<std::vector<int> > intersectionFace(const std::vector<RationalPoint>& face, Complex grid);
std::vector<std::vector<RationalPoint> > intersectionFace(const std::vector<RationalPoint>& face, const Complex& grid, bool saveFile = false);

std::vector<std::vector<RationalPoint> > intersection(const Complex& c, std::pair<RationalPoint,RationalPoint> line);
std::vector<std::vector<RationalPoint> > intersectionVertex(const Complex& c, Complex grid);
std::vector<std::vector<std::pair<RationalPoint,RationalPoint> > > intersectionEdge(const Complex& c, const Complex& grid);
//std::vector<std::vector<std::vector<int> > > intersectionFace(Complex c, Complex grid);
std::vector<std::vector<std::vector<RationalPoint> > > intersectionFace(const Complex& c, const Complex& grid);

//Function phi (h2 -> f2) and gamma (h2->g2)
std::vector<RationalPoint> mappingGrid(const std::vector<RationalPoint>& face, const Complex& grid); //find the face of grid corresponding to the given face
//std::vector<std::vector<RationalPoint> > mappingGrid(const std::vector<std::vector<RationalPoint> >& face, const Complex& grid); //find the face of grid corresponding to the given face
//Function Phi and Gamma (inverse of phi and gamma respectively)
std::vector<std::vector<RationalPoint> > inverseMapping(const std::vector<RationalPoint>& face, const Complex& grid); //find the set of faces of grid belonging to face

//find the face of grid that contain face center and return those faces (of grid)
std::vector<std::vector<RationalPoint> > digitizationFace(const std::vector<RationalPoint>& faceCenter, const Complex& grid);
std::vector<std::vector<RationalPoint> > digitizationFace(const std::vector<std::vector<RationalPoint> > & faces, const Complex& grid);
//find the points of grid that is the digitization face centers of the current complex and return those points (of grid)
std::vector<Z2i::Point> digitization(const std::vector<RationalPoint>& faceCenter, const Complex& grid);
std::vector<Z2i::Point> digitization(const std::vector<std::vector<RationalPoint> > & faces, const Complex& grid);

Board2D drawVertex(const RationalPoint vertex, Board2D& board, double size, Color c = Color::Red);
Board2D drawEdge(const std::pair<RationalPoint,RationalPoint>& edge, Board2D& board, double size, Color c = Color::Green);
Board2D drawFace(const std::vector<RationalPoint>& face, Board2D& board, double scale, Color c = Color::Blue);
Board2D drawPoint(const Z2i::Point vertex, Board2D& board, double size, Color c = Color::Black);

#endif // COMPLEX_H
