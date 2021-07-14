#include "Functions.h"
#include "Complex.h"
#include "AffineTransform.h"

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/images/ImageSelector.h"

#include "DGtal/io/colormaps/GradientColorMap.h"

using namespace std;
using namespace DGtal;

std::vector<Z2i::Point>
readImage(string filename) {
  std::vector<Z2i::Point> vecPtsTmp;
  typedef ImageSelector<Z2i::Domain, unsigned char >::Type Image;
  Image image = PGMReader<Image>::importPGM( filename );
  for ( Z2i::Domain::ConstIterator it = image.domain().begin(); it != image.domain().end(); ++it )
  {
    unsigned char val =  image (*it);
    if(val==255)
      vecPtsTmp.push_back(*it);
  }
  //Center the points
  Z2i::RealPoint cr = getBaryCenter(vecPtsTmp);
  Z2i::Point p, c(round(cr[0]),round(cr[1]));
  std::vector<Z2i::Point> vecPts;
  for(auto it=vecPtsTmp.begin(); it!=vecPtsTmp.end(); it++) {
    p = (*it) - c;
    vecPts.push_back(p);
  }
  return vecPts;
}

//Compute the set faces of CH1 belong (true) to CH2
std::vector<bool>
getBelongFace(const Complex& cH1, const Complex& cH2) {
  std::vector<bool> isBelong(cH1.face.size(),false);
  for(size_t it=0; it<cH1.face.size(); it++) {//for each face of cH1 (represented by its center)
    RationalPoint p = cH1.getFaceCenter(it);
    isBelong.at(it) = cH2.isInsideComplex(p);
  }
  return isBelong;
}
  
//Compute the set B2 of border faces of complex H
std::vector<bool>
getBorderFace(const Complex& cH) {
  std::vector<bool> isBorder(cH.face.size(), false);
  for(size_t it=0; it<cH.face.size(); it++) {//for each face h2
    isBorder.at(it) = cH.isBorderFace(it);
  }
  return isBorder;
}

std::vector<std::vector<int> > getNeigbourFace(const Complex& cH) {
  std::vector<std::vector<int> > vN = cH.getAllNeighbourFace();
  /*
  Board2D aBoard;
  cH.drawGrid(aBoard);
  cH.drawComplex(aBoard, Color::Gray);
  int id = 10;
  for(size_t it=0; it<vN.at(id).size(); it++)
    drawFace(cH.getFaceVertices(vN.at(id).at(it)), aBoard, 0.9, DGtal::Color::Green);
  drawFace(cH.getFaceVertices(id), aBoard, 0.9, DGtal::Color::Red);
  aBoard.saveSVG("../Illustration/testNeighbour.svg");
  aBoard.clear();
  */
  return vN;
}

std::vector<bool> getNeigbourFace(const Complex& gH, const std::vector<bool>& bordercH, const std::vector<bool>& iscH) {
  std::vector<bool> vecBorder (gH.face.size(),false);
  std::vector<std::vector<int> > vecNext;
  std::vector<std::vector<int> > vecNint;
  for(size_t it=0; it<bordercH.size(); it++) {
    if(bordercH.at(it)) { //if interieur border
      std::vector<int> vN = gH.getNeighbourFace(it);
      std::vector<int> vNExt;
      std::vector<int> vNInt;
      for(size_t it_bis = 0; it_bis<vN.size(); it_bis++) {
        if(!iscH.at(vN.at(it_bis))) //Not interieur neigbours
          vNExt.push_back(vN.at(it_bis));
        else
          vNInt.push_back(vN.at(it_bis));
      }
      assert(vNExt.size()!=0);
      vNExt.push_back(it);
      vecNext.push_back(vNExt);
      vecNint.push_back(vNInt);
    }
  }
  //get all int / ext border faces
  for(size_t it=0; it<vecNext.size(); it++)
    for(size_t it_bis=0; it_bis<vecNext.at(it).size(); it_bis++) {
      //if(!containElement(vecBorder, vecNext.at(it).at(it_bis)))
      //  vecBorder.push_back(vecNext.at(it).at(it_bis));
      vecBorder.at(vecNext.at(it).at(it_bis)) = true;
    }
  /*
  Board2D aBoard;
  gH.drawGrid(aBoard);
  for(size_t it=0; it<iscH.size(); it++) //dessin CH
    if(iscH.at(it))
      drawFace(gH.getFaceVertices(it), aBoard, 0.9, DGtal::Color::Green);
  for(size_t it=0; it<vecNext.size(); it++) {
    //Draw exterieur
    for(size_t it_bis=0; it_bis<vecNext.at(it).size(); it_bis++)
      drawFace(gH.getFaceVertices(vecNext.at(it).at(it_bis)), aBoard, 0.9, DGtal::Color::Red);
    //Draw interieur
    for(size_t it_bis=0; it_bis<vecNint.at(it).size(); it_bis++)
      drawFace(gH.getFaceVertices(vecNint.at(it).at(it_bis)), aBoard, 0.9, DGtal::Color::Blue);
  }
  aBoard.saveSVG("../Illustration/testBorderExt.svg");
  aBoard.clear();
  
  gH.drawGrid(aBoard);
  for(size_t it=0; it<iscH.size(); it++) //dessin CH
    if(iscH.at(it))
      drawFace(gH.getFaceVertices(it), aBoard, 0.9, DGtal::Color::Green);
  for(size_t it=0; it<vecBorder.size(); it++) { //Draw border
    if(vecBorder.at(it))
      drawFace(gH.getFaceVertices(it), aBoard, 0.9, DGtal::Color::Red);
  }
  aBoard.saveSVG("../Illustration/testBorder.svg");
  aBoard.clear();
  */
  return vecBorder;
}

//Compute the set B2 of border faces of gridH is simple wrt complex H
std::vector<bool>
getBorderFace(const Complex& gridH, const Complex& cH) {
  std::vector<bool> isBorder(gridH.face.size(), false);
  bool border = false;
  int idF = -1;
  for(size_t it=0; it<gridH.face.size(); it++) {
    RationalPoint p = gridH.getFaceCenter(it);
    border = false;
    idF = cH.isInsideComplex(p); //a face of cH
    if(idF != -1 && cH.isBorderFace(idF))//inside and and border (outside => not border by default)
      border = true;
    isBorder.at(it) = border;
  }
  return isBorder;
}
  
//Compute sigma fct indicating a pixel (its center) is majority inside (true) / outside (false) of the target complex H
bool
sigma(const RationalPoint p, const Complex& cH) {
  double a = 0.0;
  DGtal::Z2i::Point px(int(getRealValue(p.first)),int(getRealValue(p.second)));
  for(size_t it=0; it<cH.face.size(); it++) {
    RationalPoint p = cH.getFaceCenter(it);
    if(isInsidePixel(px, p))
      a += cH.getFaceArea(it);
  }
  if(a>=0.5)
    return true;
  return false;
}

std::vector<bool>
sigma(const Complex& gF, const Complex& cH) {
  std::vector<double> area(gF.face.size(), 0.0);
  for(size_t it=0; it<gF.face.size(); it++) {//for each face of gF, get all faceH inside it
    area.at(it) = sigma(gF.getFaceCenter(it), cH);
  }
  std::vector<bool> majority(gF.face.size(), false);
  for(size_t it=0; it<area.size(); it++) {
    if(area.at(it)>=0.5)
      majority.at(it) = true;
  }
  return majority;
}

//Compute zeta fct indicating fill area of pixel (its center) wrt to complexH
//insideFaces contains the set of cH face inside a pixel of gridF
double
zeta(const RationalPoint p, const Complex& cH, vector<size_t>& idInFaces) {
  DGtal::Z2i::Point px(int(getRealValue(p.first)),int(getRealValue(p.second)));
  double a = 0.0;
  for(size_t it=0; it<cH.face.size(); it++) {
    RationalPoint p = cH.getFaceCenter(it);
    if(isInsidePixel(px, p)) {
      a += cH.getFaceArea(it);
      idInFaces.push_back(it);
    }
  }
  if(a>=0.5)
    return a;
  return 1-a;
}

std::vector<double>
zeta(const Complex& gF, const Complex& cH, std::vector<vector<size_t> >& idInFaces) {
  std::vector<double> area(gF.face.size(), 0.0);
  double a = 0.0;
  for(size_t it=0; it<gF.face.size(); it++) {//for each face of gF, get all faceH inside it
    vector<size_t> inF;
    a = zeta(gF.getFaceCenter(it), cH, inF);
    idInFaces.push_back(inF);
    area.at(it) = a;
  }
  return area;
}

//Compute zeta fct indicating fill area of pixel (its center) wrt to cHinit and cHcurrent
double
zeta(const RationalPoint p, const Complex& cHinit, const Complex& cHcurrent) {
  DGtal::Z2i::Point px(int(getRealValue(p.first)),int(getRealValue(p.second)));
  double a = 0.0;
  for(size_t it=0; it<cHcurrent.face.size(); it++) {
    RationalPoint p = cHcurrent.getFaceCenter(it);
    if(isInsidePixel(px, p))
      a += cHcurrent.getFaceArea(it);
  }
  bool sigmaValue = sigma(p, cHinit);
  if(sigmaValue) //majority inside Hinit
    return a;
  return 1-a;
}

std::vector<double>
zeta(const Complex& gF, const Complex& cHinit, const Complex& cHcurrent) {
  std::vector<double> area(gF.face.size(), 0.0);
  for(size_t it=0; it<gF.face.size(); it++) //for each face of gF, compute zeta wrt complex H
    area.at(it) = zeta(gF.getFaceCenter(it), cHinit, cHcurrent);
  return area;
}

//Compute iota fct indicating whether a face h2 (represented by its center) is belong to complex H
bool
iota (RationalPoint h2, const Complex& cH) {
  bool isBelong = false;
  for(size_t it=0; it<cH.face.size() && !isBelong; it++) {
    if(isInsidePolygon(cH.getFaceVertices(it), h2))
      isBelong = true;
  }
  return isBelong;
}

//Compute iota fct indicating whether a face h2 of gridH is belong to complex H
std::vector<bool>
iota (const Complex& gridH, const Complex& cH) {
  std::vector<bool> isBelong(gridH.face.size(),false);
  for(size_t it=0; it<gridH.face.size(); it++) //for each face h2 (represented by its center)
    isBelong.at(it) = iota(gridH.getFaceCenter(it),cH);
  return isBelong;
}

//For collapse operation: whether a cell is BG / FG
std::vector<bool>
setFaceState(const Complex& cH, const Complex& gridH) {
  std::vector<bool> isBelong(gridH.face.size(),false);
  bool belong = false;
  size_t it_bis = 0;
  //std::vector<RationalPoint> f;
  //RationalPoint c;
  for(size_t it=0; it<cH.face.size(); it++) {
    belong = false;
    //f = cH.getFaceVertices(it);
    for(it_bis=0; it_bis<gridH.face.size() && !belong; it_bis++) {
      if(cH.face.at(it).size() == gridH.face.at(it_bis).size()) { //same size
        //c = gridH.getFaceCenter(it_bis);
        if (isEqual(cH.getFaceCenter(it), gridH.getFaceCenter(it_bis))) { //(isInsidePolygon(f, c)) {
          belong = true;
          isBelong.at(it_bis) = belong;
        }
      }
    }
    assert(belong==true);
  }
  /*
  std::vector<bool> isBelong(gridH.face.size(),false);
  bool belong = false;
  //std::vector<RationalPoint> f;
  RationalPoint c;
  for(size_t it=0; it<gridH.face.size(); it++) {
    belong = false;
    //f = gridH.getFaceVertices(it);
    c = gridH.getFaceCenter(it);
    for(size_t it_bis=0; it_bis<cH.face.size() && !belong; it_bis++) {
      if(gridH.face.at(it).size() == cH.face.at(it_bis).size()) { //same size
        //c = cH.getFaceCenter(it_bis);
        if(isInsidePolygon(cH.getFaceVertices(it_bis), c))//(isInsidePolygon(f, c))
          belong = true;
      }
    }
    isBelong.at(it) = belong;
  }
  */
  return isBelong;
  //return iota(gridH, cH);
}

//Compute tau fct indicating whether a face h2 of complex H is simple
std::vector<bool>
tau (const Complex& cH) {
  std::vector<bool> isSimple(cH.face.size(), false);
  for(size_t it=0; it<cH.face.size(); it++) {//for each face h2
    isSimple.at(it) = cH.isSimpleFace(it);
  }
  return isSimple;
}
//Compute tau fct indicating whether a face h2 of gridH is simple wrt complex H
std::vector<bool>
tau (const Complex& gridH, const Complex& cH) {
  std::vector<bool> isSimple(gridH.face.size(), false);
  bool simple = false;
  int idF = -1;
  for(size_t it=0; it<gridH.face.size(); it++) {
    RationalPoint p = gridH.getFaceCenter(it);
    simple = false;
    idF = cH.isInsideComplex(p);
    if(idF == -1)//outside => simple by default
      simple = true;
    else {
      if(cH.isSimpleFace(idF)) //inside and simple
        simple = true;
    }
    isSimple.at(it) = simple;
  }
  return isSimple;
}

//Find index of face f2 of gF that h2 (represented by its center) belong to
int
getIdFace(RationalPoint h2, const Complex& gF) {
  for(size_t it=0; it<gF.face.size(); it++) {//for each face h2
    if(isInsidePolygon(gF.getFaceVertices(it), h2))
      return it;
  }
  return -1;
}

//Compute priority fct (epsilon) indicating the priority of processing h2 (represented by its center)
double
epsilon (RationalPoint h2, const Complex& cHinit, const Complex& cHcurrent, const Complex& gF) {
  bool iotaValue = iota(h2, cHcurrent);
  //find index of f2 that h2 belong to compute its sigma and zeta values
  int idF = getIdFace(h2, gF);
  assert(idF != -1);
  RationalPoint cF = gF.getFaceCenter(idF);
  bool sigmaValue = sigma(cF, cHcurrent);
  double zetaV = zeta(cF, cHinit, cHcurrent);
  double iotaV = iotaValue==true ? 1 : -1;
  double sigmaV = sigmaValue==true ? 1 : -1;
  return -iotaV*sigmaV*zetaV;
}

std::vector<double>
epsilon (const Complex& cH, const Complex& cHinit, const Complex& cHcurrent, const Complex& gF) {
  std::vector<double> priority;
  double p;
  for(size_t it=0; it<cH.face.size(); it++) {//compute de priority for all cell of complex H
    p = epsilon(cH.getFaceCenter(it), cHinit, cHcurrent, gF);
    priority.push_back(p);
  }
  return priority;
}
  
size_t
maxIndexEpsilon(const std::vector<double>& epsilonValue, const std::vector<bool>& isBorder) {
  //int maxEpsilonIndex = std::max_element(epsilonValue.begin(),epsilonValue.end()) - epsilonValue.begin();
  assert(isBorder.size()==epsilonValue.size());
  double maxEpsilon = -1;
  int maxIndex = -1;
  
  for(size_t it=0; it<isBorder.size(); it++) {
    if(isBorder.at(it) && epsilonValue.at(it)>maxEpsilon) {//border cell
        maxEpsilon = epsilonValue.at(it);
        maxIndex = it;
    }
  }
  return maxIndex;
}

Complex initizeGauss(const Complex &cH, const Complex &gH) {
  std::vector<Z2i::Point> vP;
  RationalPoint p, pp;
  Z2i::Point px;
  for(size_t it=0; it<cH.face.size(); it++) {
    p = cH.getFaceCenter(it);
    px = Z2i::Point(round(getRealValue(p.first)),round(getRealValue(p.second)));
    pp = RationalPoint(px[0], px[1]);
    if(!containElement(vP, px) && isInsidePolygon(cH.getFaceVertices(it), pp))
      vP.push_back(px);
  }
  //reproject into gH
  Complex Hinit;
  Hinit.addGridLines(gH.vertical_lines, gH.horizontal_lines);
  for(size_t it=0; it<vP.size(); it++) {
    px = vP.at(it);
    for(size_t it_bis=0; it_bis<gH.face.size(); it_bis++) {
      p = gH.getFaceCenter(it_bis);
      if(isInsidePixel(px, p))
        Hinit.addFace(gH.getFaceVertices(it_bis), true, true);
    }
  }
  return Hinit;
}

void
testFunctions(int adj = 8, bool saveFile = true) {
  //Input set X
  std::vector<Z2i::Point> X = readImage("test.pgm");
  /*
  X.push_back(Z2i::Point(0,0));
  X.push_back(Z2i::Point(0,1));
  X.push_back(Z2i::Point(1,1));
  //X.push_back(Z2i::Point(2,2));
  //X.push_back(Z2i::Point(0,0));
  
  //multiple connected components
  X.push_back(Z2i::Point(-2,-2));
  X.push_back(Z2i::Point(-2,-1));
  X.push_back(Z2i::Point(-2,0));
  */
  /*
   //hole in connected component
   X.push_back(Z2i::Point(-1,0));
   X.push_back(Z2i::Point(1,0));
   X.push_back(Z2i::Point(0,1));
   X.push_back(Z2i::Point(0,-1));
   */
  assert(X.size()!=0);
  
  //Transformation
  //AffineTransform t(1,0,0,1,0,1,0,1); // Id
  //AffineTransform t(1,0.5,0.3,1.1); // scalling
  AffineTransform t(1.1,0.4,0.3,1.2,1,2,0,1); // affine
  //AffineTransform t(0.9,0.5,0.3,1.1,5,7,1,5); // affine two label in same pixel
  //AffineTransform t(1.2,0.5,0.3,1.1,5,7,1,5); // affine
  std::cout<<t<<std::endl;
  
  //Grid F
  Complex gridF;
  Z2i::Point p1 = Z2i::Point(-4,-4);
  Z2i::Point p2 = Z2i::Point(2,2);
  gridF.initGrid(p1,p2);
  
  //GridG = t(GridF)
  Complex gridG = transform(gridF,t);
  Z2i::Point tp1 = gridG.domain.first;
  Z2i::Point tp2 = gridG.domain.second;
  Complex tgridF;
  tgridF.initGrid(tp1,tp2);
  
  //Grid H = tgridF and gridG
  //vector<vector<vector<RationalPoint> > > h_faces = intersectionFace(tgridF,gridG);
  vector<vector<vector<RationalPoint> > > h_faces = intersectionFace(gridG, tgridF);
  vector<vector<RationalPoint> > h_cells; //all faces
  for(int it=0; it<h_faces.size(); it++)
    for(int it_bis=0; it_bis<h_faces.at(it).size(); it_bis++)
      h_cells.push_back(h_faces.at(it).at(it_bis));
  //Get cells of H completly embedded in tgridF
  Complex tgridH(h_cells);
  vector<vector<RationalPoint> > h_pure_cells; //all pure faces
  vector<vector<size_t> > idInFace; //all pure faces
  std::vector<double> zetaH = zeta(tgridF, tgridH, idInFace);
  std::vector<bool> sigmaH = sigma(tgridF, tgridH);
  assert(idInFace.size()==tgridF.face.size());
  for(size_t it=0; it<zetaH.size(); it++)
  if(fabs(zetaH.at(it)-1)<EPSILON && sigmaH.at(it)) {//pure and inside pixels
    //h_pure_cells.push_back(tgridF.getFaceVertices(it));
    for(size_t it_bis=0; it_bis<idInFace.at(it).size(); it_bis++) {
      size_t idF = idInFace.at(it).at(it_bis);
      h_pure_cells.push_back(tgridH.getFaceVertices(idF));
    }
  }
  //Embed into gridH
  Complex gridH(h_pure_cells);//gridH(h_cells);
  gridH.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  gridH.addGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  
  Board2D aBoard;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/gridH.svg");
    aBoard.clear();
  }
  
  cout<<"domainF = "<<gridF.domain.first<<" and "<<gridF.domain.second<<endl;
  cout<<"domainG = "<<gridG.domain.first<<" and "<<gridG.domain.second<<endl;
  cout<<"domainH = "<<gridH.domain.first<<" and "<<gridH.domain.second<<endl;
  
  //ComplexX
  Complex cX(X, adj);
  cX.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  if(saveFile) {
    cX.drawGrid(aBoard);
    cX.drawPixel(aBoard);
    ///cX.drawConectedComponent(aBoard);
    aBoard.saveSVG("../Illustration/objetX.svg");
    aBoard.clear();
  }
  
  //ComplexF
  Complex complexF;
  complexF.init(tp1,tp2);
  for(int it=0; it<X.size(); it++) {
    complexF.addCubicalFace(X.at(it));
    complexF.idConnectedComponent.push_back(cX.getFaceConnectedComponent(it));
  }
  
  if(saveFile) {
    complexF.drawGrid(aBoard);
    complexF.drawComplex(aBoard, Color::Gray);
    ///complexF.drawConectedComponent(aBoard);
    aBoard.saveSVG("../Illustration/complexF.svg");
    aBoard.clear();
    complexF.drawGrid(aBoard);
    cX.drawCenter(aBoard);
    aBoard.saveSVG("../Illustration/DigitX.svg");
    aBoard.clear();
  }
  
  //ComplexG = t(ComplexF)
  Complex complexG = transform(complexF,t);
  if(saveFile) {
    complexG.drawGrid(aBoard);
    complexG.drawComplex(aBoard, Color::Gray);
    ///complexG.drawConectedComponent(aBoard);
    aBoard.saveSVG("../Illustration/complexG.svg");
    aBoard.clear();
    
    complexF.drawComplex(aBoard, Color::Red);
    complexG.drawComplex(aBoard, Color::Blue);
    complexF.drawGrid(aBoard);
    complexG.drawGrid(aBoard);
    aBoard.saveSVG("../Illustration/complexFG.svg");
    aBoard.clear();
  }
  
  //Complex H :refine from gridF and gridG
  //Decomposition of R(X) --given by ComplexG-- by GridF => output: cells
  vector<vector<vector<RationalPoint> > > intersection_cells = intersectionFace(complexG, tgridF);
  vector<vector<RationalPoint> > cells; //all faces
  for(int it=0; it<intersection_cells.size(); it++)
    for(int it_bis=0; it_bis<intersection_cells.at(it).size(); it_bis++)
      cells.push_back(intersection_cells.at(it).at(it_bis));
  
  //Build the Complex H (decomposition of ComplexG(=t(ComplexF)) by GridF) from cells
  Complex complexH(cells);
  complexH.initGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  for(size_t it=0; it<complexH.face.size(); it++) {
    vector<RationalPoint> f = complexH.getFaceVertices(it);
    int idCC = complexG.getIndexConnectedComponent(f);
    complexH.idConnectedComponent.push_back(idCC);
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    gridG.drawGrid(aBoard);
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexH.svg");
    aBoard.saveSVG("../Illustration/complexMajorityH.svg");
    aBoard.clear();
  }
  
  Complex complexGauss = initizeGauss(complexH, gridH);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    complexGauss.drawComplex(aBoard, Color::Gray);
    complexH.drawComplex(aBoard, Color::Red);
    aBoard.saveSVG("../Illustration/complexGaussH.svg");
    aBoard.clear();
  }
  
  //Test getBelongFace fct
  std::vector<bool> solH = getBelongFace(gridH,complexH);
  assert(solH.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(solH.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/belongFaces.svg");
    aBoard.clear();
  }
  
  //Test getBorderFace fct
  /*
  std::vector<bool> borderH = getBorderFace(complexH);
  assert(borderH.size()==complexH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<complexH.face.size(); it++) {
      if(borderH.at(it))
        complexH.drawFace(aBoard, it, DGtal::Color::Red);
      else
        complexH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/borderFace.svg");
    aBoard.clear();
  }
  */
  /*
  std::vector<bool> borderH = getBorderFace(gridH);
  assert(borderH.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    //complexH.drawComplex(aBoard, Color::Gray);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(borderH.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/borderFace.svg");
    aBoard.clear();
  }
  */
  std::vector<bool> borderH = getBorderFace(gridH, complexH);
  assert(borderH.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(borderH.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/borderFace.svg");
    aBoard.clear();
  }
  
  //Test sigma fct
  std::vector<bool> sigmaValue = sigma(tgridF, complexH);
  assert(sigmaValue.size()==tgridF.face.size());
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    //gridG.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(sigmaValue.at(it))
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else
        tgridF.drawFace(aBoard, it, DGtal::Color::Red);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/sigma.svg");
    aBoard.clear();
  }
  
  //Test zeta fct
  vector<vector<size_t> > idInFaceH; //all pure faces
  std::vector<double> zetaValue = zeta(tgridF, complexH, idInFaceH);
  assert(zetaValue.size()==tgridF.face.size());
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(fabs(zetaValue.at(it)-1)<EPSILON) //pures pixels
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else { //non pures pixels
        if(zetaValue.at(it)>=0.5) //majority fill area
          tgridF.drawFace(aBoard, it, DGtal::Color::Blue);
        else
          assert(false);
      }
    }
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/zeta.svg");
    aBoard.clear();
  }
  
  //Test iota fct
  std::vector<bool> iotaValue = iota(gridH, complexH);
  assert(iotaValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it)) //cell belongs to solution
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/iota.svg");
    aBoard.clear();
  }
  
  //Test tau fct
  /*
  std::vector<bool> tauValue = tau(complexH);
  assert(tauValue.size()==complexH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<complexH.face.size(); it++) {
      if(tauValue.at(it)) //if simple
        complexH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        complexH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/tau.svg");
    aBoard.clear();
  }
  */
  std::vector<bool> tauValue = tau(gridH, complexH);
  assert(tauValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it)) //if simple
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/tau.svg");
    aBoard.clear();
  }
  
  //Test epsilon fct
  std::vector<double> epsilonValue = epsilon(gridH, complexH, complexH, tgridF);
  assert(epsilonValue.size()==gridH.face.size());
  if(saveFile) {
    // Creating colormap.
    GradientColorMap<int> cmap_grad( 0, 10 );
    cmap_grad.addColor( Color( 255, 255, 0 ) ); //from Yellow
    cmap_grad.addColor( Color( 255, 0, 0 ) ); //to Red
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(epsilonValue.at(it)>0) { //most priority and positif
        Color c = cmap_grad( int(epsilonValue.at(it)*10) );
        gridH.drawFace(aBoard, it, c);//DGtal::Color::Red
      }
      else {//negative
        if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
        else //less priority and negative
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
      }
    }
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/epsilon.svg");
    aBoard.clear();
  }
}

std::vector<Z2i::Point>
illustrationMajority(std::vector<Z2i::Point> X, AffineTransform t, int adj = 8, bool saveFile = true) {
  std::vector<Z2i::Point> Y;
  std::pair<Z2i::Point,Z2i::Point> bb = getBoundingBox(X);
  
  //Grid F
  Complex gridF;
  Z2i::Point p1 = bb.first - Z2i::Point(2,2);
  Z2i::Point p2 = bb.second + Z2i::Point(2,2);
  gridF.initGrid(p1,p2);
  
  //GridG = t(GridF)
  Complex gridG = transform(gridF,t);
  Z2i::Point tp1 = gridG.domain.first;
  Z2i::Point tp2 = gridG.domain.second;
  Complex tgridF;
  tgridF.initGrid(tp1,tp2);
  
  //Grid H = tgridF and gridG
  vector<vector<vector<RationalPoint> > > h_faces = intersectionFace(gridG, tgridF);
  vector<vector<RationalPoint> > h_cells; //all faces
  for(int it=0; it<h_faces.size(); it++)
    for(int it_bis=0; it_bis<h_faces.at(it).size(); it_bis++)
      h_cells.push_back(h_faces.at(it).at(it_bis));
  //Get cells of H completly embedded in tgridF
  Complex tgridH(h_cells);
  vector<vector<RationalPoint> > h_pure_cells; //all pure faces
  vector<vector<size_t> > idInFace; //all pure faces
  std::vector<double> zetaH = zeta(tgridF, tgridH, idInFace);
  std::vector<bool> sigmaH = sigma(tgridF, tgridH);
  assert(idInFace.size()==tgridF.face.size());
  for(size_t it=0; it<zetaH.size(); it++)
  if(fabs(zetaH.at(it)-1)<EPSILON && sigmaH.at(it)) {//pure and inside pixels
    //h_pure_cells.push_back(tgridF.getFaceVertices(it));
    for(size_t it_bis=0; it_bis<idInFace.at(it).size(); it_bis++) {
      size_t idF = idInFace.at(it).at(it_bis);
      h_pure_cells.push_back(tgridH.getFaceVertices(idF));
    }
  }
  //Embed into gridH
  Complex gridH(h_pure_cells);//gridH(h_cells);
  gridH.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  gridH.addGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  
  Board2D aBoard;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/gridH.svg");
    aBoard.clear();
  }
  
  Complex cX(X, adj);
  cX.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  if(saveFile) {
    cX.drawGrid(aBoard);
    cX.drawPixel(aBoard);
    aBoard.saveSVG("../Illustration/objetX.svg");
    aBoard.clear();
  }
  
  //ComplexF
  Complex complexF;
  complexF.init(tp1,tp2);
  for(int it=0; it<X.size(); it++) {
    complexF.addCubicalFace(X.at(it));
    complexF.idConnectedComponent.push_back(cX.getFaceConnectedComponent(it));
  }
  
  if(saveFile) {
    complexF.drawGrid(aBoard);
    complexF.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexF.svg");
    aBoard.clear();
  }
  
  //ComplexG = t(ComplexF)
  Complex complexG = transform(complexF,t);
  if(saveFile) {
    complexG.drawGrid(aBoard);
    complexG.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexG.svg");
    aBoard.clear();
    
    complexF.drawComplex(aBoard, Color::Red);
    complexG.drawComplex(aBoard, Color::Blue);
    complexF.drawGrid(aBoard);
    complexG.drawGrid(aBoard);
    aBoard.saveSVG("../Illustration/complexFG.svg");
    aBoard.clear();
  }
  
  //Complex H :refine from gridF and gridG
  //Decomposition of R(X) --given by ComplexG-- by GridF => output: cells
  vector<vector<vector<RationalPoint> > > intersection_cells = intersectionFace(complexG, tgridF);
  vector<vector<RationalPoint> > cells; //all faces
  for(int it=0; it<intersection_cells.size(); it++)
    for(int it_bis=0; it_bis<intersection_cells.at(it).size(); it_bis++)
      cells.push_back(intersection_cells.at(it).at(it_bis));
  
  //Build the Complex H (decomposition of ComplexG(=t(ComplexF)) by GridF) from cells
  Complex complexH(cells);
  complexH.initGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  for(size_t it=0; it<complexH.face.size(); it++) {
    vector<RationalPoint> f = complexH.getFaceVertices(it);
    int idCC = complexG.getIndexConnectedComponent(f);
    complexH.idConnectedComponent.push_back(idCC);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexH.svg");
    aBoard.clear();
  }
  
  //Line 1
  Complex Hinit = complexH;
  std::vector<bool> stateFaceHinit = setFaceState(Hinit, gridH);
  assert(stateFaceHinit.size()==gridH.face.size());
  Complex Hcurrent = complexH;
  std::vector<bool> stateFaceHcurrent = stateFaceHinit;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(stateFaceHinit.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/stageFace.svg");
    aBoard.clear();
  }
  
  //Optimization of data structure
  //Compute h2 area of Hinit (wrt Hgrid)
  std::vector<double> areaFace(gridH.face.size(),0.0);
  for(size_t it=0; it<gridH.face.size(); it++) {
    if(stateFaceHinit.at(it))
      areaFace.at(it) = gridH.getFaceArea(it);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(areaFace.at(it)>0) //face of Hinit
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/areaFace.svg");
    aBoard.clear();
  }
  //Compute fill area f2 of Hinit (wrt Hgrid)
  std::vector<double> fillAreaFace(gridH.face.size(),0.0);
  std::vector<bool> fillAreaDone(gridH.face.size(),false);
  std::pair<Z2i::Point,Z2i::Point> d = tgridF.domain;
  for(int x=d.first[0]; x<=d.second[0]; x++) {
    for(int y=d.first[1]; y<=d.second[1]; y++) {
      Z2i::Point px = Z2i::Point(x,y); //pixel center
      double a = 0;
      std::vector<size_t> faces;
      for(size_t it=0; it<gridH.face.size(); it++) { //find all face in the same pixel
        if(!fillAreaDone.at(it) && isInsidePixel(px, gridH.getFaceCenter(it))) {
          if(stateFaceHinit.at(it)) // a face of Hinit
            a += gridH.getFaceArea(it);
          faces.push_back(it);
        }
        for(size_t it=0; it<faces.size(); it++) { //update area of all faces in the pixel
          fillAreaFace.at(faces.at(it)) = a;
          fillAreaDone.at(faces.at(it)) = true;
        }
      }
    }
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(fillAreaFace.at(it)>0.5) //majority interieur
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else {
        if(fillAreaFace.at(it)<EPSILON || fabs(fillAreaFace.at(it)-1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        else //majority exterieur
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
    Hcurrent.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/fillAreaFace.svg");
    aBoard.clear();
  }
  
  //Line 2-4
  //Sigma fct on F2
  std::vector<bool> sigmaValue = sigma(tgridF, Hcurrent);
  assert(sigmaValue.size()==tgridF.face.size());
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(sigmaValue.at(it))
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else
        tgridF.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hcurrent.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/sigma.svg");
    aBoard.clear();
  }
  //Zeta fct on F2
  vector<vector<size_t> > idInFaceH; //all pure faces
  std::vector<double> zetaValue = zeta(tgridF, complexH, idInFaceH);
  assert(zetaValue.size()==tgridF.face.size());
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(fabs(zetaValue.at(it)-1)<EPSILON) //pures pixels
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else { //non pures pixels
        if(zetaValue.at(it)>=0.5) //majority fill area
          tgridF.drawFace(aBoard, it, DGtal::Color::Blue);
        else
          assert(false);
      }
    }
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/zeta.svg");
    aBoard.clear();
  }
  
  //Line 5-8
  //Iota fct on h2
  std::vector<bool> iotaValue = iota(gridH, Hcurrent);
  assert(iotaValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/iota.svg");
    aBoard.clear();
  }
  
  //Epsilon fct on h2
  std::vector<double> epsilonValue = epsilon(gridH, Hinit, Hcurrent, tgridF);
  assert(epsilonValue.size()==gridH.face.size());
  if(saveFile) {
    // Creating colormap.
    GradientColorMap<int> cmap_grad( 0, 10 );
    cmap_grad.addColor( Color( 255, 255, 0 ) ); //from Yellow
    cmap_grad.addColor( Color( 255, 0, 0 ) ); //to Red
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(epsilonValue.at(it)>0) { //most priority and positif
        Color c = cmap_grad( int(epsilonValue.at(it)*10) );
        gridH.drawFace(aBoard, it, c);//DGtal::Color::Red
      }
      else {//negative
        if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
        else //less priority and negative
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
      }
    }
    //complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/epsilon.svg");
    aBoard.clear();
  }
  
  //Tau fct on h2
  std::vector<bool> tauValue = tau(gridH, Hcurrent);
  assert(tauValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/tau.svg");
    aBoard.clear();
  }
  
  //Line 9: Compute B2 of Hcurrent (of same index as gridH)
  std::vector<bool> borderH = getBorderFace(gridH, Hcurrent);
  assert(borderH.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(borderH.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/borderFace.svg");
    aBoard.clear();
  }
  
  //Line 10-11: find border face of max priority
  //int maxEpsilonIndex = std::max_element(epsilonValue.begin(),epsilonValue.end()) - epsilonValue.begin();
  double maxEpsilon = -1;
  int maxEpsilonIndex = -1;
  
  for(size_t it=0; it<borderH.size(); it++) {
    if(borderH.at(it)) {//border cell
      if(epsilonValue.at(it)>maxEpsilon) {
        maxEpsilon = epsilonValue.at(it);
        maxEpsilonIndex = it;
      }
    }
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, DGtal::Color::Green);
    gridH.drawFace(aBoard, maxEpsilonIndex, DGtal::Color::Red);
    aBoard.saveSVG("../Illustration/maxEpsilon.svg");
    aBoard.clear();
  }
  
  return Y;
}

std::vector<Z2i::Point>
illustrationMajority2(std::vector<Z2i::Point> X, AffineTransform t, int adj = 8, bool saveFile = true) {
  std::vector<Z2i::Point> Y;
  std::pair<Z2i::Point,Z2i::Point> bb = getBoundingBox(X);
  string filename;
  
  //Grid F
  Complex gridF;
  Z2i::Point p1 = bb.first - Z2i::Point(2,2);
  Z2i::Point p2 = bb.second + Z2i::Point(2,2);
  gridF.initGrid(p1,p2);
  
  //GridG = t(GridF)
  Complex gridG = transform(gridF,t);
  Z2i::Point tp1 = gridG.domain.first;
  Z2i::Point tp2 = gridG.domain.second;
  Complex tgridF;
  tgridF.initGrid(tp1,tp2);
  
  //Grid H = tgridF and gridG
  vector<vector<vector<RationalPoint> > > h_faces = intersectionFace(gridG, tgridF);
  vector<vector<RationalPoint> > h_cells; //all faces
  for(int it=0; it<h_faces.size(); it++)
    for(int it_bis=0; it_bis<h_faces.at(it).size(); it_bis++)
      h_cells.push_back(h_faces.at(it).at(it_bis));
  //Get cells of H completly embedded in tgridF
  Complex tgridH(h_cells);
  vector<vector<RationalPoint> > h_pure_cells; //all pure faces
  vector<vector<size_t> > idInFace; //all pure faces
  std::vector<double> zetaH = zeta(tgridF, tgridH, idInFace);
  std::vector<bool> sigmaH = sigma(tgridF, tgridH);
  assert(idInFace.size()==tgridF.face.size());
  for(size_t it=0; it<zetaH.size(); it++)
  if(fabs(zetaH.at(it)-1)<EPSILON && sigmaH.at(it)) {//pure and inside pixels
    //h_pure_cells.push_back(tgridF.getFaceVertices(it));
    for(size_t it_bis=0; it_bis<idInFace.at(it).size(); it_bis++) {
      size_t idF = idInFace.at(it).at(it_bis);
      h_pure_cells.push_back(tgridH.getFaceVertices(idF));
    }
  }
  //Embed into gridH
  Complex gridH(h_pure_cells);//gridH(h_cells);
  gridH.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  gridH.addGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  
  Board2D aBoard;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/gridH.svg");
    aBoard.clear();
  }
  
  Complex cX(X, adj);
  cX.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  if(saveFile) {
    cX.drawGrid(aBoard);
    cX.drawPixel(aBoard);
    aBoard.saveSVG("../Illustration/objetX.svg");
    aBoard.clear();
  }
  
  //ComplexF
  Complex complexF;
  complexF.init(tp1,tp2);
  for(int it=0; it<X.size(); it++) {
    complexF.addCubicalFace(X.at(it));
    complexF.idConnectedComponent.push_back(cX.getFaceConnectedComponent(it));
  }
  
  if(saveFile) {
    complexF.drawGrid(aBoard);
    complexF.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexF.svg");
    aBoard.clear();
  }
  
  //ComplexG = t(ComplexF)
  Complex complexG = transform(complexF,t);
  if(saveFile) {
    complexG.drawGrid(aBoard);
    complexG.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexG.svg");
    aBoard.clear();
    
    complexF.drawComplex(aBoard, Color::Red);
    complexG.drawComplex(aBoard, Color::Blue);
    complexF.drawGrid(aBoard);
    complexG.drawGrid(aBoard);
    aBoard.saveSVG("../Illustration/complexFG.svg");
    aBoard.clear();
  }
  
  //Complex H :refine from gridF and gridG
  //Decomposition of R(X) --given by ComplexG-- by GridF => output: cells
  vector<vector<vector<RationalPoint> > > intersection_cells = intersectionFace(complexG, tgridF);
  vector<vector<RationalPoint> > cells; //all faces
  for(int it=0; it<intersection_cells.size(); it++)
    for(int it_bis=0; it_bis<intersection_cells.at(it).size(); it_bis++)
      cells.push_back(intersection_cells.at(it).at(it_bis));
  
  //Build the Complex H (decomposition of ComplexG(=t(ComplexF)) by GridF) from cells
  Complex complexH(cells);
  complexH.initGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  for(size_t it=0; it<complexH.face.size(); it++) {
    vector<RationalPoint> f = complexH.getFaceVertices(it);
    int idCC = complexG.getIndexConnectedComponent(f);
    complexH.idConnectedComponent.push_back(idCC);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexH.svg");
    aBoard.saveSVG("../Illustration/complexMajorityH.svg");
    aBoard.clear();
  }
  
  //Line 1
  Complex Hinit = complexH;
  std::vector<bool> stateFaceHinit = setFaceState(Hinit, gridH);
  int sum = 0;
  for(size_t it = 0; it<stateFaceHinit.size(); it++)
    if(stateFaceHinit.at(it))
      sum++;
  assert(sum==Hinit.face.size());
  assert(stateFaceHinit.size()==gridH.face.size());
  std::vector<bool> stateFaceHcurrent = stateFaceHinit;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(stateFaceHinit.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/stageFaceMajority.svg");
    aBoard.clear();
  }
  
  //Optimization of data structure
  //Compute h2 area of Hinit (wrt Hgrid)
  std::vector<double> areaFace(gridH.face.size(),0.0);
  for(size_t it=0; it<gridH.face.size(); it++) {
    if(stateFaceHinit.at(it))
      areaFace.at(it) = gridH.getFaceArea(it);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(areaFace.at(it)>0) //face of Hinit
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/areaFaceMajority.svg");
    aBoard.clear();
  }
  //Compute fill area f2 of Hinit (wrt Hgrid)
  std::vector<double> fillAreaFace(gridH.face.size(),0.0);
  std::vector<bool> fillAreaDone(gridH.face.size(),false);
  std::vector<std::vector<int> > Phi; //Phi (f2) = {h2}
  std::vector<int> phi (gridH.face.size(), -1); //phi (h2) = f2
  std::pair<Z2i::Point,Z2i::Point> d = tgridF.domain;
  for(size_t it=0; it<tgridF.face.size(); it++) { //for each pixel of gridF
    RationalPoint p = tgridF.getFaceCenter(it);
    Z2i::Point px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second))); //pixel center
    double a = 0;
    std::vector<int> insideFaces;
    for(size_t it_bis=0; it_bis<gridH.face.size(); it_bis++) { //find all face h2 in the same pixel
      if(!fillAreaDone.at(it_bis) && isInsidePixel(px, gridH.getFaceCenter(it_bis))) {
        if(stateFaceHinit.at(it_bis)) // a face of Hinit
          a += gridH.getFaceArea(it_bis);
        phi.at(it_bis) = it;
        insideFaces.push_back(it_bis);
      }
      for(size_t it_bis=0; it_bis<insideFaces.size(); it_bis++) { //update area of all faces in the pixel
        fillAreaFace.at(insideFaces.at(it_bis)) = a;
        fillAreaDone.at(insideFaces.at(it_bis)) = true;
      }
    }
    Phi.push_back(insideFaces);
  }
  assert(Phi.size()==tgridF.face.size());
  sum = 0;
  for(size_t it=0; it<Phi.size(); it++)
    sum += Phi.at(it).size();
  assert(sum==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(fillAreaFace.at(it)>0.5) //majority interieur
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else {
        if(fillAreaFace.at(it)<EPSILON || fabs(fillAreaFace.at(it)-1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        else //majority exterieur
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/fillAreaFaceMajority.svg");
    aBoard.clear();
  }
  
  //Line 2-4
  //Sigma fct on F2
  std::vector<bool> sigmaValue(tgridF.face.size(), false);
  for(size_t it=0; it<tgridF.face.size(); it++) {
    if(Phi.at(it).size()!=0 && fillAreaFace.at(Phi.at(it).front())>=0.5) //contains Hinit cell and in majority
      sigmaValue.at(it) = true;
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(sigmaValue.at(it))
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else
        tgridF.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/sigmaMajority.svg");
    aBoard.clear();
  }
  
  //Zeta fct on F2
  std::vector<double> zetaValue(tgridF.face.size(), 0.0);
  double z = 0.0, s = -1.0;
  for(size_t it=0; it<tgridF.face.size(); it++) {
    s = sigmaValue.at(it) ? 1.0 : -1.0;
    z = (1 - s) / 2.0;
    if(Phi.at(it).size()!=0) //contains Hinit cell
       z += s*fillAreaFace.at(Phi.at(it).front());
    zetaValue.at(it) = z;
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(fabs(zetaValue.at(it)-1)<EPSILON) //pures pixels
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else { //non pures pixels
        if(zetaValue.at(it)>=0.5) //majority fill area
          tgridF.drawFace(aBoard, it, DGtal::Color::Blue);
        else
          assert(false);
      }
    }
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/zetaMajority.svg");
    aBoard.clear();
  }
  
  //Line 5-8
  //Iota fct on h2
  std::vector<bool> iotaValue = stateFaceHcurrent;
  assert(iotaValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/iotaMajority.svg");
    aBoard.clear();
  }
  
  //Epsilon fct on h2
  std::vector<double> epsilonValue(gridH.face.size(), 0.0);
  double s1, s2;
  for(size_t it=0; it<gridH.face.size(); it++) {
    s1 = iotaValue.at(it) ? 1.0 : -1.0;
    s2 = sigmaValue.at(phi.at(it)) ? 1.0 : -1.0;
    epsilonValue.at(it) = -s1*s2*zetaValue.at(phi.at(it));
  }
  if(saveFile) {
    // Creating colormap.
    GradientColorMap<int> cmap_grad( 0, 10 );
    cmap_grad.addColor( Color( 255, 255, 0 ) ); //from Yellow
    cmap_grad.addColor( Color( 255, 0, 0 ) ); //to Red
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(epsilonValue.at(it)>0) { //most priority and positif
        Color c = cmap_grad( int(epsilonValue.at(it)*10) );
        gridH.drawFace(aBoard, it, c);//DGtal::Color::Red
      }
      else {//negative
        if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
        else //less priority and negative
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
      }
    }
    //complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/epsilonMajority.svg");
    aBoard.clear();
  }
  
  //Tau fct on h2
  std::vector<bool> tauValue = tau(gridH, Hinit);
  assert(tauValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/tauMajority.svg");
    aBoard.clear();
  }
  
  //Line 9: Compute B2 of Hcurrent (of same index as gridH)
  std::vector<bool> borderH = getBorderFace(gridH, Hinit);
  std::vector<bool> borderHcurrent = getNeigbourFace(gridH, borderH, iotaValue);
  assert(borderHcurrent.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(borderHcurrent.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/borderFaceMajority.svg");
    aBoard.clear();
  }
  //Update tau value of border faces
  for(size_t it=0; it<borderHcurrent.size(); it++) {
    if(borderHcurrent.at(it) && !iotaValue.at(it)) { //border face but not belong to Hcurrent
      Complex tmp = Hinit;
      tmp.addFace(gridH.getFaceVertices(it), true, true);
      //Verify new face simplicity
      tauValue.at(it) = tmp.isSimpleFace(tmp.face.size()-1);
    }
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/tauMajority2.svg");
    aBoard.clear();
  }
  
  //Line 10-11: find border face of max priority
  size_t maxIndex = maxIndexEpsilon(epsilonValue, borderHcurrent);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, DGtal::Color::Green);
    gridH.drawFace(aBoard, maxIndex, DGtal::Color::Red);
    aBoard.saveSVG("../Illustration/maxEpsilonMajority.svg");
    aBoard.clear();
  }
  //Line 12-26: Loop
  int nbIter=0;
  double a, h2 = maxIndex, f2 = phi.at(h2);
  while(epsilonValue.at(h2)>=0) {
    std::cout<<"nbIter="<<nbIter<<": h2="<<h2<<", f2="<<f2;
    //Line 13
    Complex Hcurrent;
    //Hcurrent.initGridLines(gridH.vertical_lines, gridH.horizontal_lines);
    Hcurrent.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
    std::vector<int> inFace = Phi.at(f2);
    a = fillAreaFace.at(h2);
    if(tauValue.at(h2)) {//if simple cell
      //Line 14-15: Update iota of h2
      //std::cout<<"Before: iotaValue.at(h2)="<<iotaValue.at(h2)<<std::endl;
      iotaValue.at(h2) = !iotaValue.at(h2); //Inverse state face
      //std::cout<<"After: iotaValue.at(h2)="<<iotaValue.at(h2)<<std::endl;
      for(size_t it=0; it<iotaValue.size(); it++)
        if(iotaValue.at(it))
          Hcurrent.addFace(gridH.getFaceVertices(it), true, true);
      std::cout<<", Hcurrent="<<Hcurrent.face.size()<<", cH="<<gridH.face.size()<<std::endl;
      if(saveFile) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(iotaValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        Hcurrent.drawComplex(aBoard, Color::Gray);
        filename = "../Illustration/HcurrentMajority_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
      
      //Line 16 => out of condition
      //std::cout<<"Before: fillAreaFace.at(h2)="<<fillAreaFace.at(maxIndex)<<std::endl;
      if(iotaValue.at(h2))
        a += areaFace.at(h2);
      else
        a -= areaFace.at(h2);
      for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) //update fill area for all inside faces
        fillAreaFace.at(inFace.at(it_bis)) = a;
      //std::cout<<"After: fillAreaFace.at(h2)="<<fillAreaFace.at(maxIndex)<<std::endl;
      //Line 17-19
      std::vector<bool> tauHcurrent = tau(gridH, Hcurrent);
      for(size_t it=0; it<gridH.face.size(); it++)
        tauValue.at(it) = tauHcurrent.at(it);
      /*
      if(saveFile) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(tauValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        complexH.drawComplex(aBoard, Color::Gray);
        filename = "../Illustration/tau_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
      */
      //Line 19 ????
    }
    else { //Line 20-21
      sigmaValue.at(f2) = !sigmaValue.at(f2);
      for(size_t it=0; it<iotaValue.size(); it++)
        if(iotaValue.at(it))
          Hcurrent.addFace(gridH.getFaceVertices(it), true, true);
      std::cout<<", Hcurrent="<<Hcurrent.face.size()<<", cH="<<gridH.face.size()<<std::endl;
      if(saveFile) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(iotaValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        Hcurrent.drawComplex(aBoard, Color::Gray);
        filename = "../Illustration/HcurrentGauss_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
    }
    //Line 16/22: Update zeta of f2
    //std::cout<<"Before: zetaValue.at(f2)="<<zetaValue.at(f2)<<std::endl;
    s = sigmaValue.at(f2) ? 1 : -1;
    zetaValue.at(f2) = (1.0-s)/2.0 + s*a;
    //std::cout<<"After: zetaValue.at(f2)="<<zetaValue.at(f2)<<std::endl;
    //Line 23-25
    //std::vector<int> inFace = Phi.at(f2);
    s2 = sigmaValue.at(f2) ? 1 : -1;
    a = zetaValue.at(f2);
    for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) {
      //std::cout<<"Before: epsilonValue.at("<<inFace.at(it_bis)<<")="<<epsilonValue.at(inFace.at(it_bis))<<std::endl;
      s1 = iotaValue.at(inFace.at(it_bis)) ? 1 : -1;
      epsilonValue.at(inFace.at(it_bis)) = -s1*s2*a;
      //std::cout<<"After: epsilonValue.at("<<inFace.at(it_bis)<<")="<<epsilonValue.at(inFace.at(it_bis))<<std::endl;
      //Line 25 out of loop
    }
    if(saveFile) {
      gridH.drawGrid(aBoard);
      for(size_t it=0; it<gridH.face.size(); it++) {
        if(epsilonValue.at(it)>0) { //most priority and positif
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        else {//negative
          if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
          else //less priority and negative
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
        }
      }
      Hcurrent.drawComplex(aBoard, Color::Gray);
      filename = "../Illustration/epsilonMajority_s" + std::to_string(nbIter) + ".svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
    
    //Line 25: Update border faces
    std::vector<bool> borderH = getBorderFace(gridH, Hcurrent);
    std::vector<bool> borderHTmp = getNeigbourFace(gridH, borderH, iotaValue);
    for(size_t it_bis=0; it_bis<borderHTmp.size(); it_bis++)
      borderHcurrent.at(it_bis) = borderHTmp.at(it_bis);
    
    for(size_t it=0; it<borderHcurrent.size(); it++) {
      if(borderHcurrent.at(it) && !iotaValue.at(it)) { //border face but not belong to Hcurrent
        Complex tmp = Hcurrent;
        tmp.addFace(gridH.getFaceVertices(it), true, true);
        //Verify new face simplicity
        tauValue.at(it) = tmp.isSimpleFace(tmp.face.size()-1);
      }
    }
    if(saveFile) {
      gridH.drawGrid(aBoard);
      for(size_t it=0; it<gridH.face.size(); it++) {
        if(tauValue.at(it))
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
        else
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
      }
      Hcurrent.drawComplex(aBoard, Color::Gray);
      filename = "../Illustration/tauGauss_s" + std::to_string(nbIter) + ".svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
    
    //Line 26
    h2 = maxIndexEpsilon(epsilonValue, borderHcurrent);
    f2 = phi.at(h2);
    nbIter++;
  }
  
  Complex Hopt;
  Hopt.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  for(size_t it=0; it<iotaValue.size(); it++)
    if(iotaValue.at(it))
      Hopt.addFace(gridH.getFaceVertices(it), true, true);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hopt.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/HoptMajority.svg");
    aBoard.clear();
  }
  //Retreive integer points
  RationalPoint p;
  Z2i::Point px;
  for(size_t it=0; it<iotaValue.size(); it++) {
    if(iotaValue.at(it)) {
      p = tgridF.getFaceCenter(phi.at(it));
      px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second)));
      if(!containElement(Y, px))
        Y.push_back(px);
    }
  }
  if(saveFile) {
    Complex cY(Y, adj);
    cY.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
    if(saveFile) {
      cY.drawGrid(aBoard);
      cY.drawPixel(aBoard);
      aBoard.saveSVG("../Illustration/objetYMajority.svg");
      aBoard.clear();
    }
  }
  return Y;
}

std::vector<Z2i::Point>
illustrationMajority3(std::vector<Z2i::Point> X, AffineTransform t, string file = "test", int border = 2, int adj = 8, bool saveExpe = true, bool saveFile = false, bool useSaveFile = false) {
  std::vector<Z2i::Point> Y;
  std::pair<Z2i::Point,Z2i::Point> bb = getBoundingBox(X);
  
  //Grid F
  Complex gridF;
  Z2i::Point p1 = bb.first - Z2i::Point(border,border);
  Z2i::Point p2 = bb.second + Z2i::Point(border,border);
  gridF.initGrid(p1,p2);
  
  //GridG = t(GridF)
  Complex gridG = transform(gridF,t);
  Z2i::Point tp1 = gridG.domain.first;
  Z2i::Point tp2 = gridG.domain.second;
  Complex tgridF;
  tgridF.initGrid(tp1,tp2);
  
  //Grid H = tgridF and gridG
  vector<vector<vector<RationalPoint> > > h_faces = intersectionFace(gridG, tgridF);
  vector<vector<RationalPoint> > h_cells; //all faces
  for(int it=0; it<h_faces.size(); it++)
    for(int it_bis=0; it_bis<h_faces.at(it).size(); it_bis++)
      h_cells.push_back(h_faces.at(it).at(it_bis));
  //Get cells of H completly embedded in tgridF
  Complex tgridH(h_cells);
  vector<vector<RationalPoint> > h_pure_cells; //all pure faces
  vector<vector<size_t> > idInFace; //all pure faces
  std::vector<double> zetaH = zeta(tgridF, tgridH, idInFace);
  std::vector<bool> sigmaH = sigma(tgridF, tgridH);
  assert(idInFace.size()==tgridF.face.size());
  for(size_t it=0; it<zetaH.size(); it++)
  if(fabs(zetaH.at(it)-1)<EPSILON && sigmaH.at(it)) {//pure and inside pixels
    //h_pure_cells.push_back(tgridF.getFaceVertices(it));
    for(size_t it_bis=0; it_bis<idInFace.at(it).size(); it_bis++) {
      size_t idF = idInFace.at(it).at(it_bis);
      h_pure_cells.push_back(tgridH.getFaceVertices(idF));
    }
  }
  //Embed into gridH
  Complex gridH(h_pure_cells);//gridH(h_cells);
  gridH.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  gridH.addGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  
  Board2D aBoard;
  string filename;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, Color::Gray);
    filename = file + "_gridH.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  Complex cX(X, adj);
  cX.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  if(saveFile || saveExpe) {
    cX.drawGrid(aBoard);
    cX.drawPixel(aBoard);
    filename = file + "_objetX.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //ComplexF
  Complex complexF;
  complexF.init(tp1,tp2);
  for(int it=0; it<X.size(); it++) {
    complexF.addCubicalFace(X.at(it));
    complexF.idConnectedComponent.push_back(cX.getFaceConnectedComponent(it));
  }
  
  if(saveFile) {
    complexF.drawGrid(aBoard);
    complexF.drawComplex(aBoard, Color::Gray);
    filename = file + "_complexF.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //ComplexG = t(ComplexF)
  Complex complexG = transform(complexF,t);
  if(saveFile) {
    complexG.drawGrid(aBoard);
    complexG.drawComplex(aBoard, Color::Gray);
    filename = file + "_complexG.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
    
    complexF.drawComplex(aBoard, Color::Red);
    complexG.drawComplex(aBoard, Color::Blue);
    complexF.drawGrid(aBoard);
    complexG.drawGrid(aBoard);
    filename = file + "_complexFG.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Complex H :refine from gridF and gridG
  //Decomposition of R(X) --given by ComplexG-- by GridF => output: cells
  vector<vector<vector<RationalPoint> > > intersection_cells = intersectionFace(complexG, tgridF);
  vector<vector<RationalPoint> > cells; //all faces
  for(int it=0; it<intersection_cells.size(); it++)
    for(int it_bis=0; it_bis<intersection_cells.at(it).size(); it_bis++)
      cells.push_back(intersection_cells.at(it).at(it_bis));
  
  //Build the Complex H (decomposition of ComplexG(=t(ComplexF)) by GridF) from cells
  Complex complexH(cells);
  complexH.initGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  for(size_t it=0; it<complexH.face.size(); it++) {
    vector<RationalPoint> f = complexH.getFaceVertices(it);
    int idCC = complexG.getIndexConnectedComponent(f);
    complexH.idConnectedComponent.push_back(idCC);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_complexH.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  if(saveExpe) {
    gridH.drawGrid(aBoard);
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_HinitMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 1
  Complex Hinit = complexH;
  std::vector<bool> stateFaceHinit = setFaceState(Hinit, gridH);
  int sum = 0;
  for(size_t it = 0; it<stateFaceHinit.size(); it++)
    if(stateFaceHinit.at(it))
      sum++;
  assert(sum==Hinit.face.size());
  assert(stateFaceHinit.size()==gridH.face.size());
  std::vector<bool> stateFaceHcurrent = stateFaceHinit;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(stateFaceHinit.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_stageFaceMajority.svg";
    aBoard.clear();
  }
  
  //Optimization of data structure
  //Compute h2 area of Hinit (wrt Hgrid)
  std::vector<double> areaFace(gridH.face.size(),0.0);
  for(size_t it=0; it<gridH.face.size(); it++) {
    if(stateFaceHinit.at(it))
      areaFace.at(it) = gridH.getFaceArea(it);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(areaFace.at(it)>0) //face of Hinit
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_areaFaceMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  //Compute fill area f2 of Hinit (wrt Hgrid)
  std::vector<double> fillAreaFace(gridH.face.size(),0.0);
  std::vector<bool> fillAreaDone(gridH.face.size(),false);
  std::vector<std::vector<int> > Phi; //Phi (f2) = {h2}
  std::vector<int> phi (gridH.face.size(), -1); //phi (h2) = f2
  std::pair<Z2i::Point,Z2i::Point> d = tgridF.domain;
  for(size_t it=0; it<tgridF.face.size(); it++) { //for each pixel of gridF
    RationalPoint p = tgridF.getFaceCenter(it);
    Z2i::Point px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second))); //pixel center
    double a = 0;
    std::vector<int> insideFaces;
    for(size_t it_bis=0; it_bis<gridH.face.size(); it_bis++) { //find all face h2 in the same pixel
      if(!fillAreaDone.at(it_bis) && isInsidePixel(px, gridH.getFaceCenter(it_bis))) {
        if(stateFaceHinit.at(it_bis)) // a face of Hinit
          a += gridH.getFaceArea(it_bis);
        phi.at(it_bis) = it;
        insideFaces.push_back(it_bis);
      }
      for(size_t it_bis=0; it_bis<insideFaces.size(); it_bis++) { //update area of all faces in the pixel
        fillAreaFace.at(insideFaces.at(it_bis)) = a;
        fillAreaDone.at(insideFaces.at(it_bis)) = true;
      }
    }
    Phi.push_back(insideFaces);
  }
  assert(Phi.size()==tgridF.face.size());
  sum = 0;
  for(size_t it=0; it<Phi.size(); it++)
    sum += Phi.at(it).size();
  assert(sum==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(fillAreaFace.at(it)>0.5) //majority interieur
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else {
        if(fillAreaFace.at(it)<EPSILON || fabs(fillAreaFace.at(it)-1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        else //majority exterieur
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    filename = file + "_fillAreaFaceMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 2-4
  //Sigma fct on F2
  std::vector<bool> sigmaValue(tgridF.face.size(), false);
  for(size_t it=0; it<tgridF.face.size(); it++) {
    if(Phi.at(it).size()!=0 && fillAreaFace.at(Phi.at(it).front())>=0.5) //contains Hinit cell and in majority
      sigmaValue.at(it) = true;
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(sigmaValue.at(it))
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else
        tgridF.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    filename = file + "_sigmaMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Zeta fct on F2
  std::vector<double> zetaValue(tgridF.face.size(), 0.0);
  double z = 0.0, s = -1.0;
  for(size_t it=0; it<tgridF.face.size(); it++) {
    s = sigmaValue.at(it) ? 1.0 : -1.0;
    z = (1 - s) / 2.0;
    if(Phi.at(it).size()!=0) //contains Hinit cell
       z += s*fillAreaFace.at(Phi.at(it).front());
    zetaValue.at(it) = z;
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(fabs(zetaValue.at(it)-1)<EPSILON) //pures pixels
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else { //non pures pixels
        if(zetaValue.at(it)>=0.5) //majority fill area
          tgridF.drawFace(aBoard, it, DGtal::Color::Blue);
        else
          assert(false);
      }
    }
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_zetaMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 5-8
  //Iota fct on h2
  std::vector<bool> iotaValue = stateFaceHcurrent;
  assert(iotaValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    filename = file + "_iotaMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Epsilon fct on h2
  std::vector<double> epsilonValue(gridH.face.size(), 0.0);
  double s1, s2;
  for(size_t it=0; it<gridH.face.size(); it++) {
    s1 = iotaValue.at(it) ? 1.0 : -1.0;
    s2 = sigmaValue.at(phi.at(it)) ? 1.0 : -1.0;
    epsilonValue.at(it) = -s1*s2*zetaValue.at(phi.at(it));
  }
  if(saveFile) {
    // Creating colormap.
    GradientColorMap<int> cmap_grad( 0, 10 );
    cmap_grad.addColor( Color( 255, 255, 0 ) ); //from Yellow
    cmap_grad.addColor( Color( 255, 0, 0 ) ); //to Red
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(epsilonValue.at(it)>0) { //most priority and positif
        Color c = cmap_grad( int(epsilonValue.at(it)*10) );
        gridH.drawFace(aBoard, it, c);//DGtal::Color::Red
      }
      else {//negative
        if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
        else //less priority and negative
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
      }
    }
    //complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_epsilonMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Tau fct on h2
  std::vector<bool> tauValue = tau(gridH, Hinit);
  assert(tauValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_tauMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 9: Compute B2 of Hcurrent (of same index as gridH)
  std::vector<bool> borderHcurrent = getBorderFace(gridH, Hinit);
  assert(borderHcurrent.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(borderHcurrent.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    filename = file + "_borderFaceMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 10-11: find border face of max priority
  size_t maxIndex = maxIndexEpsilon(epsilonValue, borderHcurrent);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, DGtal::Color::Green);
    gridH.drawFace(aBoard, maxIndex, DGtal::Color::Red);
    filename = file + "_maxEpsilonMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  //Line 12-26: Loop
  int nbIter=0;
  double a, h2 = maxIndex, f2 = phi.at(h2);
  while(epsilonValue.at(h2)>=0) {
    std::cout<<"nbIter="<<nbIter<<": h2="<<h2<<", f2="<<f2;
    //Line 13
    Complex Hcurrent;
    //Hcurrent.initGridLines(gridH.vertical_lines, gridH.horizontal_lines);
    Hcurrent.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
    std::vector<int> inFace = Phi.at(f2);
    if(tauValue.at(h2)) {//if simple cell
      //Line 14-15: Update iota of h2
      //std::cout<<"Before: iotaValue.at(h2)="<<iotaValue.at(h2)<<std::endl;
      iotaValue.at(h2) = !iotaValue.at(h2); //Inverse state face
      //std::cout<<"After: iotaValue.at(h2)="<<iotaValue.at(h2)<<std::endl;
      for(size_t it=0; it<iotaValue.size(); it++)
        if(iotaValue.at(it))
          Hcurrent.addFace(gridH.getFaceVertices(it), true, true);
      std::cout<<", Hcurrent="<<Hcurrent.face.size()<<", cH="<<gridH.face.size()<<std::endl;
      if(saveFile || saveExpe) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(iotaValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        Hcurrent.drawComplex(aBoard, Color::Gray);
        filename = file + "_HcurrentMajority_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
      
      //Line 16 => out of condition
      a = fillAreaFace.at(h2);
      //std::cout<<"Before: fillAreaFace.at(h2)="<<fillAreaFace.at(maxIndex)<<std::endl;
      if(iotaValue.at(h2))
        a += areaFace.at(h2);
      else
        a -= areaFace.at(h2);
      for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) //update fill area for all inside faces
        fillAreaFace.at(inFace.at(it_bis)) = a;
      //std::cout<<"After: fillAreaFace.at(h2)="<<fillAreaFace.at(maxIndex)<<std::endl;
      //Line 17-19
      std::vector<bool> tauHcurrent = tau(gridH, Hcurrent);
      for(size_t it=0; it<gridH.face.size(); it++)
        tauValue.at(it) = tauHcurrent.at(it);
      //Line 19 ????
    }
    else //Line 20-21
      sigmaValue.at(f2) = -sigmaValue.at(f2);
    //Line 16/22: Update zeta of f2
    //std::cout<<"Before: zetaValue.at(f2)="<<zetaValue.at(f2)<<std::endl;
    s = sigmaValue.at(f2) ? 1 : -1;
    zetaValue.at(f2) = (1.0-s)/2.0 + s*a;
    //std::cout<<"After: zetaValue.at(f2)="<<zetaValue.at(f2)<<std::endl;
    //Line 23-25
    //std::vector<int> inFace = Phi.at(f2);
    s2 = sigmaValue.at(f2) ? 1 : -1;
    a = zetaValue.at(f2);
    for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) {
      //std::cout<<"Before: epsilonValue.at("<<inFace.at(it_bis)<<")="<<epsilonValue.at(inFace.at(it_bis))<<std::endl;
      s1 = iotaValue.at(inFace.at(it_bis)) ? 1 : -1;
      epsilonValue.at(inFace.at(it_bis)) = -s1*s2*a;
      //std::cout<<"After: epsilonValue.at("<<inFace.at(it_bis)<<")="<<epsilonValue.at(inFace.at(it_bis))<<std::endl;
      //Line 25 out of loop
    }
    
    if(saveFile || saveExpe) {
      gridH.drawGrid(aBoard);
      for(size_t it=0; it<gridH.face.size(); it++) {
        if(epsilonValue.at(it)>0) { //most priority and positif
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        else {//negative
          if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
          else //less priority and negative
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
        }
      }
      Hcurrent.drawComplex(aBoard, Color::Gray);
      filename = file + "_epsilonMajority_s" + std::to_string(nbIter) + ".svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
    
    //Line 25
    std::vector<bool> borderH = getBorderFace(gridH, Hcurrent);
    for(size_t it_bis=0; it_bis<borderH.size(); it_bis++)
      borderHcurrent.at(it_bis) = borderH.at(it_bis);
    //Line 26
    h2 = maxIndexEpsilon(epsilonValue, borderHcurrent);
    f2 = phi.at(h2);
    nbIter++;
  }
  
  Complex Hopt;
  Hopt.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  for(size_t it=0; it<iotaValue.size(); it++)
    if(iotaValue.at(it))
      Hopt.addFace(gridH.getFaceVertices(it), true, true);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hopt.drawComplex(aBoard, Color::Gray);
    filename = file + "_HoptMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  //Retreive integer points
  RationalPoint p;
  Z2i::Point px;
  for(size_t it=0; it<iotaValue.size(); it++) {
    if(iotaValue.at(it)) {
      p = tgridF.getFaceCenter(phi.at(it));
      px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second)));
      if(!containElement(Y, px))
        Y.push_back(px);
    }
  }
  if(saveFile || saveExpe) {
    Complex cY(Y, adj);
    cY.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
    if(saveFile) {
      cY.drawGrid(aBoard);
      cY.drawPixel(aBoard);
      filename = file + "_objetYMajority.svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
  }
  return Y;
}

std::vector<Z2i::Point>
illustrationGauss2(std::vector<Z2i::Point> X, AffineTransform t, int adj = 8, bool saveFile = true) {
  std::vector<Z2i::Point> Y;
  std::pair<Z2i::Point,Z2i::Point> bb = getBoundingBox(X);
  string filename;
  
  //Grid F
  Complex gridF;
  Z2i::Point p1 = bb.first - Z2i::Point(2,2);
  Z2i::Point p2 = bb.second + Z2i::Point(2,2);
  gridF.initGrid(p1,p2);
  
  //GridG = t(GridF)
  Complex gridG = transform(gridF,t);
  Z2i::Point tp1 = gridG.domain.first;
  Z2i::Point tp2 = gridG.domain.second;
  Complex tgridF;
  tgridF.initGrid(tp1,tp2);
  
  //Grid H = tgridF and gridG
  vector<vector<vector<RationalPoint> > > h_faces = intersectionFace(gridG, tgridF);
  vector<vector<RationalPoint> > h_cells; //all faces
  for(int it=0; it<h_faces.size(); it++)
    for(int it_bis=0; it_bis<h_faces.at(it).size(); it_bis++)
      h_cells.push_back(h_faces.at(it).at(it_bis));
  //Get cells of H completly embedded in tgridF
  Complex tgridH(h_cells);
  vector<vector<RationalPoint> > h_pure_cells; //all pure faces
  vector<vector<size_t> > idInFace; //all pure faces
  std::vector<double> zetaH = zeta(tgridF, tgridH, idInFace);
  std::vector<bool> sigmaH = sigma(tgridF, tgridH);
  assert(idInFace.size()==tgridF.face.size());
  for(size_t it=0; it<zetaH.size(); it++)
  if(fabs(zetaH.at(it)-1)<EPSILON && sigmaH.at(it)) {//pure and inside pixels
    //h_pure_cells.push_back(tgridF.getFaceVertices(it));
    for(size_t it_bis=0; it_bis<idInFace.at(it).size(); it_bis++) {
      size_t idF = idInFace.at(it).at(it_bis);
      h_pure_cells.push_back(tgridH.getFaceVertices(idF));
    }
  }
  //Embed into gridH
  Complex gridH(h_pure_cells);//gridH(h_cells);
  gridH.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  gridH.addGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  
  Board2D aBoard;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/gridH.svg");
    aBoard.clear();
  }
  
  Complex cX(X, adj);
  cX.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  if(saveFile) {
    cX.drawGrid(aBoard);
    cX.drawPixel(aBoard);
    aBoard.saveSVG("../Illustration/objetX.svg");
    aBoard.clear();
  }
  
  //ComplexF
  Complex complexF;
  complexF.init(tp1,tp2);
  for(int it=0; it<X.size(); it++) {
    complexF.addCubicalFace(X.at(it));
    complexF.idConnectedComponent.push_back(cX.getFaceConnectedComponent(it));
  }
  
  if(saveFile) {
    complexF.drawGrid(aBoard);
    complexF.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexF.svg");
    aBoard.clear();
  }
  
  //ComplexG = t(ComplexF)
  Complex complexG = transform(complexF,t);
  if(saveFile) {
    complexG.drawGrid(aBoard);
    complexG.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexG.svg");
    aBoard.clear();
    
    complexF.drawComplex(aBoard, Color::Red);
    complexG.drawComplex(aBoard, Color::Blue);
    complexF.drawGrid(aBoard);
    complexG.drawGrid(aBoard);
    aBoard.saveSVG("../Illustration/complexFG.svg");
    aBoard.clear();
  }
  
  //Complex H :refine from gridF and gridG
  //Decomposition of R(X) --given by ComplexG-- by GridF => output: cells
  vector<vector<vector<RationalPoint> > > intersection_cells = intersectionFace(complexG, tgridF);
  vector<vector<RationalPoint> > cells; //all faces
  for(int it=0; it<intersection_cells.size(); it++)
    for(int it_bis=0; it_bis<intersection_cells.at(it).size(); it_bis++)
      cells.push_back(intersection_cells.at(it).at(it_bis));
  
  //Build the Complex H (decomposition of ComplexG(=t(ComplexF)) by GridF) from cells
  Complex complexH(cells);
  complexH.initGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  for(size_t it=0; it<complexH.face.size(); it++) {
    vector<RationalPoint> f = complexH.getFaceVertices(it);
    int idCC = complexG.getIndexConnectedComponent(f);
    complexH.idConnectedComponent.push_back(idCC);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    complexH.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/complexH.svg");
    aBoard.clear();
  }
  
  /*
  getNeigbourFace(complexH);
  std::vector<bool> borderH = getBorderFace(gridH, complexH);
  std::vector<bool> stateH = setFaceState(complexH, gridH);
  getNeigbourFace(gridH, borderH, stateH);
  */
  
  Complex complexGauss = initizeGauss(complexH, gridH);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    complexGauss.drawComplex(aBoard, Color::Gray);
    complexH.drawComplex(aBoard, Color::Red);
    aBoard.saveSVG("../Illustration/complexGaussH.svg");
    aBoard.clear();
  }
  
  //Line 1
  Complex Hinit = complexH;
  std::vector<bool> stateFaceHinit = setFaceState(Hinit, gridH);
  int sum = 0;
  for(size_t it = 0; it<stateFaceHinit.size(); it++)
    if(stateFaceHinit.at(it))
      sum++;
  assert(sum==Hinit.face.size());
  assert(stateFaceHinit.size()==gridH.face.size());
  //std::vector<bool> stateFaceHcurrent = stateFaceHinit;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(stateFaceHinit.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/stateFaceInitGauss.svg");
    aBoard.clear();
  }
  
  Complex Htarget = complexGauss;
  std::vector<bool> stateFaceHtarget = setFaceState(Htarget, gridH);
  sum = 0;
  for(size_t it = 0; it<stateFaceHtarget.size(); it++)
    if(stateFaceHtarget.at(it))
      sum++;
  assert(sum==Htarget.face.size());
  assert(stateFaceHtarget.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(stateFaceHtarget.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/stateFaceTargetGauss.svg");
    aBoard.clear();
  }
  
  //Optimization of data structure
  //Compute h2 area of Hinit (wrt Hgrid)
  std::vector<double> areaFace(gridH.face.size(),0.0);
  for(size_t it=0; it<gridH.face.size(); it++) {
    if(stateFaceHinit.at(it))
      areaFace.at(it) = gridH.getFaceArea(it);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(areaFace.at(it)>0) //face of Hinit
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/areaFaceGauss.svg");
    aBoard.clear();
  }
  //Compute fill area f2 of Hinit (wrt Hgrid)
  std::vector<double> fillAreaHinitFace(gridH.face.size(),0.0);
  std::vector<bool> fillAreaDone(gridH.face.size(),false);
  std::vector<std::vector<int> > Phi; //Phi (f2) = {h2}
  std::vector<int> phi (gridH.face.size(), -1); //phi (h2) = f2
  std::pair<Z2i::Point,Z2i::Point> d = tgridF.domain;
  for(size_t it=0; it<tgridF.face.size(); it++) { //for each pixel of gridF
    RationalPoint p = tgridF.getFaceCenter(it);
    Z2i::Point px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second))); //pixel center
    double a = 0;
    std::vector<int> insideFaces;
    for(size_t it_bis=0; it_bis<gridH.face.size(); it_bis++) { //find all face h2 in the same pixel
      if(!fillAreaDone.at(it_bis) && isInsidePixel(px, gridH.getFaceCenter(it_bis))) {
        if(stateFaceHinit.at(it_bis)) // a face of Hinit
          a += gridH.getFaceArea(it_bis);
        phi.at(it_bis) = it;
        insideFaces.push_back(it_bis);
      }
      for(size_t it_bis=0; it_bis<insideFaces.size(); it_bis++) { //update area of all faces in the pixel
        fillAreaHinitFace.at(insideFaces.at(it_bis)) = a;
        fillAreaDone.at(insideFaces.at(it_bis)) = true;
      }
    }
    Phi.push_back(insideFaces);
  }
  assert(Phi.size()==tgridF.face.size());
  sum = 0;
  for(size_t it=0; it<Phi.size(); it++)
    sum += Phi.at(it).size();
  assert(sum==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(fillAreaHinitFace.at(it)>0.5) //majority interieur
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else {
        if(fillAreaHinitFace.at(it)<EPSILON || fabs(fillAreaHinitFace.at(it)-1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        else //majority exterieur
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/fillAreaFaceHinitGauss.svg");
    aBoard.clear();
  }
  
  //Compute fill area f2 of Hinit (wrt Hgrid)
  std::vector<double> fillAreaHtargetFace(gridH.face.size(),0.0);
  std::vector<bool> fillAreaHtargetDone(gridH.face.size(),false);
  for(size_t it=0; it<tgridF.face.size(); it++) { //for each pixel of gridF
    RationalPoint p = tgridF.getFaceCenter(it);
    Z2i::Point px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second))); //pixel center
    double a = 0;
    std::vector<int> insideFaces;
    for(size_t it_bis=0; it_bis<gridH.face.size(); it_bis++) { //find all face h2 in the same pixel
      if(!fillAreaHtargetDone.at(it_bis) && isInsidePixel(px, gridH.getFaceCenter(it_bis))) {
        if(stateFaceHtarget.at(it_bis)) // a face of Hinit
          a += gridH.getFaceArea(it_bis);
        insideFaces.push_back(it_bis);
      }
      for(size_t it_bis=0; it_bis<insideFaces.size(); it_bis++) { //update area of all faces in the pixel
        fillAreaHtargetFace.at(insideFaces.at(it_bis)) = a;
        fillAreaHtargetDone.at(insideFaces.at(it_bis)) = true;
      }
    }
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(fillAreaHtargetFace.at(it)>0.5) //majority interieur
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else {
        if(fillAreaHtargetFace.at(it)<EPSILON || fabs(fillAreaHtargetFace.at(it)-1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        else //majority exterieur
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
    Htarget.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/fillAreaFaceHtargetGauss.svg");
    aBoard.clear();
  }
  
  //Line 2-4
  //Sigma fct on F2
  std::vector<bool> sigmaValue(tgridF.face.size(), false);
  for(size_t it=0; it<tgridF.face.size(); it++) {
    if(Phi.at(it).size()!=0 && fillAreaHtargetFace.at(Phi.at(it).front())>=0.5) //contains Hinit cell and in majority
      sigmaValue.at(it) = true;
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(sigmaValue.at(it))
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else
        tgridF.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Htarget.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/sigmaGauss.svg");
    aBoard.clear();
  }
  
  //Zeta fct on F2
  std::vector<double> zetaValue(tgridF.face.size(), 0.0);
  double z = 0.0, s = -1.0;
  for(size_t it=0; it<tgridF.face.size(); it++) {
    s = sigmaValue.at(it) ? 1.0 : -1.0;
    z = (1 - s) / 2.0;
    if(Phi.at(it).size()!=0) //contains Hinit cell
       z += s*fillAreaHinitFace.at(Phi.at(it).front());
    zetaValue.at(it) = z;
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(fabs(zetaValue.at(it)-1)<EPSILON || fabs(zetaValue.at(it))<EPSILON) //pure and (in)correct pixels
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else { //non pures pixels
        if(zetaValue.at(it)>=0.5) //majority fill area and correct
          tgridF.drawFace(aBoard, it, DGtal::Color::Blue);
        else //majority fill area and incorrect
          tgridF.drawFace(aBoard, it, DGtal::Color::Red);
      }
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    Htarget.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/zetaGauss.svg");
    aBoard.clear();
  }
  
  //Line 5-8
  //Iota fct on h2
  std::vector<bool> iotaValue = stateFaceHinit;
  assert(iotaValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    aBoard.saveSVG("../Illustration/iotaGauss.svg");
    aBoard.clear();
  }
  
  //Epsilon fct on h2
  std::vector<double> epsilonValue(gridH.face.size(), 0.0);
  double s1, s2;
  for(size_t it=0; it<gridH.face.size(); it++) {
    s1 = iotaValue.at(it) ? 1.0 : -1.0;
    s2 = sigmaValue.at(phi.at(it)) ? 1.0 : -1.0;
    epsilonValue.at(it) = -s1*s2*zetaValue.at(phi.at(it));
  }
  if(saveFile) {
    // Creating colormap.
    GradientColorMap<int> cmap_grad( 0, 10 );
    cmap_grad.addColor( Color( 255, 255, 0 ) ); //from Yellow
    cmap_grad.addColor( Color( 255, 0, 0 ) ); //to Red
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(epsilonValue.at(it)>0) { //most priority and positif
        Color c = cmap_grad( int(epsilonValue.at(it)*10) );
        gridH.drawFace(aBoard, it, c);//DGtal::Color::Red
      }
      else {//negative
        if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
        else //less priority and negative
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
      }
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    Htarget.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/epsilonGauss.svg");
    aBoard.clear();
  }
  
  //Tau fct on h2
  std::vector<bool> tauValue = tau(gridH, Hinit);
  assert(tauValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/tauGauss.svg");
    aBoard.clear();
  }
  
  //Line 9: Compute B2 of Hcurrent (of same index as gridH)
  std::vector<bool> borderH = getBorderFace(gridH, Hinit);
  std::vector<bool> borderHcurrent = getNeigbourFace(gridH, borderH, iotaValue);
  assert(borderHcurrent.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(borderHcurrent.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    aBoard.saveSVG("../Illustration/borderFaceGauss.svg");
    aBoard.clear();
  }
  //Update tau value of border faces
  for(size_t it=0; it<borderHcurrent.size(); it++) {
    if(borderHcurrent.at(it) && !iotaValue.at(it)) { //border face but not belong to Hcurrent
      Complex tmp = Hinit;
      tmp.addFace(gridH.getFaceVertices(it), true, true);
      //Verify new face simplicity
      tauValue.at(it) = tmp.isSimpleFace(tmp.face.size()-1);
    }
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/tauGauss2.svg");
    aBoard.clear();
  }
  
  //Line 10-11: find border face of max priority
  size_t maxIndex = maxIndexEpsilon(epsilonValue, borderHcurrent);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, DGtal::Color::Green);
    gridH.drawFace(aBoard, maxIndex, DGtal::Color::Red);
    aBoard.saveSVG("../Illustration/maxEpsilonGauss.svg");
    aBoard.clear();
  }
  //Line 12-26: Loop
  int nbIter=0;
  double a, h2 = maxIndex, f2 = phi.at(h2);
  while(epsilonValue.at(h2)>=0) {
    std::cout<<"nbIter="<<nbIter<<": h2="<<h2<<", f2="<<f2;
    //Line 13
    Complex Hcurrent;
    //Hcurrent.initGridLines(gridH.vertical_lines, gridH.horizontal_lines);
    Hcurrent.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
    std::vector<int> inFace = Phi.at(f2);
    a = fillAreaHinitFace.at(h2);
    if(tauValue.at(h2)) {//if simple cell
      //Line 14-15: Update iota of h2
      //std::cout<<"Before: iotaValue.at(h2)="<<iotaValue.at(h2)<<std::endl;
      iotaValue.at(h2) = !iotaValue.at(h2); //Inverse state face
      //std::cout<<"After: iotaValue.at(h2)="<<iotaValue.at(h2)<<std::endl;
      for(size_t it=0; it<iotaValue.size(); it++)
        if(iotaValue.at(it))
          Hcurrent.addFace(gridH.getFaceVertices(it), true, true);
      std::cout<<", Hcurrent="<<Hcurrent.face.size()<<", cH="<<gridH.face.size()<<std::endl;
      if(saveFile) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(iotaValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        Hcurrent.drawComplex(aBoard, Color::Gray);
        filename = "../Illustration/HcurrentGauss_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
      
      //Line 16 => out of condition
      //std::cout<<"Before: fillAreaFace.at(h2)="<<fillAreaFace.at(maxIndex)<<std::endl;
      if(iotaValue.at(h2))
        a += areaFace.at(h2);
      else
        a -= areaFace.at(h2);
      for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) //update fill area for all inside faces
        fillAreaHinitFace.at(inFace.at(it_bis)) = a;
      //std::cout<<"After: fillAreaFace.at(h2)="<<fillAreaFace.at(maxIndex)<<std::endl;
      //Line 17-19
      std::vector<bool> tauHcurrent = tau(gridH, Hcurrent);
      for(size_t it=0; it<gridH.face.size(); it++)
        tauValue.at(it) = tauHcurrent.at(it);
      /*
      if(saveFile) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(tauValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        complexH.drawComplex(aBoard, Color::Gray);
        filename = "../Illustration/tau_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
      */
      //Line 19 ????
    }
    else { //Line 20-21
      sigmaValue.at(f2) = !sigmaValue.at(f2);
      for(size_t it=0; it<iotaValue.size(); it++)
        if(iotaValue.at(it))
          Hcurrent.addFace(gridH.getFaceVertices(it), true, true);
      std::cout<<", Hcurrent="<<Hcurrent.face.size()<<", cH="<<gridH.face.size()<<std::endl;
      if(saveFile) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(iotaValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        Hcurrent.drawComplex(aBoard, Color::Gray);
        filename = "../Illustration/HcurrentGauss_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
    }
    //Line 16/22: Update zeta of f2
    //std::cout<<"Before: zetaValue.at(f2)="<<zetaValue.at(f2)<<std::endl;
    s = sigmaValue.at(f2) ? 1 : -1;
    zetaValue.at(f2) = (1.0-s)/2.0 + s*a;
    //std::cout<<"After: zetaValue.at(f2)="<<zetaValue.at(f2)<<std::endl;
    //Line 23-25
    //std::vector<int> inFace = Phi.at(f2);
    s2 = sigmaValue.at(f2) ? 1 : -1;
    a = zetaValue.at(f2);
    for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) {
      //std::cout<<"Before: epsilonValue.at("<<inFace.at(it_bis)<<")="<<epsilonValue.at(inFace.at(it_bis))<<std::endl;
      s1 = iotaValue.at(inFace.at(it_bis)) ? 1 : -1;
      epsilonValue.at(inFace.at(it_bis)) = -s1*s2*a;
      //std::cout<<"After: epsilonValue.at("<<inFace.at(it_bis)<<")="<<epsilonValue.at(inFace.at(it_bis))<<std::endl;
      //Line 25 out of loop
    }
    if(saveFile) {
      gridH.drawGrid(aBoard);
      for(size_t it=0; it<gridH.face.size(); it++) {
        if(epsilonValue.at(it)>0) { //most priority and positif
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        else {//negative
          if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
          else //less priority and negative
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
        }
      }
      Hcurrent.drawComplex(aBoard, Color::Gray);
      filename = "../Illustration/epsilonGauss_s" + std::to_string(nbIter) + ".svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
    
    //Line 25: Update border faces
    std::vector<bool> borderH = getBorderFace(gridH, Hcurrent);
    std::vector<bool> borderHTmp = getNeigbourFace(gridH, borderH, iotaValue);
    for(size_t it_bis=0; it_bis<borderHTmp.size(); it_bis++)
      borderHcurrent.at(it_bis) = borderHTmp.at(it_bis);
    
    for(size_t it=0; it<borderHcurrent.size(); it++) {
      if(borderHcurrent.at(it) && !iotaValue.at(it)) { //border face but not belong to Hcurrent
        Complex tmp = Hcurrent;
        tmp.addFace(gridH.getFaceVertices(it), true, true);
        //Verify new face simplicity
        tauValue.at(it) = tmp.isSimpleFace(tmp.face.size()-1);
      }
    }
    if(saveFile) {
      gridH.drawGrid(aBoard);
      for(size_t it=0; it<gridH.face.size(); it++) {
        if(tauValue.at(it))
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
        else
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
      }
      Hcurrent.drawComplex(aBoard, Color::Gray);
      filename = "../Illustration/tauGauss_s" + std::to_string(nbIter) + ".svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
    
    //Line 26
    h2 = maxIndexEpsilon(epsilonValue, borderHcurrent);
    f2 = phi.at(h2);
    nbIter++;
  }
  
  Complex Hopt;
  Hopt.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  for(size_t it=0; it<iotaValue.size(); it++)
    if(iotaValue.at(it))
      Hopt.addFace(gridH.getFaceVertices(it), true, true);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hopt.drawComplex(aBoard, Color::Gray);
    aBoard.saveSVG("../Illustration/HoptGauss.svg");
    aBoard.clear();
  }
  //Retreive integer points
  RationalPoint p;
  Z2i::Point px;
  for(size_t it=0; it<iotaValue.size(); it++) {
    if(iotaValue.at(it)) {
      p = tgridF.getFaceCenter(phi.at(it));
      px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second)));
      if(!containElement(Y, px))
        Y.push_back(px);
    }
  }
  if(saveFile) {
    Complex cY(Y, adj);
    cY.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
    if(saveFile) {
      cY.drawGrid(aBoard);
      cY.drawPixel(aBoard);
      aBoard.saveSVG("../Illustration/objetYGauss.svg");
      aBoard.clear();
    }
  }
  return Y;
}

std::vector<Z2i::Point>
illustrationGauss3(std::vector<Z2i::Point> X, AffineTransform t, string file = "test", int border = 2, int adj = 8, bool saveExpe = true, bool saveFile = false, bool useSaveFile = false) {
  std::vector<Z2i::Point> Y;
  std::pair<Z2i::Point,Z2i::Point> bb = getBoundingBox(X);
  
  //Grid F
  Complex gridF;
  Z2i::Point p1 = bb.first - Z2i::Point(border,border);
  Z2i::Point p2 = bb.second + Z2i::Point(border,border);
  gridF.initGrid(p1,p2);
  
  //GridG = t(GridF)
  Complex gridG = transform(gridF,t);
  Z2i::Point tp1 = gridG.domain.first;
  Z2i::Point tp2 = gridG.domain.second;
  Complex tgridF;
  tgridF.initGrid(tp1,tp2);
  
  //Grid H = tgridF and gridG
  vector<vector<vector<RationalPoint> > > h_faces = intersectionFace(gridG, tgridF);
  vector<vector<RationalPoint> > h_cells; //all faces
  for(int it=0; it<h_faces.size(); it++)
    for(int it_bis=0; it_bis<h_faces.at(it).size(); it_bis++)
      h_cells.push_back(h_faces.at(it).at(it_bis));
  //Get cells of H completly embedded in tgridF
  Complex tgridH(h_cells);
  vector<vector<RationalPoint> > h_pure_cells; //all pure faces
  vector<vector<size_t> > idInFace; //all pure faces
  std::vector<double> zetaH = zeta(tgridF, tgridH, idInFace);
  std::vector<bool> sigmaH = sigma(tgridF, tgridH);
  assert(idInFace.size()==tgridF.face.size());
  for(size_t it=0; it<zetaH.size(); it++)
  if(fabs(zetaH.at(it)-1)<EPSILON && sigmaH.at(it)) {//pure and inside pixels
    //h_pure_cells.push_back(tgridF.getFaceVertices(it));
    for(size_t it_bis=0; it_bis<idInFace.at(it).size(); it_bis++) {
      size_t idF = idInFace.at(it).at(it_bis);
      h_pure_cells.push_back(tgridH.getFaceVertices(idF));
    }
  }
  //Embed into gridH
  Complex gridH(h_pure_cells);//gridH(h_cells);
  gridH.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  gridH.addGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  
  Board2D aBoard;
  string filename;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, Color::Gray);
    filename = file + "_gridH.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  Complex cX(X, adj);
  cX.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  if(saveFile || saveExpe) {
    cX.drawGrid(aBoard);
    cX.drawPixel(aBoard);
    filename = file + "_objetX.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //ComplexF
  Complex complexF;
  complexF.init(tp1,tp2);
  for(int it=0; it<X.size(); it++) {
    complexF.addCubicalFace(X.at(it));
    complexF.idConnectedComponent.push_back(cX.getFaceConnectedComponent(it));
  }
  
  if(saveFile) {
    complexF.drawGrid(aBoard);
    complexF.drawComplex(aBoard, Color::Gray);
    filename = file + "_complexF.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //ComplexG = t(ComplexF)
  Complex complexG = transform(complexF,t);
  if(saveFile) {
    complexG.drawGrid(aBoard);
    complexG.drawComplex(aBoard, Color::Gray);
    filename = file + "_complexG.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
    
    complexF.drawComplex(aBoard, Color::Red);
    complexG.drawComplex(aBoard, Color::Blue);
    complexF.drawGrid(aBoard);
    complexG.drawGrid(aBoard);
    filename = file + "_complexFG.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Complex H :refine from gridF and gridG
  //Decomposition of R(X) --given by ComplexG-- by GridF => output: cells
  vector<vector<vector<RationalPoint> > > intersection_cells = intersectionFace(complexG, tgridF);
  vector<vector<RationalPoint> > cells; //all faces
  for(int it=0; it<intersection_cells.size(); it++)
    for(int it_bis=0; it_bis<intersection_cells.at(it).size(); it_bis++)
      cells.push_back(intersection_cells.at(it).at(it_bis));
  
  //Build the Complex H (decomposition of ComplexG(=t(ComplexF)) by GridF) from cells
  Complex complexH(cells);
  complexH.initGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  for(size_t it=0; it<complexH.face.size(); it++) {
    vector<RationalPoint> f = complexH.getFaceVertices(it);
    int idCC = complexG.getIndexConnectedComponent(f);
    complexH.idConnectedComponent.push_back(idCC);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_complexH.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  Complex complexGauss = initizeGauss(complexH, gridH);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    complexGauss.drawComplex(aBoard, Color::Gray);
    complexH.drawComplex(aBoard, Color::Red);
    filename = file + "_HinitGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 1
  Complex Hinit = complexH;
  std::vector<bool> stateFaceHinit = setFaceState(Hinit, gridH);
  int sum = 0;
  for(size_t it = 0; it<stateFaceHinit.size(); it++)
    if(stateFaceHinit.at(it))
      sum++;
  assert(sum==Hinit.face.size());
  assert(stateFaceHinit.size()==gridH.face.size());
  std::vector<bool> stateFaceHcurrent = stateFaceHinit;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(stateFaceHinit.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_stageFaceGauss.svg";
    aBoard.clear();
  }
  
  //Optimization of data structure
  //Compute h2 area of Hinit (wrt Hgrid)
  std::vector<double> areaFace(gridH.face.size(),0.0);
  for(size_t it=0; it<gridH.face.size(); it++) {
    if(stateFaceHinit.at(it))
      areaFace.at(it) = gridH.getFaceArea(it);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(areaFace.at(it)>0) //face of Hinit
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_areaFaceGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  //Compute fill area f2 of Hinit (wrt Hgrid)
  std::vector<double> fillAreaFace(gridH.face.size(),0.0);
  std::vector<bool> fillAreaDone(gridH.face.size(),false);
  std::vector<std::vector<int> > Phi; //Phi (f2) = {h2}
  std::vector<int> phi (gridH.face.size(), -1); //phi (h2) = f2
  std::pair<Z2i::Point,Z2i::Point> d = tgridF.domain;
  for(size_t it=0; it<tgridF.face.size(); it++) { //for each pixel of gridF
    RationalPoint p = tgridF.getFaceCenter(it);
    Z2i::Point px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second))); //pixel center
    double a = 0;
    std::vector<int> insideFaces;
    for(size_t it_bis=0; it_bis<gridH.face.size(); it_bis++) { //find all face h2 in the same pixel
      if(!fillAreaDone.at(it_bis) && isInsidePixel(px, gridH.getFaceCenter(it_bis))) {
        if(stateFaceHinit.at(it_bis)) // a face of Hinit
          a += gridH.getFaceArea(it_bis);
        phi.at(it_bis) = it;
        insideFaces.push_back(it_bis);
      }
      for(size_t it_bis=0; it_bis<insideFaces.size(); it_bis++) { //update area of all faces in the pixel
        fillAreaFace.at(insideFaces.at(it_bis)) = a;
        fillAreaDone.at(insideFaces.at(it_bis)) = true;
      }
    }
    Phi.push_back(insideFaces);
  }
  assert(Phi.size()==tgridF.face.size());
  sum = 0;
  for(size_t it=0; it<Phi.size(); it++)
    sum += Phi.at(it).size();
  assert(sum==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(fillAreaFace.at(it)>0.5) //majority interieur
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else {
        if(fillAreaFace.at(it)<EPSILON || fabs(fillAreaFace.at(it)-1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        else //majority exterieur
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    filename = file + "_fillAreaFaceGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 2-4
  //Sigma fct on F2
  std::vector<bool> sigmaValue(tgridF.face.size(), false);
  for(size_t it=0; it<tgridF.face.size(); it++) {
    if(Phi.at(it).size()!=0 && fillAreaFace.at(Phi.at(it).front())>=0.5) //contains Hinit cell and in majority
      sigmaValue.at(it) = true;
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(sigmaValue.at(it))
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else
        tgridF.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hinit.drawComplex(aBoard, Color::Gray);
    filename = file + "_sigmaGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Zeta fct on F2
  std::vector<double> zetaValue(tgridF.face.size(), 0.0);
  double z = 0.0, s = -1.0;
  for(size_t it=0; it<tgridF.face.size(); it++) {
    s = sigmaValue.at(it) ? 1.0 : -1.0;
    z = (1 - s) / 2.0;
    if(Phi.at(it).size()!=0) //contains Hinit cell
       z += s*fillAreaFace.at(Phi.at(it).front());
    zetaValue.at(it) = z;
  }
  if(saveFile) {
    tgridF.drawGrid(aBoard);
    for(size_t it=0; it<tgridF.face.size(); it++) {
      if(fabs(zetaValue.at(it)-1)<EPSILON) //pures pixels
        tgridF.drawFace(aBoard, it, DGtal::Color::Green);
      else { //non pures pixels
        if(zetaValue.at(it)>=0.5) //majority fill area
          tgridF.drawFace(aBoard, it, DGtal::Color::Blue);
        else
          assert(false);
      }
    }
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_zetaGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 5-8
  //Iota fct on h2
  std::vector<bool> iotaValue = stateFaceHcurrent;
  assert(iotaValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    filename = file + "_iotaGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Epsilon fct on h2
  std::vector<double> epsilonValue(gridH.face.size(), 0.0);
  double s1, s2;
  for(size_t it=0; it<gridH.face.size(); it++) {
    s1 = iotaValue.at(it) ? 1.0 : -1.0;
    s2 = sigmaValue.at(phi.at(it)) ? 1.0 : -1.0;
    epsilonValue.at(it) = -s1*s2*zetaValue.at(phi.at(it));
  }
  if(saveFile) {
    // Creating colormap.
    GradientColorMap<int> cmap_grad( 0, 10 );
    cmap_grad.addColor( Color( 255, 255, 0 ) ); //from Yellow
    cmap_grad.addColor( Color( 255, 0, 0 ) ); //to Red
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(epsilonValue.at(it)>0) { //most priority and positif
        Color c = cmap_grad( int(epsilonValue.at(it)*10) );
        gridH.drawFace(aBoard, it, c);//DGtal::Color::Red
      }
      else {//negative
        if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
        else //less priority and negative
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
      }
    }
    //complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_epsilonGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Tau fct on h2
  std::vector<bool> tauValue = tau(gridH, Hinit);
  assert(tauValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    complexH.drawComplex(aBoard, Color::Gray);
    filename = file + "_tauGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 9: Compute B2 of Hcurrent (of same index as gridH)
  std::vector<bool> borderHcurrent = getBorderFace(gridH, Hinit);
  assert(borderHcurrent.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(borderHcurrent.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    filename = file + "_borderFaceGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 10-11: find border face of max priority
  size_t maxIndex = maxIndexEpsilon(epsilonValue, borderHcurrent);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, DGtal::Color::Green);
    gridH.drawFace(aBoard, maxIndex, DGtal::Color::Red);
    filename = file + "_maxEpsilonGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  //Line 12-26: Loop
  int nbIter=0;
  double a, h2 = maxIndex, f2 = phi.at(h2);
  while(epsilonValue.at(h2)>=0) {
    std::cout<<"nbIter="<<nbIter<<": h2="<<h2<<", f2="<<f2;
    //Line 13
    Complex Hcurrent;
    //Hcurrent.initGridLines(gridH.vertical_lines, gridH.horizontal_lines);
    Hcurrent.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
    std::vector<int> inFace = Phi.at(f2);
    if(tauValue.at(h2)) {//if simple cell
      //Line 14-15: Update iota of h2
      //std::cout<<"Before: iotaValue.at(h2)="<<iotaValue.at(h2)<<std::endl;
      iotaValue.at(h2) = !iotaValue.at(h2); //Inverse state face
      //std::cout<<"After: iotaValue.at(h2)="<<iotaValue.at(h2)<<std::endl;
      for(size_t it=0; it<iotaValue.size(); it++)
        if(iotaValue.at(it))
          Hcurrent.addFace(gridH.getFaceVertices(it), true, true);
      std::cout<<", Hcurrent="<<Hcurrent.face.size()<<", cH="<<gridH.face.size()<<std::endl;
      if(saveFile || saveExpe) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(iotaValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        Hcurrent.drawComplex(aBoard, Color::Gray);
        filename = file + "_HcurrentGauss_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
      
      //Line 16 => out of condition
      a = fillAreaFace.at(h2);
      //std::cout<<"Before: fillAreaFace.at(h2)="<<fillAreaFace.at(maxIndex)<<std::endl;
      if(iotaValue.at(h2))
        a += areaFace.at(h2);
      else
        a -= areaFace.at(h2);
      for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) //update fill area for all inside faces
        fillAreaFace.at(inFace.at(it_bis)) = a;
      //std::cout<<"After: fillAreaFace.at(h2)="<<fillAreaFace.at(maxIndex)<<std::endl;
      //Line 17-19
      std::vector<bool> tauHcurrent = tau(gridH, Hcurrent);
      for(size_t it=0; it<gridH.face.size(); it++)
        tauValue.at(it) = tauHcurrent.at(it);
      //Line 19 ????
    }
    else //Line 20-21
      sigmaValue.at(f2) = -sigmaValue.at(f2);
    //Line 16/22: Update zeta of f2
    //std::cout<<"Before: zetaValue.at(f2)="<<zetaValue.at(f2)<<std::endl;
    s = sigmaValue.at(f2) ? 1 : -1;
    zetaValue.at(f2) = (1.0-s)/2.0 + s*a;
    //std::cout<<"After: zetaValue.at(f2)="<<zetaValue.at(f2)<<std::endl;
    //Line 23-25
    //std::vector<int> inFace = Phi.at(f2);
    s2 = sigmaValue.at(f2) ? 1 : -1;
    a = zetaValue.at(f2);
    for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) {
      //std::cout<<"Before: epsilonValue.at("<<inFace.at(it_bis)<<")="<<epsilonValue.at(inFace.at(it_bis))<<std::endl;
      s1 = iotaValue.at(inFace.at(it_bis)) ? 1 : -1;
      epsilonValue.at(inFace.at(it_bis)) = -s1*s2*a;
      //std::cout<<"After: epsilonValue.at("<<inFace.at(it_bis)<<")="<<epsilonValue.at(inFace.at(it_bis))<<std::endl;
      //Line 25 out of loop
    }
    
    if(saveFile || saveExpe) {
      gridH.drawGrid(aBoard);
      for(size_t it=0; it<gridH.face.size(); it++) {
        if(epsilonValue.at(it)>0) { //most priority and positif
          gridH.drawFace(aBoard, it, DGtal::Color::Red);
        }
        else {//negative
          if(fabs(epsilonValue.at(it)+1)<EPSILON) //pure pixels
            gridH.drawFace(aBoard, it, DGtal::Color::Blue);
          else //less priority and negative
            gridH.drawFace(aBoard, it, DGtal::Color::Green);
        }
      }
      Hcurrent.drawComplex(aBoard, Color::Gray);
      filename = file + "_epsilonGauss_s" + std::to_string(nbIter) + ".svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
    
    //Line 25
    std::vector<bool> borderH = getBorderFace(gridH, Hcurrent);
    for(size_t it_bis=0; it_bis<borderH.size(); it_bis++)
      borderHcurrent.at(it_bis) = borderH.at(it_bis);
    //Line 26
    h2 = maxIndexEpsilon(epsilonValue, borderHcurrent);
    f2 = phi.at(h2);
    nbIter++;
  }
  
  Complex Hopt;
  Hopt.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  for(size_t it=0; it<iotaValue.size(); it++)
    if(iotaValue.at(it))
      Hopt.addFace(gridH.getFaceVertices(it), true, true);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(iotaValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    Hopt.drawComplex(aBoard, Color::Gray);
    filename = file + "_HoptGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  //Retreive integer points
  RationalPoint p;
  Z2i::Point px;
  for(size_t it=0; it<iotaValue.size(); it++) {
    if(iotaValue.at(it)) {
      p = tgridF.getFaceCenter(phi.at(it));
      px = Z2i::Point(int(getRealValue(p.first)), int(getRealValue(p.second)));
      if(!containElement(Y, px))
        Y.push_back(px);
    }
  }
  if(saveFile || saveExpe) {
    Complex cY(Y, adj);
    cY.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
    if(saveFile) {
      cY.drawGrid(aBoard);
      cY.drawPixel(aBoard);
      filename = file + "_objetYGauss.svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
  }
  return Y;
}

int main()
{
  //testFunctions();
  //return 1;
  
  std::vector<Z2i::Point> X1 = readImage("test.pgm");
  assert(X1.size()!=0);
  AffineTransform t1(1.1,0.4,0.3,1.5,1,2,0,1);
  illustrationMajority2(X1,t1);
  AffineTransform t2(1.1,0.4,0.3,1.2,1,2,0,1);
  illustrationGauss2(X1,t1);
  return 1;
  
  bool saveFile = true;
  bool saveExpe = true;
  bool useSaveFile = true;
  int adj = 8;
  int border = 2; //5 ellipse10-6
  
  string intputDir = "../Samples/";
  string outputDir = "../Results/";//ball_r5//ellipse_10-6
  string filename = "jmiv";//ball_r5 dgmm3 rings5 smile4 ellipse10-6 ball_r5 jmiv
  string extenstion = ".pgm";
  string inputFile = intputDir + filename + extenstion;
  string outputFile = outputDir + filename;
  
  std::vector<Z2i::Point> X = readImage(inputFile);
  assert(X.size()!=0);
  std::cout<<"X.size="<<X.size()<<std::endl;
  
  //AffineTransform t(1.2,0.1,0.7,1.5,1,3,2,5); //dgmm3
  //AffineTransform t(1.2,-0.5,-0.1,1.1,1,2,0,1); //dgmm3 2
  //AffineTransform t(1.5,-0.2,-0.5,1.2,1,2,0,1); //rings5
  //AffineTransform t(1.2,-1.5,0.2,1.5,1,3,2,5); //smile4
  //AffineTransform t(1.1,0.8,0.3,1.7,1,3,2,5); // ellipse10-6
  //AffineTransform t(1.5,0.2,0.5,1.2,1,5,1,5); // ball_r5
  //AffineTransform t(0.9,0.3,0.2,1.1,1,5,1,5); // ellipse_10-6 2
  //AffineTransform t(1.2,1.5,0.2,1.5,1,3,2,5); //jmiv
  AffineTransform t(1.1,0.4,0.3,1.2,1,2,0,1); // affine
  
  illustrationMajority3(X, t, outputFile, border, adj, saveExpe, saveFile);
  //illustrationGauss3(X, t, outputFile, border, adj, saveExpe, saveFile);
  
  return 1;
}
