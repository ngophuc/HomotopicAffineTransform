#include "Functions.h"
#include "Complex.h"
#include "AffineTransform.h"

#include "CLI11.hpp"

/************************************/
/* Homotopic transform */
/************************************/
template <typename T>
std::vector<size_t> sort_indexes_decrease(const std::vector<T> &v) {
  
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
  return idx;
}

//Read image and return pixels in the image
std::vector<Z2i::Point>
readImage(std::string filename) {
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

//Compute border faces of gridH is simple wrt complex H
std::vector<bool>
getBorderFace(const Complex& gridH, const std::vector<bool>& cH) {
  std::vector<bool> isBorder(gridH.face.size(), false);
  for(size_t it=0; it<gridH.edge.size(); it++) {
    std::pair<int, int> f = gridH.star_edge.at(it);
    int f1 = f.first;
    int f2 = f.second;
    if(f1!=-1 && f2!=-1){ //must have 2 faces
      bool idF1 = cH.at(f1); //is a face of cH
      bool idF2 = cH.at(f2); //is a face of cH
      if((idF1 && !idF2) || (!idF1 && idF2)) { //one face inside and one outside
        isBorder.at(f1) = true;
        isBorder.at(f2) = true;
      }
    }
  }
  return isBorder;
}

//Compute tau fct indicating whether a face h2 of gridH is simple wrt complex H
bool isSimpleFace(const Complex& gridH, const std::vector<bool>& cH, int idF) {
  //Compute vetex border
  int nbVertexBorder = 0;
  std::vector<int> cell_vertex = gridH.cell_face_vertex.at(idF); //cell of faces (=set of vertices by id)
  for(size_t it=0; it<cell_vertex.size(); it++) { //For each vertex of idF, verify if it is a border vertex
    std::vector<int> star_vertex = gridH.star_vertex_face.at(cell_vertex.at(it)); // Get star_face of each vertex
    int nbFace = 0; //compute nb of faces it belongs to
    for(size_t it_bis=0; it_bis<star_vertex.size(); it_bis++) {
      if(cH.at(star_vertex.at(it_bis)))
        nbFace++;
    }
    if (nbFace==1) //vertex is border if it belonsg to one face only
      nbVertexBorder++;
  }
  //Compute edge border
  int nbEdgeBorder = 0;
  std::vector<int> cell_edge= gridH.cell_face_edge.at(idF); //cell of faces (=set of edges by id)
  for(size_t it=0; it<cell_edge.size(); it++) { //For each edge of idF, verify if it is a border edge
    std::pair<int,int> star_edge = gridH.star_edge.at(cell_edge.at(it)); // Get star of each edge
    int f1 = star_edge.first;
    int f2 = star_edge.second;
    if(f1!=-1 && f2!=-1){ //must have 2 faces
      bool idF1 = cH.at(f1); //is a face of cH
      bool idF2 = cH.at(f2); //is a face of cH
      if((idF1 && !idF2) || (!idF1 && idF2)) { //one face inside and one outside
        nbEdgeBorder++;
      }
    }
    else //one face in, one face out
      nbEdgeBorder++;
  }
  return (nbEdgeBorder-nbVertexBorder)==1;
}

//Compute tau fct indicating whether a face h2 of complex H is simple
std::vector<bool>
tau (const Complex& gridH, const std::vector<bool>& cH) {
  std::vector<bool> isSimple(gridH.face.size(), false);
  std::vector<bool> borderH = getBorderFace(gridH, cH);
  bool simple = false;
  for(size_t it=0; it<gridH.face.size(); it++) {
    simple = false;
    if(cH.at(it)==true) { //inside face
      simple = isSimpleFace(gridH,cH,it); //verify the face simplicity
    }
    else { //outside face
      if(borderH.at(it)) { //and a border face
        std::vector<bool> cH_tmp = cH;
        cH_tmp.at(it) = true;
        simple = isSimpleFace(gridH,cH_tmp,it); //verify new face simplicity
      }
      else //not simple by default
        simple = false;
    }
    isSimple.at(it) = simple;
  }
  return isSimple;
}

//Count number of complex H in the current grid
int countNbFaceComplex(const std::vector<bool>& cH) {
  int nbFace = 0;
  for(size_t it=0; it<cH.size(); it++)
    if(cH.at(it))
      nbFace++;
  return nbFace;
}

//Initialize the gauss solution
std::vector<bool> initizeGauss(const Complex &gridH, const std::vector<bool>& cH) {
  std::vector<Z2i::Point> vP;
  RationalPoint p, pp;
  Z2i::Point px;
  for(size_t it=0; it<cH.size(); it++) {
    if(cH.at(it)==true) { //a face of CH
      p = gridH.getFaceCenter(it);
      px = Z2i::Point(round(getRealValue(p.first)),round(getRealValue(p.second)));
      pp = RationalPoint(px[0], px[1]);
      if(!containElement(vP, px) && isInsidePolygon(gridH.getFaceVertices(it), pp))
        vP.push_back(px);
    }
  }
  //reproject into gH
  std::vector<bool> Hinit (gridH.face.size(), false);
  for(size_t it=0; it<vP.size(); it++) {
    px = vP.at(it);
    for(size_t it_bis=0; it_bis<gridH.face.size(); it_bis++) {
      p = gridH.getFaceCenter(it_bis);
      if(isInsidePixel(px, p))
        Hinit.at(it_bis) = true;
    }
  }
  return Hinit;
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

//Find index of complex H face of max epsilon value
size_t
maxIndexEpsilon(const std::vector<double>& epsilonValue, const std::vector<bool>& isBorder) {
  assert(isBorder.size()==epsilonValue.size());
  double maxEpsilon = -1;
  std::vector<size_t> epsilonValueIndex = sort_indexes_decrease(epsilonValue);
  int index = 0;
  int maxIndex = epsilonValueIndex.at(index);//-1;
  while (!isBorder.at(maxIndex)) {
    index++;
    maxIndex = epsilonValueIndex.at(index);
  }
  return maxIndex;
}

//Find a random index of complex H face
size_t
randomIndexEpsilon(const std::vector<bool>& isBorder, const std::vector<bool>& isSimple) {
  int index = std::rand() % isBorder.size();
  while(!isSimple.at(index) || !isBorder.at(index))
    index = std::rand() % isBorder.size();
  return index;
}

//Find a index of complex H face according to a stats function
size_t
statIndexEpsilon(const std::vector<double>& epsilonValue, const std::vector<bool>& isBorder, const std::vector<bool>& isSimple) {
  double coefficient = 1000000.0;
  std::vector<size_t> epsilonValueIndex = sort_indexes_decrease(epsilonValue);
  double x = double(rand()) / double(RAND_MAX);
  double p = log(1.0+coefficient*x)/log(coefficient + 1.0);//log(1.0+x)/log(2.0);
  p = p*epsilonValue.size();
  int indexStat = floor(p);
  int index = epsilonValueIndex.at(indexStat);
  while(!isSimple.at(index) || !isBorder.at(index)) {
    x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    p = log(1.0+coefficient*x)/log(coefficient + 1.0);//log(1.0+x)/log(2.0);
    p = p*epsilonValue.size();
    indexStat = floor(p);
    index = epsilonValueIndex.at(indexStat);
  }
  std::cout<<"statIndexEpsilon, p="<<p<<" and index="<<indexStat<<std::endl;
  return index;
}

std::vector<Z2i::Point>
homotopicAffineMajority(std::vector<Z2i::Point> X, AffineTransform t, string file, int border = 2, bool saveFile = false) {
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
  vector<RationalPoint> borderGridG;
  borderGridG.push_back(gridG.vertical_lines.front().second);
  borderGridG.push_back(gridG.vertical_lines.front().first);
  borderGridG.push_back(gridG.vertical_lines.back().first);
  borderGridG.push_back(gridG.vertical_lines.back().second);
  borderGridG.push_back(gridG.vertical_lines.front().second);
  vector<vector<RationalPoint> > borderGridG_cells;
  borderGridG_cells.push_back(borderGridG);
  Complex tgridG(borderGridG_cells);
  
  //Embed into gridH
  Complex gridH(h_cells);//gridH(h_pure_cells);//SIMPLIFIED
  gridH.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  gridH.addGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  gridH.updateAllCell();
  gridH.updateAllStar();
  
  Board2D aBoard;
  string filename;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, Color::Gray);
    filename = file + "_gridH.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  Complex cX(X);
  cX.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  cX.drawGrid(aBoard);
  cX.drawPixel(aBoard);
  filename = file + "_objetX.svg";
  aBoard.saveSVG(filename.c_str());
  aBoard.clear();
  
  
  //ComplexF
  Complex complexF;
  complexF.init(tp1,tp2);
  for(int it=0; it<X.size(); it++) {
    complexF.addCubicalFace(X.at(it));
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
  
  //Embed complex G into the grid H
  std::vector<bool> belongVertex, belongEdge, belongFace;
  embedComplexInGrid(complexG,gridH,belongFace);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<belongFace.size(); it++) {
      if(belongFace.at(it)==true)
        gridH.drawFace(aBoard, it, Color::Gray);
    }
    filename = file + "_complexH.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Prepare Htarget structure
  std::vector<bool> stateFaceHtarget = belongFace;
  
  //Compute h2 area of Htarget (wrt Hgrid)
  std::vector<double> areaFace(gridH.face.size(),0.0);
  for(size_t it=0; it<gridH.face.size(); it++) {
    if(stateFaceHtarget.at(it))
      areaFace.at(it) = gridH.getFaceArea(it);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(areaFace.at(it)>0) //face of Hinit
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
    filename = file + "_areaFaceMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Compute fill area f2 of Htarget (wrt Hgrid)
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
        if(stateFaceHtarget.at(it_bis)) // a face of Hinit
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
 
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(fillAreaFace.at(it)<EPSILON || fabs(fillAreaFace.at(it)-1)<EPSILON) //pure pixels
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else {
        if(fillAreaFace.at(it)>0.5) //majority interieur
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
        else //majority exterieur
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
    filename = file + "_fillAreaFaceMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 1
  std::vector<bool> stateFaceHcurrent = stateFaceHtarget;
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
    // Creating colormap
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
    filename = file + "_epsilonMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Tau fct on h2
  std::vector<bool> tauValue = tau(gridH, stateFaceHcurrent);
  assert(tauValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    filename = file + "_tauMajority.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 9: Compute B2 of Hcurrent (of same index as gridH)
  std::vector<bool> borderHcurrent = getBorderFace(gridH, stateFaceHcurrent);
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
    std::cout<<"nbIter="<<nbIter<<": h2="<<h2<<", f2="<<f2<<", prior="<<epsilonValue.at(h2);
    //Line 13
    std::vector<int> inFace = Phi.at(f2);
    a = fillAreaFace.at(h2);
    if(tauValue.at(h2)) {//Line 13: if simple cell
      //Line 14
      if(iotaValue.at(h2)) //remove face
        stateFaceHcurrent.at(h2) = false;
      else //add face
        stateFaceHcurrent.at(h2) = true;
      //Line 15: Update iota of h2
      iotaValue.at(h2) = !iotaValue.at(h2); //Inverse state face
      
      std::cout<<", Hcurrent="<<countNbFaceComplex(stateFaceHcurrent)<<", cH="<<gridH.face.size()<<std::endl;
      if(saveFile) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(iotaValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Black);//Green
          else
            gridH.drawFace(aBoard, it, DGtal::Color::White);//Red
        }
        filename = file + "_HcurrentMajority_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
      
      //Line 16 => out of condition
      if(iotaValue.at(h2))
        a += areaFace.at(h2);
      else
        a -= areaFace.at(h2);
      for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) //update fill area for all inside faces
        fillAreaFace.at(inFace.at(it_bis)) = a;
      //Line 17-19
      std::vector<bool> tauHcurrent = tau(gridH, stateFaceHcurrent);
      for(size_t it=0; it<gridH.face.size(); it++)
        tauValue.at(it) = tauHcurrent.at(it);
    }
    else { //Line 20-21
      sigmaValue.at(f2) = !sigmaValue.at(f2);
      std::cout<<", Hcurrent="<<countNbFaceComplex(stateFaceHcurrent)<<", cH="<<gridH.face.size()<<std::endl;
      //f2 non simple
      gridH.drawGrid(aBoard);
      for(size_t it=0; it<gridH.face.size(); it++) {
        if(iotaValue.at(it))
          gridH.drawFace(aBoard, it, DGtal::Color::Black);//Green
        else
          gridH.drawFace(aBoard, it, DGtal::Color::White);//Red
      }
      gridH.drawFace(aBoard, h2, DGtal::Color::Red);
      filename = file + "_HcurrentMajority_s" + std::to_string(nbIter) + ".svg";
      aBoard.saveSVG(filename.c_str());
      aBoard.clear();
    }
    
    //Line 16/22: Update zeta of f2
    s = sigmaValue.at(f2) ? 1 : -1;
    zetaValue.at(f2) = (1.0-s)/2.0 + s*a;
    //Line 23-25
    s2 = sigmaValue.at(f2) ? 1 : -1;
    a = zetaValue.at(f2);
    for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) {
      s1 = iotaValue.at(inFace.at(it_bis)) ? 1 : -1;
      epsilonValue.at(inFace.at(it_bis)) = -s1*s2*a;
    }
    //Line 25: Update border faces
    borderHcurrent = getBorderFace(gridH, stateFaceHcurrent);
    //Line 26
    h2 = maxIndexEpsilon(epsilonValue, borderHcurrent);
    //h2 = randomIndexEpsilon(borderHcurrent, tauValue);
    //h2 = statIndexEpsilon(epsilonValue, borderHcurrent, tauValue);
    //Line 26
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

  Complex cY(Y);
  cY.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  cY.drawGrid(aBoard);
  cY.drawPixel(aBoard);
  filename = file + "_objetYMajority.svg";
  aBoard.saveSVG(filename.c_str());
  aBoard.clear();
  
  return Y;
}

std::vector<Z2i::Point>
homotopicAffineGauss(std::vector<Z2i::Point> X, AffineTransform t, string file, int border = 2, bool saveFile = false) {
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
  vector<RationalPoint> borderGridG;
  borderGridG.push_back(gridG.vertical_lines.front().second);
  borderGridG.push_back(gridG.vertical_lines.front().first);
  borderGridG.push_back(gridG.vertical_lines.back().first);
  borderGridG.push_back(gridG.vertical_lines.back().second);
  borderGridG.push_back(gridG.vertical_lines.front().second);
  vector<vector<RationalPoint> > borderGridG_cells;
  borderGridG_cells.push_back(borderGridG);
  Complex tgridG(borderGridG_cells);
  
  //Embed into gridH
  Complex gridH(h_cells);//gridH(h_pure_cells); //SIMPLIFIED
  gridH.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  gridH.addGridLines(gridG.vertical_lines, gridG.horizontal_lines);
  gridH.updateAllCell();
  gridH.updateAllStar();
  
  Board2D aBoard;
  string filename;
  if(saveFile) {
    gridH.drawGrid(aBoard);
    gridH.drawComplex(aBoard, Color::Gray);
    filename = file + "_gridH.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  Complex cX(X);
  cX.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  cX.drawGrid(aBoard);
  cX.drawPixel(aBoard);
  filename = file + "_objetX.svg";
  aBoard.saveSVG(filename.c_str());
  aBoard.clear();
  
  //ComplexF
  Complex complexF;
  complexF.init(tp1,tp2);
  for(int it=0; it<X.size(); it++) {
    complexF.addCubicalFace(X.at(it));
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
  
  //Embed complex G into the grid H
  std::vector<bool> belongVertex, belongEdge, belongFace;
  embedComplexInGrid(complexG,gridH,belongFace);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<belongFace.size(); it++) {
      if(belongFace.at(it)==true)
        gridH.drawFace(aBoard, it, Color::Gray);
    }
    filename = file + "_complexH.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Prepare Htarget structure
  std::vector<bool> stateFaceHtarget = initizeGauss(gridH, belongFace);
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<stateFaceHtarget.size(); it++) {
      if(stateFaceHtarget.at(it)==true)
        gridH.drawFace(aBoard, it, Color::Gray);
    }
    filename = file + "_HinitGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  int sum = 0;
  for(size_t it = 0; it<stateFaceHtarget.size(); it++)
    if(stateFaceHtarget.at(it))
      sum++;
  assert(stateFaceHtarget.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(stateFaceHtarget.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    filename = file + "_stateFaceTargetGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Compute fill area f2 of Htarget (wrt Hgrid)
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
      if(fillAreaHtargetFace.at(it)<EPSILON || fabs(fillAreaHtargetFace.at(it)-1)<EPSILON) //pure pixels
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else {
        if(fillAreaHtargetFace.at(it)>0.5) //majority interieur
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
        else //majority exterieur
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
    filename = file + "_fillAreaFaceHtargetGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 1
  std::vector<bool> stateFaceHcurrent = belongFace;
  sum = 0;
  for(size_t it = 0; it<stateFaceHcurrent.size(); it++)
    if(stateFaceHcurrent.at(it))
      sum++;
  assert(stateFaceHcurrent.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(stateFaceHcurrent.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    filename = file + "_stageFaceGauss.svg";
    aBoard.clear();
  }
  //Prepare Hcurrent structure
  //Compute h2 area of Hinit (wrt Hgrid)
  std::vector<double> areaFace(gridH.face.size(),0.0);
  for(size_t it=0; it<gridH.face.size(); it++) {
    if(stateFaceHcurrent.at(it))
      areaFace.at(it) = gridH.getFaceArea(it);
  }
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(areaFace.at(it)>0) //face of Hinit
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
    }
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
        if(stateFaceHcurrent.at(it_bis)) // a face of Hinit
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
  
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(fillAreaFace.at(it)<EPSILON || fabs(fillAreaFace.at(it)-1)<EPSILON) //pure pixels
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
      else {
        if(fillAreaFace.at(it)>0.5) //majority interieur
          gridH.drawFace(aBoard, it, DGtal::Color::Green);
        else //majority exterieur
          gridH.drawFace(aBoard, it, DGtal::Color::Blue);
      }
    }
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
    // Creating colormap
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
    filename = file + "_epsilonGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Tau fct on h2
  std::vector<bool> tauValue = tau(gridH, stateFaceHcurrent);
  assert(tauValue.size()==gridH.face.size());
  if(saveFile) {
    gridH.drawGrid(aBoard);
    for(size_t it=0; it<gridH.face.size(); it++) {
      if(tauValue.at(it))
        gridH.drawFace(aBoard, it, DGtal::Color::Green);
      else
        gridH.drawFace(aBoard, it, DGtal::Color::Red);
    }
    filename = file + "_tauGauss.svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
  }
  
  //Line 9: Compute B2 of Hcurrent (of same index as gridH)
  std::vector<bool> borderHcurrent = getBorderFace(gridH, stateFaceHcurrent);
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
    std::cout<<"nbIter="<<nbIter<<": h2="<<h2<<", f2="<<f2<<", prior="<<epsilonValue.at(h2);
    std::vector<int> inFace = Phi.at(f2);
    if(tauValue.at(h2)) {//Line 13: if simple cell
      //Line 14
      if(iotaValue.at(h2)) //remove face
        stateFaceHcurrent.at(h2) = false;
      else //add face
        stateFaceHcurrent.at(h2) = true;
      //Line 15: Update iota of h2
      iotaValue.at(h2) = !iotaValue.at(h2); //Inverse state face
      
      std::cout<<", Hcurrent="<<countNbFaceComplex(stateFaceHcurrent)<<", cH="<<gridH.face.size()<<std::endl;
      if(saveFile) {
        gridH.drawGrid(aBoard);
        for(size_t it=0; it<gridH.face.size(); it++) {
          if(iotaValue.at(it))
            gridH.drawFace(aBoard, it, DGtal::Color::Black);//Green);
          else
            gridH.drawFace(aBoard, it, DGtal::Color::White);//Red);
        }
        filename = file + "_HcurrentGauss_s" + std::to_string(nbIter) + ".svg";
        aBoard.saveSVG(filename.c_str());
        aBoard.clear();
      }
      
      //Line 16 => out of condition
      a = fillAreaFace.at(h2);
      if(iotaValue.at(h2))
        a += areaFace.at(h2);
      else
        a -= areaFace.at(h2);
      for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) //update fill area for all inside faces
        fillAreaFace.at(inFace.at(it_bis)) = a;
      //Line 17-19
      std::vector<bool> tauHcurrent = tau(gridH, stateFaceHcurrent);
      for(size_t it=0; it<gridH.face.size(); it++)
        tauValue.at(it) = tauHcurrent.at(it);
    }
    else { //Line 20-21
      sigmaValue.at(f2) = !sigmaValue.at(f2);
      std::cout<<", Hcurrent="<<countNbFaceComplex(stateFaceHcurrent)<<", cH="<<gridH.face.size()<<std::endl;
    }
    //Line 16/22: Update zeta of f2
    s = sigmaValue.at(f2) ? 1 : -1;
    zetaValue.at(f2) = (1.0-s)/2.0 + s*a;
    //Line 23-25
    //std::vector<int> inFace = Phi.at(f2);
    s2 = sigmaValue.at(f2) ? 1 : -1;
    a = zetaValue.at(f2);
    for(size_t it_bis=0; it_bis<inFace.size(); it_bis++) {
      s1 = iotaValue.at(inFace.at(it_bis)) ? 1 : -1;
      epsilonValue.at(inFace.at(it_bis)) = -s1*s2*a;
      //Line 25 out of loop
    }
    
    //Line 25: Update border faces
    borderHcurrent = getBorderFace(gridH, stateFaceHcurrent);
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
        gridH.drawFace(aBoard, it, DGtal::Color::Black);//Green
      else
        gridH.drawFace(aBoard, it, DGtal::Color::White);//Red
    }
    filename = file + "_HcurrentMajority_s" + std::to_string(nbIter) + ".svg";
    aBoard.saveSVG(filename.c_str());
    aBoard.clear();
    
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

  Complex cY(Y);
  cY.initGridLines(tgridF.vertical_lines, tgridF.horizontal_lines);
  cY.drawGrid(aBoard);
  cY.drawPixel(aBoard);
  filename = file + "_objetYGauss.svg";
  aBoard.saveSVG(filename.c_str());
  aBoard.clear();
  
  return Y;
}

int main(int argc, char**argv) {
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFile = "../Samples/ball_r3.pgm";
  string outputDir = "./";
  int model {0}; //0 : majority vote and 1 : gauss
  bool saveFile {false}; //save intermediate files
  //affine paramters
  double a11 {1.0};
  double a12 {0.0};
  double a21 {0.0};
  double a22 {1.0};
  double tx {0.0};
  double ty {0.0};
  
  app.description("Apply homotopic affine transforma on a given image.\n Example:\n \t HomotopicAffineTransform --input <PgmFileName> --outputDir <OutputDir> -s --a11 1.5 --a12 0.2 --a21 0.5 --a22 1.2 --tx 0 --ty 0 -m <0|1>\n");
  app.add_option("-i,--input,1",inputFile,"Input file.")->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output",outputDir,"Output directory (default ./).");
  app.add_option("-m,--model",model,"Transformation model: Majority vote (0, default), Gauss (1).");
  auto saveF = app.add_flag("-s,--save",saveFile,"Save all steps (defaut no).");
  if(saveF->count() > 0)
    saveFile = true;
  app.add_option("--a11",a11,"affine transform parameter (default 1.0)",true);
  app.add_option("--a12",a12,"affine transform parameter (default 0.0)",true);
  app.add_option("--a21",a21,"affine transform parameter (default 0.0)",true);
  app.add_option("--a22",a22,"affine transform parameter (default 1.0)",true);
  app.add_option("--tx",tx,"X component of translation vector (default 0.0)",true);
  app.add_option("--ty",ty,"Y component of translation vector (default 0.0)",true);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  int border = 2;
  
  std::string extenstion = inputFile.substr(inputFile.find_last_of(".")+1);
  if(extenstion != "pgm") {
    std::cout<<"Error: input image must be in pgm format !"<<std::endl;
    exit(-1);
  }
  string filename = inputFile.substr(inputFile.find_last_of("/")+1);
  filename = filename.substr(0, filename.find_last_of("."));
  string outputFile = outputDir + "/" + filename;
  
  std::vector<Z2i::Point> X = readImage(inputFile);
  assert(X.size()!=0);
  std::cout<<"X.size="<<X.size()<<std::endl;
  
  AffineTransform t(a11,a12,a21,a22,tx,ty);
  
  if(model==1)
    homotopicAffineGauss(X, t, outputFile, border, saveFile);
  else
    homotopicAffineMajority(X, t, outputFile, border, saveFile);
  
  return 1;
}
