#ifndef AFFINETRANSFORM_H
#define AFFINETRANSFORM_H

#include "Functions.h"
#include "Complex.h"

class AffineTransform
{
    public :
        Rational a11;
        Rational a12;
        Rational a21;
        Rational a22;
        Rational tx;
        Rational ty;
    
    AffineTransform () : a11(1), a12(1), a21(1), a22(1), tx(0), ty(0) {}
    AffineTransform (Rational v11, Rational v12, Rational v21, Rational v22, Rational x, Rational y);
    AffineTransform (double v11, double v12, double v21, double v22);
    AffineTransform (double v11, double v12, double v21, double v22, double x, double y);
    AffineTransform (double v11, double v12, double v21, double v22, int t1x, int t1y, int t2x, int t2y);
    friend std::ostream& operator<<(std::ostream& os, const AffineTransform& t);
    
    AffineTransform setCoeff(Rational v11, Rational v12, Rational v21, Rational v22);
    AffineTransform setCoeff(int v11, int v12, int v21, int v22);
    AffineTransform setCoeff(double v11, double v12, double v21, double v22);
    AffineTransform setTranslation(int x, int y);
    AffineTransform setTranslation(double x, double y);
    AffineTransform setTranslation(Rational x, Rational y);
    AffineTransform setTranslation(int t1x, int t1y, int t2x, int t2y);
    
    RationalPoint transform(RationalPoint p);
    RationalPoint transform(Z2i::Point p);
    std::pair<RationalPoint,RationalPoint> transform(RationalPoint p1, RationalPoint p2);
    std::pair<RationalPoint,RationalPoint> transform(std::pair<RationalPoint,RationalPoint> p);
    //Line transform(Line l);
};

RationalPoint transformPoint(RationalPoint p, AffineTransform t);
std::pair<RationalPoint,RationalPoint> transformSegment(RationalPoint p1, RationalPoint p2, AffineTransform t);
std::pair<RationalPoint,RationalPoint> transformSegment(std::pair<RationalPoint,RationalPoint> p, AffineTransform t);
Complex transform(Complex c, AffineTransform t);

#endif // AFFINETRANSFORM_H
