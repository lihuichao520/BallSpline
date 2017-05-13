#include "nurbs.h"

namespace PLib {

// double specialization

template <>
int intersectLine(const Point_nD<double,2>& p1, const Point_nD<double,2>& t1, const Point_nD<double,2>& p2, const Point_nD<double,2>& t2, Point_nD<double,2>& p){
  cout << "PLEASE, DEFINE THIS FUNCTION\n" ; 

  return 1 ;
}

template <>
Point_nD<double,2> NurbsCurve<double,2>::normal(double u, const Point_nD<double,2>& v) const{
  cerr << "YOU CAN'T COMPUTE THE NORMAL in 2D of a 2D vector!\n" ; 
  return firstDn(u) ;
}

void NurbsCurve<double,2>::makeCircle(const Point_nD<double,2>& O, double r, double as, double ae){
  makeCircle(O,Point_nD<double,2>(1,0),Point_nD<double,2>(0,1),r,as,ae) ;
}

template <>
int NurbsCurve<double,2>::writeVRML(const char* filename,double radius,int K, const Color& color,int Nu,int Nv, double u_s, double u_e) const{
  NurbsCurve<double,3> C ;
  to3D(*this,C) ; 
  return C.writeVRML(filename,radius,K,color,Nu,Nv,u_s,u_e) ;
}

template <>
int NurbsCurve<double,2>::writeVRML97(const char* filename,double radius,int K, const Color& color,int Nu,int Nv, double u_s, double u_e) const{
  NurbsCurve<double,3> C ;
  to3D(*this,C) ; 
  return C.writeVRML97(filename,radius,K,color,Nu,Nv,u_s,u_e) ;
}

template <>
int NurbsCurve<double,2>::writeVRML(ostream& fout,double radius,int K, const Color& color,int Nu,int Nv, double u_s, double u_e) const{
  NurbsCurve<double,3> C ;
  to3D(*this,C) ; 
  return C.writeVRML(fout,radius,K,color,Nu,Nv,u_s,u_e) ;
}

template <>
int NurbsCurve<double,2>::writeVRML97(ostream& fout,double radius,int K, const Color& color,int Nu,int Nv, double u_s, double u_e) const{
  NurbsCurve<double,3> C ;
  to3D(*this,C) ; 
  return C.writeVRML97(fout,radius,K,color,Nu,Nv,u_s,u_e) ;
}

template <>
void NurbsCurve<double,2>::drawAaImg(Image_Color& Img, const Color& color, int precision, int alpha){
  NurbsCurve<double,3> C ;
  to3D(*this,C) ; 
  C.drawAaImg(Img,color,precision,alpha) ;
}



}

