#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

struct Point { double x; double y; };

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
static inline int isLeft(Point P0, Point P1, Point P2)
{
  return ( (P1.x - P0.x) * (P2.y - P0.y)
           - (P2.x -  P0.x) * (P1.y - P0.y) );
}

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
int cn_PnPoly(Point P, Point* V, int n)
{
  int    cn = 0;    // the  crossing number counter
    
  // loop through all edges of the polygon
    for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
      if (((V[i].y <= P.y) && (V[i+1].y > P.y))     // an upward crossing
            || ((V[i].y > P.y) && (V[i+1].y <=  P.y))) { // a downward crossing
        // compute  the actual edge-ray intersect x-coordinate
        float vt = (float)(P.y  - V[i].y) / (V[i+1].y - V[i].y);
        if (P.x <  V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
          ++cn;   // a valid crossing of y=P.y right of P.x
      }
    }
    return (cn&1);    // 0 if even (out), and 1 if  odd (in)
}
//===================================================================


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
int
  wn_PnPoly( Point P, Point* V, int n )
  {
    int    wn = 0;    // the  winding number counter
    
    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
      if (V[i].y <= P.y) {          // start y <= P.y
        if (V[i+1].y  > P.y)      // an upward crossing
          if (isLeft( V[i], V[i+1], P) > 0)  // P left of  edge
            ++wn;            // have  a valid up intersect
      }
      else {                        // start y > P.y (no test needed)
        if (V[i+1].y  <= P.y)     // a downward crossing
          if (isLeft( V[i], V[i+1], P) < 0)  // P right of  edge
            --wn;            // have  a valid down intersect
      }
    }
    return wn;
  }
//===================================================================


// [[Rcpp::export(".point_in_polygon_cpp")]]
LogicalVector point_in_polygon(NumericVector x, NumericVector y, NumericVector poly_x, NumericVector poly_y)
{
  std::vector<Point> polygon;
  for (int i = 0; i != poly_x.size(); i++) {
    polygon.emplace_back(Point{poly_x[i], poly_y[i]});
  }
  polygon.emplace_back(Point{poly_x[0], poly_y[0]});
  LogicalVector inside_polygon(x.size());
  for (int i = 0; i != x.size(); i++) {
    auto wn = wn_PnPoly(Point{x[i], y[i]}, polygon.data(), polygon.size()-1);
    inside_polygon[i] = wn != 0;
  }
  return inside_polygon;
}
