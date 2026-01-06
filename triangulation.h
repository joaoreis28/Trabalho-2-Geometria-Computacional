#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "opencv2/opencv.hpp"
#include "CGAL/Surface_mesh.h"

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
using Halfedge = CGAL::Surface_mesh<Point>::Halfedge_index;
using Vertex = CGAL::Surface_mesh<Point>::Vertex_index;
using Face = CGAL::Surface_mesh<Point>::Face_index;

namespace triangulation {

struct Location {
    enum Type { FACE, LINE, VERTEX, OUTSIDE };
    Type t;
    Face face;
    Halfedge halfedge;
    Vertex vertex;

    Location(Face face);
    Location(Halfedge halfedge);
    Location(Vertex vertex);
    Location();
};

class Delaunay {
   private:
    CGAL::Surface_mesh<Point> mesh;
    void add_point_to_face(Point p, Face f);
    void add_point_to_edge(Point p, Halfedge he);
    void flip(Halfedge he);
   public:
    Delaunay(int img_w, int img_h);
    void draw(cv::Mat& img);
    void add_point(Point p);
    Location locate(Point p);
    Point point(Vertex v);
    std::array<Point, 2> points(Halfedge l);
    std::array<Point, 3> points(Face f);
};

}  
