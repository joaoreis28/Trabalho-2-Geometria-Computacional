#include "triangulation.h"

#include <CGAL/draw_surface_mesh.h>

#include <filesystem>
#include <ostream>

using namespace triangulation;

Location::Location(Vertex vertex) : t(Type::VERTEX), vertex(vertex) {}
Location::Location(Halfedge halfedge) : t(Type::LINE), halfedge(halfedge) {}
Location::Location(Face face) : t(Type::FACE), face(face) {}
Location::Location() : t(Type::OUTSIDE) {}

bool is_in_line(Point p1, Point p2, Point p)
{
    if (!CGAL::collinear(p1, p2, p)) return false;
    auto between = [](int a, int b, int c) -> bool
    {
        return c > std::min(a, b) && c < std::max(a, b);
    };

    if (p1.x() != p2.x())
        return between(p1.x(), p2.x(), p.x());
    else
        return between(p1.y(), p2.y(), p.y());
    
}

bool is_in_triangle(Point p1, Point p2, Point p3, Point p)
{
    bool left = CGAL::left_turn(p1, p2, p) && CGAL::left_turn(p2, p3, p) &&
                CGAL::left_turn(p3, p1, p);
    bool right = CGAL::right_turn(p1, p2, p) && CGAL::right_turn(p2, p3, p) &&
                 CGAL::right_turn(p3, p1, p);
    return left || right;
}

bool is_convex_quadrilateral(Point p1, Point p2, Point p3, Point p4)
{
    return !is_in_triangle(p1, p2, p3, p4);
}

void Delaunay::flip(Halfedge he) {
    if (mesh.is_border(he) || mesh.is_border(mesh.opposite(he))) return;
    Point p1 = mesh.point(mesh.target(he));
    Point p2 = mesh.point(mesh.target(mesh.next(he)));
    Point p3 = mesh.point(mesh.source(he));
    Point p4 = mesh.point(mesh.target(mesh.next(mesh.opposite(he))));

    if (!is_convex_quadrilateral(p1, p2, p3, p4)) return;
    Halfedge h1 = mesh.next(he);
    Halfedge h2 = mesh.next(h1);
    Halfedge h3 = mesh.next(mesh.opposite(he));
    Halfedge h4 = mesh.next(h3);
    if (CGAL::side_of_oriented_circle(p1, p3, p4, p2) ==
        CGAL::ON_NEGATIVE_SIDE) {
        std::cout << "is in circle" << std::endl;
        CGAL::Euler::flip_edge(he, mesh);
        flip(h1);
        flip(h2);
        flip(h3);
        flip(h4);
    }
}


Point Delaunay::point(Vertex v) {
    return mesh.point(v);
}

std::array<Point, 2> Delaunay::points(Halfedge l) {
    return {mesh.point(mesh.source(l)), mesh.point(mesh.target(l))};
}

std::array<Point, 3> Delaunay::points(Face f) {
    Halfedge he = mesh.halfedge(f);
    return {
        mesh.point(mesh.source(he)),
        mesh.point(mesh.target(he)),
        mesh.point(mesh.target(mesh.next(he)))
    };
}

void Delaunay::add_point_to_face(Point p, Face f) {
    std::cout << "adding point to face" << std::endl;
    Halfedge h = mesh.halfedge(f);
    Halfedge new_he = CGAL::Euler::add_center_vertex(h, mesh);
    Vertex vp = mesh.target(new_he);
    mesh.point(vp) = p;
    Halfedge h1 = mesh.next(mesh.opposite(new_he));
    Halfedge h2 = mesh.next(mesh.opposite(mesh.next(h1)));
    Halfedge h3 = mesh.next(mesh.opposite(mesh.next(h2)));
    flip(h1);
    flip(h2);
    flip(h3);
}

void Delaunay::add_point_to_edge(Point p, Halfedge he) {
    std::cout << "adding point to edge" << std::endl;
    Halfedge new_he = CGAL::Euler::split_edge(he, mesh);
    Vertex vp = mesh.target(new_he);
    mesh.point(vp) = p;

    Halfedge new_he_next = mesh.next(new_he);
    if (!mesh.is_border(new_he_next)) {
        CGAL::Euler::split_face(new_he, mesh.next(new_he_next), mesh);
    }
    Halfedge new_he_opp = mesh.opposite(new_he);
    if (!mesh.is_border(new_he_opp)) {
        CGAL::Euler::split_face(mesh.opposite(new_he_next), mesh.next(new_he_opp), mesh);
    }

    Halfedge h1 = mesh.next(new_he_opp);
    Halfedge h2 = mesh.prev(mesh.opposite(new_he_next));
    Halfedge h3 = mesh.next(new_he_next);
    Halfedge h4 = mesh.prev(new_he);
    flip(h1);
    flip(h2);
    flip(h3);
    flip(h4);
}

void Delaunay::add_point(Point p) {
    Location l = locate(p);
    using namespace triangulation;
    switch (l.t) {
        case Location::Type::OUTSIDE:
            std::cout << "OUTSIDE" << std::endl;
            return;
        case Location::Type::VERTEX:
            std::cout << "VERTEX" << std::endl;
            return;
        case Location::Type::FACE:
            std::cout << "FACE" << std::endl;
            add_point_to_face(p, l.face);
            break;
        case Location::Type::LINE:
            std::cout << "LINE" << std::endl;
            add_point_to_edge(p, l.halfedge);
            break;
    }
}

Location Delaunay::locate(Point p) {
    Halfedge h = *mesh.halfedges_begin();
    if (mesh.is_border(h)) {
        h = mesh.opposite(h);
    }

    int path_size = 0;

    while (path_size < mesh.num_halfedges()) {
        path_size++;

        Halfedge h_start = h;
        bool found_face = true;

        for (int i = 0; i < 3; i++) {
            Point p1 = mesh.point(mesh.source(h));
            Point p2 = mesh.point(mesh.target(h));

            if (p == p1) return Location(mesh.source(h));
            if (p == p2) return Location(mesh.target(h));
            if (is_in_line(p1, p2, p)) return Location(h);

            if (CGAL::left_turn(p1, p2, p)) {
                Halfedge opp = mesh.opposite(h);
                if (mesh.is_border(opp)) {
                    return Location();  // Outside mesh
                }
                h = mesh.next(opp);
                found_face = false;
                break;
            }

            h = mesh.next(h);
        }

        if (found_face) {
            return Location(mesh.face(h_start));
        }
    }

    return Location();
}

Delaunay::Delaunay(int img_w, int img_h) {
    Vertex p1 = mesh.add_vertex(Point(0, 0));
    Vertex p2 = mesh.add_vertex(Point(0, img_h - 1));
    Vertex p3 = mesh.add_vertex(Point(img_w - 1, img_h - 1));
    Vertex p4 = mesh.add_vertex(Point(img_w - 1, 0));
    mesh.add_face(p1, p2, p3);
    mesh.add_face(p1, p3, p4);
}

void Delaunay::draw(cv::Mat& img) {
    for (Halfedge he : mesh.halfedges()) {
        if (he < mesh.opposite(he)) {
            Point p1 = mesh.point(mesh.source(he));
            Point p2 = mesh.point(mesh.target(he));
            cv::line(img, cv::Point((int)p1.x(), (int)p1.y()),
                     cv::Point((int)p2.x(), (int)p2.y()),
                     cv::Scalar(255, 255, 255), 2);
        }
    }
    for (Vertex v : mesh.vertices()) {
        Point p1 = mesh.point(v);
        cv::circle(img, cv::Point((int)p1.y(), (int)p1.x()), 2,
                   cv::Scalar(0, 0, 255), -1);
    }
}
