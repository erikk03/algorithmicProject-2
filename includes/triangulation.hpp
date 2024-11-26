#ifndef TRIANGULATION_HPP
#define TRIANGULATION_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polygon_2.h>

#include <sstream>
#include <cmath>

// Define Kernel
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;

// Define a custom vertex base to hold the `info` field (the index)
typedef CGAL::Triangulation_vertex_base_with_info_2<int, Kernel> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

using Vector = Kernel::Vector_2;

// Custom CDT class with new insert methods
template <class Gt, class Tds = CGAL::Default, class Itag = CGAL::Default>

class Custom_Constrained_Delaunay_triangulation_2

    : public CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>
{

public:
    using Base = CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>;

    using typename Base::Face_handle;

    using typename Base::Point;

    using typename Base::Vertex_handle;

    using typename Base::Locate_type;

    // Constructors

    Custom_Constrained_Delaunay_triangulation_2(const Gt &gt = Gt())

        : Base(gt)
    {
    }

    Custom_Constrained_Delaunay_triangulation_2(typename Base::List_constraints &lc, const Gt &gt = Gt())

        : Base(lc, gt)
    {
    }

    template <class InputIterator>

    Custom_Constrained_Delaunay_triangulation_2(InputIterator it, InputIterator last, const Gt &gt = Gt())

        : Base(it, last, gt)
    {
    }

    // New insert method without flips

    Vertex_handle insert_no_flip(const Point &a, Face_handle start = Face_handle())
    {

        // Call Ctr::insert without flip_around

        Vertex_handle va = this->Base::Ctr::insert(a, start); // Directly call Ctr::insert from the base

        return va;
    }

    // Another insert method with known location

    Vertex_handle insert_no_flip(const Point &a, Locate_type lt, Face_handle loc, int li)
    {

        Vertex_handle va = this->Base::Ctr::insert(a, lt, loc, li); // Directly call Ctr::insert from the base

        return va;
    }

    void remove_no_flip(Vertex_handle v)
    {

        this->Base::Ctr::remove(v);
    }
};

typedef Custom_Constrained_Delaunay_triangulation_2<Kernel, Tds> CDT; // Defines the CDT type after class

template <typename TPoint>
bool isObtuse(const TPoint &p1, const TPoint &p2, const TPoint &p3, const Polygon_2 &regionPolygon)
{
    // Lambda to compute squared distances between two points and convert to double
    auto sq_dist = [](const TPoint &a, const TPoint &b) -> double
    {
        return CGAL::to_double(CGAL::squared_distance(a, b));
    };

    // Calculate the squared distances between the points and convert them to double
    double a2 = sq_dist(p1, p2); // squared distance between p1 and p2
    double b2 = sq_dist(p2, p3); // squared distance between p2 and p3
    double c2 = sq_dist(p3, p1); // squared distance between p3 and p1

    // Return whether the triangle is obtuse based on the squared distances
    double epsilon = 1e-8;
    bool is_obtuse = a2 + b2 < c2 - epsilon || b2 + c2 < a2 - epsilon || c2 + a2 < b2 - epsilon;

    if (!is_obtuse)
    {
        return false;
    }

    TPoint centroid = CGAL::centroid(p1, p2, p3);

    CGAL::Bounded_side location = regionPolygon.bounded_side(centroid);
    bool is_inside_boundary = location != CGAL::ON_UNBOUNDED_SIDE;

    return is_obtuse && is_inside_boundary;
}

// Helper function to check if a point is inside the triangle using orientation tests
template <typename TPoint>
bool isPointInTriangle(const TPoint &p, const TPoint &a, const TPoint &b, const TPoint &c)
{
    // Compute the orientations
    CGAL::Orientation o1 = CGAL::orientation(a, b, p);
    CGAL::Orientation o2 = CGAL::orientation(b, c, p);
    CGAL::Orientation o3 = CGAL::orientation(c, a, p);

    // Point is inside if all orientations are the same (either all counterclockwise or all clockwise)
    return (o1 == o2) && (o2 == o3);
}

// Function to get the centroid of a triangle
template <typename TPoint>
TPoint getCentroid(const TPoint &p1, const TPoint &p2, const TPoint &p3)
{
    return CGAL::centroid(p1, p2, p3);
}

// Function to get the projection of a triangle
template <typename TPoint>
TPoint getProjection(const TPoint &p1, const TPoint &p2, const TPoint &p3)
{
    TPoint obtuseVertex, a, b;

    // Determine which vertex forms the obtuse angle
    if (CGAL::angle(p2, p1, p3) == CGAL::OBTUSE)
    {
        obtuseVertex = p1;
        a = p2;
        b = p3;
    }
    else if (CGAL::angle(p1, p2, p3) == CGAL::OBTUSE)
    {
        obtuseVertex = p2;
        a = p1;
        b = p3;
    }
    else
    {
        obtuseVertex = p3;
        a = p1;
        b = p2;
    }

    // Project the obtuse vertex onto the opposite side (line segment a-b)
    Vector ab = b - a;            // Vector from a to b
    Vector ao = obtuseVertex - a; // Vector from a to obtuseVertex

    // Scalar projection of ao onto ab (dot product based projection)
    auto dotProduct_ao_ab = ao * ab; // Dot product of ao and ab
    auto dotProduct_ab_ab = ab * ab; // Dot product of ab with itself

    double t = CGAL::to_double(dotProduct_ao_ab) / CGAL::to_double(dotProduct_ab_ab);

    // Clamp t between [0, 1] to ensure the projection falls on the segment
    t = std::max(0.0, std::min(1.0, t));

    // Calculate the projection point
    TPoint projection = a + t * ab; // Compute the projected point
    return projection;
}

// Calculates the midpoint of the longest edge of a triangle.
template <typename TPoint>
TPoint getMidpointOfLongestEdge(const TPoint &p1, const TPoint &p2, const TPoint &p3)
{
    double d12 = CGAL::to_double(CGAL::squared_distance(p1, p2)); // distance between p1 and p2
    double d23 = CGAL::to_double(CGAL::squared_distance(p2, p3)); // distance between p2 and p3
    double d31 = CGAL::to_double(CGAL::squared_distance(p3, p1)); // distance between p3 and p1

    // Compare the distances to find the longest edge and return its midpoint
    if (d12 >= d23 && d12 >= d31)
    {
        return CGAL::midpoint(p1, p2); // Longest edge between p1 and p2
    }
    else if (d23 >= d12 && d23 >= d31)
    {
        return CGAL::midpoint(p2, p3); // Longest edge between p2 and p3
    }
    else
    {
        return CGAL::midpoint(p3, p1); // Longest edge between p3 and p1
    }
}

// Counts the number of obtuse triangles in the CDT.
template <typename TCDT>
int countObtuseTriangles(TCDT &cdt, const Polygon_2 &regionPolygon)
{
    int obtuseCount = 0;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();
        if (isObtuse<Point>(p1, p2, p3, regionPolygon))
        {
            obtuseCount++;
        }
    }
    return obtuseCount;
}

void processCluster(CDT &cdt, const std::vector<CDT::Face_handle> &cluster)
{
    std::cerr << "Entering processCluster with " << cluster.size() << " faces.\n";
    std::set<CDT::Edge> sharedEdges;
    int skipEdge = false;

    // Step 1: Identify shared edges in the cluster
    for (const auto &face : cluster)
    {
        if (!cdt.is_infinite(face))
        { // Ensure the face is valid
            for (int i = 0; i < 3; ++i)
            {
                CDT::Edge edge(face, i);
                CDT::Face_handle neighbor = face->neighbor(i);

                // Check if the edge is shared within the cluster
                if (std::find(cluster.begin(), cluster.end(), neighbor) != cluster.end())
                {
                    // Only add non-constrained edges to sharedEdges
                    if (!cdt.is_constrained(edge))
                    {
                        sharedEdges.insert(edge);
                    }
                }
            }
        }
    }
    std::cerr << "Shared edges identified: " << sharedEdges.size() << ".\n";

    // Step 2: Process each shared edge only once
    std::vector<CDT::Edge> unconstrainedEdgesToRemove; // Track unconstrained edges for later removal

    for (const auto &edge : sharedEdges)
    {
        // Retrieve vertices associated with the shared edge
        CDT::Vertex_handle vh1 = edge.first->vertex(cdt.cw(edge.second));
        CDT::Vertex_handle vh2 = edge.first->vertex(cdt.ccw(edge.second));

        // Ensure vertices are valid before proceeding
        if (cdt.is_infinite(vh1) || cdt.is_infinite(vh2))
        {
            std::cerr << "Skipping infinite vertex handle in shared edge processing.\n";
            continue;
        }

        Point p1 = vh1->point();
        Point p2 = vh2->point();
        std::cerr << "Processing shared edge between points " << p1 << " and " << p2 << ".\n";

        // Step 3: Collect unconstrained edges incident to both vertices
        std::vector<CDT::Edge> unconstrainedEdges;
        int maxIncidentEdges = 50; // Limit on incident edges to avoid infinite loops

        auto CollectUnconstrainedEdges = [&](CDT::Vertex_handle vh)
        {
            auto edgeIt = cdt.incident_edges(vh);
            auto start = edgeIt;
            int edgeCount = 0;

            do
            {
                if (!cdt.is_infinite(*edgeIt))
                {
                    Point p1_inc = edgeIt->first->vertex(cdt.cw(edgeIt->second))->point();
                    Point p2_inc = edgeIt->first->vertex(cdt.ccw(edgeIt->second))->point();

                    if (cdt.is_constrained(*edgeIt))
                    {
                        skipEdge = true;
                        return;
                    }
                    else if (!((p1_inc == p1 && p2_inc == p2) || (p1_inc == p2 && p2_inc == p1)))
                    {
                        unconstrainedEdges.push_back(*edgeIt);
                    }
                    else
                    {
                        std::cerr << "Skipping shared edge in incident edge collection.\n";
                    }
                }
                // if (edgeCount++ >= maxIncidentEdges) {
                //     std::cerr << "Max incident edges limit reached for vertex.\n";
                //     break;
                // }
                ++edgeIt;
            } while (edgeIt != start);
        };

        CollectUnconstrainedEdges(vh1);
        CollectUnconstrainedEdges(vh2);
        std::cerr << "Collected " << unconstrainedEdges.size() << " unconstrained edges.\n";

        // Step 4: Insert constraints for unconstrained edges, if needed
        for (const auto &e : unconstrainedEdges)
        {
            Point source = e.first->vertex(cdt.cw(e.second))->point();
            Point target = e.first->vertex(cdt.ccw(e.second))->point();

            if (skipEdge)
            {
                std::cerr << "Skipping edge due to constraint.\n";
                return;
            }

            for (const auto &edge : unconstrainedEdges)
            {
                Point source = edge.first->vertex(cdt.cw(edge.second))->point();
                Point target = edge.first->vertex(cdt.ccw(edge.second))->point();
                cdt.insert_constraint(source, target);
            }

            Point edgeV1 = vh1->point();
            Point edgeV2 = vh2->point();
            cdt.remove(cdt.insert_no_flip(vh1->point()));
            cdt.remove(cdt.insert_no_flip(vh2->point()));

            CDT::Vertex_handle newVh1 = cdt.insert_no_flip(edgeV1);
            CDT::Vertex_handle newVh2 = cdt.insert_no_flip(edgeV2);

            for (const auto &edge : unconstrainedEdges)
            {
                cdt.remove_constraint(edge.first, edge.second);
            }

            return;
        }
    }
}

bool checkForCircumcenter(CDT &cdt, typename CDT::Face_handle face)
{
    // std::cout << "draw1" << std::endl;
    // CGAL::draw(cdt);
    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();
    // std::cout << "triangle:" << p1 << "," << p2 << "," << p3 << std::endl;

    CDT::Vertex_handle v1, v2;
    CDT::Edge exceedingEdge;
    if (CGAL::angle(p2, p1, p3) == CGAL::OBTUSE)
    {
        // p1 is the obtuse vertex; opposite edge is between p2 and p3
        exceedingEdge = CDT::Edge(face, 0);
        v1 = face->vertex(1);
        v2 = face->vertex(2);
    }
    else if (CGAL::angle(p1, p2, p3) == CGAL::OBTUSE)
    {
        // p2 is the obtuse vertex; opposite edge is between p1 and p3
        exceedingEdge = CDT::Edge(face, 1);
        v1 = face->vertex(0);
        v2 = face->vertex(2);
    }
    else if (CGAL::angle(p1, p3, p2) == CGAL::OBTUSE)
    {
        // p3 is the obtuse vertex; opposite edge is between p1 and p2
        exceedingEdge = CDT::Edge(face, 2);
        v1 = face->vertex(0);
        v2 = face->vertex(1);
    }
    else
    {
        return false;
    }

    // Print the exceeding edge for reference
    Point p1_exceeding = exceedingEdge.first->vertex(cdt.cw(exceedingEdge.second))->point();
    Point p2_exceeding = exceedingEdge.first->vertex(cdt.ccw(exceedingEdge.second))->point();
    // std::cout << "Exceeding edge between (" << p1_exceeding.x() << ", " << p1_exceeding.y() << ") and ("
    //           << p2_exceeding.x() << ", " << p2_exceeding.y() << ")" << std::endl;
    // std::cout << "v1:" << v1->point() << "v2:" << v2->point() << std::endl;

    if (cdt.is_constrained(exceedingEdge))
    {
        return false;
    }

    bool skipEdge = false;
    std::vector<CDT::Edge> incidentEdges;

    auto CollectEdges = [&](CDT::Vertex_handle vh)
    {
        auto edgeIt = cdt.incident_edges(vh);
        auto start = edgeIt;

        // std::cout << "Collecting edges for vertex at ("
        //           << vh->point().x() << ", " << vh->point().y() << "):" << std::endl;

        do
        {
            if (!cdt.is_infinite(*edgeIt))
            {
                Point p1 = edgeIt->first->vertex(cdt.cw(edgeIt->second))->point();
                Point p2 = edgeIt->first->vertex(cdt.ccw(edgeIt->second))->point();
                Point excP1 = exceedingEdge.first->vertex(cdt.cw(exceedingEdge.second))->point();
                Point excP2 = exceedingEdge.first->vertex(cdt.ccw(exceedingEdge.second))->point();

                // std::cout << "  Edge between (" << p1.x() << ", " << p1.y() << ") and ("
                //           << p2.x() << ", " << p2.y() << ")";

                if (cdt.is_constrained(*edgeIt))
                {
                    // std::cout << " - Constrained Edge, skipping." << std::endl;
                    skipEdge = true;
                    return;
                }
                else if (!((p1 == excP1 && p2 == excP2) || (p1 == excP2 && p2 == excP1)))
                {
                    // std::cout << " - Unconstrained Edge, adding to list." << std::endl;
                    incidentEdges.push_back(*edgeIt);
                }
                else
                {
                    // std::cout << " - Exceeding edge, not adding." << std::endl;
                }
            }
            ++edgeIt;
        } while (edgeIt != start);
    };

    CollectEdges(v1);
    CollectEdges(v2);

    // std::cout << "Incident edges collected:" << std::endl;
    for (const auto &edge : incidentEdges)
    {
        Point pp1 = edge.first->vertex(cdt.cw(edge.second))->point();
        Point pp2 = edge.first->vertex(cdt.ccw(edge.second))->point();
        // std::cout << "  Edge between (" << pp1.x() << ", " << pp1.y() << ") and ("
        //           << pp2.x() << ", " << pp2.y() << ")" << std::endl;
    }
    // std::cout << "End of incident edges list." << std::endl;

    if (skipEdge)
    {
        return false;
    }

    for (const auto &edge : incidentEdges)
    {
        Point source = edge.first->vertex(cdt.cw(edge.second))->point();
        Point target = edge.first->vertex(cdt.ccw(edge.second))->point();
        cdt.insert_constraint(source, target);
    }
    // std::cout << "draw2" << std::endl;
    // CGAL::draw(cdt);

    Point edgeV1 = v1->point();
    Point edgeV2 = v2->point();
    // std::cout << "edgeV1:" << edgeV1 << "edgeV2:" << edgeV2 << std::endl;
    cdt.remove(cdt.insert_no_flip(v1->point()));
    cdt.remove(cdt.insert_no_flip(v2->point()));
    // cdt.remove(v1);
    // cdt.remove(v2);
    // cdt.remove_no_flip(v1);
    // cdt.remove_no_flip(v2);

    // std::cout << "draw3" << std::endl;
    // CGAL::draw(cdt);

    CDT::Vertex_handle newVh1 = cdt.insert_no_flip(edgeV1);
    CDT::Vertex_handle newVh2 = cdt.insert_no_flip(edgeV2);
    // std::cout << "draw4" << std::endl;
    // CGAL::draw(cdt);

    for (const auto &edge : incidentEdges)
    {
        cdt.remove_constraint(edge.first, edge.second);
        // std::cout << "Removed constraint between (" << edge.first->vertex(cdt.cw(edge.second))->point().x() << ", "
        //           << edge.first->vertex(cdt.cw(edge.second))->point().y() << ") and ("
        //           << edge.first->vertex(cdt.ccw(edge.second))->point().x() << ", "
        //           << edge.first->vertex(cdt.ccw(edge.second))->point().y() << ")" << std::endl;
    }
    // CGAL::draw(cdt);
    return true;
}

// Function to insert a point and count the number of obtuse triangles
template <typename TCDT, typename TPoint>
int tryPointInsertion(TCDT &cdt, const TPoint &test_point, const Polygon_2 &regionPolygon, bool merge = false, bool circ = false, std::optional<std::vector<CDT::Face_handle>> obtuseCluster = std::nullopt, std::optional<CDT::Face_handle> face = std::nullopt)
{
    TCDT temp_cdt = cdt; // Make a copy of the current triangulation

    if (merge)
    {
        processCluster(temp_cdt, *obtuseCluster);
    }
    if (circ)
    {
        // std::cout << "Checking for circumcenter" << std::endl;
        if (checkForCircumcenter(temp_cdt, *face) == false)
        {
            return std::numeric_limits<int>::max();
        }
    }

    temp_cdt.insert_no_flip(test_point); // Insert the test point into the copy

    return countObtuseTriangles<CDT>(temp_cdt, regionPolygon); // Count how many obtuse triangles remain
}

// Function to check if two triangles share an edge
template <typename TFaceHandle>
bool areNeighbours(TFaceHandle face1, TFaceHandle face2)
{
    int sharedVertices = 0;
    for (int i = 0; i < 3; ++i)
    {
        Point p1 = face1->vertex(i)->point();
        for (int j = 0; j < 3; ++j)
        {
            Point p2 = face2->vertex(j)->point();
            if (p1 == p2)
            {
                sharedVertices++;
            }
        }
    }
    return (sharedVertices == 2); // True if two triangles share exactly 2 vertices
}

// Function to collect neighboring obtuse triangles that form a convex hull
template <typename TCDT, typename TFaceHandle>
std::vector<CDT::Face_handle> collectNeighbouringObtuseTriangles(TCDT &cdt, TFaceHandle face, const Polygon_2 &regionPolygon)
{
    std::vector<TFaceHandle> obtuseCluster;
    std::vector<TFaceHandle> facesToCheck;
    std::set<TFaceHandle> visitedFaces;

    facesToCheck.push_back(face);
    visitedFaces.insert(face);

    while (!facesToCheck.empty())
    {
        TFaceHandle currentFace = facesToCheck.front();
        facesToCheck.erase(facesToCheck.begin());

        obtuseCluster.push_back(currentFace);

        // Check neighbors
        for (int i = 0; i < 3; ++i)
        {
            TFaceHandle neighborFace = currentFace->neighbor(i);

            Point p1 = neighborFace->vertex(0)->point();
            Point p2 = neighborFace->vertex(1)->point();
            Point p3 = neighborFace->vertex(2)->point();

            // Check if the neighbor is a valid face, not already visited, and is also obtuse
            if (neighborFace != TFaceHandle() && visitedFaces.find(neighborFace) == visitedFaces.end() && isObtuse<Point>(p1, p2, p3, regionPolygon))
            {
                if (areNeighbours<CDT::Face_handle>(currentFace, neighborFace))
                {
                    facesToCheck.push_back(neighborFace);
                    visitedFaces.insert(neighborFace);
                }
            }
        }
    }

    return obtuseCluster;
}

// Function to find the centroid of a convex polygon from the merged triangles
template <typename TPoint>
TPoint getConvexHullCentroid(const std::vector<TPoint> &convexHullPoints)
{
    double xSum = 0.0, ySum = 0.0;
    for (const auto &point : convexHullPoints)
    {
        xSum += CGAL::to_double(point.x());
        ySum += CGAL::to_double(point.y());
    }

    // Compute the average x and y coordinates
    double xCentroid = xSum / convexHullPoints.size();
    double yCentroid = ySum / convexHullPoints.size();

    TPoint centroid = TPoint(xCentroid, yCentroid);

    // CGAL::Polygon_2<typename TPoint::R> polygon(convexHullPoints.begin(), convexHullPoints.end());

    // // Check if the centroid is inside the convex hull
    // if(CGAL::bounded_side_2(polygon.vertices_begin(), polygon.vertices_end(), centroid, typename TPoint::R()) == CGAL::ON_BOUNDED_SIDE) {
    //     return centroid;
    // }
    // else{
    //     return std::nullopt;
    // }
    return centroid;
}

// Function to try merging obtuse triangles, form the convex hull, and insert a Steiner point
template <typename TCDT, typename TFaceHandle>
std::optional<Point> tryMergingObtuseTriangles(TCDT &cdt, const std::vector<TFaceHandle> &obtuseCluster)
{
    // Make a temporary copy of the triangulation
    TCDT temp_cdt = cdt;

    // Collect points from the obtuse triangles
    std::set<Point> uniquePoints;
    for (const auto &face : obtuseCluster)
    {
        for (int i = 0; i < 3; ++i)
        {
            uniquePoints.insert(face->vertex(i)->point());
        }
    }

    // Form the convex hull of the unique points
    std::vector<Point> convexHullPoints(uniquePoints.begin(), uniquePoints.end());
    Point centroid = getConvexHullCentroid<Point>(convexHullPoints);

    CGAL::Polygon_2<typename Point::R> ConvexHullPolygon(convexHullPoints.begin(), convexHullPoints.end());

    if (CGAL::bounded_side_2(ConvexHullPolygon.vertices_begin(), ConvexHullPolygon.vertices_end(), centroid, typename Point::R()) != CGAL::ON_BOUNDED_SIDE)
    {
        return centroid;
    }
    else
    {
        return std::nullopt;
    }
}

// Function to add Steiner points to remove obtuse triangles
template <typename TCDT, typename TPoint>
void addSteinerPoints(TCDT &cdt, std::vector<TPoint> &steiner_points, const Polygon_2 &regionPolygon)
{
    int numVerticesBefore = cdt.number_of_vertices(); // Get current vertex count
    int steinerPointsNum = 0;
    bool steinerPointInserted = false;

    // Loop until no more Steiner points are inserted
    do
    {
        steinerPointInserted = false; // Reset the flag for Steiner point insertion
        int globalMinObtuseTriangles = std::numeric_limits<int>::max();
        TPoint bestPointToInsert;
        CDT::Face_handle bestTriangle;
        bool edgeRemoval = false;
        bool mergeEdgeRemoval = false;
        std::vector<CDT::Face_handle> obtuseCluster;

        // Iterate over all triangles
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
        {
            TPoint p1 = face->vertex(0)->point();
            TPoint p2 = face->vertex(1)->point();
            TPoint p3 = face->vertex(2)->point();

            // Check if the triangle is obtuse
            if (isObtuse<Point>(p1, p2, p3, regionPolygon))
            {
                // Try circumcenter method
                TPoint circumcenter = CGAL::circumcenter(p1, p2, p3);
                int obtuseAfterCircumcenter = std::numeric_limits<int>::max();
                // if (regionPolygon.bounded_side(circumcenter) != CGAL::ON_UNBOUNDED_SIDE)
                // {
                //     int obtuseAfterCircumcenter = tryPointInsertion<CDT, Point>(cdt, circumcenter, regionPolygon, false, true, std::nullopt, face);
                // }

                // Try midpoint of the longest edge
                TPoint midpoint = getMidpointOfLongestEdge<Point>(p1, p2, p3);
                int obtuseAfterMidpoint = tryPointInsertion<CDT, Point>(cdt, midpoint, regionPolygon);

                // Try centroid of the triangle
                TPoint centroid = getCentroid<Point>(p1, p2, p3);
                int obtuseAfterCentroid = tryPointInsertion<CDT, Point>(cdt, centroid, regionPolygon);

                // Try projection of the triange
                TPoint projection = getProjection<Point>(p1, p2, p3);
                int obtuseAfterProjection = tryPointInsertion<CDT, Point>(cdt, projection, regionPolygon);

                // Try merging obtuse triangles
                obtuseCluster = collectNeighbouringObtuseTriangles<CDT, CDT::Face_handle>(cdt, face, regionPolygon);
                int obtuseAfterMerge = std::numeric_limits<int>::max(); // Initialize to maximum value

                TPoint mergeCentroid;
                if (obtuseCluster.size() > 1)
                {
                    std::optional<Point> mergeCentroid = tryMergingObtuseTriangles<CDT, CDT::Face_handle>(cdt, obtuseCluster);
                    if (mergeCentroid)
                    {
                        obtuseAfterMerge = tryPointInsertion<CDT, Point>(cdt, *mergeCentroid, regionPolygon, true, false, obtuseCluster, face);
                    }
                    else
                    {
                        obtuseAfterMerge = std::numeric_limits<int>::max();
                    }
                }

                // Find the method that reduces the number of obtuse triangles the most
                int minObtuseTrianglesForThisFace = std::min({obtuseAfterCircumcenter, obtuseAfterMidpoint, obtuseAfterCentroid, obtuseAfterProjection, obtuseAfterMerge});

                // Update the global minimum if needed
                if (minObtuseTrianglesForThisFace < globalMinObtuseTriangles)
                {
                    globalMinObtuseTriangles = minObtuseTrianglesForThisFace;

                    // Determine the best point to insert based on the method
                    if (minObtuseTrianglesForThisFace == obtuseAfterProjection)
                    {
                        bestPointToInsert = projection;
                    }
                    else if (minObtuseTrianglesForThisFace == obtuseAfterCircumcenter)
                    {
                        bestPointToInsert = circumcenter;
                    }
                    else if (minObtuseTrianglesForThisFace == obtuseAfterMidpoint)
                    {
                        bestPointToInsert = midpoint;
                    }
                    else if (minObtuseTrianglesForThisFace == obtuseAfterCentroid)
                    {
                        bestPointToInsert = centroid;
                    }
                    else if (minObtuseTrianglesForThisFace == obtuseAfterMerge)
                    {
                        bestPointToInsert = mergeCentroid;
                    }

                    bestTriangle = face;
                    edgeRemoval = (minObtuseTrianglesForThisFace == obtuseAfterCircumcenter);
                    mergeEdgeRemoval = (minObtuseTrianglesForThisFace == obtuseAfterMerge);
                }
            }
        }
        // Apply the best insertion if it reduces obtuse triangles
        if (globalMinObtuseTriangles <= countObtuseTriangles<CDT>(cdt, regionPolygon))
        {

            if (edgeRemoval)
            {
                checkForCircumcenter(cdt, bestTriangle);
            }
            if (mergeEdgeRemoval)
            {
                processCluster(cdt, obtuseCluster);
            }
            // Insert the best point
            steiner_points.push_back(bestPointToInsert);
            CDT::Vertex_handle steiner_vh = cdt.insert_no_flip(bestPointToInsert);
            steiner_vh->info() = numVerticesBefore + steiner_points.size() - 1;
            steinerPointsNum++;          // Increment the number of Steiner points added
            steinerPointInserted = true; // Mark that we inserted a point
        }
    } while (steinerPointInserted);

    std::cout << "Steiner points added: " << steinerPointsNum << std::endl;
}

template <typename TCDT, typename TFaceHandle, typename T>
bool isFlipableForNonObtuse(TCDT &cdt, TFaceHandle face, T index, const Polygon_2 &regionPolygon)
{
    // Get the neighbor face
    TFaceHandle neighbor_face = face->neighbor(index);

    // Check if the edge is infinite (on the convex hull or boundary)
    if (cdt.is_infinite(face) || cdt.is_infinite(neighbor_face))
    {
        return false; // Don't flip edges near the boundary or convex hull
    }

    // Get the vertices of the two triangles that share the edge
    auto vertex1 = face->vertex((index + 1) % 3)->point();
    auto vertex2 = face->vertex((index + 2) % 3)->point();
    auto vertex3 = face->vertex(index)->point();

    auto n_vertex1 = neighbor_face->vertex((neighbor_face->index(face) + 1) % 3)->point();
    auto n_vertex2 = neighbor_face->vertex((neighbor_face->index(face) + 2) % 3)->point();
    auto n_vertex3 = neighbor_face->vertex(neighbor_face->index(face))->point();

    // Check if either of the triangles before the flip is obtuse
    bool isObtuseBefore = isObtuse<Point>(vertex1, vertex2, vertex3, regionPolygon) || isObtuse<Point>(n_vertex1, n_vertex2, n_vertex3, regionPolygon);

    // Simulate the flip to get the new edge
    auto new_edge_start = face->vertex((index + 1) % 3)->point();
    auto new_edge_end = neighbor_face->vertex((neighbor_face->index(face) + 2) % 3)->point();

    // Check if the new edge intersects any other edges in the triangulation
    for (auto e = cdt.edges_begin(); e != cdt.edges_end(); ++e)
    {
        // Get the edge vertices
        auto edge_vertex1 = e->first->vertex((e->second + 1) % 3)->point();
        auto edge_vertex2 = e->first->vertex((e->second + 2) % 3)->point();

        // Skip checking the current edge
        if ((edge_vertex1 == vertex1 && edge_vertex2 == vertex2) || (edge_vertex1 == vertex2 && edge_vertex2 == vertex1))
        {
            continue;
        }

        // Check if the new edge would intersect an existing edge
        if (CGAL::do_intersect(CGAL::Segment_2<Kernel>(new_edge_start, new_edge_end), CGAL::Segment_2<Kernel>(edge_vertex1, edge_vertex2)))
        {
            return false; // The new edge would intersect an existing edge, so don't flip
        }
    }

    // Check if either of the triangles after the flip would be obtuse
    auto new_vertex1 = vertex1;
    auto new_vertex2 = vertex3;
    auto new_vertex3 = n_vertex2;

    auto n_new_vertex1 = n_vertex1;
    auto n_new_vertex2 = vertex3;
    auto n_new_vertex3 = vertex2;

    bool isObtuseAfter = isObtuse<Point>(new_vertex1, new_vertex2, new_vertex3, regionPolygon) || isObtuse<Point>(n_new_vertex1, n_new_vertex2, n_new_vertex3, regionPolygon);

    // Flip the edge only if it reduces the number of obtuse triangles
    return (isObtuseBefore && !isObtuseAfter);
}

// Function to try flipping edges to remove obtuse triangles
template <typename TCDT>
void flipEdgesToMakeNonObtuse(TCDT &cdt, const Polygon_2 &regionPolygon)
{
    for (auto edge = cdt.edges_begin(); edge != cdt.edges_end(); ++edge)
    {
        // Get the face handle and the index of the edge
        CDT::Face_handle face = edge->first;
        int index = edge->second;

        // Check if the edge is flippable
        if (isFlipableForNonObtuse<CDT, CDT::Face_handle, int>(cdt, face, index, regionPolygon))
        {
            // Get the triangles before the flip
            Point p1 = face->vertex((index + 1) % 3)->point();
            Point p2 = face->vertex((index + 2) % 3)->point();
            Point p3 = face->vertex(index)->point();

            CDT::Face_handle neighbor_face = face->neighbor(index);
            int neighbor_index = neighbor_face->index(face);
            Point np1 = neighbor_face->vertex((neighbor_index + 1) % 3)->point();
            Point np2 = neighbor_face->vertex((neighbor_index + 2) % 3)->point();
            Point np3 = neighbor_face->vertex(neighbor_index)->point();

            bool isObtuseBefore = isObtuse<Point>(p1, p2, p3, regionPolygon) || isObtuse<Point>(np1, np2, np3, regionPolygon);

            // Perform the edge flip
            cdt.flip(face, index);

            // Get the triangles after the flip
            Point p1_new = face->vertex((index + 1) % 3)->point();
            Point p2_new = face->vertex((index + 2) % 3)->point();
            Point p3_new = face->vertex(index)->point();

            Point np1_new = neighbor_face->vertex((neighbor_index + 1) % 3)->point();
            Point np2_new = neighbor_face->vertex((neighbor_index + 2) % 3)->point();
            Point np3_new = neighbor_face->vertex(neighbor_index)->point();

            bool isObtuseAfter = isObtuse<Point>(p1_new, p2_new, p3_new, regionPolygon) || isObtuse<Point>(np1_new, np2_new, np3_new, regionPolygon);

            // Revert the flip if it didn't reduce obtuseness
            if (isObtuseAfter && !isObtuseBefore)
            {
                std::cout << "Reverting flip for edge: " << index << std::endl;
                cdt.flip(face, index); // Revert the flip if it doesn't help
            }
            else
            {
                std::cout << "Flip successful for edge: " << index << std::endl;
            }
        }
    }
}

// Primary template definition for general types (including `double`)
#include <regex>
template <typename T>
std::string toFraction(const T &value)
{
    // Convert the value to double using CGAL::to_double
    double doubleValue = CGAL::to_double(value);

    double intPart;
    double fracPart = modf(doubleValue, &intPart);

    if (fracPart == 0.0)
    {
        return std::to_string(static_cast<int>(intPart)); // Integer number
    }

    int denominator = 10000; // Large denominator for accuracy
    int numerator = static_cast<int>(round(fracPart * denominator));

    std::string json_string_array = std::to_string(static_cast<int>(intPart * denominator + numerator)) + "/" + std::to_string(denominator);
    json_string_array = std::regex_replace(json_string_array, std::regex(R"(\\/)"), "/");
    return json_string_array;
}

#endif // TRIANGULATION_HPP