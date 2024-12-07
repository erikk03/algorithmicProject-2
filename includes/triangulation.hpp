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
typedef CGAL::Segment_2<Kernel> Segment;

using Vector = Kernel::Vector_2;

////////////////////////////////////////////
//// CUSTOM CONSTRAINED DELAUNAY CLASS /////
////////////////////////////////////////////

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

///////////////////////////////////////////
//// HELPER FUNCTIONS FOR TRIANGULATION ///
///////////////////////////////////////////

// Helper function to check if a triangle is obtuse
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

////////////////////////
// For MergeTriangles //
////////////////////////

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

    return centroid;
}

// Function to try merging obtuse triangles, form the convex hull, and insert a Steiner point
template <typename TCDT, typename TFaceHandle>
Point tryMergingObtuseTriangles(TCDT &cdt, const std::vector<TFaceHandle> &obtuseCluster)
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

    return centroid;
}

// Function to process a cluster of triangles and refine the triangulation
bool processCluster(CDT &cdt, const std::vector<CDT::Face_handle> &cluster)
{
    // std::cerr << "Processing cluster with " << cluster.size() << " obtuse faces.\n";

    // Step 1: Collect shared edges within the cluster
    std::set<CDT::Edge> sharedEdges;
    for (const auto &face : cluster)
    {
        if (!cdt.is_infinite(face))
        {
            for (int i = 0; i < 3; ++i)
            {
                CDT::Edge edge(face, i);
                CDT::Face_handle neighbor = face->neighbor(i);

                // If the neighbor is also in the cluster, add the edge to the shared edges
                if (std::find(cluster.begin(), cluster.end(), neighbor) != cluster.end())
                {
                    if (cdt.is_constrained(edge))
                    {
                        // std::cerr << "Shared edge is constrained. Aborting.\n";
                        return false; // Return false if the shared edge is constrained
                    }
                    sharedEdges.insert(edge);
                }
            }
        }
    }

    // Step 2: Collect unique points from the cluster
    std::vector<Point> uniquePoints;
    for (const auto &face : cluster)
    {
        for (int i = 0; i < 3; ++i)
        {
            CDT::Vertex_handle vh = face->vertex(i);
            if (!cdt.is_infinite(vh))
            {
                uniquePoints.push_back(vh->point());
            }
        }
    }

    // Erase duplicates and sort the points
    std::sort(uniquePoints.begin(), uniquePoints.end());
    uniquePoints.erase(std::unique(uniquePoints.begin(), uniquePoints.end()), uniquePoints.end());

    // std::cerr << "Collected " << uniquePoints.size() << " unique points.\n";

    // Step 3: Check if the points form a convex hull
    std::vector<Point> convexHullPoints(uniquePoints.begin(), uniquePoints.end());
    CGAL::Polygon_2<typename Point::R> ConvexHullPolygon(convexHullPoints.begin(), convexHullPoints.end());

    if (!ConvexHullPolygon.is_simple() || !ConvexHullPolygon.is_convex())
    {
        // std::cerr << "The polygon is not convex. Aborting.\n";
        return false; // Return false if the polygon is not convex
    }

    // Centroid of the convex hull
    Point centroid = getConvexHullCentroid<Point>(convexHullPoints);

    // Check if centroid is inside the convex hull
    if (CGAL::bounded_side_2(ConvexHullPolygon.vertices_begin(), ConvexHullPolygon.vertices_end(),
                             centroid, typename Point::R()) != CGAL::ON_BOUNDED_SIDE)
    {
        // std::cerr << "Centroid is outside the convex polygon. Aborting.\n";
        return false;
    }

    // Step 4: Remove non constraint edges
    for (const auto &edge : sharedEdges)
    {
        if (!cdt.is_constrained(edge))
        {
            cdt.remove_constraint(edge.first, edge.second);
        }
    }
    // std::cerr << "Removed non-constrained edges within the cluster.\n";

    // Step 5: Insert the steiner point
    CDT::Vertex_handle steinerHandle = cdt.insert(centroid);
    if (steinerHandle == nullptr)
    {
        // std::cerr << " Failed to insert Steiner point. Aborting.\n";
        return false;
    }
    // std::cerr << "Inserted Steiner point at " << centroid << ".\n";

    // Step6: Insert constraints again for the convex hull edges
    for (size_t i = 0; i < convexHullPoints.size(); ++i)
    {
        Point p1 = convexHullPoints[i];
        Point p2 = convexHullPoints[(i + 1) % convexHullPoints.size()];
        cdt.insert_constraint(p1, p2);
    }
    // std::cerr << "Inserted constraints for the convex hull edges.\n";

    // std::cerr << "Cluster processed successfully.\n";
    return true;
}

//////////////////////
// For Circumcenter //
//////////////////////

bool checkForCircumcenter(CDT &cdt, typename CDT::Face_handle face, const Polygon_2 &regionPolygon)
{

    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();

    CDT::Edge exceedingEdge;
    Point obtuseVertex;
    Point other1, other2;

    // Identify obtuse vertex and the opposite edge
    if (CGAL::angle(p2, p1, p3) == CGAL::OBTUSE)
    {
        // p1 is the obtuse vertex; opposite edge is between p2 and p3
        exceedingEdge = CDT::Edge(face, 0);
        obtuseVertex = p1;
        other1 = p2;
        other2 = p3;
    }
    else if (CGAL::angle(p1, p2, p3) == CGAL::OBTUSE)
    {
        // p2 is the obtuse vertex; opposite edge is between p1 and p3
        exceedingEdge = CDT::Edge(face, 1);
        obtuseVertex = p2;
        other1 = p1;
        other2 = p3;
    }
    else if (CGAL::angle(p1, p3, p2) == CGAL::OBTUSE)
    {
        // p3 is the obtuse vertex; opposite edge is between p1 and p2
        exceedingEdge = CDT::Edge(face, 2);
        obtuseVertex = p3;
        other1 = p1;
        other2 = p2;
    }
    else
    {
        // std::cerr << "No obtuse angle found in the triangle.\n";
        return false; // No obtuse angle, skip
    }

    // Compute circumcenter
    Point circumcenter = CGAL::circumcenter(p1, p2, p3);

    // Step 1: Check if circumcenter lies within the region boundary
    if (!regionPolygon.bounded_side(circumcenter) == CGAL::ON_BOUNDED_SIDE)
    {
        // std::cerr << "Circumcenter lies outside the region boundary.\n";
        return false;
    }

    // Step 2: Create a segment from obtuseVertex to circumcenter
    Segment segmentToCircumcenter(obtuseVertex, circumcenter);

    // Check intersections with other segments
    std::vector<CDT::Edge> intersectingEdges;
    for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit)
    {
        Point p1 = eit->first->vertex(cdt.cw(eit->second))->point();
        Point p2 = eit->first->vertex(cdt.ccw(eit->second))->point();
        Segment edgeSegment(p1, p2);
        if (CGAL::do_intersect(segmentToCircumcenter, edgeSegment))
        {
            // Calculate the intersection point
            auto result = CGAL::intersection(segmentToCircumcenter, edgeSegment);
            if (const Point *intersectionPoint = boost::get<Point>(&*result))
            {
                // Check if the intersection point is not the same as the obtuse vertex
                if (*intersectionPoint != obtuseVertex)
                {
                    intersectingEdges.push_back(*eit);
                }
            }
        }
    }

    // Step 3: If segment intersects more than one edge or if it intersects a constrained edge, cancel
    if (intersectingEdges.size() > 1)
    {
        // std::cerr << "Circumcenter intersects more than one edge.\n";
        return false;
    }
    else if (intersectingEdges.size() == 1 && cdt.is_constrained(intersectingEdges[0]))
    {
        // std::cerr << "Circumcenter intersects a constrained edge.\n";
        return false;
    }

    // Track original constraints
    std::set<std::pair<Point, Point>> originalConstraints;

    for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit)
    {
        if (cdt.is_constrained(*eit))
        {
            Point p1 = eit->first->vertex(cdt.cw(eit->second))->point();
            Point p2 = eit->first->vertex(cdt.ccw(eit->second))->point();
            originalConstraints.insert({std::min(p1, p2), std::max(p1, p2)});
        }
    }

    // Step 4: Get vertices of the intersecting edge
    CDT::Vertex_handle v1 = intersectingEdges[0].first->vertex(cdt.cw(intersectingEdges[0].second));
    CDT::Vertex_handle v2 = intersectingEdges[0].first->vertex(cdt.ccw(intersectingEdges[0].second));

    // Step 5: Check if these vertices touch constrained edges, remove constraints temporarily
    std::vector<CDT::Edge> edgesToRestore;
    auto RemoveEdgeConstraints = [&](CDT::Vertex_handle vh)
    {
        auto edgeIt = cdt.incident_edges(vh);
        auto start = edgeIt;
        int i = 0;
        do
        {
            if (!cdt.is_infinite(*edgeIt) && cdt.is_constrained(*edgeIt))
            {
                edgesToRestore.push_back(*edgeIt);
                cdt.remove_constraint(edgeIt->first, edgeIt->second);
            }
            ++edgeIt;
            i++;
            if (i > 100)
            {
                break;
            }
        } while (edgeIt != start);
    };

    RemoveEdgeConstraints(v1);
    RemoveEdgeConstraints(v2);

    // Step 6: Remove the vertices themselves, storing any other edges removed
    std::vector<Segment> unconstrainedSegments;
    Point v1Point = v1->point();
    Point v2Point = v2->point();

    auto CollectUnconstrainedEdges = [&](CDT::Vertex_handle vh)
    {
        auto edgeIt = cdt.incident_edges(vh);
        auto start = edgeIt;

        do
        {
            if (!cdt.is_infinite(*edgeIt) && !cdt.is_constrained(*edgeIt))
            {
                Segment edgeSegment = cdt.segment(*edgeIt);
                unconstrainedSegments.push_back(edgeSegment);
            }
            ++edgeIt;
        } while (edgeIt != start);
    };

    CollectUnconstrainedEdges(v1);
    CollectUnconstrainedEdges(v2);

    cdt.remove(v1);
    cdt.remove(v2);

    // Step 7: Restore the edges
    for (const auto &seg : unconstrainedSegments)
    {
        if (!(seg.source() == exceedingEdge.first->vertex(cdt.cw(exceedingEdge.second))->point() &&
              seg.target() == exceedingEdge.first->vertex(cdt.ccw(exceedingEdge.second))->point()) &&
            !(seg.source() == exceedingEdge.first->vertex(cdt.ccw(exceedingEdge.second))->point() &&
              seg.target() == exceedingEdge.first->vertex(cdt.cw(exceedingEdge.second))->point()))
        {
            cdt.insert_constraint(seg.source(), seg.target());
        }
    }

    // Step 8: Insert the circumcenter
    CDT::Vertex_handle circumcenterHandle = cdt.insert(circumcenter);

    // Step 9: Remove edges that were not constrained originally
    for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit)
    {
        if (cdt.is_constrained(*eit))
        {
            // Get the vertices of the edge
            Point p1 = eit->first->vertex(cdt.cw(eit->second))->point();
            Point p2 = eit->first->vertex(cdt.ccw(eit->second))->point();

            // Create a pair for easy comparison
            std::pair<Point, Point> edgePair = {std::min(p1, p2), std::max(p1, p2)};

            // If the edge is constrained but wasn't part of the original constraints, remove it
            if (originalConstraints.find(edgePair) == originalConstraints.end())
            {
                // We need to get the face and the edge index for the constrained edge
                CDT::Face_handle f = eit->first;
                int index = eit->second;

                cdt.remove_constraint(f, index); // Remove the constraint using the correct parameters
            }
        }
    }

    return true;
}

//////////////////////////////////
// For Temporary Steiner Points //
//////////////////////////////////

// Function to insert a point and count the number of obtuse triangles
template <typename TCDT, typename TPoint>
int tryPointInsertion(TCDT &cdt, const TPoint &test_point, const Polygon_2 &regionPolygon, bool merge = false, bool circ = false, std::optional<std::vector<CDT::Face_handle>> obtuseCluster = std::nullopt, std::optional<CDT::Face_handle> face = std::nullopt)
{
    TCDT temp_cdt = cdt; // Make a copy of the current triangulation

    if (merge)
    {
        if (processCluster(temp_cdt, *obtuseCluster) == false)
        {
            return std::numeric_limits<int>::max();
        }
    }
    if (circ)
    {
        if (checkForCircumcenter(temp_cdt, *face, regionPolygon) == false)
        {
            return std::numeric_limits<int>::max();
        }
    }

    temp_cdt.insert_no_flip(test_point); // Insert the test point into the copy

    return countObtuseTriangles<CDT>(temp_cdt, regionPolygon); // Count how many obtuse triangles remain
}

///////////////////////////////////////////////////
//// THREE METHODS TO ADD STEINER POINTS //////////
///////////////////////////////////////////////////

// Method 1: Insert Steiner points using local search optimization
template <typename TCDT, typename TPoint>
void localSearchOptimization(TCDT &cdt, std::vector<TPoint> &steiner_points, const Polygon_2 &regionPolygon, int L)
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
                if (regionPolygon.bounded_side(circumcenter) != CGAL::ON_UNBOUNDED_SIDE)
                {
                    int obtuseAfterCircumcenter = tryPointInsertion<CDT, Point>(cdt, circumcenter, regionPolygon, false, true, std::nullopt, face);
                }

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
                int obtuseAfterMerge = std::numeric_limits<int>::max(); // Initialize to maximum value
                obtuseCluster = collectNeighbouringObtuseTriangles<CDT, CDT::Face_handle>(cdt, face, regionPolygon);
                TPoint mergeCentroid = tryMergingObtuseTriangles<CDT, CDT::Face_handle>(cdt, obtuseCluster);

                if (obtuseCluster.size() > 1)
                {
                    obtuseAfterMerge = tryPointInsertion<CDT, Point>(cdt, mergeCentroid, regionPolygon, true, false, obtuseCluster, face);
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
                checkForCircumcenter(cdt, bestTriangle, regionPolygon);
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

// Method 2: Insert Steiner points using simulated annealing optimization
template <typename TCDT, typename TPoint>
void simulatedAnnealingOptimization(TCDT &cdt, std::vector<TPoint> &steiner_points, const Polygon_2 &regionPolygon, double alpha, double beta, int L)
{
    // Initialize temperature and compute initial energy
    double T = 1.0; // Initial temperature
    int obtuseTriangles = countObtuseTriangles(cdt, regionPolygon);
    int steinerPoints = steiner_points.size();
    double currentEnergy = alpha * obtuseTriangles + beta * steinerPoints;
    double R = 0.05; // Acceptance probability threshold

    std::cout << "Initial Energy: " << currentEnergy << std::endl;

    int iteration = 0;
    // Main simulated annealing loop
    while (T > 0 && iteration < L)
    {
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
        {
            TPoint p1 = face->vertex(0)->point();
            TPoint p2 = face->vertex(1)->point();
            TPoint p3 = face->vertex(2)->point();

            // Check if the triangle is obtuse
            if (isObtuse<Point>(p1, p2, p3, regionPolygon))
            {
                // Randomly select a method (1 to 5)
                int method = 1 + rand() % 5;
                int newObtuseTriangles = std::numeric_limits<int>::max();
                TPoint selectedPoint;
                std::vector<CDT::Face_handle> obtuseCluster;
                bool processOriginal = false;

                // Make a temporary CDT copy to evaluate the new configuration
                // TCDT temp_cdt = cdt;

                switch (method)
                {
                case 1: // Circumcenter
                {
                    // TPoint circumcenter = CGAL::circumcenter(p1, p2, p3);
                    // if (regionPolygon.bounded_side(circumcenter) != CGAL::ON_UNBOUNDED_SIDE)
                    // {
                    //         newObtuseTriangles = tryPointInsertion<CDT, Point>(cdt, circumcenter, regionPolygon, false, true, std::nullopt, face);
                    //         selectedPoint = circumcenter;
                    //         processOriginal = true;
                    // }
                    // break;
                }
                case 2: // Midpoint
                {
                    TPoint midpoint = getMidpointOfLongestEdge<Point>(p1, p2, p3);
                    newObtuseTriangles = tryPointInsertion<CDT, Point>(cdt, midpoint, regionPolygon);
                    selectedPoint = midpoint;
                    break;
                }
                case 3: // Centroid
                {
                    TPoint centroid = getCentroid<Point>(p1, p2, p3);
                    newObtuseTriangles = tryPointInsertion<CDT, Point>(cdt, centroid, regionPolygon);
                    selectedPoint = centroid;
                    break;
                }
                case 4: // Projection
                {
                    TPoint projection = getProjection<Point>(p1, p2, p3);
                    newObtuseTriangles = tryPointInsertion<CDT, Point>(cdt, projection, regionPolygon);
                    selectedPoint = projection;
                    break;
                }
                case 5: // Merge
                {
                    // obtuseCluster = collectNeighbouringObtuseTriangles<CDT, CDT::Face_handle>(temp_cdt, face, regionPolygon);
                    // if (obtuseCluster.size() > 1)
                    // {
                    //     auto mergeResult = tryMergingObtuseTriangles<CDT, CDT::Face_handle>(temp_cdt, obtuseCluster);
                    //     if (mergeResult.has_value())
                    //     {
                    //         selectedPoint = *mergeResult;
                    //         newObtuseTriangles = tryPointInsertion<CDT, Point>(temp_cdt, selectedPoint, regionPolygon, true, false, obtuseCluster, face);
                    //         processOriginal = true;
                    //     }
                    // }
                    // break;
                }
                default:
                    break;
                }

                // Calculate the energy change
                int newSteinerPoints = steiner_points.size() + 1; // Assume one more Steiner point
                double newEnergy = alpha * newObtuseTriangles + beta * newSteinerPoints;
                double deltaE = newEnergy - currentEnergy;

                // Accept or reject the new configuration
                if (deltaE < 0 || exp(-deltaE / T) >= R)
                {
                    // Apply to original CDT
                    // if (method == 1) // Circumcenter
                    // {
                    //     checkForCircumcenter(cdt, face, regionPolygon);
                    // }
                    // if (method == 5) // Merge
                    // {
                    //     processCluster(cdt, obtuseCluster);
                    // }

                    cdt.insert_no_flip(selectedPoint);
                    steiner_points.push_back(selectedPoint);
                    currentEnergy = newEnergy;

                    std::cout << "Accepted new configuration using method " << method
                              << " with Energy: " << currentEnergy << std::endl;
                }
                else
                {
                    std::cout << "Rejected new configuration using method " << method
                              << " with Energy: " << newEnergy << std::endl;
                }
            }
        }

        // Decrease the temperature
        T -= 1.0 / L;
        iteration++;
        std::cout << "Temperature reduced to: " << T << std::endl;
    }

    std::cout << "Final Energy: " << currentEnergy << std::endl;
    std::cout << "Total Steiner Points: " << steiner_points.size() << std::endl;
}

template <typename TPoint>
double calculateRadiusToHeight(const TPoint &p1, const TPoint &p2, const TPoint &p3)
{
    // Step 1: Calculate circumradius (R) using the circumcenter
    TPoint circumcenter = CGAL::circumcenter(p1, p2, p3);
    double R = std::sqrt(CGAL::to_double(CGAL::squared_distance(circumcenter, p1)));

    // Step 2: Find the longest side (base) and calculate the height (h)
    double d1 = CGAL::to_double(CGAL::squared_distance(p1, p2));
    double d2 = CGAL::to_double(CGAL::squared_distance(p2, p3));
    double d3 = CGAL::to_double(CGAL::squared_distance(p3, p1));

    // Find the longest side
    double maxDistance = std::max({d1, d2, d3});
    double baseLength = std::sqrt(maxDistance);

    // Determine the height corresponding to the longest side
    double area = std::abs(CGAL::to_double(CGAL::area(p1, p2, p3)));
    double height = (2 * area) / baseLength;

    // Step 3: Compute and return the Radius-to-Height ratio (ρ)
    if (height == 0)
    {
        throw std::runtime_error("Degenerate triangle with zero height.");
    }
    return R / height;
}

// Method 3: Insert Steiner points using ant colony optimization

template <typename TCDT, typename TPoint>
void antColonyOptimization(TCDT &cdt, std::vector<TPoint> &steiner_points, const Polygon_2 &regionPolygon, double alpha, double beta, double xi, double psi, double lambda, int kappa, int L)
{
    int n = cdt.number_of_vertices(); // Number of input points
    int K = std::max(1, n / 4);       // Number of ants (at least n/4)

    // Initialize pheromone values for all methods (1 to 4)
    std::map<int, double> pheromoneTrails = {
        {1, 1.0}, // Projection
        {2, 1.0}, // Circumcenter
        {3, 1.0}, // Midpoint
        {4, 1.0}  // Merge
    };

    TCDT bestTriangulation = cdt; // Store the best overall triangulation
    double bestEnergy = alpha * countObtuseTriangles(cdt, regionPolygon) + beta * steiner_points.size();

    // Main ACO cycle
    for (int cycle = 1; cycle <= L; ++cycle)
    {
        std::cout << "Cycle " << cycle << "/" << L << std::endl;

        // Store temporary triangulations and energies for all ants
        std::vector<TCDT> antTriangulations(K, cdt);
        std::vector<double> antEnergies(K, bestEnergy);

        // Track pheromone reinforcements for each method
        std::map<int, double> pheromoneReinforcement = {{1, 0.0}, {2, 0.0}, {3, 0.0}, {4, 0.0}};

        for (int k = 0; k < K; ++k)
        {
            TCDT &antTriangulation = antTriangulations[k];
            double antEnergy = antEnergies[k];
            bool processedTriangle = false; // Track if the ant has processed a triangle

            for (auto face = antTriangulation.finite_faces_begin(); face != antTriangulation.finite_faces_end(); ++face)
            {
                TPoint p1 = face->vertex(0)->point();
                TPoint p2 = face->vertex(1)->point();
                TPoint p3 = face->vertex(2)->point();

                if (isObtuse<TPoint>(p1, p2, p3, regionPolygon))
                {
                    double rho = calculateRadiusToHeight(p1, p2, p3);

                    // Compute heuristic values ηsp
                    double etaProjection = std::max(0.0, (rho - 1) / rho);
                    double etaCircumcenter = rho / (2 + rho);
                    double etaMidpoint = std::max(0.0, (3 - 2 * rho) / 3);
                    double etaMerge = 1.0;

                    // Calculate probabilities Psp(k)
                    std::vector<double> probabilities = {
                        pheromoneTrails[1] * std::pow(etaProjection, psi),
                        pheromoneTrails[2] * std::pow(etaCircumcenter, psi),
                        pheromoneTrails[3] * std::pow(etaMidpoint, psi),
                        pheromoneTrails[4] * std::pow(etaMerge, psi)};

                    // Normalize probabilities
                    double sumProbabilities = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
                    for (auto &prob : probabilities)
                    {
                        prob /= sumProbabilities;
                    }

                    // Select a method based on probabilities
                    double randomValue = ((double)rand() / RAND_MAX);
                    int chosenMethod = std::distance(probabilities.begin(), std::lower_bound(probabilities.begin(), probabilities.end(), randomValue));

                    TPoint selectedPoint;
                    switch (chosenMethod)
                    {
                    case 0: // Projection
                        selectedPoint = getProjection<TPoint>(p1, p2, p3);
                        break;
                    case 1: // Circumcenter
                        // selectedPoint = CGAL::circumcenter(p1, p2, p3);
                        selectedPoint = getProjection<TPoint>(p1, p2, p3);
                        break;
                    case 2: // Midpoint
                        selectedPoint = getMidpointOfLongestEdge<TPoint>(p1, p2, p3);
                        break;
                    case 3: // Merge
                    {
                        selectedPoint = getProjection<TPoint>(p1, p2, p3);
                        // std::vector<CDT::Face_handle> obtuseCluster = collectNeighbouringObtuseTriangles<TCDT, CDT::Face_handle>(antTriangulation, face, regionPolygon);
                        // if (obtuseCluster.size() > 1)
                        // {
                        //     auto mergeResult = tryMergingObtuseTriangles<TCDT, CDT::Face_handle>(antTriangulation, obtuseCluster);
                        //     if (mergeResult.has_value())
                        //     {
                        //         selectedPoint = *mergeResult;
                        //     }
                        // }
                        break;
                    }
                    default:
                        continue;
                    }

                    // if (!selectedPoint.is_valid())
                    // {
                    //     continue; // Skip invalid selections
                    // }

                    // Insert the Steiner point and update energy
                    antTriangulation.insert(selectedPoint);
                    steiner_points.push_back(selectedPoint);

                    int newObtuseTriangles = countObtuseTriangles(antTriangulation, regionPolygon);
                    antEnergy = alpha * newObtuseTriangles + beta * steiner_points.size();

                    if (antEnergy < antEnergies[k])
                    {
                        antEnergies[k] = antEnergy;

                        // Pheromone reinforcement for this method
                        pheromoneReinforcement[chosenMethod + 1] += 1.0 / (1.0 + alpha * newObtuseTriangles + beta * steiner_points.size());
                    }

                    // processedTriangle = true; // Mark this triangle as processed
                    // break;                   // Exit loop after processing one triangle
                }
            }

            // // Skip to the next ant if no triangle was processed
            // if (!processedTriangle)
            // {
            //     continue;
            // }
        }

        // Update pheromone trails after processing all ants
        for (auto &pheromone : pheromoneTrails)
        {
            int method = pheromone.first;
            pheromone.second = (1 - lambda) * pheromone.second + pheromoneReinforcement[method];
        }

        // Resolve conflicts and save the best triangulation of this cycle
        for (int k = 0; k < K; ++k)
        {
            if (antEnergies[k] < bestEnergy)
            {
                bestEnergy = antEnergies[k];
                bestTriangulation = antTriangulations[k];
            }
        }

        std::cout << "Best energy for cycle " << cycle << ": " << bestEnergy << std::endl;
    }

    cdt = bestTriangulation; // Update the original CDT with the best triangulation
    std::cout << "Final best energy: " << bestEnergy << std::endl;
}

////////////////////////////////////////
//////////// FLIP FUNCTIONS ////////////
////////////////////////////////////////

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
