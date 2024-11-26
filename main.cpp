#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
// #include <CGAL/Polygon_2.h>
#include "includes/triangulation.hpp" // Keep this for the CDT definition and for helper functions

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = CDT::Point; // Use the Point type from the CDT defined in triangulation.hpp
namespace pt = boost::property_tree;
// typedef CGAL::Polygon_2<Kernel> Polygon_2;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Please provide the input JSON file as a command-line argument." << std::endl;
        return 1;
    }

    std::ifstream inputFile(argv[1]);
    if (!inputFile.is_open())
    {
        std::cerr << "Failed to open the input file." << std::endl;
        return 1;
    }

    pt::ptree inputData;

    // Read the input JSON file
    try
    {
        pt::read_json(inputFile, inputData);
    }
    catch (const pt::json_parser_error &e)
    {
        std::cerr << "Error parsing JSON file: " << e.what() << std::endl;
        return 1;
    }

    std::string instanceUid = inputData.get<std::string>("instance_uid");
    int numPoints = inputData.get<int>("num_points");

    std::vector<int> pointsX, pointsY, regionBoundary;
    for (const auto &item : inputData.get_child("points_x"))
    {
        pointsX.push_back(item.second.get_value<int>());
    }
    for (const auto &item : inputData.get_child("points_y"))
    {
        pointsY.push_back(item.second.get_value<int>());
    }
    for (const auto &item : inputData.get_child("region_boundary"))
    {
        regionBoundary.push_back(item.second.get_value<int>());
    }

    std::vector<std::vector<int>> additionalConstraints;
    for (const auto &constraint : inputData.get_child("additional_constraints"))
    {
        std::vector<int> constraintPair;
        for (const auto &point : constraint.second)
        {
            constraintPair.push_back(point.second.get_value<int>());
        }
        additionalConstraints.push_back(constraintPair);
    }

    Polygon_2 regionPolygon;
    for (int idx : regionBoundary)
    {
        regionPolygon.push_back(Point(pointsX[idx], pointsY[idx]));
    }

    // std::cout <<"Region boundary:\n";
    // for(auto it = regionPolygon.vertices_begin(); it != regionPolygon.vertices_end(); ++it) {
    //     std::cout << it->x() << " " << it->y() << std::endl;
    // }

    // Use the CDT defined in triangulation.hpp
    CDT cdt;

    // Insert points into the triangulation and assign index info
    for (int i = 0; i < numPoints; ++i)
    {
        CDT::Vertex_handle vh = cdt.insert(Point(pointsX[i], pointsY[i]));
        vh->info() = i; // Assign index to the vertex
    }

    // Insert constraints for the region boundary
    for (size_t i = 0; i < regionBoundary.size(); ++i)
    {
        int nextIndex = (i + 1) % regionBoundary.size();
        cdt.insert_constraint(Point(pointsX[regionBoundary[i]], pointsY[regionBoundary[i]]),
                              Point(pointsX[regionBoundary[nextIndex]], pointsY[regionBoundary[nextIndex]]));
    }

    // Insert additional constraints
    for (const auto &constraint : additionalConstraints)
    {
        cdt.insert_constraint(Point(pointsX[constraint[0]], pointsY[constraint[0]]),
                              Point(pointsX[constraint[1]], pointsY[constraint[1]]));
    }
    std::cout << "Triangulation done, adding steiner points and flipping edges..." << std::endl;
    std::cout << "Number of all triangles at start: " << cdt.number_of_faces() << std::endl;
    std::cout << "Number of obtuse triangles at start: " << countObtuseTriangles<CDT>(cdt, regionPolygon) << std::endl;

    std::vector<Point> steiner_points;

    // Add Steiner points for obtuse triangles
    addSteinerPoints<CDT, Point>(cdt, steiner_points, regionPolygon);

    // Try to flip edges to make non-obtuse triangles
    flipEdgesToMakeNonObtuse<CDT>(cdt, regionPolygon);

    std::cout << ">>FINAL RESULTS<<" << std::endl;
    // Count the number of triangles in the triangulation
    int allTriangles = cdt.number_of_faces();
    std::cout << "All triangles: " << allTriangles << std::endl;

    // Count the number of obtuse triangles remaining
    int obtuseCount = countObtuseTriangles<CDT>(cdt, regionPolygon);
    std::cout << "Obtuse triangles remaining: " << obtuseCount << std::endl;

    std::cout << "Outputting results..." << std::endl;

    // Prepare output JSON using Boost PropertyTree
    pt::ptree outputData;
    outputData.put("content_type", "CG_SHOP_2025_Solution");
    outputData.put("instance_uid", instanceUid);

    // Store steiner points
    pt::ptree steinerPointsXNode, steinerPointsYNode;
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex)
    {
        if (vertex->point().x() != pointsX[vertex->info()] || vertex->point().y() != pointsY[vertex->info()])
        {
            pt::ptree xNode, yNode;
            xNode.put("", toFraction(vertex->point().x())); // Use fraction format
            yNode.put("", toFraction(vertex->point().y())); // Use fraction format
            steinerPointsXNode.push_back(std::make_pair("", xNode));
            steinerPointsYNode.push_back(std::make_pair("", yNode));
        }
    }
    outputData.add_child("steiner_points_x", steinerPointsXNode);
    outputData.add_child("steiner_points_y", steinerPointsYNode);

    // Store edges
    pt::ptree edgesNode;
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge)
    {
        int index1 = edge->first->vertex((edge->second + 1) % 3)->info();
        int index2 = edge->first->vertex((edge->second + 2) % 3)->info();
        pt::ptree edgeNode;
        edgeNode.push_back(std::make_pair("", pt::ptree(std::to_string(index1))));
        edgeNode.push_back(std::make_pair("", pt::ptree(std::to_string(index2))));
        edgesNode.push_back(std::make_pair("", edgeNode));
    }
    outputData.add_child("edges", edgesNode);

    // Write output to JSON file
    std::ofstream outputFile("../output/output.json");
    if (outputFile.is_open())
    {
        pt::write_json(outputFile, outputData);
        outputFile.close();
    }
    else
    {
        std::cerr << "Failed to open output file." << std::endl;
        return 1;
    }

    // Draw the triangulation using CGAL's draw function
    CGAL::draw(cdt);

    return 0;
}
