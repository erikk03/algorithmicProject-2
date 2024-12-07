#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "includes/triangulation.hpp" // Keep this for the CDT definition and for helper functions

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = CDT::Point; // Use the Point type from the CDT defined in triangulation.hpp
namespace pt = boost::property_tree;

int main(int argc, char *argv[])
{
    if (argc < 5 || std::string(argv[1]) != "-i" || std::string(argv[3]) != "-o")
    {
        std::cerr << "Usage: ./opt_triangulation -i /path/to/input.json -o /path/to/output.json" << std::endl;
        return 1;
    }

    std::string inputPath = argv[2];
    std::string outputPath = argv[4];

    std::ifstream inputFile(inputPath);
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

    // Parse the input data
    std::string instanceUid = inputData.get<std::string>("instance_uid");
    int numPoints = inputData.get<int>("num_points");
    bool delaunay = inputData.get<bool>("delaunay", true);
    std::string method = inputData.get<std::string>("method", "ant");

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

    // Parse algorithm-specific parameters
    pt::ptree parameters;
    if (auto params = inputData.get_child_optional("parameters"))
    {
        parameters = *params;
    }
    int L = parameters.get<int>("L", 100);

    // Create the region boundary polygon
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

    if (delaunay)
    {
        std::cout << "Ensuring Delaunay triangulation..." << std::endl;
        // cdt.make_delaunay();
    }

    std::cout << "Triangulation done. Starting optimization with method: " << method << "..." << std::endl;
    std::cout << "Number of all triangles at start: " << cdt.number_of_faces() << std::endl;
    std::cout << "Number of obtuse triangles at start: " << countObtuseTriangles<CDT>(cdt, regionPolygon) << std::endl;

    std::vector<Point> steiner_points;
    if (method == "local")
    {
        localSearchOptimization<CDT>(cdt, steiner_points, regionPolygon, L);
    }
    else if (method == "sa")
    {
        // double alpha = parameters.get<double>("alpha", 1.0);
        // double beta = parameters.get<double>("beta", 1.0);
        double alpha = 5.0;
        double beta = 0.8;
        int L = 300;

        simulatedAnnealingOptimization<CDT>(cdt, steiner_points, regionPolygon, alpha, beta, L);
    }
    else if (method == "ant")
    {
        // double alpha = parameters.get<double>("alpha", 1.0);
        // double beta = parameters.get<double>("beta", 1.0);
        // double xi = parameters.get<double>("xi", 1.0);
        // double psi = parameters.get<double>("psi", 1.0);
        // double lambda = parameters.get<double>("lambda", 0.5);
        // int kappa = parameters.get<int>("kappa", 10);
        double alpha = 4.0;
        double beta = 0.5;
        int xi = 1;
        int psi = 4;
        double lambda = 0.5;
        int kappa = 10;
        int L = 50;


        antColonyOptimization<CDT>(cdt, steiner_points, regionPolygon, alpha, beta, xi, psi, lambda, kappa, L);
    }
    else
    {
        std::cerr << "Invalid method specified: " << method << std::endl;
        return 1;
    }

    //////////// Terminal Results ////////////
    // Count the number of triangles in the triangulation
    int allTriangles = cdt.number_of_faces();
    std::cout << "Number of all triangles at the end: " << allTriangles << std::endl;

    // Count the number of obtuse triangles remaining
    int obtuseCount = countObtuseTriangles<CDT>(cdt, regionPolygon);
    std::cout << "Number of obtuse triangles remaining: " << obtuseCount << std::endl;
    //////////////////////////////////////////

    // Output results
    std::cout << "Optimization completed. Preparing results..." << std::endl;

    pt::ptree outputData;
    outputData.put("content_type", "CG_SHOP_2025_Solution");
    outputData.put("instance_uid", instanceUid);

    // Store steiner points
    pt::ptree steinerPointsXNode, steinerPointsYNode;
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex)
    {
        if (cdt.is_infinite(vertex))
        {
            continue; // Skip infinite vertices
        }

        if (vertex->info() < 0 )
        {
            std::cerr << "Invalid vertex info: " << vertex->info() << std::endl;
            continue; // Skip vertices with invalid or uninitialized info
        }

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

    outputData.put("obtuse_count", obtuseCount);
    outputData.put("method", method);
    outputData.add_child("parameters", parameters);

    // Write output to JSON file
    std::ofstream outputFile(outputPath);
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
