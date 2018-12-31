
#include <fstream>
#include <iomanip>
#include <iostream>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <functional>
#include "cgalDistMesh.h"
#include "meshutils.h"
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Implicit_mesh_domain_3.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::Vector_3 Vector_3;
typedef GT::FT FT;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef FT (Function)(const K::Point_3&);
//typedef CGAL::Implicit_mesh_domain_3<Function,K> Mesh_domain;

typedef FT (*dist_fx_point)(Point_3);
//typedef std::function<FT (Point_3)> dist_fx_point; // doesn't compile on ab-initio (g++ version 4.7)

typedef CGAL::Implicit_surface_3<GT, dist_fx_point> Surface_3;

dist_fx df;

// convert my dist_fx to the proper signature as required by CGAL
FT returnDist(Point_3 p) {
    return df(p.x(), p.y(), p.z());
}

//FT spherical_sizing_field(const Point_3 &p, const int, const Mesh_domain::Index&) {
    //FT sq_d_to_origin = CGAL::squared_distance(p, Point_3(CGAL::ORIGIN));
    //return 0.001 + CGAL::sqrt(sq_d_to_origin);
//}

// maxLength = max edge length
// numPoints = number of random points to find the max. radius
void cgalDistMesh_single(dist_fx dist_function, double bound, mesh_inputs *mi, mesh_outputs *mo) {
    df = dist_function;

    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
    Surface_3 surface(returnDist, Sphere_3(CGAL::ORIGIN, bound)); // define surface 
    double maxEdgeLength = mi->maxEdgeLength;
    
    // inputs are angular bound, radius bound, then distance bound
    //CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., spherical_sizing_field, spherical_sizing_field); 
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., 3.*maxEdgeLength, maxEdgeLength); 
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

    std::ofstream out(mi->offFile);
    out << std::setprecision(16);
    CGAL::output_surface_facets_to_off (out, c2t3, 0);

    unsigned numPanels = c2t3.number_of_facets();
    double *panelAreas = new double[numPanels];
    double *centroids = new double[3*numPanels];

    GT::Orientation_3 orientation;
    C2t3::Facet_iterator fi;
    C2t3::Cell_handle cell;
    Vector_3 a, b, c, p0, p1, p2, pc; // pc = centroid
    Point_3 zero(CGAL::ORIGIN);
    int index;
    double sign;
    unsigned cnt=0;
    double surfArea = 0;
    double volume = 0;
    for(fi = c2t3.facets_begin(); fi != c2t3.facets_end(); ++fi, ++cnt) {
        cell = fi->first;
        index = fi->second;
        p0 = Vector_3(zero, cell->vertex(tr.vertex_triple_index(index, 0))->point());
        p1 = Vector_3(zero, cell->vertex(tr.vertex_triple_index(index, 1))->point());
        p2 = Vector_3(zero, cell->vertex(tr.vertex_triple_index(index, 2))->point());
        a = p1 - p0;
        b = p2 - p0;
        c = cross_product(a,b);
        panelAreas[cnt] = 0.5 * sqrt( c.squared_length() );
        surfArea += panelAreas[cnt];

        pc = (p0+p1+p2)/3.;
        sign = 1;
        if(c*pc < 0)
            sign = -1;
        volume += sign * (1./6.) * ( p0 * cross_product(p1,p2) );

        centroids[3*cnt] = pc[0];
        centroids[3*cnt+1] = pc[1];
        centroids[3*cnt+2] = pc[2];
    }

    if(mi->verbose>0) {
        std::cout << "numPanels: " << numPanels << std::endl;
        std::cout << "surfaceArea: " << surfArea << std::endl;
        std::cout << "volume: " << volume << std::endl;
    }

    if(mi->convertToGMSH)
        ConvertOFFtoGMSH(mi->offFile, mi->gmshFile, mi->verbose);

    mo->numPanels = numPanels;
    mo->panelAreas = panelAreas;
    mo->centroids = centroids;
    mo->surfArea = surfArea;
    mo->volume = volume;
}

void cgalDistMesh(dist_fx dist_function, double bound, mesh_inputs *mi, mesh_outputs *mo) {
    unsigned numPanels;
    double avgNumPanels = 0.5 * (mi->minNumPanels + mi->maxNumPanels);
    //std::cout << "min number of panels: " << mi->minNumPanels << std::endl;
    //std::cout << "max number of panels: " << mi->maxNumPanels << std::endl;
    //std::cout << "ideal number of panels: " << avgNumPanels << std::endl;
    mesh_inputs mi_copy = {mi->minNumPanels, mi->maxNumPanels, mi->verbose, mi->maxEdgeLength, mi->convertToGMSH, mi->offFile, mi->gmshFile};

    unsigned maxRetry = 20;
    unsigned cnt = 0;
    while(cnt++ < maxRetry) {
        std::cout << "maxEdgeLength: " << mi_copy.maxEdgeLength << std::endl;
        cgalDistMesh_single(dist_function, bound, &mi_copy, mo);
        if ((mo->numPanels > mi->maxNumPanels) || (mo->numPanels < mi->minNumPanels))
            mi_copy.maxEdgeLength *= sqrt(mo->numPanels / avgNumPanels);
        else
            break;
    }
}
