// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_SURFACE_MESH_CUSTOM_CRITERIA_3_H
#define CGAL_SURFACE_MESH_CUSTOM_CRITERIA_3_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Surface_mesher/Standard_criteria.h>
#include <iostream>

namespace CGAL {

template <class Tr>
class Surface_mesh_custom_criteria_3
{
  typedef Surface_mesher::Refine_criterion<Tr> Criterion;
  typedef Surface_mesher::Standard_criteria<Criterion> Criteria;

public:
  typedef Tr Triangulation;
  typedef typename Tr::Geom_traits::FT FT;

  typedef typename Criteria::Quality Quality;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Point Point;
  
  Surface_mesh_custom_criteria_3(const FT angle_bound,
				  const FT radius_bound1,
				  const FT distance_bound1,
                  const FT radius_bound2,
				  const FT distance_bound2, double d_min)
    : curvature_size_criterion1(distance_bound1),
      uniform_size_criterion1(radius_bound1),
      curvature_size_criterion2(distance_bound2),
      uniform_size_criterion2(radius_bound2),
      aspect_ratio_criterion(angle_bound),
      Dmin(d_min)
      
  {
    criterion_vector1.reserve(4);
    
    criterion_vector1.push_back (&aspect_ratio_criterion);
    criterion_vector1.push_back (&uniform_size_criterion1);
    criterion_vector1.push_back (&curvature_size_criterion1);
    
    criteria1.set_criteria(criterion_vector1);

    criterion_vector2.reserve(4);
    
    criterion_vector2.push_back (&aspect_ratio_criterion);
    criterion_vector2.push_back (&uniform_size_criterion2);
    criterion_vector2.push_back (&curvature_size_criterion2);
    criteria2.set_criteria(criterion_vector2);
  }

  bool is_bad (const Facet& f, Quality& q) const
  {
    Point p1 = f.first->vertex ((f.second+1)&3)->point();
    Point p2 = f.first->vertex ((f.second+2)&3)->point();
    Point p3 = f.first->vertex ((f.second+3)&3)->point();
    double dtemp = sqrt((p1.x()+p2.x()+p3.x())*(p1.x()+p2.x()+p3.x())+(p1.y()+p2.y()+p3.y())*(p1.y()+p2.y()+p3.y())+(p1.z()+p2.z()+p3.z())*(p1.z()+p2.z()+p3.z()))/3;
    if (dtemp<Dmin)
        return criteria1.is_bad(f, q);
    else
        return criteria2.is_bad(f, q);
  }
private:
  Surface_mesher::Curvature_size_criterion<Tr> curvature_size_criterion1;
  Surface_mesher::Curvature_size_criterion<Tr> curvature_size_criterion2;
  // bound on Hausdorff distance does not play any role if bigger than
  // the square of the Uniform_size_criterion

  Surface_mesher::Uniform_size_criterion<Tr> uniform_size_criterion1;
  Surface_mesher::Uniform_size_criterion<Tr> uniform_size_criterion2;
  // bound on radii of surface Delaunay balls
  
  Surface_mesher::Aspect_ratio_criterion<Tr> aspect_ratio_criterion;
  // lower bound on minimum angle in degrees

  std::vector<Criterion*> criterion_vector1,criterion_vector2;
  Criteria criteria1, criteria2;
  double Dmin;
}; // end class Surface_mesh_custom_criteria_3

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_CUSTOME_CRITERIA_3_H
