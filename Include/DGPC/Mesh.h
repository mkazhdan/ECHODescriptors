/************************************************************
 * This file is part of the DGPC library. The library computes
 * Discrete Geodesic Polar Coordinates on a polygonal mesh.
 *
 * More info:
 *   http://folk.uio.no/eivindlm/dgpc/
 *
 * Authors: Eivind Lyche Melvær and Martin Reimers
 * Centre of Mathematics and Department of Informatics
 * University of Oslo, Norway, 2012
 ************************************************************/
#ifndef DGPC_MESH_H
#define DGPC_MESH_H

#include <limits>
#include <fstream>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Misha/Geometry.h>
#include "Vector3.h"

namespace DGPC {

  template<class P>
  struct OpenMeshTraits : OpenMesh::DefaultTraits {
    typedef P Point;
  };

  template<class Point>
    class MeshOM : public OpenMesh::PolyMesh_ArrayKernelT< OpenMeshTraits<Point> > {
  public:
    typedef Point point_type;

    bool openOBJ(const char* filename) {

      std::ifstream in(filename);
      std::string line;

      int maxpoly = 0;

      while ( std::getline(in, line) ) {

        std::istringstream i(line.c_str ());

        if(line[0] == 'v' && isblank(line[1])) {
          char t;
          typename point_type::value_type v0, v1, v2;
          i >> t;
          i >> v0;
          i >> v1;
          i >> v2;

          this->add_vertex(Point(v0, v1, v2));
        } else if(line[0] == 'f') {
          char t;
          i >> t;
 
          std::vector< typename MeshOM<Point>::VertexHandle > face;

          while (!i.eof()) {
            int v;
            i >> v;	    
            v -= 1; //OBJ is 1-based, we are 0-based
            face.push_back( this->vertex_handle(v) );

	    char next = i.peek();	  
	    while(!isspace(next) && !i.eof()) {
	      //Ignore normals and texture coordinates
	      i >> next;
	      next = i.peek();
	      if(next == '\r') i >> next; //Raise eof. Fix for obj files with '\r\n' line endings
	    }
	    
          }
          this->add_face(face);

          int poly = face.size();
          if(poly > maxpoly) maxpoly = poly;

        }
      }
      in.close();

      const int nv = this->n_vertices();
      const int nf = this->n_faces();
      printf("%s: %d vertices and %d faces (most complex face has %d nodes)\n", filename, nv, nf, maxpoly);

      return true;
    }

      void fromTriangles( const std::vector< Point3D<float> >& vertices, const std::vector< TriangleIndex >& triangles) 
      {


         for (int l = 0; l < vertices.size (); l++)
         {
            typename point_type::value_type v0 = vertices[l][0], v1 = vertices[l][1], v2 = vertices[l][2];
 
            this->add_vertex(Point(v0, v1, v2));
         }

         for (int l = 0; l < triangles.size (); l++)
         {
            std::vector<  typename MeshOM<Point>::VertexHandle > face;

            for (int j = 0; j < 3; j++)
            {
               face.push_back ( this->vertex_handle(triangles[l][j]) );
            }

            this->add_face (face);
         }


      }




      void fromTriangles( const std::vector< Point3D<double> >& vertices, const std::vector< TriangleIndex >& triangles) 
      {


         for (int l = 0; l < vertices.size (); l++)
         {
            typename point_type::value_type v0 = vertices[l][0], v1 = vertices[l][1], v2 = vertices[l][2];
 
            this->add_vertex(Point(v0, v1, v2));
         }

         for (int l = 0; l < triangles.size (); l++)
         {
            std::vector<  typename MeshOM<Point>::VertexHandle > face;

            for (int j = 0; j < 3; j++)
            {
               face.push_back ( this->vertex_handle(triangles[l][j]) );
            }

            this->add_face (face);
         }


      }


  };

}

#endif
