/*
Copyright (c) 2020, Michael Kazhdan and Thomas Mitchel
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef TRIMESH_INCLUDED
#define TRIMESH_INCLUDED

#include <vector>
#include <omp.h>
#include <limits.h>
#include <Misha/Geometry.h>
#include <Misha/FEM.h>
#include <Eigen/Eigen>
#include <Misha/Spectrum.h>


#ifndef M_PI
#define M_PI 3.1415926535897932384
#endif

// === Geodesics ===
#include "DGPC/Generator.h"

typedef DGPC::Vector3< double > PointOM;
typedef DGPC::MeshOM<PointOM> repOM;
typedef DGPC::Generator<repOM> DGPCgenerator;


// This class represents a set of regular samples of a real-valued, periodic, function in 3D
template< class Real=float >
class TriMesh
{
protected:
	const std::vector< Point3D< Real > > &_vertices;
	const std::vector< TriangleIndex > &_triangles;

	// Computed Quantities
	std::vector< Point3D< Real > > _vertexNormals , _triangleNormals;

	std::vector< Real > _vertexAreas , _triangleAreas;

	std::vector< std::vector< TriangleIndex > > _vertexStars;
	std::vector< std::vector< int > > _starTriangles;
	std::vector< std::vector< int > > _starVertices;

   std::vector< Point3D< double > > _triangleMetrics; // g = [a, b; b, c], (a, b, c)

	Real _maxEL , _minEL , _meanEL;

	// === Geodesics ===
	repOM _surfaceMesh;
	double _logEpsilon = 1.0e-12;

	Spectrum< double > _spectrum;

	// Initalizes class from triangle mesh
	void _init( int threads=1 );

public:

	TriMesh( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles , int threads=1 );

	//  Utility clamp function
	Real clamp( Real v , Real lo , Real hi );

	// Compute closest point on triangle to target point
	static Point3D< Real > ClosestPointOnTri( Point3D< Real > P , Point3D< Real > v1 , Point3D< Real > v2 , Point3D< Real > v3 );

	// Get barycentric coordinates for point on mesh 
	std::pair< int , Point3D< double > > getBarycentric( Point3D< Real > P ) const;
	Point3D< double > getBarycentric( Point3D< Real > P , Point3D< Real > v1 , Point3D< Real > v2 , Point3D< Real > v3 ) const; 
	Point3D< double > getBarycentric( Point2D< double > P, Point2D< double > v1 , Point2D< double > v2 , Point2D< double > v3 ) const ;

   
	static Point3D< double > ComposeCoordinates( Point3D< double > P , Point3D< double > v1 , Point3D< double > v2 , Point3D< double > v3 ); // P - barycoords for triangle with verts v1, v2, v3, vi - vertices passed in terms of barycentric coords

	Point3D<Real> asPoint ( std::pair<int, Point3D< double >> baryC ) const;


	// Project point onto triangle plane
	Point3D< float > projectToTri ( Point3D<Real> P , int triIndex) const;

	// Return vertex and triangle vectors respectively
	const std::vector< Point3D<Real> > &vertices( void ) const { return _vertices; }
	const std::vector< TriangleIndex> &triangles( void ) const { return _triangles; }

	// Returns indexed normal
	Point3D< Real >   vertexNormal( int l ) const;
	Point3D< Real > triangleNormal( int l ) const;

	// Returns indexed area
	Real   vertexArea( int l ) const;
	Real triangleArea( int l ) const;
	Real    totalArea( void )  const;

	double triangleArea ( const Point3D< double >& c1, const Point3D< double >& c2, const Point3D< double >& c3 ) const;

	double triangleArea ( const Point2D< double >& c1, const Point2D< double >& c2, const Point2D< double >& c3) const;

	double subTriangleArea ( int l, const Point3D< double >& c1, const Point3D< double >& c2, const Point3D< double >& c3) const;

	// Returns indexed vertex star
	const std::vector< TriangleIndex > &vertexStar( int l ) const;
	const std::vector< int > &vertexStarList ( int l ) const;

	// Returns edge length info
	Real meanEdgeLength( void ) const;
	Real  maxEdgeLength( void ) const;
	Real  minEdgeLength( void ) const;

   // Riemannian Metrics (Euclidean distance)
   void initEuclideanMetrics( void );

   // Gets value of gradient of implicit function in tangent space at triangles w.r.t. to defined metric
    Point2D< double > getMetricGradient( int l, const double phi1, const double phi2, const double phi3) const;
	Point2D< double > getMetricGradient( int l, const std::vector< double >& Implicit) const;

	// Computes gradient of implicit function at all triangles
   void metricGradient( const std::vector< double >& Implciit , std::vector<Point2D< double > > &triGrads , int threads=1) const;

   // Metric dot
   double metricDot( int l, const Point2D< double >& w1, const Point2D< double >& w2 ) const;

   double metricSquareNorm( int l , const Point2D< double >& w ) const;

   Point2D< double > metricRotate90 ( int l, const Point2D< double >& w) const;


	// Curvature
	void computeFundamentalMatricesAndShapeOperator( int l, SquareMatrix<Real, 2>& firstForm, SquareMatrix<Real, 2>& secondForm, SquareMatrix<Real, 2>& shapeOperator) const;

	// Computes Gaussian Curvature at vertices or triangles
	void vertexGaussCurvature ( std::vector< double >& vertGauss, int threads = 1) const;
	void triangleGaussCurvature ( std::vector< double >& triGauss, int threads = 1) const;

	// Computes mean curvature at vertices or triangles
	void vertexMeanCurvature (std::vector< double >& vertMean, int threads = 1) const;
	void triangleMeanCurvature (std::vector< double >& triMean, int threads = 1) const;

	// Smooth vertex signal
	void returnMassStiffness( Real diffTime , SparseMatrix< double , int > &mass , SparseMatrix< double , int >& stiffness );
	void smoothVertexSignal( std::vector< double >& Implicit , Real diffTime );


	// =================
	// === Geodesics ===
	// =================

	void setGeodesicEpsilon( double epsilon );

	void initGeodesicCalc( void ); // Initalizes mesh for geodesic calculations

	std::vector< double > computeGeodesicsAbout( int nodeIndex, float rho ) const;

	std::vector< double > computeGeodesicsAbout( std::pair< int , Point3D< double > > P , float rho ) const;

	std::pair<std::vector< double >, std::vector< double >> computeLogarithmAbout( int nodeIndex , float rho ) const;

	std::pair<std::vector< double >, std::vector< double >> computeLogarithmAbout( std::pair<int, Point3D< double >> P, float rho ) const;

	// ===============================================
	// === Laplace-Beltrami spectral decomposition ===
	// ===============================================

	// Compute spectral decomposition of mesh
	void setSpectralDecomposition( int specDimension=200 , Real specOffset=(Real)100. , bool specLump=false );
	void readSpectralDecomposition( std::string fileName );
	const Spectrum< double > &spectrum( void ) const { return _spectrum; }

	// Heat Kernel Signature
	void vertexHKS( std::vector< double > &vertHKS , double diffTime=0.1 , int threads=1 ) const;

	// Biharmonic distance
	double biharmonicDist( int nodeX , int nodeY ) const;
	double biharmonicDist( std::pair< int , Point3D< double > > P , int node ) const;

	std::vector< double > computeBiharmonicsAbout( int nodeIndex, double rho = std::numeric_limits< double >::max () ) const;

	std::vector< double > computeBiharmonicsAbout ( std::pair<int, Point3D< double >> & P, double rho = std::numeric_limits< double >::max () ) const;

   double getAreaBiharmonic (int l) const;
   void initAreaBiharmonic (int threads = 1);

   void initMetricsBiharmonic (int threads = 1);


	// Diffusion distance
	double diffusionDist( int nodeX , int nodeY , double t0=0.2 ) const;
	double diffusionDist( std::pair< int , Point3D< double > > P , int node , double t0=0.2 ) const;

	std::vector< double > computeDiffusionsAbout( int nodeIndex , double t0=0.2 , double rho=std::numeric_limits< double >::max() ) const;
	std::vector< double > computeDiffusionsAbout( std::pair< int , Point3D< double > > &P , double t0=0.2 , double rho=std::numeric_limits< double >::max () ) const;

	void initMetricsDiffusion (double t = 0.1, int threads = 1);

   double getAreaDiffusion (int l, double t = 0.1) const;
   void initAreaDiffusion (double t = 0.1, int threads = 1);
};

#include "TriMesh.inl"
#endif // TRIMESH_INCLUDED
