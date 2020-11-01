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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <Misha/Solver.h>


template< class Real >
TriMesh< Real >::TriMesh( const std::vector< Point3D< Real > > &vertices , const std::vector< TriangleIndex > &triangles ) : _minEL(0) , _maxEL(0) , _meanEL(0) , _vertices(vertices) , _triangles(triangles)
{  
    _init();
}

template< class Real >
void TriMesh< Real >::_init( void )
{
    // Clear all
    _vertexNormals.resize (0); _triangleNormals.resize (0);
    _vertexAreas.resize (0); _triangleAreas.resize (0);
    _vertexStars.resize (0); _starTriangles.resize (0);

    // Compute normals and areas
    _vertexAreas.resize( _vertices.size() , 0 );
    _triangleAreas.resize( _triangles.size() );

    _vertexNormals.resize( _vertices.size() , Point3D< Real >() );
    _triangleNormals.resize( _triangles.size() );

#pragma omp parallel for
    for( int l=0 ; l<_triangles.size() ; l++ )
    {
        Point3D< Real > v1 = _vertices[ _triangles[l][1] ] - _vertices[ _triangles[l][0] ];
        Point3D< Real > v2 = _vertices[ _triangles[l][2] ] - _vertices[ _triangles[l][0] ];

        _triangleNormals[l] = Point3D<Real>::CrossProduct (v1, v2);

        _triangleAreas[l] = Real (0.5) * std::sqrt( _triangleNormals[l].squareNorm() );

        if( _triangleAreas[l] < ( 16.0 * FLT_EPSILON ) ) _triangleNormals[l] = Point3D< Real >();
        else _triangleNormals[l] /= std::sqrt (_triangleNormals[l].squareNorm ());

    }

    _minEL = std::numeric_limits< Real >::max();
    int eC = 0;
    for( int l=0 ; l<_triangles.size() ; l++ )
    {
        std::vector< Real > EL { std::sqrt( ( _vertices[ _triangles[l][1] ] - _vertices[ _triangles[l][0]] ).squareNorm () ) , std::sqrt( ( _vertices[ _triangles[l][2] ] - _vertices[ _triangles[l][0]]).squareNorm () ) , std::sqrt( ( _vertices[ _triangles[l][2] ] - _vertices[ _triangles[l][1] ] ).squareNorm () ) };

        for (int i = 0; i < 3; i++)
        {
            _vertexNormals[ _triangles[l][i] ] += _triangleNormals[l] * _triangleAreas[l];
            _vertexAreas[ _triangles[l][i] ] += (Real)( _triangleAreas[l] / 3.0 );

            if (EL[i] < _minEL) { _minEL = EL[i]; }
            if (EL[i] > _maxEL) { _maxEL = EL[i]; }
            _meanEL += EL[i];
            eC ++;
        }
    }

    _meanEL /= eC;

#pragma omp parallel for
    for( int l=0 ; l<_vertices.size() ; l++ )
    {
        if( std::sqrt( _vertexNormals[l].squareNorm () )< ( 16.0 * FLT_EPSILON ) ) _vertexNormals[l] = Point3D< Real >();
        else _vertexNormals[l] /= std::sqrt ( _vertexNormals[l].squareNorm () );
    }


    // Compute stars
    _vertexStars.resize( _vertices.size() );
    _starTriangles.resize ( _vertices.size() );

    for( int l=0 ; l<_triangles.size() ; l++ )
    {
        for (int i = 0; i < 3; i++)
        {
            _vertexStars[ _triangles[l][i] ].push_back( _triangles[l] );
            _starTriangles[ _triangles[l][i] ].push_back (l);
        }
    }

    _starVertices.resize( _vertices.size() );

    for( int l=0 ; l<_vertices.size() ; l++ )
        for( int i=0 ; i<_starTriangles[l].size() ; i++ )
            for( int j=0 ; j<3 ; j++ )
                if( _triangles[ _starTriangles[l][i] ][j]!=l )
                    if ( std::find ( _starVertices[l].begin() , _starVertices[l].end() , _triangles[ _starTriangles[l][i] ][j] )==_starVertices[l].end() )
                        _starVertices[l].push_back( _triangles[ _starTriangles[l][i] ][j] );

    // Init Euclidean metric
    initEuclideanMetrics ();

}

// Adapted from:
// https://www.gamedev.net/forums/topic/552906-closest-point-on-triangle/
template< class Real >
Real TriMesh< Real >::clamp( Real v, Real lo, Real hi )
{
    assert( !(hi < lo) );
    return (v < lo) ? lo : (hi < v) ? hi : v;
}

template< class Real >
Point3D<Real> TriMesh< Real >::ClosestPointOnTri ( Point3D<Real> P, Point3D<Real> v1, Point3D<Real> v2, Point3D<Real> v3 )
{
    Point3D<Real> edge0 = v2 - v1;
    Point3D<Real> edge1 = v3 - v1;
    Point3D<Real> v0 = v1 - P;

    Real a = edge0.squareNorm ();
    Real b = Point3D<Real>::Dot (edge0, edge1);
    Real c = edge1.squareNorm ();
    Real d = Point3D<Real>::Dot ( edge0, v0);
    Real e = Point3D<Real>::Dot ( edge1, v0);

    Real det = a*c - b*b;
    Real s = b*e - c*d;
    Real t = b*d - a*e;

    if ( s + t < det )
    {
        if ( s < Real(0.0) )
        {
            if ( t < Real(0.0) )
            {
                if ( d < Real(0.0) )
                {
                    s = clamp( -d/a, Real(0.0), Real(1.0) );
                    t = Real(0.0);
                }
                else
                {
                    s = Real(0.0);
                    t = clamp( -e/c, Real(0.0), Real(1.0) );
                }
            }
            else
            {
                s = Real(0.0);
                t = clamp( -e/c, Real(0.0), Real(1.0) );
            }
        }
        else if ( t < Real(0.0) )
        {
            s = clamp( -d/a, Real(0.0), Real(1.0) );
            t = Real(0.0);
        }
        else
        {
            Real invDet = Real(1.0) / det;
            s *= invDet;
            t *= invDet;
        }
    }
    else
    {
        if ( s < Real(0.0) )
        {
            Real tmp0 = b+d;
            Real tmp1 = c+e;
            if ( tmp1 > tmp0 )
            {
                Real numer = tmp1 - tmp0;
                Real denom = a-2*b+c;
                s = clamp( numer/denom, Real(0.0), Real(1.0) );
                t = 1-s;
            }
            else
            {
                t = clamp( -e/c, Real(0.0), Real(1.0) );
                s = Real(0.0);
            }
        }
        else if ( t < Real(0.0) )
        {
            if ( a+d > b+e )
            {
                Real numer = c+e-b-d;
                Real denom = a-2*b+c;
                s = clamp( numer/denom, Real(0.0), Real(1.0) );
                t = 1-s;
            }
            else
            {
                s = clamp( -e/c, Real(0.0), Real(1.0) );
                t = Real(0.0);
            }
        }
        else
        {
            Real numer = c+e-b-d;
            Real denom = a-2*b+c;
            s = clamp( numer/denom, Real(0.0), Real(1.0) );
            t = 1 - s;
        }
    }

    return v1 + s * edge0 + t * edge1;

}

// Get barycentric coordinates for point on mesh
// Coordinate calculation from https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates

template< class Real >
std::pair<int, Point3D< double > > TriMesh< Real >::getBarycentric (Point3D<Real> P) const
{
    // Map P to closest point on mesh 
    int tIdx = -1;
    float d;
    Point3D<Real> Q;
    for( int i=0 ; i<_triangles.size() ; i++ )
    {
        Point3D< float > q = ClosestPointOnTri( P , _vertices[ _triangles[i][0] ] , _vertices[ _triangles[i][1] ] , _vertices[ _triangles[i][2] ] );

        float _d = ( P - q ).squareNorm();

        if( tIdx==-1 || _d<d ) 
        {
            tIdx = i , d = _d, Q = q; 

            if ( std::sqrt (d) < 16 * FLT_EPSILON )
            {
                break;
            }
        }

    }

    P = Q;

    Point3D< Real > u0 = _vertices[ _triangles[tIdx][1] ] - _vertices[ _triangles[tIdx][0] ];
    Point3D< Real > u1 = _vertices[ _triangles[tIdx][2] ] - _vertices[ _triangles[tIdx][0] ];
    Point3D< Real > u2 = P - _vertices[ _triangles[tIdx][0] ];

    Real d00 = u0.squareNorm ();
    Real d01 = Point3D<Real>::Dot (u0, u1);
    Real d11 = u1.squareNorm ();
    Real d20 = Point3D<Real>::Dot (u2, u0);
    Real d21 = Point3D<Real>::Dot (u2, u1);

    Real det = d00 * d11 - d01 * d01;

    Point3D< double > coords;

    coords[1] = (double) (d11 * d20 - d01 * d21 ) / det;
    coords[2] = (double) (d00 * d21 - d01 * d20 ) / det;
    coords[0] = (double) (1 - coords[1] - coords[2]);

    double coordSum = 0.0;
    for (int i = 0; i < 3; i++)
    {
        if (coords[i] < 0)
        {
            coords[i] = 0.0;
        }

        coordSum += coords[i];
    }

    coords /= coordSum;


    std::pair<int, Point3D< double >> coordInfo;

    coordInfo.second = coords;
    coordInfo.first = tIdx;

    return coordInfo;
}

template< class Real >
Point3D< double > TriMesh< Real >::getBarycentric (Point3D<Real> P, Point3D<Real> v0, Point3D<Real> v1, Point3D<Real> v2) const
{
    // Map P to closest point on triangle 
    P = ClosestPointOnTri ( P, v0, v1, v2);

    Point3D<Real> u0 = v1 - v0;
    Point3D<Real> u1 = v2 - v1;
    Point3D<Real> u2 = P - v0;

    Real d00 = u0.squareNorm ();
    Real d01 = Point3D<Real>::Dot (u0, u1);
    Real d11 = u1.squareNorm ();
    Real d20 = Point3D<Real>::Dot (u2, u0);
    Real d21 = Point3D<Real>::Dot (u2, u1);

    Real det = d00 * d11 - d01 * d01;

    Point3D< double > coords;

    coords[1] = (double) (d11 * d20 - d01 * d21 ) / det;
    coords[2] = (double) (d00 * d21 - d01 * d20 ) / det;
    coords[0] = (double) (1 - coords[1] - coords[2]);

    double coordSum = 0.0;
    for (int i = 0; i < 3; i++)
    {
        if (coords[i] < 0)
        {
            coords[i] = 0.0;
        }

        coordSum += coords[i];
    }

    coords /= coordSum;

    return coords;
}

template< class Real >
Point3D< double > TriMesh< Real >::getBarycentric (Point2D< double > P, Point2D< double > v1, Point2D< double > v2, Point2D< double > v3) const
{
    double detT = (v2[1] - v3[1]) * (v1[0] - v3[0]) + (v3[0] - v2[0]) * ( v1[1] - v3[1]);

    Point3D< double > coords;

    coords[0] = ( (v2[1] - v3[1]) * (P[0] - v3[0]) + (v3[0] - v2[0]) * ( P[1] - v3[1]) ) / detT;
    coords[1] = ( (v3[1] - v1[1]) * (P[0] - v3[0]) + (v1[0] - v3[0]) * ( P[1] - v3[1]) ) / detT;
    coords[2] = 1.0 - coords[0] - coords[1];

    return coords;
}

template< class Real >
Point3D< double > TriMesh< Real >::ComposeCoordinates ( Point3D< double > P, Point3D< double > v1, Point3D< double > v2, Point3D< double > v3)
{
    Point3D< double > coords;

    coords[0] = v1[0] * P[0] + v2[0] * P[1] + v3[0] * P[2];
    coords[1] = v1[1] * P[0] + v2[1] * P[1] + v3[1] * P[2];
    coords[2] = v1[2] * P[0] + v2[2] * P[1] + v3[2] * P[2];

    return coords;
}

template< class Real >
Point3D<Real> TriMesh< Real >::asPoint ( std::pair<int, Point3D< double >> baryC) const
{
    return Real( baryC.second[0] ) * _vertices[ _triangles[ baryC.first][0] ] + Real( baryC.second[1] ) * _vertices[ _triangles[baryC.first][1] ] + Real( baryC.second[2] ) * _vertices[ _triangles[baryC.first][2] ];
}

template< class Real >
Point3D<Real> TriMesh< Real >::vertexNormal( int l ) const
{
    return _vertexNormals[l];
}

template< class Real >
Point3D<Real> TriMesh< Real >::triangleNormal( int l ) const
{
    return _triangleNormals[l];
}

template< class Real >
Real TriMesh< Real >::vertexArea( int l ) const
{
    return _vertexAreas[l];
}

template< class Real >
Real TriMesh< Real >::triangleArea( int l ) const
{
    return _triangleAreas[l];
}

template< class Real >
Real TriMesh< Real >::totalArea( void ) const
{
    Real A = 0.0;

    for (int l = 0; l < _triangleAreas.size (); l++)
    {
        A += _triangleAreas[l];
    }

    return A;
}

template< class Real >
double TriMesh< Real >::triangleArea( const Point3D< double >& c1, const Point3D< double >& c2, const Point3D< double >& c3) const
{
    return 0.5 * std::sqrt( Point3D< double >::CrossProduct (c2 - c1, c3 - c1).squareNorm () );
}

template< class Real >
double TriMesh< Real >::triangleArea ( const Point2D< double >& c1, const Point2D< double >& c2, const Point2D< double >& c3) const
{
    return 0.5 * std::sqrt( Point3D< double >::CrossProduct (Point3D< double > (c2[0] - c1[0], c2[1] - c1[1], 0.0), Point3D< double > (c3[0] - c1[0], c3[1] - c1[1], 0.0) ).squareNorm () );
}


template< class Real >
double TriMesh< Real >::subTriangleArea ( int l , const Point3D< double > &c1 , const Point3D< double > &c2 , const Point3D< double > &c3 ) const
{
    Point3D< double > P1 = _vertices[ _triangles[l][0] ] * c1[0] + _vertices[ _triangles[l][1] ] * c1[1] + _vertices[ _triangles[l][2] ] *c1[2];
    Point3D< double > P2 = _vertices[ _triangles[l][0] ] * c2[0] + _vertices[ _triangles[l][1] ] * c2[1] + _vertices[ _triangles[l][2] ] *c2[2];
    Point3D< double > P3 = _vertices[ _triangles[l][0] ] * c3[0] + _vertices[ _triangles[l][1] ] * c3[1] + _vertices[ _triangles[l][2] ] *c3[2];

    Point3D< double > v1 = P2 - P1;
    Point3D< double > v2 = P3 - P1;

    return 0.5 * std::sqrt( Point3D< double >::CrossProduct (v1, v2).squareNorm () );
}

template< class Real >
const std::vector< TriangleIndex > &TriMesh< Real >::vertexStar ( int l ) const
{
    return _vertexStars[l];
}

template< class Real >
const std::vector< int > &TriMesh< Real >::vertexStarList ( int l ) const
{
    return _starTriangles[l];
}

template< class Real >
Real TriMesh< Real >::meanEdgeLength (void) const
{
    return _meanEL;
}

template< class Real >
Real TriMesh< Real >::maxEdgeLength (void) const
{
    return _maxEL;
}

template< class Real >
Real TriMesh< Real >::minEdgeLength (void) const
{
    return _minEL;
}

// Riemannian Metric (Euclidean distance)
template< class Real >
void TriMesh< Real >::initEuclideanMetrics (void)
{
    _triangleMetrics.resize( _triangles.size() );

    for( int i=0 ; i<_triangles.size() ; i++ )
    {
        double L12 = ( _vertices[ _triangles[i][0] ] - _vertices[ _triangles[i][1] ] ).squareNorm ();
        double L13 = ( _vertices[ _triangles[i][0] ] - _vertices[ _triangles[i][2] ] ).squareNorm ();
        double L23 = ( _vertices[ _triangles[i][1] ] - _vertices[ _triangles[i][2] ] ).squareNorm ();

        _triangleMetrics[i][0] = L12;
        _triangleMetrics[i][1] = (L12 + L13 - L23) / 2;
        _triangleMetrics[i][2] = L13;

    }
}


// Gradients

template< class Real >
Point2D< double > TriMesh< Real >::getMetricGradient( int l , const double phi1 , const double phi2 , const double phi3 ) const
{
    
    Point2D< double > grad(0, 0);

    double a = _triangleMetrics[l][0];
    double b = _triangleMetrics[l][1];
    double c = _triangleMetrics[l][2];

    
    double mDet = a * c - b * b;

    if (mDet > 0)
    {
       double v1 = phi2 - phi1;
       double v2 = phi3 - phi1;

       grad[0] = ( c * v1 - b * v2 ) / mDet;
       grad[1] = ( a * v2 - b * v1 ) / mDet;

    }
 
    return grad;
    
    
}

template< class Real >
Point2D< double > TriMesh< Real >::getMetricGradient( int l , const std::vector< double > &Implicit ) const
{
    
    Point2D< double > grad(0, 0);

    double a = _triangleMetrics[l][0];
    double b = _triangleMetrics[l][1];
    double c = _triangleMetrics[l][2];

    
    double mDet = a * c - b * b;

    if (mDet > 0)
    {
        double v1 = Implicit[ _triangles[l][1] ] - Implicit[ _triangles[l][0] ];
        double v2 = Implicit[ _triangles[l][2] ] - Implicit[ _triangles[l][0] ];

       grad[0] = ( c * v1 - b * v2 ) / mDet;
       grad[1] = ( a * v2 - b * v1 ) / mDet;

    }
 
    return grad;
}

template< class Real >
void TriMesh< Real >::metricGradient ( const std::vector< double >& Implicit , std::vector< Point2D< double > > & triGrads ) const
{
    triGrads.resize( _triangles.size() );

#pragma omp parallel for
    for( int l=0 ; l<_triangles.size() ; l++ ) triGrads[l] = getMetricGradient( l , Implicit );
}


// Metric dot products
template< class Real >
double TriMesh< Real >::metricDot( int l , const Point2D< double >& w1, const Point2D< double >& w2) const
{
    double a = _triangleMetrics[l][0];
    double b = _triangleMetrics[l][1];
    double c = _triangleMetrics[l][2];


    return (a * w1[0] + b * w1[1]) * w2[0] + (b * w1[0] + c * w1[1]) * w2[1];
  
}

template< class Real >
double TriMesh< Real >::metricSquareNorm( int l, const Point2D< double >& w) const
{
   return metricDot (l, w, w);
}

template< class Real >
Point2D< double > TriMesh< Real >::metricRotate90 ( int l , const Point2D< double >& w) const
{

    double a = _triangleMetrics[l][0];
    double b = _triangleMetrics[l][1];
    double c = _triangleMetrics[l][2];

    double detRoot = std::sqrt(a * c - b * b);

    if (detRoot > 0)
    {

       double alphaR = -(b * w[0] + c * w[1]) / detRoot;

       double betaR = (a * w[0] + b * w[1]) / detRoot;

       return Point2D< double >(alphaR, betaR);
    }
    else
    {
       return Point2D< double > (0, 0);
    }
}
    



// Curvatures
template< class Real >
void TriMesh< Real >::computeFundamentalMatricesAndShapeOperator( int l , SquareMatrix< Real , 2 > &firstForm , SquareMatrix< Real , 2 > &secondForm , SquareMatrix< Real , 2 > &shapeOperator) const
{
    // f(u,v) = v0 + u(v01) + v(v02) ==> fu = v01, fv = v02
    // N(u,v) = n0 + u(n1 - n0) + v(n2 - n0) ==> Nu = n1 - n0, Nv = n2 - n0
    Point3D<Real> v01 = _vertices[ _triangles[l][1] ] - _vertices[ _triangles[l][0] ];
    Point3D<Real> v02 = _vertices[ _triangles[l][2] ] - _vertices[ _triangles[l][0] ];
    firstForm(0, 0) = Point3D<Real>::Dot(v01, v01); // <fu, fu>
    firstForm(0, 1) = Point3D<Real>::Dot(v01, v02); // <fu, fv>
    firstForm(1, 0) = Point3D<Real>::Dot(v02, v01); // <fv, fu>
    firstForm(1, 1) = Point3D<Real>::Dot(v02, v02); // <fv, fv>
    Point3D<Real> n01 = _vertexNormals[ _triangles[l][1] ] - _vertexNormals[ _triangles[l][0] ];
    Point3D<Real> n02 = _vertexNormals[ _triangles[l][2] ] - _vertexNormals[ _triangles[l][0] ];
    secondForm(0, 0) = Point3D<Real>::Dot(n01, v01); // <Nu, fu>
    secondForm(0, 1) = Point3D<Real>::Dot(n02, v01); // <Nv, fu>
    secondForm(1, 0) = Point3D<Real>::Dot(n01, v02); // <Nu, fv>
    secondForm(1, 1) = Point3D<Real>::Dot(n02, v02); // <Nv, fv>
    shapeOperator = firstForm.inverse() * secondForm;
}

// Gaussian curvature
template< class Real >
void TriMesh< Real >::triangleGaussCurvature (std::vector< double >& triGauss ) const
{
    triGauss.resize( _triangles.size() , 0.0);

#pragma omp parallel for
    for( int l=0 ; l<_triangles.size() ; l++ )
    {
        SquareMatrix<Real, 2> firstForm , secondForm , shapeOperator;

        computeFundamentalMatricesAndShapeOperator( l , firstForm , secondForm , shapeOperator );

        triGauss[l] = (double) shapeOperator.determinant();
    }
}

template< class Real >
void TriMesh< Real >::vertexGaussCurvature( std::vector< double > &vertGauss ) const
{
    vertGauss.resize( _vertices.size(), 0.0);

    std::vector< double > triGauss;
    triangleGaussCurvature( triGauss );

    for( int l=0 ; l<_triangles.size() ; l++ ) for( int i=0 ; i<3 ; i++ ) vertGauss[ _triangles[l][i] ] += triGauss[l] * _triangleAreas[l] / 3.0;

#pragma omp parallel for
    for( int l=0 ; l<_vertices.size () ; l++ ) if( _vertexAreas[l]>16.0 * FLT_EPSILON ) vertGauss[l] /= _vertexAreas[l];
}


// Mean curvature
template< class Real >
void TriMesh< Real >::triangleMeanCurvature( std::vector< double > &triMean ) const
{
    triMean.resize(0);
    triMean.resize( _triangles.size() , 0.0 );

#pragma omp parallel for
    for( int l=0 ; l<_triangles.size () ; l++ )
    {
        SquareMatrix<Real, 2> firstForm, secondForm, shapeOperator;

        computeFundamentalMatricesAndShapeOperator( l , firstForm , secondForm , shapeOperator );

        triMean[l] = (double)shapeOperator.trace ();
    }
}

template< class Real >
void TriMesh< Real >::vertexMeanCurvature( std::vector< double > &vertMean ) const
{
    vertMean.resize(0);
    vertMean.resize( _vertices.size() , 0.0 );

    std::vector< double > triMean;
    triangleMeanCurvature( triMean );

    for( int l=0 ; l<_triangles.size() ; l++ ) for( int i=0 ; i<3 ; i++ )
        vertMean[ _triangles[l][i] ] += triMean[l] * _triangleAreas[l] / 3.0;
 
#pragma omp parallel for
    for( int l=0 ; l<_vertices.size() ; l++) if( _vertexAreas[l]>16.0*FLT_EPSILON ) vertMean[l] /= _vertexAreas[l];
}

template< class Real >
void TriMesh< Real >::returnMassStiffness( Real diffTime , SparseMatrix< double , int >& mass , SparseMatrix< double , int >& stiffness )
{
    FEM::RiemannianMesh< double > mesh( GetPointer( _triangles ) , _triangles.size() );
    mesh.setMetricFromEmbedding<3>( [&]( unsigned int idx ){ return _vertices[idx]; } );

    mesh.makeUnitArea();

    mass = mesh.template massMatrix< FEM::BASIS_0_WHITNEY >() * diffTime;
    stiffness = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >();
}

template< class Real >
void TriMesh< Real >::smoothVertexSignal( std::vector< double > &Implicit , Real diffTime )
{
    // Normalize mesh to have unit area
    FEM::RiemannianMesh< double > mesh( GetPointer( _triangles ) , _triangles.size() );
    mesh.template setMetricFromEmbedding< 3 >( [&]( unsigned int idx ){ return Point3D< double >( _vertices[idx] ); } );
    double area = mesh.area();
    mesh.makeUnitArea();

    // Compute and solve the system

    SparseMatrix< double , int > m = mesh.template massMatrix< FEM::BASIS_0_WHITNEY >() * diffTime;
    SparseMatrix< double , int > s = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >();

    SparseMatrix< double , int > M = m + s;

    EigenSolverCholeskyLLt< double > llt( M , false );

    llt.update( M );


    Pointer( double ) x = AllocPointer< double >( _vertices.size() );
    Pointer( double ) b = AllocPointer< double >( _vertices.size() );

    for( int i=0 ; i<_vertices.size() ; i++ ) 
    { 
        if ( !std::isnan( Implicit[i] ) ) x[i] = Implicit[i];
        else                              x[i] = 0.0;
    }

    m.Multiply( x , b );

    for( int i=0 ; i<_vertices.size () ; i++ ) x[i] *= 0.0;

    s.Multiply (x, b, MULTIPLY_ADD);

    llt.solve( b , x );

    for( int i=0 ; i<_vertices.size() ; i++ ) Implicit[i] = x[i];

    FreePointer (x);
    FreePointer (b);
}


// Adapted from https://math.stackexchange.com/q/100766
template< class Real >
Point3D<float> TriMesh< Real >::projectToTri ( Point3D<Real> P, int triIndex) const
{
    Point3D<Real> triNormal = _triangleNormals[triIndex];
    Point3D<Real> q = vertices[triangles[triIndex][0]];

    Real dist = Point3D<Real>::Dot (triNormal, q - P) / triNormal.squareNorm ();

    return P + (dist * triNormal);
}



// ==============================
// ======== Geodesics ===========
// ==============================

template< class Real >
void TriMesh< Real >::setGeodesicEpsilon( double epsilon )
{
    _logEpsilon = epsilon;
}

template< class Real >
void TriMesh< Real >::initGeodesicCalc( void )
{
    _surfaceMesh.fromTriangles( _vertices , _triangles );
}

template< class Real >
std::vector< double > TriMesh< Real >::computeGeodesicsAbout( int nodeIndex , float rho ) const
{
    DGPCgenerator dgpc( _surfaceMesh );
    dgpc.setEps( _logEpsilon );
    dgpc.setStopDist( rho );
    dgpc.setNodeSource( nodeIndex );
    dgpc.run ();

    return dgpc.getDistances ();
}

template< class Real >
std::vector< double > TriMesh< Real >::computeGeodesicsAbout( std::pair< int , Point3D< double > > P , float rho ) const
{
    DGPCgenerator dgpc( _surfaceMesh );
    dgpc.setEps( _logEpsilon );
    dgpc.setStopDist( rho );
    Point3D< Real > pt = _vertices[ _triangles[P.first][0] ] * (Real)P.second[0] + _vertices[ _triangles[P.first][1] ] * (Real)P.second[1] + _vertices[ _triangles[P.first][2] ] * (Real)P.second[2];
    dgpc.setSource( PointOM( pt[0] , pt[1] , pt[2] ) , P.first );
    dgpc.run();

    return dgpc.getDistances ();
}


template< class Real >
std::pair<std::vector< double >, std::vector< double >> TriMesh< Real >::computeLogarithmAbout( int nodeIndex , float rho ) const
{
    DGPCgenerator dgpc( _surfaceMesh );
    dgpc.setEps( _logEpsilon );
    dgpc.setStopDist( rho );
    dgpc.setNodeSource( nodeIndex );
    dgpc.run ();

    return std::pair< std::vector< double >, std::vector< double > >( dgpc.getDistances () , dgpc.getAngles () );

}

template< class Real >
std::pair<std::vector< double >, std::vector< double >> TriMesh< Real >::computeLogarithmAbout( std::pair< int , Point3D< double > > P , float rho ) const
{
    DGPCgenerator dgpc( _surfaceMesh );
    dgpc.setEps( _logEpsilon );
    dgpc.setStopDist( rho );
    Point3D<Real> pt = vertices[triangles[P.first][0]] * P.second[0] +  vertices[triangles[P.first][1]] * P.second[1] + vertices[triangles[P.first][2]] * P.second[2];
    dgpc.setSource( PointOM( pt[0] , pt[1] , pt[2] ) , P.first );
    dgpc.run();

    return std::pair< std::vector< double >, std::vector< double > > ( dgpc.getDistances (), dgpc.getAngles () );
}

// =============================
// === Heat Kernel Signature ===
// =============================

template< class Real >
void TriMesh< Real >::vertexHKS( const Spectrum< double > &spectrum , std::vector< double > &vertHKS , double diffTime ) const
{
    // Compute HKS
    vertHKS.resize( _vertices.size() , 0.0 );

#pragma omp parallel for
    for( int i=0 ; i<_vertices.size() ; i++ ) for( int j=0 ; j<spectrum.size() ; j++ )
        vertHKS[i] += (double) exp( - spectrum.eValue(j) * diffTime ) * spectrum.eVector(j)[i] * spectrum.eVector(j)[i];
}

// ==========================
// === Spectral Distances ===
// ==========================
template< class Real >
double TriMesh< Real >::spectralDist( const Spectrum< double > &spectrum , int nodeX , int nodeY ) const
{
    return spectrum.spectralDistance( nodeX , nodeY , 1 , std::numeric_limits< int >::max() );
}

template< class Real >
double TriMesh< Real >::spectralDist( const Spectrum< double > &spectrum , std::pair< int , Point3D< double > > P , int node ) const
{
    return spectrum.spectralDistance( node , _triangles[ P.first ] , P.second , 1 , std::numeric_limits< int >::max() );
}

template< class Real >
std::vector< double > TriMesh< Real >::computeSpectralDistancesAbout( const Spectrum< double > &spectrum , int nodeIndex, double rho ) const
{
    std::vector< double > sDistances; 
    sDistances.resize( _vertices.size () , std::numeric_limits< double >::max() );

    std::vector< bool > processed;
    processed.resize( _vertices.size() , false );

    std::vector< int > Q { nodeIndex };
    processed[ nodeIndex ] = true;
    sDistances[ nodeIndex ] = 0.0;

    while( Q.size() )
    {
        int q = Q.back();
        Q.pop_back();

        for( int l=0 ; l<_starVertices[q].size() ; l++ )
        {
            int r = _starVertices[q][l];
            if( !processed[r] )
            {
                sDistances[r] = spectralDist( spectrum , nodeIndex , r );
                processed[r] = true;
                if( sDistances[r]<=rho ) Q.push_back ( r );
            }
        }
    }

    return sDistances;
}

template< class Real >
std::vector< double > TriMesh< Real >::computeSpectralDistancesAbout( const Spectrum< double > &spectrum , std::pair<int, Point3D< double >> & P, double rho ) const
{
    std::vector< double > sDistances; 
    sDistances.resize( _vertices.size() , std::numeric_limits< double >::max() );

    std::vector<bool> processed;
    processed.resize( _vertices.size () , false );

    std::vector< int > Q { (int)_triangles[P.first][0] , (int)_triangles[P.first][1] , (int)_triangles[P.first][2] };
    processed[ _triangles[P.first][0] ] = processed[ _triangles[P.first][1] ] = processed[ _triangles[P.first][2] ] = true;

    sDistances[ _triangles[P.first][0] ] = spectralDist( spectrum , P , _triangles[P.first][0] );
    sDistances[ _triangles[P.first][1] ] = spectralDist( spectrum , P , _triangles[P.first][1] );
    sDistances[ _triangles[P.first][2] ] = spectralDist( spectrum , P , _triangles[P.first][2] );

    while( Q.size () > 0 )
    {
        int q = Q.back ();
        Q.pop_back ();

        for (int l = 0; l < _starVertices[q].size (); l++)
        {
            int r = _starVertices[q][l];
            if ( !processed[r] )
            {
                sDistances[r] = spectralDist( spectrum , P , r );
                processed[r] = true;
                if( sDistances[r]<=rho ) Q.push_back ( r );
            }
        }
    }

    return sDistances;
}

template< class Real >
double TriMesh< Real >::getSpectralArea( const Spectrum< double > &spectrum , int l ) const
{
    double d1 = spectralDist( spectrum , _triangles[l][0] , _triangles[l][1] );
    double d2 = spectralDist( spectrum , _triangles[l][0] , _triangles[l][2] );
    double d3 = spectralDist( spectrum , _triangles[l][1] , _triangles[l][2] );

    double s = (d1 + d2 + d3) / 2;

    return std::sqrt ( s * (s - d1) * (s - d2) * (s - d3) );
}

template< class Real >
void TriMesh< Real >::initSpectralArea( const Spectrum< double > &spectrum )
{
#pragma omp parallel for
    for( int l=0 ; l<_triangles.size() ; l++ ) _triangleAreas[l] = (Real)getSpectralArea( spectrum , l );
}

template< class Real >
void TriMesh< Real >::initSpectralMetrics( const Spectrum< double > &spectrum )
{
    _triangleMetrics.resize( _triangles.size() );

#pragma omp parallel for
    for( int i=0 ; i<_triangles.size() ; i++ )
    {
        double l12 = spectralDist( spectrum , _triangles[i][0] , _triangles[i][1] );
        double l13 = spectralDist( spectrum , _triangles[i][0] , _triangles[i][2] );
        double l23 = spectralDist( spectrum , _triangles[i][1] , _triangles[i][2] );

        _triangleMetrics[i][0] = l12 * l12;
        _triangleMetrics[i][1] = (l12 * l12 + l13 * l13 - l23 * l23) / 2;
        _triangleMetrics[i][2] = l13 * l13;
    }
}
