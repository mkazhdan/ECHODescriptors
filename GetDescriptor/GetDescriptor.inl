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

#ifndef ECHO
#define ECHO

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif 

#include <iostream>
#include <Misha/Ply.h>
#include <Misha/Geometry.h>
#include <Misha/RegularGrid.h>
#include <Misha/TriMesh.h>
#include <Misha/RightTriangleQuadrature.h>

enum
{
	DISTANCE_BIHARMONIC ,
	DISTANCE_GEODESIC ,
	DISTANCE_DIFFUSION ,
    DISTANCE_COUNT
};
const static std::string DistanceNames[] = { "biharmonic" , "geodesic" , "diffusion" };


#ifndef T_DIFF
#define T_DIFF 0.1
#endif

struct QuadratureSample
{
    QuadratureSample( Point2D< double > p=Point2D< double >() , double i=0 ) : integral(i) , position(p) {}
    Point2D< double > position;
    double integral;
};


void splatBoundingBox ( Point2D<int>& rllc, Point2D<int>& rurc, const Point2D<double>& v1, const Point2D<double>& v2, const Point2D<double>& v3, int descriptorRadius, int win )
{
    rllc[0] = std::max( -descriptorRadius, (int) std::round ( std::min ( v1[0], std::min (v2[0], v3[0]) ) - win ));
    rllc[1] = std::max( -descriptorRadius, (int) std::round ( std::min ( v1[1], std::min (v2[1], v3[1]) ) - win) );

    rurc[0] = std::min( descriptorRadius, (int) std::round ( std::max ( v1[0], std::max (v2[0], v3[0]) ) + win ) );
    rurc[1] = std::min( descriptorRadius, (int) std::round ( std::max ( v1[1], std::max (v2[1], v3[1]) ) + win ) );
}


template< unsigned int nSamples >
unsigned int computeTriangleWeights( const Point2D< double > v[] , const double f[] , double rho , QuadratureSample quadrature[] )
{
    unsigned int count = 0;
    typedef TriangleIntegrator< nSamples > Integrator;
    for( int l=0 ; l<nSamples ; l++ )
    {
        const double w[] = { Integrator::Positions[l][0] , Integrator::Positions[l][1] , 1. - Integrator::Positions[l][0] - Integrator::Positions[l][1] };
        Point2D< double > q = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
        if( q.squareNorm()<=rho*rho ) quadrature[count++] = QuadratureSample( q , ( f[0] * w[0] + f[1] * w[1] + f[2] * w[2] ) * Integrator::Weights[l] );
    }
    return count;
}

template< class Real >
void tangentsAtVertices( const TriMesh<Real> & tMesh,  const std::vector<Point2D<double>>& triangleGradients, std::vector<double>& sDistances, std::vector<Point2D<double> >& logTan, std::vector<double>& mags)
{
    std::vector< Point2D< double > > triDir;
    const std::vector< Point3D< Real > > &vertices = tMesh.vertices();
    const std::vector< TriangleIndex > &triangles = tMesh.triangles();
    triDir.resize ( triangles.size() );

    for( int l=0 ; l<triangles.size(); l++)
    {
        bool degenerate = false;
        for( int i=0 ; i<3 ; i++ ) { if( sDistances[ triangles[l][i] ] >= std::numeric_limits<double>::max() ) { degenerate = true; break; } }

        if( degenerate ) triDir[l][0] = triDir[l][1] = std::numeric_limits<Real>::max();
        else
        {
            Point2D<double> tangentX = triangleGradients[l] / std::sqrt (tMesh.metricSquareNorm (l, triangleGradients[l]));

            Point2D<double> tangentY = tMesh.metricRotate90 (l, tangentX);

            // Compute logarithm direction
            Point2D< double > geoGrad = tMesh.getMetricGradient( l , sDistances[ triangles[l][0] ] , sDistances[ triangles[l][1] ] , sDistances[ triangles[l][2]] );

            geoGrad /= std::sqrt(tMesh.metricSquareNorm(l, geoGrad)) * -1;


            double detR = tangentX[0] * tangentY[1] - tangentX[1] * tangentY[0];

            triDir[l][0] = (tangentY[1] * geoGrad[0] - tangentY[0] * geoGrad[1]) / detR;
            triDir[l][1] = (tangentX[0] * geoGrad[1] - tangentX[1] * geoGrad[0]) / detR;

         
        }
    }

    logTan.resize ( vertices.size() , Point2D< double >() );
    mags.resize( vertices.size() , 0 );

    for( int l=0 ; l<vertices.size() ; l++ )
    {
        if (sDistances[l] >= std::numeric_limits<double>::max ()) logTan[l] = Point2D<double> ( std::numeric_limits<double>::max (), std::numeric_limits<double>::max () );
        else
        {
            std::vector<int> starTris = tMesh.vertexStarList(l);
            double totalArea = 0.0;
            int tCount = 0;
            for (int k = 0; k < starTris.size (); k++)
            {
                if (triDir[starTris[k]][0] < std::numeric_limits<Real>::max ())
                {
                    
                    tCount ++;
                    double tArea = tMesh.triangleArea (starTris[k]);

                    logTan[l] += triDir[starTris[k]] * tArea;

                    mags[l] += std::sqrt( tMesh.metricSquareNorm (starTris[k], triangleGradients[starTris[k]]) ) * tArea;
                    totalArea += tArea;
                } 
            }

            if (tCount > 0)
            {
                logTan[l] *= sDistances[l] / std::sqrt(logTan[l].squareNorm ());
                mags[l] /= totalArea;
            }
            else
            {
                logTan[l] = Point2D<double> ( std::numeric_limits<double>::max (), std::numeric_limits<double>::max () );
            }
        }
    }
}




template< class Real, unsigned int nSamples>
RegularGrid< Real , 2 > echo( const TriMesh< Real > &tMesh , const std::vector< Point2D< double > >& logTan , std::vector< double > &mags , double supportRadius, unsigned int descriptorRadius )
{
   auto Mod = []( int i , unsigned int modulus )
   {
      if( i<0 ) return modulus - ( (-i) % modulus );
      else      return i % modulus;
   };

   const std::vector< TriangleIndex > &triangles = tMesh.triangles();
   const std::vector< Point3D< Real > > &vertices = tMesh.vertices();

   RegularGrid< Real , 2 > F;
   unsigned int res = 2 * descriptorRadius + 1;
   F.resize( res , res );
   for( int i=0 ; i<F.resolution() ; i++ ) F[i] = 0;

   // The mesh-to-descriptor scale factor
   Real surface2descriptor = (Real)( descriptorRadius / supportRadius );

   // Gaussian fall-off
   Real splatRad = (Real)( 1.3 * descriptorRadius / 5.0 );
   Real sigma2 = (Real)( splatRad * splatRad / std::log (0.05) ); 
   int win = (int) std::ceil (splatRad);

   double rho = surface2descriptor * supportRadius;

   for( int l=0 ; l<triangles.size() ; l++ )
   { 
      // Check if triangle is on boundary of ball
      bool degenerate = false;
      int nOut = 0;
      for (int i = 0; i < 3; i++)
      {
          double n2T = logTan[ triangles[l][i] ].squareNorm ();

         if ( n2T >= std::numeric_limits< double >::max() ) 
         {
             degenerate = true;
             break;
         }
         
         if (n2T >= supportRadius * supportRadius ) 
         {
            nOut ++;
         }
      }


      if (!degenerate && nOut<3)
      { 
          Point2D< double > nodes[] = { logTan[ triangles[l][0] ] * surface2descriptor , logTan[ triangles[l][1] ] * surface2descriptor , logTan[ triangles[l][2] ] * surface2descriptor };

          double phi[] = { mags[ triangles[l][0] ] , mags[ triangles[l][1] ], mags[ triangles[l][2] ] };

         QuadratureSample q[nSamples];

         unsigned int qCount = computeTriangleWeights< nSamples >( nodes , phi , rho , q );

         Point2D<int> rllc, rurc;
         splatBoundingBox( rllc, rurc, nodes[0] , nodes[1] , nodes[2] , descriptorRadius , win );

         for (int x0 = rllc[0]; x0 <= rurc[0]; x0++ )
         {
            for (int y0 = rllc[1]; y0 <= rurc[1]; y0++)
            {
               if ( ( x0 * x0 + y0 * y0 ) <= ( ( descriptorRadius + 0.5) * (descriptorRadius + 0.5) ) )
               {

                  double value = 0;
                  for( int k=0 ; k<(int)qCount ; k++ )
                  { 
                     value += std::exp ( ( Point2D< double >( x0 , y0 ) - q[k].position ).squareNorm () / sigma2 ) * q[k].integral;
 
                     F ( Mod( x0 + descriptorRadius, F.res(0) ) , Mod (y0 + descriptorRadius, F.res(1) ) ) += (Real) ( value * tMesh.triangleArea (l) );
                  }

          	   }
                
            }
            
         }
      
      } // end if (!degenerate)
      
   } // end for (int l = 0; l < tMesh.triangles (); l++)
   
   return F;
}



template< class Real , unsigned int nSamples = 7>
RegularGrid< Real , 2 > echo( const TriMesh< Real > &tMesh , const std::vector< Point2D< double > > &triangleGradients , int nodeIndex , double supportRadius, unsigned int descriptorRadius, int distFlag = DISTANCE_BIHARMONIC)
{

   // Compute surface distances
	std::vector<double> sDistances;

   if (distFlag == DISTANCE_GEODESIC)
   {
      sDistances = tMesh.computeGeodesicsAbout( nodeIndex , (float)supportRadius);
   }
   else if ( distFlag == DISTANCE_DIFFUSION)
   {
      sDistances = tMesh.computeDiffusionsAbout( nodeIndex , T_DIFF , supportRadius );
   }
   else
   {
      sDistances = tMesh.computeBiharmonicsAbout (nodeIndex, supportRadius  );
   }

   // Pre-compute tangent space coordinates + signal
   std::vector< Point2D< double > > logTan;
   std::vector< double > mags;

   tangentsAtVertices( tMesh , triangleGradients , sDistances , logTan , mags );

   // Compute descriptor
   return echo<Real, nSamples> (tMesh, logTan, mags, supportRadius, descriptorRadius);

}

template< class Real, unsigned int nSamples = 7>
RegularGrid< Real , 2 > echo( const TriMesh< Real > &tMesh , const std::vector< Point2D< double > > &triangleGradients , std::pair< int , Point3D< double > > P , double supportRadius, unsigned int descriptorRadius, int distFlag = DISTANCE_BIHARMONIC )
{
    // Compute surface distances
	std::vector<double> sDistances;

   if (distFlag == DISTANCE_GEODESIC)
   {
      sDistances = tMesh.computeGeodesicsAbout ( P, (float)supportRadius );
   }
   else if ( distFlag == DISTANCE_DIFFUSION)
   {
      sDistances = tMesh.computeDiffusionsAbout (P, T_DIFF, supportRadius );
   }
   else
   {
      sDistances = tMesh.computeBiharmonicsAbout (P, supportRadius);
   }

   // Pre-compute tangent space coordinates + signal
   std::vector< Point2D< double > > logTan;
   std::vector< double > mags;

   tangentsAtVertices( tMesh , triangleGradients , sDistances , logTan , mags );

   // Compute descriptor
   return echo<Real, nSamples> (tMesh, logTan, mags, supportRadius, descriptorRadius);
}





#endif

