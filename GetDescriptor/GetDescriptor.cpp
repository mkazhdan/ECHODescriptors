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

#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <Misha/Ply.h>
#include <Misha/CmdLineParser.h>
#include <Misha/Geometry.h>
#include <Misha/Image.h>
#include <Misha/TriMesh.h>
#include <Misha/Miscellany.h>
#include <Misha/PlyVertexData.h>
#include <Misha/Miscellany.h>
#include <Misha/Polynomial.h>

#undef DEBUG_DESCRIPTOR
#include "GetDescriptor.inl"



Misha::CmdLineParameter< std::string > In( "in" ) , Spec( "spec" ) , Out( "out" ); 
Misha::CmdLineParameter< int > SourceNode( "vertex" ) , HistogramRadius( "hRadius" , 5 ) , OutResolution( "resolution" ) , SourceFace( "tri") , DistanceType( "distance" , DISTANCE_BIHARMONIC );
Misha::CmdLineParameter< float > ThreshFactor( "tau" , 0.08f ) , Deviation( "dev" , -1. ) , DiffusionTime( "diffusion" , 0.1f );
Misha::CmdLineParameterArray< float , 3 > BC( "bc" );
Misha::CmdLineReadable  Verbose( "verbose" ) , NoDiskSupport( "noDisk" ) , ExactGaussian( "exactG" );

Misha::CmdLineReadable* params[] =
{
    &In ,
    &Spec ,
    &Out ,
    &SourceNode ,
    &HistogramRadius ,
    &ThreshFactor ,
    &Verbose ,
    &Deviation ,
    &OutResolution ,
    &NoDiskSupport ,
    &SourceFace ,
    &BC ,
    &DistanceType ,
    &DiffusionTime ,
    &ExactGaussian ,
    NULL
};

void ShowUsage( const char* ex )
{
    printf( "Usage %s:\n" , ex );
    printf( "\t --%s <input mesh>\n" , In.name.c_str() );
    printf( "\t --%s <source vertex index>\n" , SourceNode.name.c_str() );
    printf( "\t --%s <source face index>\n" , SourceFace.name.c_str() );
    printf( "\t --%s <barycentric coordinate 1> <barycentric coordinate 2> <barycentric coordinate 3>\n" , BC.name.c_str() );
    printf( "\t[--%s <spectral decomposition>]\n" , Spec.name.c_str() );
    printf( "\t[--%s <output ECHO descriptor>]\n" , Out.name.c_str() );
    printf( "\t[--%s <distance type>=%d]\n" , DistanceType.name.c_str() , DistanceType.value );
    for( int i=0 ; i<DISTANCE_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i , DistanceNames[i].c_str() );
    printf( "\t[--%s <diffusion time>=%f]\n" , DiffusionTime.name.c_str() , DiffusionTime.value );
    printf( "\t[--%s <mesh area to support radius scale>=%.2f]\n" , ThreshFactor.name.c_str() , ThreshFactor.value );
    printf( "\t[--%s <histogram radius (in bin units)>=%d, size of histogram will be (2 * n + 1)^2]\n" , HistogramRadius.name.c_str() , HistogramRadius.value );
    printf( "\t[--%s <target deviation (for color output)>=%f]\n" , Deviation.name.c_str() , Deviation.value );
    printf( "\t[--%s <resampled output resolution>=<histogram radius>*2+1]\n" , OutResolution.name.c_str() );
    printf( "\t[--%s]\n" , NoDiskSupport.name.c_str() );
    printf( "\t[--%s]\n" , ExactGaussian.name.c_str() );
    printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

template< typename Real > using PlyVertexFactory = VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , 3 > >;
template< typename Real > using PlyVertex = typename PlyVertexFactory< Real >::VertexType;

RegularGrid< float , 2 > TransposeSignal( const RegularGrid< float , 2 > &in )
{
    RegularGrid< float , 2 > out;
    out.resize( in.res(1) , in.res(0) );
    for( unsigned int i=0 ; i<in.res(0) ; i++ ) for( unsigned int j=0 ; j<in.res(1) ; j++ ) out(j,i) = in(i,j);
    return out;
}

RegularGrid< float , 2 > ResampleSignal( const RegularGrid< float , 2 > &in , unsigned int resX , unsigned int resY )
{
    RegularGrid< float , 2 > out;
    out.resize( resX , resY );

    for( unsigned int i=0 ; i<resX ; i++ ) for( unsigned int j=0 ; j<resY ; j++ )
    {
        float x = (float)( ((double)i)/(resX-1) * (in.res(0)-1) );
        float y = (float)( ((double)j)/(resY-1) * (in.res(1)-1) );
        out(i,j) = in(x,y);
    }
    return out;
}

RegularGrid< float , 2 > ResampleSignalDisk( const RegularGrid< float , 2 > &in , unsigned int resX , unsigned int resY )
{
    RegularGrid< float , 2 > out;
    out.resize( resX , resY );

    int diskCenter = (in.res(0) - 1) / 2;

    for( unsigned int i=0 ; i<resX ; i++ ) for( unsigned int j=0 ; j<resY ; j++ )
    {
        float x = (float)( ((double)i)/(resX-1) * (in.res(0)-1) );
        float y = (float)( ((double)j)/(resY-1) * (in.res(1)-1) );
        if ( (x - diskCenter) * (x - diskCenter) + (y - diskCenter) * (y - diskCenter) <= diskCenter * diskCenter )
        {        
           out(i,j) = in(x,y);
        }
        else
        {
           out(i, j) = std::numeric_limits< float >::infinity();
        }

    }
    return out;
}

template< unsigned int Derivatives >
Polynomial< 1 , 2*Derivatives+1 > FitDerivativeValuesToPolynomial( double startX , double endX , double startDerivatives[] , double endDerivatives[] )
{
    Polynomial< 1 , 2*Derivatives+1 > monomials[2*Derivatives+2];
    for( int i=0 ; i<=2*Derivatives+1 ; i++ ) monomials[i][i] = 1;
    // P(x) = a_0 + a_1 x + ... a_{2d+1} x^{2d+1}
    SquareMatrix< double , 2*Derivatives+2 > A;
    Point< double , 2*Derivatives+2 > b , x;
    for( int d=0 ; d<=Derivatives ; d++ )
    {
        for( int j=0 ; j<2*Derivatives+2 ; j++ )
        {
            A(j,2*d+0) = monomials[j]( startX );
            A(j,2*d+1) = monomials[j](   endX );
            monomials[j] = monomials[j].d();
        }
        b[2*d+0] = startDerivatives[d];
        b[2*d+1] =   endDerivatives[d];
    }
    x = A.inverse() * b;
    Polynomial< 1 , 2*Derivatives+1 > P;
    for( int i=0 ; i<2*Derivatives+2 ; i++ ) P[i] = x[i];
    return P;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// MAIN
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
template< class Real >
void run( void )
{
    static const int NSamples = 7;
    Miscellany::Timer timer;
    std::vector< double > signalValues;
    std::vector< double > signal;
    std::vector< std::pair< Point2D< double > , Point2D< double > > > dualFrameField;
    std::vector< double > hks;
    Spectrum< double > spectrum;
    std::function< double ( double ) > spectralFunction = SpectralFunction( DistanceType.value , DiffusionTime.value );
    std::function< double ( double ) > weightFunction;

    int nRadialBins = HistogramRadius.value;

    //==The Weight function==
    static const int Derivatives = 1;
    if( ExactGaussian.set ) weightFunction = []( double d2 ){ return exp( -d2 ); };
    else
    {
        Polynomial< 1 , 2*Derivatives+1 > P[2];
        static const double SupportRadius = 3.;
        static const double MidPoint = SupportRadius / 2.;
        {
            double startX = 0 , endX = MidPoint;
            double startValues[ Derivatives+1 ] , endValues[ Derivatives+1 ];
            for( int d=0 ; d<=Derivatives ; d++ )
            {
                startValues[d] = (d&1) ? -exp(-startX) : exp(-startX);
                endValues[d] = (d&1) ? -exp(-endX) : exp(-endX);
            }
            P[0] = FitDerivativeValuesToPolynomial< Derivatives >( startX , endX , startValues , endValues );
        }
        {
            double startX = MidPoint , endX = SupportRadius;
            double startValues[ Derivatives+1 ] , endValues[ Derivatives+1 ];
            for( int d=0 ; d<=Derivatives ; d++ )
            {
                startValues[d] = (d&1) ? -exp(-startX) : exp(-startX);
                endValues[d] = (d&1) ? -exp(-endX) : exp(-endX);
            }
            P[1] = FitDerivativeValuesToPolynomial< Derivatives >( startX , endX , startValues , endValues );
        }
        weightFunction = [=]( double x )
        {
            if     ( x<MidPoint ) return P[0](x);
            else if( x<SupportRadius ) return P[1](x);
            else return 0.;
        };
    }
#ifdef DEBUG_DESCRIPTOR
    {
        double e=0;
        static const int COUNT = 10000;
        for( int i=0 ; i<10000 ; i++ )
        {
            double x = (i+0.5)/COUNT * SupportRadius;
            double d = exp( -x ) - weightFunction( x );
            e += d*d;
        }
        std::cout << "Error: " << e/COUNT << std::endl;
    }
#endif // DEBUG_DESCRIPTOR


    //==Load Mesh===
    std::vector< TriangleIndex > triangles;
    std::vector< Point3D< float > > vertices;
    timer.reset();
    {
        int fileType;
        PlyVertexFactory< float > vFactory;
        std::vector< PlyVertex< float > > _vertices;
        PLY::ReadTriangles( In.value , vFactory , _vertices , triangles , NULL , fileType );
        vertices.resize( _vertices.size() );
        for( int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = Point3D< float >( _vertices[i].get<0>() );
    }
    TriMesh< float > tMesh( vertices , triangles );
    if( Verbose.set ) std::cout << "\tGot mesh: " << timer.elapsed() << "(s)" << std::endl;

    //==Load Spectral Decomposition==
    timer.reset();
    if( Spec.set ) spectrum.read( Spec.value );
    else           spectrum.set( vertices , triangles , 200 , 100.f , false );
    if( Verbose.set ) std::cout << "\tGot spectrum: " << timer.elapsed() << "(s)" << std::endl;

    // Compute + smooth HKS
    timer.reset();
    hks.resize( vertices.size() );
#pragma omp parallel for
    for( int i=0 ; i<vertices.size() ; i++ ) hks[i] = spectrum.HKS( i , 0.1 );
    if( Verbose.set ) std::cout << "\tGot HKS: " << timer.elapsed() << "(s)" <<  std::endl;

    //WARN( "Why are we smoothing the HKS instead of using a larger time-step?" );
    timer.reset();
    tMesh.smoothVertexSignal( hks , 1.0e7 ); 
    if( Verbose.set ) std::cout << "\tSmoothed HKS: " << timer.elapsed() << "(s)" <<  std::endl;

    // [WARNING] We are re-scaling the eigenvectors here
    if( IsSpectral( DistanceType.value ) )
#pragma omp parallel for
        for( int i=0 ; i<spectrum.size() ; i++ )
        {
            double scale = spectralFunction( spectrum.eValue(i) );
            auto &ev = spectrum.eVector(i);
            for( int j=0 ; j<ev.size() ; j++ ) ev[j] *= scale;
        }

    timer.reset();
    if( IsSpectral( DistanceType.value ) )
    {
        auto squareEdgeLengthFunctor = [&]( TriangleIndex tri )
        {
            Point3D< double > d;
            for( int i=0 ; i<3 ; i++ ) d[i] = spectrum.spectralDistance( tri[(i+1)%3] , tri[(i+2)%3] , 1 , (int)spectrum.size() );
            for( int i=0 ; i<3 ; i++ ) d[i] *= d[i];
            return d;
        };
        tMesh.initMetricsFromSquareEdgeLengths( squareEdgeLengthFunctor );
    }
    else
    {
        auto squareEdgeLengthFunctor = [&]( TriangleIndex tri )
        {
            Point3D< double > d;
            for( int i=0 ; i<3 ; i++ ) d[i] = ( vertices[ tri[(i+1)%3] ] - vertices[ tri[(i+2)%3] ] ).squareNorm();
            return d;
        };
        tMesh.initMetricsFromSquareEdgeLengths( squareEdgeLengthFunctor );
        tMesh.initGeodesicCalc();
    }
    if( Verbose.set ) std::cout << "\tSet triangle metrics: " << timer.elapsed() << "(s)" <<  std::endl;

    // Compute the frame field and the signal
    timer.reset();
    signal.resize( triangles.size() );
    dualFrameField.resize( triangles.size() );
    {
        std::vector< Point2D< double > > triangleGradients;
        tMesh.metricGradient( hks , triangleGradients );
#pragma omp parallel for
        for( int l=0 ; l<triangles.size(); l++)
        {
            SquareMatrix< double , 2 > m = tMesh.triangleMetric( l );
            signal[l] = std::sqrt( Point2D< double >::Dot( triangleGradients[l] , m * triangleGradients[l] ) );
            dualFrameField[l].first  = triangleGradients[l] / signal[l];
            dualFrameField[l].second = tMesh.metricRotate90( l , dualFrameField[l].first );
            dualFrameField[l].first = m * dualFrameField[l].first;
            dualFrameField[l].second = m * dualFrameField[l].second;
        }
    }
    if( Verbose.set ) std::cout <<"\tGot frame field and signal: " << timer.elapsed() << "(s)" << std::endl;

    float rho = (float)( ThreshFactor.value * std::sqrt( tMesh.totalArea() / M_PI ) );

    if( SourceNode.set && SourceNode.value<0 )
    {
        timer.reset();
#ifdef DEBUG_DESCRIPTOR
        verticesInNeighborhood = 0;
#endif // DEBUG_DESCRIPTOR
        for( int i=0 ; i<-SourceNode.value ; i++ )
        {
            if( IsSpectral( DistanceType.value ) ) spectralEcho< float , NSamples >( tMesh , spectrum , signal , dualFrameField , rand() % vertices.size() , rho , nRadialBins , weightFunction );
            else                                   geodesicEcho< float , NSamples >( tMesh , signal , dualFrameField, rand() % vertices.size() , rho , nRadialBins , weightFunction );
        }
#ifdef DEBUG_DESCRIPTOR
        if( Verbose.set )
        {
            std::cout << "\t\tAverage vertices in neighborhood: " << (double)verticesInNeighborhood / (-SourceNode.value) << std::endl;
            std::cout << "\t\tTime per vertex: " << timer.elapsed() * 1000 / verticesInNeighborhood << "(ms)" << std::endl;
        }
#endif // DEBUG_DESCRIPTOR
        if( Verbose.set )
        {
            double e = timer.elapsed();
            std::cout << "\tGot ECHO descriptors: " << e << "(s) / " << (-SourceNode.value) << " = " << ( e / -SourceNode.value ) << "(s)" <<  std::endl;
        }
    }
    else
    {
        // Compute ECHO descriptor
        RegularGrid< float , 2 > F;

        timer.reset();
#ifdef DEBUG_DESCRIPTOR
        verticesInNeighborhood = 0;
#endif // DEBUG_DESCRIPTOR
        if( SourceNode.set )
        {
            if( IsSpectral( DistanceType.value ) ) F = spectralEcho< float , NSamples >( tMesh , spectrum , signal , dualFrameField , SourceNode.value , rho , nRadialBins , weightFunction );
            else                                   F = geodesicEcho< float , NSamples >( tMesh , signal , dualFrameField , SourceNode.value , rho , nRadialBins , weightFunction );
        }
        else
        {
            if( IsSpectral( DistanceType.value ) ) F = spectralEcho< float , NSamples >( tMesh , spectrum , signal , dualFrameField , std::pair< int , Point3D< double > >( SourceFace.value , Point3D< double >( BC.values[0] , BC.values[1] , BC.values[2] ) ) , rho , nRadialBins , weightFunction );
            else                                   F = geodesicEcho< float , NSamples >( tMesh , signal , dualFrameField , std::pair< int , Point3D< double > >( SourceFace.value , Point3D< double >( BC.values[0] , BC.values[1] , BC.values[2] ) ) , rho , nRadialBins , weightFunction );
        }
#ifdef DEBUG_DESCRIPTOR
        if( Verbose.set )
        {
            std::cout << "\t\tVertices in neighborhood: " << verticesInNeighborhood << std::endl;
            std::cout << "\t\tTime per vertex: " << timer.elapsed() * 1000 / verticesInNeighborhood << "(ms)" << std::endl;
        }
#endif // DEBUG_DESCRIPTOR
        if( Verbose.set ) std::cout << "\tGot ECHO descriptor: " << timer.elapsed() << "(s)" <<  std::endl;

        if( NoDiskSupport.set ) F = ResampleSignal( F , OutResolution.value , OutResolution.value );
        else                    F = ResampleSignalDisk( F , OutResolution.value , OutResolution.value );

        F = TransposeSignal( F );

        if( Out.set )
        {
            std::string ext = Misha::ToLower( Misha::GetFileExtension( Out.value ) );
            if( ext==std::string( "txt" ) )
            {
                std::ofstream echoD;
                echoD.open( Out.value );

                for( unsigned int x=0 ; x<F.res(0) ; x++ ) for( unsigned int y=0 ; y<F.res(1) ; y++ )
                {
                    if( x==F.res(0)-1 && y==F.res(1)-1 ) echoD << F(x,y) << std::endl;
                    else                                 echoD << F(x,y) << " ";
                }

                echoD.close();
            }
            else
            {
                unsigned char *pixels = new unsigned char[ F.resolution()*3 ];

                int count = 0;
                double dev = 0;
                double sumD = 0;
                for( int i=0 ; i<F.resolution() ; i++ ) if( F[i]<std::numeric_limits< float >::infinity() ) dev += F[i] * F[i] , sumD += F[i], count++;

                dev = sqrt( dev/count );
                double hue , saturation;
                if( Deviation.value<=0 )
                {
                    saturation = 0.;
                    hue = 0;
                    std::cout << "Deviation: " << dev << std::endl;
                    std::cout << "Sum: " << sumD << std::endl;
                }
                else
                {
                    saturation = 1.;
                    hue = 4.*M_PI/3. * dev / Deviation.value;
                }
                for( int i=0 ; i<F.resolution() ; i++ )
                {
                    if( F[i]==std::numeric_limits< float >::infinity() ) pixels[3*i+0] = pixels[3*i+1] = pixels[3*i+2] = 255;
                    else
                    {
                        double d = std::max( 0. , std::min( F[i] / (3.*dev) , 1. ) );
                        Point3D< double > rgb , hsv( hue , saturation , d );
                        Miscellany::HSVtoRGB( &hsv[0] , &rgb[0] );
                        for( int c=0 ; c<3 ; c++ ) pixels[3*i+c] = (unsigned char)(int)floor( rgb[c] * 255 );
                    }
                }
                ImageWriter::Write( Out.value.c_str() , pixels , OutResolution.value , OutResolution.value , 3 );
                delete[] pixels;
            }
        }
    }
}

int main( int argc , char* argv[] )
{
    Misha::CmdLineParse( argc-1 , argv+1 , params );

    if( !In.set || ( !SourceNode.set && !SourceFace.set ) || ( SourceFace.set && !BC.set ) )
    {
        ShowUsage( argv[0] );
        return EXIT_FAILURE;
    }
    if( !OutResolution.set ) OutResolution.value = HistogramRadius.value * 2 + 1;

    Miscellany::Timer timer;
    run< float >();
    if( Verbose.set )
        if( SourceNode.value<1 ) std::cout << "Got desriptors in: " << timer.elapsed() << "(s)" <<  std::endl;
        else                     std::cout << "Got desriptor in: " << timer.elapsed() << "(s)" <<  std::endl;

    return 1;
}

