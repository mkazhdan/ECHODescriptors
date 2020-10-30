/*
Copyright (c) 2018, Michael Kazhdan, Alex Baden, and Keenan Crane
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

#include <cstdlib>
#include <vector>
#include <Misha/CmdLineParser.h>
#include <Misha/Ply.h>
#include <Misha/PlyVertexData.h>
#include <Misha/CmdLineParser.h>
#include <Misha/Spectrum.h>
#include <Misha/Miscellany.h>

Misha::CmdLineParameter< std::string >
	In( "in" ) ,				// Input mesh
	Out( "out" );				// Output spectrum
Misha::CmdLineParameter< int >
	SpectralDimension( "dim" , 200 );	// Maximum dimension of the spectrum
Misha::CmdLineParameter< float >
	Offset( "off" , 100.f );	// Offset used for the "shift" part of "invert-and-shift"
Misha::CmdLineReadable
	Lump( "lump" ) ,			// Should the mass matrix be lumped?
	Verbose( "verbose" );		// Should verbose output be provided?

void Usage( const char* ex ) 
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name.c_str() );
	printf( "\t[--%s <spectral dimension>=%d]\n" , SpectralDimension.name.c_str() );
	printf( "\t[--%s <output spectrum>]\n" , Out.name.c_str() );
	printf( "\t[--%s <offset>=%f]\n" , Offset.name.c_str() , Offset.value );
	printf( "\t[--%s]\n" , Lump.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}
Misha::CmdLineReadable* params[] = { &In , &SpectralDimension , &Out , &Verbose , &Lump , &Offset , NULL };


template< typename Real > using PlyVertexFactory = VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , 3 > >;
template< typename Real > using PlyVertex = typename PlyVertexFactory< Real >::VertexType;

template< typename Real >
void _main( int argc , char *argv[] )
{
	std::vector< TriangleIndex > triangles;
	std::vector< Point3D< Real > > vertices;
	//////////////////////
	// Read in the data //
	{
		int file_type;
		PlyVertexFactory< float > vFactory;
		std::vector< PlyVertex< float > > _vertices;
		PLY::ReadTriangles( In.value , vFactory , _vertices , triangles , NULL , file_type );
		vertices.resize( _vertices.size() );
		for( int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[i].get<0>() );
	}
	// Read in the data //
	//////////////////////

	Spectrum< Real > spectrum;
	spectrum.set( vertices , triangles , SpectralDimension.value , Offset.value , Lump.set );
	if( Out.set ) spectrum.write( Out.value );
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	Miscellany::Timer tmr;
	_main< double >( argc , argv );
	if( Verbose.set ) printf( "Got spectrum: %.2f(s)\n" , tmr.elapsed() );
	return EXIT_SUCCESS;
}