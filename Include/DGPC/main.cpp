#include <iostream>

#include "Generator.h"
#include "Mesh.h"
#include <iostream>
#include <limits>

void usage_error(const char*);

int main(int argc, char** argv) {

  //Parse options
  if(argc < 3) usage_error(argv[0]);
  const char* filename = argv[1];
  const int source_idx = atoi(argv[2]);

  double stopdist = -1;
  if(argc > 3) {
    stopdist = atof(argv[3]);
  }
  double epsilon = -1;
  if(argc > 4) {
    epsilon = atof(argv[4]);
  }

  if(argc > 5) usage_error(argv[0]);

  //Declare point, mesh and generator
  typedef DGPC::Vector3<double> Point;
  typedef DGPC::MeshOM<Point> Mesh;
  typedef DGPC::Generator<Mesh> DGPCgenerator;

  //Read mesh from file
  Mesh my_mesh;
  my_mesh.openOBJ(filename);
  
  //Make a DGPC generator
  DGPCgenerator my_dgpc(my_mesh);

  //Set options
  if(stopdist > 0) 
    my_dgpc.setStopDist(stopdist);

  if(epsilon > 0)
    my_dgpc.setEps(epsilon);

  //Set source node
  my_dgpc.setNodeSource(source_idx);

  //Optionally, the library also supports to set the source point:
  //Point p = Point(0,0,0);  //R3 coordinate of source
  //int face_idx = 0;                 //Face index of mesh where point lies
  //my_dgpc.setSource(point, face_idx);

  //Compute DGPC
  int last_node = my_dgpc.run();

  //Fetch and print result
  std::cout << "Computed distances until node " << last_node << std::endl;
  std::cout << std::endl;

  std::cout << "i      r      theta" << std::endl;
  std::cout << "-------------------" << std::endl;

  for(int i = 0; i < my_mesh.n_vertices(); i++) {
    const double r = my_dgpc.getDistance(i);
    if(r < std::numeric_limits<double>::max()) {
      const double theta = my_dgpc.getAngle(i);
      std::cout << i << "    " << r << "    "<< theta << std::endl;
    }
  }

}

void usage_error(const char* progname) {
  using namespace std;
  cout << "Usage: " << progname << " filename source_idx [stopdist] [epsilon]" << endl;
  cout << endl;
  cout << "Examples: " << endl;
  cout << endl;
  cout << "  " << "Compute DGPC for an entire mesh with source at node 0:" << endl;
  cout << "  " << progname << " hat6.obj 0" << endl;
  cout << endl;
  cout << "  " << "Compute DGPC for for nodes within geodesic radius 1.5 from the source at node 0:" << endl;
  cout << "  " << progname << " hat6.obj 0 1.5" << endl;
  cout << endl;
  cout << "  " << "Compute DGPC for for nodes within geodesic radius 1.5 from the source at node 0, with epsilon 1e-8:" << endl;
  cout << "  " << progname << " hat6.obj 0 1.5 1e-8" << endl;
  cout << endl;
  exit(-1);
}
