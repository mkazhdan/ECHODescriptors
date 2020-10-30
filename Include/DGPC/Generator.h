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
#ifndef DGPC_GENERATOR_H
#define DGPC_GENERATOR_H

#include "Heap.h"
#include "Mesh.h"

#include <limits>
#include <algorithm>

namespace DGPC {

  /**
   * A class to generate DGPC: Discrete Geodesic Polar Coordinates on
   * polygonal meshes. The class is templetized with a halfedge mesh
   * datastructure, and has been tested with OpenMesh. It should
   * probably also work on a CGAL::Polyhedron_3 without too much
   * effort.
   */
  template<class Mesh>
    class Generator {

    typedef typename Mesh::point_type Point;
    typedef typename Point::value_type real;

    typedef typename Mesh::FaceHandle FaceHandle;
    typedef typename Mesh::HalfedgeHandle HalfedgeHandle;
    typedef typename Mesh::VertexHandle VertexHandle;


  public:
    /**
     * Construct a Generator for a Mesh.
     */
  Generator(const Mesh& mesh) :
    mesh_(mesh) {
      eps_ = 1e-12;
      stopdist_ = std::numeric_limits<real>::max();
      const int n = (int)mesh_.n_vertices();
      distances_.resize(n);
      angles_.resize(n);
    };

    /**
     * Set epsilon. The algorithm will skip iterations which are not
     * significant for accuracy less than the given epsilon.
     */
    void setEps(real eps) { eps_ = eps; };
    /**
     * Set stop distance, geodesic radius of patch to compute.
     */
    void setStopDist(real d) { stopdist_ = d; };

    /**
     * Set source point. The point is assumed to be either one of the
     * nodes of face_idx, or on one of the faces.
     */
    void setSource(const Point& source, int face_idx);

    /**
     * Set source point on node node_idx.
     */
    void setNodeSource(int node_idx);

    /**
     * Set source on point, which lies on the face face_idx.
     */
    void setFaceSource(const Point& point, int face_idx);

    // Set x-axis for angle calculation
    void setAxes (Point3D<double> x, Point3D<double> y);

    /**
     * Start generation of DGPC. When complete, distances and angles
     * for the nodes in a geodesic disk with radius stopdist_ will be
     * available with getDistance(ni) and getAngle(ni).
     */
    int run();

    /**
     * Get Gamma, the smallest ring of nodes connected by edges
     * surrounding the source point.  These nodes are initialized with
     * angles and distance after setSource is called.
     */
    const std::vector<int>& getGamma() { return gamma_; };

    /**
     * Get DGPC polar distance for node ni
     */
    real getDistance(int ni) { return distances_[ni]; };

    /**
     * Get DGPC polar angle for node ni
     */
    real getAngle(int ni) { return angles_[ni]; };

    /**
     * Get DGPC polar distances
     */
    const std::vector<real> getDistances() { return distances_; };

    /**
     * Get DGPC polar angles
     */
    const std::vector<real>& getAngles() { return angles_; };

    // Get index of direction of tangent vector from source node
    int getTangentIndex (void) const 
    {
       return gamma_[0];
    }


  protected:
    const Mesh& mesh_;
    real eps_;
    real stopdist_;

    DGPC::Heap<real> heap_;

    std::vector<real> distances_;
    std::vector<real> angles_;

    std::vector<int> gamma_;

    void initialize();

    real initializeGamma(const Point& point);

    bool tryComputeNodeFromEdge(int node, int edge[2]);
    real computeDistance(const Point& pt, int edge[2], real& alpha);
    real computeAngle(int node, int edge[2], real alpha);

    Point xAxis, yAxis;
    bool axisSet = false;

  };

  //Implementation of setSource
  template<class Mesh>
    void
    Generator<Mesh>::setSource(const Point& point, int face_idx)
    {

      const real proximity_threshold = 10e-5;

      //Fetch nodes of the face
      std::vector<int> nodes;
      FaceHandle face = mesh_.face_handle(face_idx);
      HalfedgeHandle heh = mesh_.halfedge_handle(face);
      HalfedgeHandle start = heh;
      do {
        VertexHandle vh = mesh_.to_vertex_handle(heh);
        nodes.push_back(vh.idx());
        heh = mesh_.next_halfedge_handle(heh);
      } while(heh != start);

      //Is the source on a node?
      for(int i = 0; i < nodes.size(); i++) {
        VertexHandle vh = mesh_.vertex_handle(nodes[i]);
        const Point& np = mesh_.point(vh);
        
        if(np.dist(point) < proximity_threshold) {
          setNodeSource(nodes[i]);
          return;
        }
      }

      //Assume the source is on the face
      setFaceSource(point, face_idx);
      return;

    }

  //Implementation of setNodeSource
  template<class Mesh>
    void
    Generator<Mesh>::setNodeSource(int node_idx)
    {
      //Clear distances, angles and gamma
      initialize();

      //Initialize source node
      distances_[node_idx] = 0;
      angles_[node_idx] = 0;

      //Find gamma, walk along the 1-ring around source
      VertexHandle source = mesh_.vertex_handle(node_idx);
      HalfedgeHandle heh = mesh_.halfedge_handle(source);
      
      if(mesh_.is_boundary(source)) {
        //Skip anticlockwise around source until heh is the last non-boundary halfedge
        HalfedgeHandle b = mesh_.opposite_halfedge_handle(heh);
        while(!mesh_.is_boundary(b)) {
          heh = mesh_.next_halfedge_handle(b);
          b = mesh_.opposite_halfedge_handle(heh);
        }
      }

      HalfedgeHandle start = heh;
      VertexHandle to;

      //Traverse all halfedges pointing into source
      do {
        heh = mesh_.next_halfedge_handle(heh);
        to = mesh_.to_vertex_handle(heh);

        //Traverse all nodes on the edge of this face, except source and first anticlockwise neighhbour
        while (to != source) {
          gamma_.push_back(to.idx());
          heh = mesh_.next_halfedge_handle(heh);
          to = mesh_.to_vertex_handle(heh);
        } 
        //heh is now pointing to source
        heh = mesh_.opposite_halfedge_handle(heh);

      } while(heh != start);

      Point source_pt = mesh_.point(source);

      //Initialize gamma with distances and angles
      real phitot = initializeGamma(source_pt);

      if(!mesh_.is_boundary(source)) {
        //Scale angles to sum to 2pi
        const real alpha = (2*M_PI)/phitot;
        const int num = (int)gamma_.size();
        for( int i = 0; i < num; i++) {
          //Store the angle for this node
          angles_[gamma_[i]] *= alpha;
        }
      }

    }

  //Implementation of setFaceSource
  template<class Mesh>
    void
    Generator<Mesh>::setFaceSource(const Point& point, int face_idx)
    {

      //Clear distances, angles and gamma
      initialize();

      //Find gamma, the nodes of this face.
      FaceHandle face = mesh_.face_handle(face_idx);
      HalfedgeHandle heh = mesh_.halfedge_handle(face);
      HalfedgeHandle start = heh;
      do {
        VertexHandle vh = mesh_.to_vertex_handle(heh);
        int ni = vh.idx();
        gamma_.push_back(ni);
        heh = mesh_.next_halfedge_handle(heh);
      } while(heh != start);

      //Initialize gamma with distances and angles
      initializeGamma(point);
    }


  //Implementation of run
  template<class Mesh>
    int
    Generator<Mesh>::run()
    {

      int last_finished = -1;

      int edges[3];
      std::vector<int> next;

      HalfedgeHandle heh, end;

      while (!heap_.empty()) {
    
        int curr = heap_.getCandidate();
        if (curr == -1) break;
    
        last_finished = curr;

        VertexHandle curr_vertex = mesh_.vertex_handle(curr);

        //Iterate halfedges pointing into current vertex (one for each
        //face adjacent to curr)
        HalfedgeHandle face_start = mesh_.halfedge_handle(curr_vertex);
        HalfedgeHandle face = face_start;

        do {
          face = mesh_.opposite_halfedge_handle(face);
          if(!mesh_.is_boundary(face)) {
            heh = mesh_.prev_halfedge_handle(face);
            end = heh;

            //For this face, we will attempt to compute DGPC from each
            //of the two edges connected to source
            edges[0] = mesh_.to_vertex_handle(heh).idx();
            heh = mesh_.next_halfedge_handle(heh);
            edges[1] = mesh_.to_vertex_handle(heh).idx();
            heh = mesh_.next_halfedge_handle(heh);
            edges[2] = mesh_.to_vertex_handle(heh).idx();
            heh = mesh_.next_halfedge_handle(heh);

            assert(edges[1] == curr);

            //We can now attempt to compute DGPC from the two edges
            // [edges[0], edges[1]] and [edges[1], edges[2]]

            //We will attempt to compute DGPC for all nodes in this
            //face (except source). Build a list of the nodes in "next".
            next.clear();
            next.push_back(edges[2]);

            while(heh != end) {
              next.push_back(mesh_.to_vertex_handle(heh).idx());
              heh = mesh_.next_halfedge_handle(heh);
            }

            next.push_back(edges[0]);
      
            for(int i = 0; i < next.size(); i++) {

              int n = next[i];

              if( n != edges[0] ) {
                //Compute distance to n over the edge [edges[0], edges[1]]
                tryComputeNodeFromEdge(n, &edges[0]);
              }

              if( n != edges[2] ) {
                //Compute distance to n over the edge [edges[1], edges[2]]
                tryComputeNodeFromEdge(n, &edges[1]);
              }
            }
          }
          face = mesh_.next_halfedge_handle(face);
        } while(face != face_start);
      }
  
      return last_finished;
    }

    template<class Mesh>
    void Generator<Mesh>::setAxes (Point3D<double> x, Point3D<double> y)
    {
       axisSet = true;
       
       double xNorm = std::sqrt (x.squareNorm ()), yNorm = std::sqrt (y.squareNorm ());

       for (int i = 0; i < 3; i++)
       {
          xAxis[i] = x[i] / xNorm;
          yAxis[i] = y[i] / yNorm;
       }

    }


  ////////////////////////////////////////////
  // Implementation of protected methods below
  ////////////////////////////////////////////

  // Implementation of initialize
  template<class Mesh>
    void
    Generator<Mesh>::initialize()
    {
      std::fill(distances_.begin(), distances_.end(), std::numeric_limits<real>::max());
      heap_.initialize(&distances_);
      gamma_.clear();
    }

  template<class Mesh>
    typename Generator<Mesh>::real
    Generator<Mesh>::initializeGamma(const Point& point)
    {

      const int num = gamma_.size();
      real phitot = 0;

      //For each node in gamma_
      // * Compute distances from point, store in distances_
      // * Compute angles spanned in point, store in angles_
      // * Insert node in heap
      // return total sum of angles spanned in point.

      if (!axisSet)
      {

         for(int i = 0; i < num; i++) 
         {
           int ni = gamma_[i];

           VertexHandle nbvh = mesh_.vertex_handle(ni);
           const Point& nb = mesh_.point(nbvh);
           real dist = (point-nb).length();
           distances_[ni] = dist;

           int ip = i+1;
           if(ip >= num) ip = 0;
           int nip = gamma_[ip];
           VertexHandle nbvhp = mesh_.vertex_handle(nip);
           const Point& nbp = mesh_.point(nbvhp);

           const Point nb_t =  (nb  - point).normalize();
           const Point nbp_t = (nbp - point).normalize();
           real cos_phi = nb_t * nbp_t;
           real phi = acos(cos_phi);

           angles_[ni] = phitot;

           heap_.push(ni);

           phitot += phi;

         }
      }
      else
      {
         phitot = 0.0;

         for (int i = 0; i < num; i++)
         {
            int ni = gamma_[i];

            VertexHandle nbvh = mesh_.vertex_handle(ni);
            const Point& nb = mesh_.point(nbvh);

            real dist = (point - nb).length ();
            distances_[ni] = dist;

            const Point vec = nb - point;

            real ang = std::atan2 (yAxis * vec, xAxis * vec);

            if (ang < 0) { ang += 2.0 * M_PI;}

            angles_[ni] = ang;

            heap_.push (ni);

            if (ang > phitot)
            {
               phitot = ang;
            }
        }

     }
        
      return phitot;
    }

  
  //Implementation of tryComputeNodeFromEdge
  template<class Mesh>
    bool
    Generator<Mesh>::tryComputeNodeFromEdge(int node, int edge[2])
    {
      real  thresh = 1.0+eps_;
      real alpha;

      VertexHandle h = mesh_.vertex_handle(node);
      const Point& pt = mesh_.point(h);

      real newdist = computeDistance(pt, edge, alpha); 

      if (distances_[node]/newdist > thresh) {
        //Store new distance, and compute angle
        distances_[node] = newdist;
        angles_[node] = computeAngle(node, edge, alpha);
	    
        if(newdist < stopdist_) {
          heap_.push(node);
        }
        return true;
      }
      return false;
    }

  //Implementation of computeNodeFromEdge
  template<class Mesh>
    typename Generator<Mesh>::real
    Generator<Mesh>::computeDistance(const Point& pt, int edge[2], real& alpha)
    {

      const Point& Nk = mesh_.point(mesh_.vertex_handle(edge[0]));
      const Point& Nj = mesh_.point(mesh_.vertex_handle(edge[1]));

      const real Uk = distances_[edge[0]];
      const real Uj = distances_[edge[1]];

      const real djptsq = Nj.dist2(pt);
      const real djpt = sqrt(djptsq);
      const real dkptsq = Nk.dist2(pt);
      const real dkpt = sqrt(dkptsq);

      const Point ekj = Nk-Nj;
      const real djksq = ekj.length2();
      const real djk = sqrt(djksq);

      //Stable evaluation of Herons formula
      using namespace std; //For max and min
      const real a = max(djk, max(Uj, Uk));
      const real c = min(djk, min(Uj, Uk));
      real b;
      if(a == djk) {
        b = max(Uj, Uk);
      } else if(a == Uj) {
        b = max(djk, Uk);
      } else if(a == Uk) {
        b = max(djk, Uj);
      }

      const real H_under_root = ( (a + (b+c)) *
                                  (c - (a-b)) *
                                  (c + (a-b)) *
                                  (a + (b-c)) );

      if(H_under_root < 0 || djk < 10e-12) {
        // Triangle inequality fails, return Dijkstra instead
        const real dijkstra_j = Uj + Nj.dist(pt);
        const real dijkstra_k = Uk + Nk.dist(pt);
        if(dijkstra_j < dijkstra_k) {
          alpha = 0;
          return dijkstra_j;
        } else {
          alpha = 1;
          return dijkstra_k;
        }
      }


      const real H = sqrt( H_under_root );

      const Point ej = Nj-pt;
      const Point ek = Nk-pt;

      const real A2 = ej.crossProd(ek).length();

      const real ej_ekj = ej * ekj;
      const real ek_ekj = ek * ekj;

      const real f31 = djksq - (Uj-Uk)*(Uj-Uk);
      const real f32 = (Uj+Uk)*(Uj+Uk) - djksq;

      const real f3 = f31*f32; // ( djksq - (Uj-Uk)*(Uj-Uk) ) * ( (Uj+Uk)*(Uj+Uk) - djksq );

      const real Ujsq = Uj*Uj;
      const real Uksq = Uk*Uk;

      const real f1_j = A2 * (djksq + Uksq - Ujsq);
      const real f1_k = A2 * (djksq + Ujsq - Uksq);

      const real xj = (f1_j + ek_ekj*H);
      const real xk = (f1_k - ej_ekj*H);

      if(xj < 0 || xk < 0) {
        // Update from outside triangle, return Dijkstra instead
        const real dijkstra_j = Uj + Nj.dist(pt);
        const real dijkstra_k = Uk + Nk.dist(pt);
        if(dijkstra_j < dijkstra_k) {
          alpha = 0;
          return dijkstra_j;
        } else {
          alpha = 1;
          return dijkstra_k;
        }
      }

      const real f4 = 2*A2*djksq;

      real Ui = sqrt(xj*xj*djptsq + 2*xj*xk*ej*ek + xk*xk*dkptsq)/f4;

      const real cos_jk = (  Ujsq+Uksq - djksq)/(2*Uj*Uk);
      const real cos_ji  = ( Ujsq+Ui*Ui - djptsq)/(2*Uj*Ui);
      alpha = acos(cos_ji) /  acos (cos_jk);

      return Ui;

    }

  template<class Mesh>
    typename Generator<Mesh>::real
    Generator<Mesh>::computeAngle(int node, int edge[2], real alpha)
    {
      real nkphi = angles_[edge[0]];
      real njphi = angles_[edge[1]];

      const real diff = fabs(njphi-nkphi);

      if(diff < eps_)
        return njphi;

      if(diff > M_PI) {
        //Make the interpolation modulo 2pi
        if(njphi < nkphi) njphi += 2*M_PI;
        else nkphi += 2*M_PI;
      }
      
      real angle = (1-alpha) * njphi  + alpha * nkphi;
      
      if(angle > 2*M_PI)
        angle -= 2*M_PI;
      
      return angle;
    }

}; //End namespace DGPC

#endif //DGPC_GENERATOR_H
