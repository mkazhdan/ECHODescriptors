<center><h2>Echo Descriptors (Version 1.0)</h2></center>
<center>
<a href="#LINKS">links</a>
<a href="#EXECUTABLES">executables</a>
<a href="#USAGE">usage</a>
<a href="#CHANGES">changes</a>
</center>
<hr>
<a name="LINKS"><b>LINKS</b></a><br>
<ul>
This code-base implements the ECHO descriptor presented in <A HREF="http://www.cs.jhu.edu/~misha/MyPapers/CGF21.pdf">ECHO: Extended Convolution Histogram of Orientations for Local Surface Description</A>.
<br>
<b>Executables: </b>
<a href="http://www.cs.jhu.edu/~misha/Code/ECHODescriptors/Version1.00/ECHODescriptors.x64.zip">Win64</a><br>
<b>Source Code:</b>
<a href="http://www.cs.jhu.edu/~misha/Code/ECHODescriptors/Version1.00/ECHODescriptors.zip">ZIP</a> <a href="https://github.com/mkazhdan/EchoDescriptors">GitHub</a><br>
<!--
<b>Older Versions:</b>
<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version1/">V1</a>
-->
</ul>
<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>
<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>GetSpectrum</b></font>:
Computes the spectrum of the Laplace-Beltrami operator by solving the associated generalized eigen-system.
</SUMMARY>
<dt><b>--in</b> &lt;<i>input triangle mesh</i>&gt;
</dt><dd> This string is the name of the file from which the triangle mesh will be read.<br>
The file is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format.
</dd>

</dd><dt>[<b>--out</b> &lt;<i>output spectrum file</i>&gt;]
</dt><dd> This string is the name of the file to which the computed spectrum will be written. 

</dd><dt>[<b>--dim</b> &lt;<i>spectral dimension</i>&gt;]
</dt><dd> This integer specifies the number of eigenvalue/eigenvector pairs to be computed.<br>
The default value for this parameter is 200.

</dd><dt>[<b>--off</b> &lt;<i>shift offset</i>&gt;]
</dt><dd> This floating point value specifies the offset to be used in the invert-and-shift implementation computing the lower-frequency part of the spectrum.<br>
The default value for this parameter is 100.

</dd><dt>[<b>--verbose</b>]
</dt><dd> Enabling this flag provides a more verbose gives the running time for computing the spectrum.

</dd>
</DETAILS>
</dl>
</ul>


<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>GetDescriptor</b></font>:
Computes the ECHO descriptor at the prescribed vertex.
</SUMMARY>
<dt><b>--in</b> &lt;<i>input triangle mesh</i>&gt;
</dt><dd> This string is the name of the file from which the triangle mesh will be read.<br>
The file is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format.

</dd><dt><b>--vertex</b> &lt;<i>vertex index</i>&gt;]
</dt><dd> This integer specifies the index of the vertex at which the ECHO descriptor should be computed.<BR>
If the prescribed value is negative, the code will compute descriptors at -(<i>vertex index</i>) random locations and will not write out the results.

</dd><dt>[<b>--spec</b> &lt;<i>spectral decomposition file</i>&gt;]
</dt><dd> This string is the name of the file from which the spectral decomposition of the triangle mesh will be read.<BR>
If no file specified, the code will first compute the spectral decomposition (using the lowest 200 eigenvector/eigenvalue pairs) and use that for processing.

</dd><dt>[<b>--out</b> &lt;<i>output ECHO descriptor</i>&gt;]
</dt><dd> This string is the name of the file to which the computed ECHO descriptor will be written.<BR>
If the extension of the filename is "<I>txt</I>" the descriptor will be written out in ASCII format. Otherwise, the descriptor will be written out as an image. (Currently, BMP, JPEG, PNG, and PBM formats are supported.)

</dd><dt>[<b>--distance</b> &lt;<i>distance type</i>&gt;]
</dt><dd> This integer specifies the type of distance used for computing the ECHO descriptor. Supported values are:
<UL>
<LI>0: geodesic
<LI>1: biharmonic
<LI>2: diffusion
<LI>3: commute
</UL>
<BR>
The default value for this parameter is 1 (i.e. biharmonic).

</dd><dt>[<b>--diffusion</b> &lt;<i>diffusion time</i>&gt;]
</dt><dd> If the diffusion distance is used as the distance, this floating point gives the diffusion time.<BR>
The default value for this parameter is 0.1.

</dd><dt>[<b>--tau</b> &lt;<i>radius scale</i>&gt;]
</dt><dd> This floating point value defines the support radius for computing the ECHO descriptor. Specifically, string is the name of the file to which the the octree and solution coefficients are to be written. Specifically, the radius of support is defined as:
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\tau\cdot\sqrt{\frac{|M|}{\pi}}" title="\Large \tau\cdot\sqrt{\frac{|M|}{\pi}}" HEIGHT="36"> where <I>|M|</I> is the area of the mesh.<BR>
The default value for this parameter is 0.08.

</dd><dt>[<b>--hRadius</b> &lt;<i>histogram radius</i>&gt;]
</dt><dd> This integer specifies the radius of the histogram used to discrtize the ECHO descriptor. (If the prescribed radius is <I>r</I> then the ECHO descriptor will be sampled on a (2<I>r</I>+1)x(2<I>r</I>+1) grid.<BR>
The default value for this parameter is 5.

</dd><dt>[<b>--resolution</b> &lt;<i>output resolution</i>&gt;]
</dt><dd> This integer specifies the resolution to which the ECHO descriptor will be resampled prior to output.<BR>
If no value is specfied, the resolution of the output will match the resolution of the histogram.

</dd><dt>[<b>--dev</b> &lt;<i>deviation for color mapping</i>&gt;]
</dt><dd> If the ECHO descriptor is written out as an image. If the prescribed deviation <I>dev</I> is negative, an ECHO value of <I>v</I> is computed to the color whose HSV representation is (0,0,<I>v</I>/(3x&sigma;)), where &sigma; is the standard deviation of ECHO values over the descriptor. If the prescribed deviation deviation is positiven, the HSV representation is (4&pi;/3x<I>dev</I>/&sigma;,1,<I>v</I>/(3x&sigma;)).<BR>
The default value for this parameter is -1.



</dd>
</DETAILS>
</dl>
</ul>


<hr>
<a name="USAGE"><b>USAGE EXAMPLES (WITH SAMPLE DATA)</b></a><br>

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>PoissonRecon / SSDRecon / SurfaceTrimmer / ChunkPly</b></font>
</SUMMARY>
For testing purposes, three point sets are provided:
<ol>

<li> <a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/horse.npts"><b>Horse</b></a>:
A set of 100,000 oriented point samples (represented in ASCII format) was obtained by sampling a virtual horse model with a sampling density proportional to curvature, giving a set of non-uniformly distributed points.<br>
The surface of the model can be reconstructed by calling the either Poisson surface reconstructor:
<blockquote><code>% PoissonRecon --in horse.npts --out horse.ply --depth 10</code></blockquote>
or the SSD surface reconstructor
<blockquote><code>% SSDRecon --in horse.npts --out horse.ply --depth 10</code></blockquote>
</li>

<li> <a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/bunny.points.ply"><b>Bunny</b></a>:
A set of 362,271 oriented point samples (represented in PLY format) was obtained by merging the data from the original Stanford Bunny
<a href="ftp://graphics.stanford.edu/pub/3Dscanrep/bunny.tar.gz">range scans</a>. The orientation of the sample points was estimated
using the connectivity information within individual range scans.<br>
The original (unscreened) Poisson reconstruction can be obtained by setting the point interpolation weight to zero:
<blockquote><code>% PoissonRecon --in bunny.points.ply --out bunny.ply --depth 10 --pointWeight 0</code></blockquote>
By default, the Poisson surface reconstructor uses degree-2 B-splines. A more efficient reconstruction can be obtained using degree-1 B-splines:
<blockquote><code>% PoissonRecon --in bunny.points.ply --out bunny.ply --depth 10 --pointWeight 0 --degree 1</code></blockquote>
(The SSD reconstructor requires B-splines of degree at least 2 since second derivatives are required to formulate the bi-Laplacian energy.)
</li>

<li> <a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/eagle.points.ply"><b>Eagle</b></a>:
A set of 796,825 oriented point samples with color (represented in PLY format) was obtained in the EPFL <a href="http://lgg.epfl.ch/statues.php">Scanning 3D Statues from Photos</a> course.<br>
A reconstruction of the eagle can be obtained by calling:
<blockquote><code>% PoissonRecon --in eagle.points.ply --out eagle.pr.ply --depth 10</code></blockquote>
(with the RGBA color properties automatically detected from the .ply header).<BR>
A reconstruction of the eagle that does not close up the holes can be obtained by first calling:
<blockquote><code>% SSDRecon --in eagle.points.ply --out eagle.ssd.ply --depth 10 --density</code></blockquote>
using the <b>--density</b> flag to indicate that density estimates should be output with the vertices of the mesh, and then calling:
<blockquote><code>% SurfaceTrimmer --in eagle.ssd.ply --out eagle.ssd.trimmed.ply --trim 7</code></blockquote>
to remove all subsets of the surface where the sampling density corresponds to a depth smaller than 7.<BR>
This reconstruction can be chunked into cubes of size 4&times;4&times;4 by calling:
<blockquote><code>% ChunkPly --in 1 eagle.ssd.trimmed.ply --out eagle.ssd.trimmed.chnks --width 4</code></blockquote>
which partitions the reconstruction into 11 pieces.

<li> <a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/torso.zip"><b>Torso</b></a>:
A set of 3,488,432 (torso.points.ply) and an envelope (torso.envelope.ply).<br>
A reconstruction of the torso that constrains the reconstruction to be contained within the envelope can be obtained by calling:
<blockquote><code>% PoissonRecon --in torso.points.ply --envelope torso.envelope.ply --out torso.pr.ply --depth 10</code></blockquote>
using the <b>--envelope</b> flag to specify the water-tight mesh constraining the reconstruction.<BR>
</li>

</ol>

</DETAILS>
</dl>
</ul>

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>PointInterpolant / AdaptiveTreeVisualization</b></font>
</SUMMARY>
For testing purposes, a pair of point-sets is provided:
<ol>

<li> <a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/quadratic.2D.fitting.samples"><b>fitting samples</b></a>:
A set of 1000 random 2D samples from within the square [-1,1,]x[-1,1] along with the evaluation of the quadratic <i>f(x,y)=x*x+y*y</i> at each sample point (represented in ASCII format).
<LI> <a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/quadratic.2D.evaluation.samples"><b>evaluation samples</b></a>:
A set of 4 2D positions at which the fit function is to be evaluated (represented in ASCII format).
</ol>

The function fitting the input samples can be by calling the point interpolant:
<blockquote><code>% PointInterpolant --inValues quadratic.2D.fitting.samples --tree quadratic.2D.tree --dim 2</code></blockquote>
Then, the reconstructed function can be evaluated at the evaluation samples by calling the adaptive tree visualization:
<blockquote><code>% AdaptiveTreeVisualization --in quadratic.2D.tree --samples quadratic.2D.evaluation.samples</code></blockquote>
This will output the evaluation positions and values:
<blockquote><CODE>0 0 1.33836e-05</CODE></blockquote>
<blockquote><CODE>0.5 0 0.25001</CODE></blockquote>
<blockquote><CODE>0.5 0.5 0.500006</CODE></blockquote>
<blockquote><CODE>2 2 nan</CODE></blockquote>
Note that because the (last) evaluation position (2,2) is outside the bounding box of the fitting samples, the function cannot be evaluated at this point and a value of "nan" is output.
</DETAILS>
</dl>
</ul>

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>ImageStitching</b></font>
</SUMMARY>
For testing purposes, two panoramas are provided: <a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Jaffa.zip"><b>Jaffa</b></a> (23794 x 9492 pixels) and <a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/OldRag.zip"><b>OldRag</b></a> (87722 x 12501 pixels).

A seamless panorama can be obtained by running:
<blockquote><code>% ImageSitching --in pixels.png labels.png --out out.png</code></blockquote>

</DETAILS>
</dl>
</ul>


<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>EDTInHeat / AdaptiveTreeVisualization</b></font>
</SUMMARY>
The Euclidean Distance Tranform of the reconstructed horse can be obtained by running:
<blockquote><code>% EDTInHeat --in horse.ply --out horse.edt --depth 9</code></blockquote>
Then, the visualization code can be used to extract iso-surfaces from the implicit function.<BR>
To obtain a visualization near the input surface, use an iso-value close to zero:
<blockquote><code>% AdaptiveTreeVisualization.exe --in horse.edt --mesh horse_0.01_.ply --iso 0.01 --flip</code></blockquote>
(By default, the surface is aligned so that the outward facing normal aligns with the negative gradient. Hence, specifying the <CODE>--flip</CODE> flag is used to re-orient the surface.)<BR>
To obtain a visualization closer to the boundary of the bounding box, use an iso-value close to zero:
<blockquote><code>% AdaptiveTreeVisualization.exe --in horse.edt --mesh horse_0.25_.ply --iso 0.25 --flip</code></blockquote>
(Since the default <CODE>--scale</CODE> is 2, a value of 0.25 should still give a surface that is contained within the bounding box.)<BR>
To obtain a sampling of the implicit function over a regular grid:
<blockquote><code>% AdaptiveTreeVisualization.exe --in horse.edt --grid horse.grid</code></blockquote>

</DETAILS>
</dl>
</ul>


<hr>
<DETAILS>
<SUMMARY>
<A NAME="CHANGES"><font size="+1"><b><B>HISTORY OF CHANGES</B></b></font></A>
</SUMMARY>
<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version3/">Version 3</a>:
<ol>
<li> The implementation of the <b>--samplesPerNode</b> parameter has been modified so that a value of "1" more closely corresponds to a distribution with one sample per leaf node.
</li><li> The code has been modified to support compilation under MSVC 2010 and the associated solution and project files are now provided. (Due to a bug in the Visual Studios compiler, this required modifying the implementation of some of the bit-shifting operators.)
</li></ol>
<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version4/">Version 4</a>:
<ol>
<li> The code supports screened reconstruction, with interpolation weight specified through the <b>--pointWeight</b> parameter.
</li><li> The code has been implemented to support parallel processing, with the number of threads used for parallelization specified by the <b>--threads</b> parameter.
</li><li> The input point set can now also be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, and the file-type is determined by the extension, so that the <b>--binary</b> flag is now obsolete.
</li><li> At depths coarser than the one specified by the value <b>--minDepth</b> the octree is no longer adaptive but rather complete, simplifying the prolongation operator.
</li></ol>
<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version4.5/">Version 4.5</a>:
<ol>
<li> The algorithmic complexity of the solver was reduced from log-linear to linear.
</li></ol>
<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version4.5/">Version 4.51</a>:
<ol>
<li> Smart pointers were added to ensure that memory accesses were in bounds.
</li></ol>
<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5/">Version 5</a>:
<ol>
<li> The <b>--density</b> flag was added to the reconstructor to output the estimated depth of the iso-vertices.
</li><li> The <i>SurfaceTrimmer</i> executable was added to support trimming off the subset of the reconstructed surface that are far away from the input samples, thereby allowing for the generation of non-water-tight surface.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.1/">Version 5.1</a>:
<ol>
<li> Minor bug-fix to address incorrect neighborhood estimation in the octree finalization.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.5a/">Version 5.5a</a>:
<ol>
<li> Modified to support depths greater than 14. (Should work up to 18 or 19 now.)
</li><li> Improved speed and memory performance by removing the construction of integral and value tables.
</li><li> Fixed a bug in Version 5.5 that used memory and took more time without doing anything useful.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.6/">Version 5.6</a>:
<ol>
<li> Added the <b>--normalWeight</b> flag to support setting a point's interpolation weight in proportion to the magnitude of its normal.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.7/">Version 5.7</a>:
<ol>
<li> Modified the setting of the constraints, replacing the map/reduce implementation with OpenMP atomics to reduce memory usage.
</li><li> Fixed bugs that caused numerical overflow when processing large point clouds on multi-core machines.
</li><li> Improved efficiency of the iso-surface extraction phse.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version5.71/">Version 5.71</a>:
<ol>
<li> Added the function <i>GetSolutionValue</i> to support the evaluation of the implicit function at a specific point.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6/">Version 6</a>:
<ol>
<li> Modified the solver to use Gauss-Seidel relaxation instead of conjugate-gradients at finer resolution.
</li><li> Re-ordered the implementation of the solver so that only a windowed subset of the matrix is in memory at any time, thereby reducing the memory usage during the solver phase.
</li><li> Separated the storage of the data associated with the octree nodes from the topology.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.1/">Version 6.1</a>:
<ol>
<li> Re-ordered the implementation of the iso-surface extraction so that only a windowed subset of the octree is in memory at any time, thereby reducing the memory usage during the extracted phase.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.11/">Version 6.11</a>:
<ol>
<li> Fixed a bug that created a crash in the evaluation phase when <b>--pointWeight</b> is set zero.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.12/">Version 6.12</a>:
<ol>
<li> Removed the OpenMP <i>firstprivate</i> directive as it seemed to cause trouble under Linux compilations.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.13/">Version 6.13</a>:
<ol>
<li> Added a <b>MemoryPointStream</b> class in <i>PointStream.inl</i> to support in-memory point clouds.
</li><li> Modified the signature of <u>Octree::SetTree</u> in <i>MultiGridOctreeData.h</i> to take in a pointer to an object of type <b>PointStream</b> rather than a file-name.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version6.13a/">Version 6.13a</a>:
<ol>
<li> Modified the signature of <u>Octree::SetIsoSurface</u> to rerun a <i>void</i>. [<a href="http://www.danielgm.net/cc/">cloudcompare</a>]
</li><li> Added a definition of <u>SetIsoVertexValue</u> supporting double precision vertices. [<a href="http://www.danielgm.net/cc/">cloudcompare</a>]
</li><li> Removed <i>Time.[h/cpp]</i> from the repository. [<a href="http://www.danielgm.net/cc/">cloudcompare</a>/<a href="http://asmaloney.com/">asmaloney</a>]
</li><li> Fixed assignment bug in <u>Octree::SetSliceIsoVertices</u>. [<a href="http://asmaloney.com/">asmaloney</a>]
</li><li> Fixed initialization bug in <u>SortedTreeNodes::SliceTableData</u> and <u>SortedTreeNodes::XSliceTableData</u>. [<a href="http://asmaloney.com/">asmaloney</a>]
</li><li> Included <i>stdlib.h</i> in <i>Geometry.h</i>. [<a href="http://asmaloney.com/">asmaloney</a>]
</li><li> Fixed default value bug in declaration of <u>Octree::SetTree</u>. [<a href="http://asmaloney.com/">asmaloney</a>]
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version7.0/">Version 7.0</a>:
<ol>
<li> Added functionality to support color extrapolation if present in the input.
</li><li> Modified a bug with the way in which sample contributions were scaled.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version8.0/">Version 8.0</a>:
<ol>
<li> Added support for different degree B-splines.
(Note that as the B-spline degree is a template parameter, only degree 1 through 4 are supported.
If higher order degrees are desired, additional template parameters can be easily added in the body of the <u>Execute</u> function inside of <i>PoissonRecon.cpp</i>.
Similarly, to reduce compilation times, support for specific degrees can be removed.)
</li><li> Added the <b>--primalGrid</b> flag to support to extraction of a grid using primal sampling.
</li><li> Changed the implementation of the grid sampling so that computation is now linear, rather than log-linear, in the number of samples.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version9.0/">Version 9.0</a>:
<ol>
<li> Added support for free boundary conditions.
</li><li> Extended the solver to support more general linear systems. This makes it possible to use the same framework to implement the <a href="http://mesh.brown.edu/ssd/">Smoothed Signed Distance Reconstruction</a> of Calakli and Taubin (2011).
</li><li> Modified the implementation of density estimation and input representation. This tends to define a slightly larger system. On its own, this results in slightly increased running-time/footprint for full-res reconstructions, but provides a substantially faster implementation when the output complexity is smaller than the input.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version9.01/">Version 9.01</a>:
<ol>
<li> Reverted the density estimation to behave as in Version 8.0.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version9.01/">Version 9.011</a>:
<ol>
<li> Added a parameter for specifying the temporary directory.
</li></ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.00/">Version 10.00</a>:
<ol>
<li> The code has been reworked to support arbitrary dimensions, finite elements of arbitrary degree, generally SPD systems in the evaluated/integrated values and derivatives of the functions, etc.</LI>
<LI> For the reconstruction code, added the <B>--width</B> flag which allows the system to compute the depth of the octree given a target depth for the finest resolution nodes.</LI>
<LI> For the reconstruction code, fixed a bug in the handling of the confidence encoded in the lengths of the normals. In addition, added the flags <B>--confidence</B> and <B>--confidenceBias</B> which allow the user more control of how confidence is used to affect the contribution of a sample.</LI>
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.01/">Version 10.01</a>:
<ol>
<li> Modified the reconstruction code to facilitate interpolation of other input-sample quantities, in addition to color.</LI>
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.02/">Version 10.02</a>:
<ol>
<li> Set the default value for <b>--degree</B> in PoissonRecon to 1 and change the definitiion of <I>DATA_DEGREE</I> to 0 for sharper color interpolation.</LI>
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.03/">Version 10.03</a>:
<ol>
<li> Cleaned up memory leaks and fixed a bug causing ImageStitching and EDTInHeat to SEGFAULT on Linux.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.04/">Version 10.04</a>:
<ol>
<li> Replaced the ply I/O code with an object-oriented implementation.
<LI> Updated the code to support compilation under gcc version 4.8.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.05/">Version 10.05</a>:
<ol>
<LI> Added cleaner support for warning and error handling.
<LI> Minor bug fixes.
<LI> Added a <B>--inCore</B> flag that enables keeping the pointset in memory instead of streaming it in from disk.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.06/">Version 10.06</a>:
<ol>
<LI> Improved performance.
<LI> Modified <CODE>PoissonRecon</CODE> and <CODE>SSDRecon</CODE> to support processing of 2D point sets.
<LI> Modified the 2D implementations of <CODE>PoissonRecon</CODE>, <CODE>SSDRecon</CODE>, and <CODE>AdaptiveTreeVisualization</CODE> to support ouput to <CODE>.jpg</CODE> and <CODE>.png</CODE> image files.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version10.07/">Version 10.07</a>:
<ol>
<LI> Removed a bug that would cause memory access errors when some slices were empty.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version11.00/">Version 11.00</a>:
<ol>
<LI> Added support for processing point-sets so large that 32-bit indices for octrees are not sufficient. (Enabled by defining the preprocessor variable <B>BIG_DATA</B> in the file <I>PreProcessor.h</I>.
<LI> Added C++11 parallelism for compilers that do not support OpenMP.
<LI> Added the code for <I>ChunkPly</I> which breaks up large meshes and/or point-sets into chunks.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version11.01/">Version 11.01</a>:
<ol>
<LI> Fixed bug with <I>_mktemp</I> that caused the code to crash on Windows machine with more than 26 cores.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version11.02/">Version 11.02</a>:
<ol>
<LI> Added error handling for numerical imprecision issues arrising when too many samples fall into a leaf node.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version12.00/">Version 12.00</a>:
<ol>
<LI> Added functionality enabling <I>AdaptiveTreeVisualization</I> to output the values of a function at prescribed sample positions.
<LI> Added the implementation of <I>PointInterpolant</I> that fits a function to a discrete set of sample values.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version13.00/">Version 13.00</a>:
<ol>
<LI> Enabled passing in a constraint envelope to <I>PoissonRecon</I>, allowing one to define a region that is known to be outside the surface.
<LI> Updated <I>ChunkPLY</I> to support processing of input points in either ASCII or binary format.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version13.50/">Version 13.50</a>:
<ol>
<LI> Enabled support for automatically detecting attirbutes of input point cloud (in addition to positions and normals) when provided in .ply format.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version13.60/">Version 13.60</a>:
<ol>
<LI> Modified the implementation of <I>PointInterpolant</I> to support separately prescribing value and gradient constraints.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version13.61/">Version 13.61</a>:
<ol>
<LI> Bug fix addressing the problem that the memory for a <CODE>DynamicFactory</CODE> object is dynamically allocated and not only known at construction time.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version13.70/">Version 13.70</a>:
<ol>
<LI> Using the updated <A HREF="https://www.cc.gatech.edu/~turk/ply.tar.gz">PLY libraray</A> with the less restrictive BSD license.
</ol>

<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version13.70/">Version 13.71</a>:
<ol>
<LI> Fixed a bug that resulted in incorrect point weighting when the samples were too sparse.
</ol>

</DETAILS>


<hr>
<a name="SUPPORT"><b>SUPPORT</b></a><br>
This work genersouly supported by NSF grants #0746039 and #1422325.
