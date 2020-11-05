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
</dt><dd> Enabling this flag gives the running time for computing the spectrum.

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
</dt><dd> This integer specifies the radius of the histogram used to discrtize the ECHO descriptor. (If the prescribed radius is <I>r</I> then the ECHO descriptor will be sampled on a (2<I>r</I>+1)&times;(2<I>r</I>+1) grid.<BR>
The default value for this parameter is 5.

</dd><dt>[<b>--resolution</b> &lt;<i>output resolution</i>&gt;]
</dt><dd> This integer specifies the resolution to which the ECHO descriptor will be resampled prior to output.<BR>
If no value is specfied, the resolution of the output will match the resolution of the histogram.

</dd><dt>[<b>--dev</b> &lt;<i>deviation for color mapping</i>&gt;]
</dt><dd> If the ECHO descriptor is written out as an image. If the prescribed deviation <I>dev</I> is negative, an ECHO value of <I>v</I> is computed to the color whose HSV representation is (0,0,<I>v</I>/(3x&sigma;)), where &sigma; is the standard deviation of ECHO values over the descriptor. If the prescribed deviation deviation is positiven, the HSV representation is (4&pi;/3&times;<I>dev</I>/&sigma;,1,<I>v</I>/(3&times;&sigma;)).<BR>
The default value for this parameter is -1.

</dd><dt>[<b>--verbose</b>]
</dt><dd> Enabling this flag provides a break-down of the running times for the different steps of computing the ECHO descriptor.


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
<font size="+1"><b>GetSpectrum / GetDescriptor</b></font>
</SUMMARY>
For testing purposes, we provide the <A HREF="http://www.cs.jhu.edu/~misha/Code/ECHODescriptors/wolf0.ply"><I>wolf0</I></A> model from the publicly available <A HREF="http://tosca.cs.technion.ac.il/book/resources_data.html">TOSCA dataset</A>.

The spectrum of the Laplace-Beltrami operator can be computed by calling:
<blockquote><code>% GetSpectrum --in wolf0.ply --out wolf0.spec --verbose</code></blockquote>
This will output the spectrum to the file <I>wolf0.spec</I> and provide the running time for computing the spectrum.<BR>

A visualization of the descriptor at vertex 1000 of the mesh can be obtained by calling:
<blockquote><code>% GetDescriptor --in wolf0.ply --spect wolf0.spec --out wolf0.1000.jpeg --hRadius 10 --resolution 1024 --dev 0.005--verbose</code></blockquote>
This produces a 1024&times;1024 JPEG (color) image visualizing the ECHO descriptor computed over a histogram of size 21&times;21. Running times for the individual steps of the computation are written out to the command prompt. (Note that as the spectrum is provided as input the time for obtaining the spectrum is just the time required to read it from disk.)<BR>

</DETAILS>
</dl>
</ul>

<hr>
<DETAILS>
<SUMMARY>
<A NAME="CHANGES"><font size="+1"><b><B>HISTORY OF CHANGES</B></b></font></A>
</SUMMARY>
<a href="http://www.cs.jhu.edu/~misha/Code/ECHODescriptors/Version1.00/">Version 1.00</a>:
<ol>
<li> The original release of the source code.
</li></ol>

</DETAILS>


<hr>
<a name="SUPPORT"><b>SUPPORT</b></a><br>
This work genersouly supported by NSF grants #0746039 and #1422325.
