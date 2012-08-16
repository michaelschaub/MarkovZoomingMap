###############################################################################
Copyright (C) 2012, M. Schaub, M. Barahona

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

###############################################################################

-----------------------------------------------------------------------------
Community Detection using Markov Zooming Infomap
-----------------------------------------------------------------------------

This code implements a MATLAB/C++ version of the Markov Zooming Map equation 
method as discussed in the article:

"Encoding dynamics for multiscale community detection: Markov time sweeping for 
the Map equation", M.T.Schaub, R. Lambiotte, M. Barahona, Physical Review E, 
to appear; see also arXiv:1109.6642

The code is based on the implementation of the original map equation method
which can be found under http://www.tp.umu.se/~rosvall/code.html
See also http://www.mapequation.org for a dynamic visualisation of the map 
equation. The original Map equation method is described in the following 
publication:

"Maps of information flow reveal community structure in complex networks"
Martin Rosvall and Carl T. Bergstrom, PNAS 105, 1118 (2008); 
see also arXiv:0707.0609

Please cite the above mentioned articles if you make use of this code.

-----------------------------------------------------------------------------
INSTALLATION & USE
-----------------------------------------------------------------------------
1. Download the files and unzip the folders it neccessary
2. in a terminal, go into the infomap_dir folder and run "make" to compile 
   the C code.
3. In Matlab you can now use the function "MarkovZoomingMap" to run your 
   analysis

EXAMPLE:

     % in the MarkovZoomingMap folder, load the provided example graph 
     % (a ring of rings) into an adjacency matrix
     A = convertPajekToAdjMatrix('ring_of_rings.net');
     
     % assign an output filename
     filename = 'test';

     % specify a time interval for the analysis
     time =logspace(-1,2,100);

     % run the actual analysis; note that the time argument is optional, if not
     % provided the time interval is set to the one also used in this example
     MarkovZoomingMap(A,filename,time)

This should create a folder named "{filename}ZoomingMap" containing the results
of the analysis in a file named Map_clustering.mat, plus a (directed) pajek .net
file of the original graph.

To plot the results of the analysis you may want to load the .mat file and use 
the provided matlab script "script_plot_Map_results.m".

-----------------------------------------------------------------------------
Some potential code changes / improvements
-----------------------------------------------------------------------------
This released code has not been thoroughly optimized for speed, however there
are two obvious ways to speed up the performance.

1) At the moment for each timepoint, a temporary pajek file is created and 
written to disk, and subseqently read from the C code. This I/O is essentially
not necessary. Instead on could pass the network data straight to the C 
implementation via e.g. a mex file from Matlab (However if you are interested 
in the "intermediate" graphs you may not want to do this, and instead change 
the script as to not delete the temporary files but keep them).

2) The stationary distribution of the graph is in this code computed by a
"matrix power method". Since the actual code as it stands deals with undirected
graphs, this is however not necessary. Instead one can easily compute the 
stationary distribution analytically. However to make the code more easily 
adaptable to the directed case, this computationally more expensive variant has 
been left in  the C code.

Further you may want to stop all the output from being displayed while the 
optimization is running, this is also easily achievable by changing the C code
accordingly.


-----------------------------------------------------------------------------
Author   : M. Schaub
Email    : michael.schaub09@imperial.ac.uk 
-----------------------------------------------------------------------------
