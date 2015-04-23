/* 
Copyright (c) 2015, Daniel C. Moura, Universidade do Porto (UP), Instituto de Telecomunicações (IT)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Universidade do Porto/Instituto de Telecomunicações
      nor the names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Daniel C. Moura, Universidade do Porto (UP) or 
Instituto de Telecomunicações (IT) BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
DAMAGE.
*/

/*
  3DHEB - 3D Histogram-based Edge Bundling
  Author: Daniel C. Moura (daniel.c.moura@gmail.com)

  If you use this code please cite:
  
  @ARTICLE{moura2015_3dheb,
     author = {{Moura}, Daniel C.},
      title = "{{3D} Density Histograms for Criteria-driven Edge Bundling}",
    journal = {ArXiv e-prints},
     eprint = {1504.02687},
       year = 2015,
      month = apr
  }

  Thank you.
*/

/*
  Node: Data structure for storing nodes of a graph.
*/

class Node {
  float xOrig, yOrig;
  float x, y;
  float xNorm, yNorm;
  String id;
  int idx;
  Graph gparent;
  
  
  
  Node(float xOrig, float yOrig, String id, int idx, Graph gparent) {
    this.xOrig = xOrig;
    this.yOrig = yOrig;
    this.id = id;
    this.idx = idx;
    this.gparent = gparent;
    this.x=-1;
    this.y=-1;    
  }
  
  
  void calcCoords(float minX, float maxX, float minY, float maxY) {
    xNorm = (xOrig-minX) / (maxX-minX);
    yNorm = (yOrig-minY) / (maxY-minY);    
    x =  ((Parameters.mapWidth-6) * xNorm)+3; //always ensure a 3px margin
    y =  ((Parameters.mapHeight-6) * yNorm)+3;
  }
  
  
  
  float dist2Node(Node a) {
    return dist(x,y,a.x,a.y);
  }
  
  
}
