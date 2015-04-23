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
  Edge: Data structure that connects two nodes of a graph.
*/

class Edge {
  int sourceIdx, targetIdx;   //nodes indexes in the graph array
  Node source, target;        //nodes objects
  float weightOrig, weight;   //original weight, and normalized weight
  String id;
  Graph gparent;

  //coordinates of the edge samples
  float[] x;
  float[] y;
  
  //destinace between origin and destination
  float distanceOrig, distanceNorm;

  //instead of usig a bin for all edge it is used one bin per line segment between samples
  int[] bins;
  
  boolean firstResample=true;
  boolean firstBining = true;

  int clusterIdx=0; //in case we use clustering at the edge level, this defines the bin for all line segments of this edge

  Edge(Node source, Node target, float weightOrig, String id, Graph gparent) {
    this.source = source;
    this.target = target;
    sourceIdx = source.idx;
    targetIdx = target.idx;
    this.weightOrig = weightOrig;    
    this.id = id;    
    this.gparent = gparent;   
    this.distanceOrig = dist(source.xOrig, source.yOrig, target.xOrig, target.yOrig); //!!        
  }

  void normalize(float minW, float maxW, float minD, float maxD) {    
    //weight = (weightOrig-minW) / (maxW-minW);
    weight = (weightOrig) / (maxW);
    distanceNorm = (distanceOrig-minD) / (maxD-minD);
  }


  void initControlPoints() {    
    final float minSpace = (Parameters.removeDist*1 + Parameters.splitDist*1)/2;     
    final int n = (int) max(2, dist(source.x, source.y, target.x, target.y) / minSpace);
    final float incX = (target.x-source.x)*1.0f/n;
    final float incY = (target.y-source.y)*1.0f/n;

    //println("N CPs: " + n);

    x = new float[n+1];
    x[0] = source.x;    
    for (int i=1; i<n; i++)
      x[i] = x[i-1] + incX;
    x[n] = target.x;

    y = new float[n+1];      
    y[0] = source.y;  
    for (int i=1; i<n; i++)
      y[i] = y[i-1] + incY;        
    y[n] = target.y;

    bins = new int[n]; //even if not needed it is set to zero
  }

  //if we want, we may separate edges based on their orientation
  void translateControlPointsByNormal() {
    float dx = target.x-source.x;
    float dy = target.y-source.y;
    float dm = mag(dx, dy);
    if (dm == 0.0f)
      return;

    dx /= dm;
    dy /= dm;

    final float nx = -dy * Parameters.separationFactor;
    final float ny =  dx * Parameters.separationFactor;
    for (int i=1; i<x.length-1; i++) 
      x[i] = gparent.kde[0].correctX(x[i]+nx);       
    for (int i=1; i<y.length-1; i++)  
      y[i] = gparent.kde[0].correctY(y[i]+ny);

    //println("" + nx + " --- " + ny + "   " + y.length);
  }


  /// Sampling extrapolate a new point from two point if thoose are too far from each other.
  /// But also can suppress a point if he is to close to another.
  /// Sampling is done only one vertices from the same trajectory.
  void resampleControlPoints()
  {
    if (firstResample) {
      firstResample=false; //in the first call Control Points are good
      initControlPoints();

      if (Parameters.nBins > 1 && Parameters.separationFactor > 1e-8)                
        translateControlPointsByNormal();

      return;
    }

    final int initialCapacity= x.length*2;
    int n=1;
    
    float[] tmpX = new float[initialCapacity];
    float[] tmpY = new float[initialCapacity];

    tmpX[0] = x[0]; 
    tmpY[0] = y[0];

    boolean hasChanged = false;

    for (int j = 0; j < x.length - 1; j++)
    {
      final float curX = x[j];
      final float curY = y[j];
      final float nextX = x[j+1];
      final float nextY = y[j+1];
      final float ddist = dist(curX, curY, nextX, nextY);
      //final float dist = abs(curX-nextX) + abs(curY-nextY);

      if (ddist > Parameters.splitDist) {
        tmpX[n]   = (curX+nextX)/2.0f;
        tmpY[n++] = (curY+nextY)/2.0f;
        hasChanged = true;
      }

      if (!(ddist < Parameters.removeDist) || j == x.length - 2) {
        tmpX[n]   = nextX;
        tmpY[n++] = nextY;
      }
      
    }

    if (hasChanged) 
    { //actually, if there are only removals there is no update...      
      if (x.length != n) {
        bins = new int[n-1]; //even if not needed it is set to zero      
        x = new float[n];
        y = new float[n];
      }
      //x = Arrays.copyOf(tmpX, n);
      //y = Arrays.copyOf(tmpY, n);
      // } else {
      arrayCopy( tmpX, 0, x, 0, n );
      arrayCopy( tmpY, 0, y, 0, n );
      //}
    }
  }


  void applyGradientToControlPoints() { 
    final float halfBins = Parameters.nBins * 0.5f;    

    final int idxMargin = 1; 
    for (int i=idxMargin; i<x.length-idxMargin; i++) { 
      final KdeCpu map = gparent.kde[bins[min(i, bins.length-1)]];
      final int curX = (int) (x[i]); //gradient coulb be improvede by using the real coords
      final int curY = (int) (y[i]); 

      float[] g = map.getNormalizedGradientAtPoint(curX, curY);

      float h = gparent.h; //step magnitude
      float stepY = g[0] * h;
      float stepX = g[1] * h;
      boolean update = true;

      if (!Parameters.blindStep) {
        if (map.getDensityAtPoint(curX, curY) <= -MAX_FLOAT)
          Log.log("Alert! Point @ " + curX + " , " + curY);
        
        final float originalDensity = map.getDensityAtPoint(curX, curY);
        while ( map.getDensityAtPoint ( (int)(x[i]+stepX), (int)(y[i]+stepY)) < originalDensity) {         
          stepY /= 2.0f;
          stepX /= 2.0f;  
          
          if (stepX < 1f && stepY <1f) {  
            update = false;
            break;
          }
        }
      }

      if (update) {  
        y[i] += stepY;
        x[i] += stepX;        
      }

      //println("" + curX + " , " + curY + "  -->  "  + x[i] + " , " + y[i]);
    }
  }

  void calculateBins() {     
    if (bins.length != x.length-1)
      bins = new int[x.length-1];

    if (Parameters.binType == Parameters.BIN_ORIENTATION) {
      for (int j = 0; j < x.length - 1; j++) {          
        final float angle = atan2(y[j+1]-y[j], x[j+1]-x[j]);          
        bins[j] = round((angle + PI)*Parameters.nBins/TWO_PI) 
          % Parameters.nBins;
      }
      return;
    }

    if (firstBining) {   
      firstBining = false;  
      switch (Parameters.binType) {          
      case Parameters.BIN_SOURCE_MATRIX:
      case Parameters.BIN_TARGET_MATRIX:      
        if (Parameters.nBins==1)
          this.clusterIdx = 0;
        else {
          Node refNode = (Parameters.binType==Parameters.BIN_TARGET_MATRIX ? target : source); 
          int hw = round(sqrt(Parameters.nBins));
          int tmpx = floor((refNode.x*1.0f/Parameters.mapWidth)*hw);
          if (tmpx >= hw) tmpx = hw-1;
          int tmpy = floor((refNode.y*1.0f/Parameters.mapHeight)*hw);
          if (tmpy >= hw) tmpy = hw-1;
          int bin = tmpx + tmpy * hw;         
          this.clusterIdx = bin;
        }
        break;
      }
    }

    //Arrays.fill(bins, this.clusterIdx);
    for (int i=0; i<bins.length; i++)
      bins[i]= this.clusterIdx;
  }


  void smoothControlPoints(int n) {
    final float a = Parameters.smoothingAlpha / 2.0f;
    for (int j=0; j<n; j++) {
      float prevx = x[0], prevy = y[0];     
      for (int i=1; i<x.length-1; i++) {
        final float incX = a*(prevx-x[i]+x[i+1]-x[i]);
        prevx = x[i];
        x[i] += incX;

        final float incY = a*(prevy-y[i]+y[i+1]-y[i]);                
        prevy = y[i];
        y[i] += incY;
      }
    }
  }


  float[] getNodePosition(Node node) {    
    if (node == source) {
      final float[] pos = {
        y[0], x[0]
      }; 
      return pos;
    }
    if (node == target) {
      final float[] pos = {
        y[y.length-1], x[x.length-1]
      }; 
      return pos;
    }
    return null;
  }

  void updateNodePosition(Node node) {
    if (node == source) {
      x[0] = node.x;
      y[0] = node.y;
    } else if (node == target) {
      x[x.length-1] = node.x;
      y[y.length-1] = node.y;
    }
  }
}

