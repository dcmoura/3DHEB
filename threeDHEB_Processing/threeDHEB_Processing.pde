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
  Main file
*/



Graph graph=null;             //the graph data structure and bundling functions         

//auxiliar variables
int curIt = 0;                //current iteration
boolean needsReset = false;   //used by the JavaScript app



void setup() {  
  Parameters.validate();

  if (Parameters.isJavaScript)
    size(800, 800, P2D);
    
  frameRate(Parameters.vizFrameRate);
  colorMode(HSB, 1.0, 1.0, 1.0, 1.0);

  loadGraph(Parameters.dataFilename);

  if (Parameters.isJavaScript)
    noLoop();
}


void draw() { 
  if (!Parameters.isJavaScript && curIt<10) 
    iteration(1);
    
  background(0);
  
  PGraphics pg = draw2pg(Parameters.screenFactor); 
  image(pg, 0,0);    
}



public void loadGraph(String graphLocation) {   
  if (Parameters.isJavaScript || graph == null) {  
    graph = new Graph();                //creates a new Graph object 
    graph.loadGraph(graphLocation);     //loads graph data from a file
    graph.makeDensityMap();             //creates the initial density map (required for rendering)
  }

  //size of the sketch maches the bounding box of the graph
  //screenFactor parameter allows to render higher resolution images
  size(round(Parameters.mapWidth*Parameters.screenFactor), 
       round(Parameters.mapHeight*Parameters.screenFactor), 
       (Parameters.isJavaScript ? P2D : P3D)); //I found problems with P3D in javascript...

  Log.log("Res: " + width + " x " + height);
  Log.log("Hist: " + Parameters.mapWidth + " x " + Parameters.mapHeight);
  
  Parameters.dataFilename = graphLocation;
  curIt = 0;
  
  redraw();  //draws the graph
}


//draws to a PGraphics object for later rendering onScreen or saving to image file
//receives a screen factor as parameter for controling the final resolution
PGraphics draw2pg(int sf) { 
  PGraphics pg = createGraphics(width*sf, height*sf, (Parameters.isJavaScript ? P2D : P3D));
  pg.colorMode(HSB, 1.0, 1.0, 1.0, 1.0);
  pg.beginDraw();
  pg.colorMode(HSB, 1.0, 1.0, 1.0, 1.0);
  pg.background(0);
  if (curIt==0)
    graph.drawGraph(pg, 0.25, sf, false, true); //draw raw (original) graph
  else 
    graph.drawGraph(pg, Parameters.edgeAlpha, sf, Parameters.useStrokeWeight, !true); //draw current bundling

  if (Parameters.invertColors) 
    pg.filter(INVERT);  
  pg.endDraw();
  
  return pg;
}


public void iteration(int n) {
  if (needsReset && curIt > 0) { 
    reset();
    needsReset = false;
  }

  for (int i=0; i<n; i++) {
    curIt++;
    graph.makeDensityMap(); //calculates density maps
    graph.updateEdges();    //update edges trajectories based on the density maps
  }

  Log.log("H: " + graph.h + "  Hrel: " + graph.hrel  + "  Sigma: " + graph.sigma);

  if (Parameters.isJavaScript)
    redraw();
}


// Functions for interacting with the online app

public void reset() {
  needsReset=false;
  loadGraph(Parameters.dataFilename);
}


public void setCriteria(int val) {
  if (val==0) {    //Default criteria: edges are bundled based on proximity 
    Parameters.nBins = 1;
    Parameters.histBlendMode = Parameters.NO_BLEND;
    Parameters.colorMode = Parameters.CM_HUE;
    Parameters.binType = Parameters.BIN_TARGET_MATRIX;
  } else if (val==1) {     //Flow separation based on edge orientation with a modal blending function
    Parameters.nBins = 8;
    Parameters.histBlendMode = Parameters.BLEND_MODAL;
    Parameters.colorMode = Parameters.CM_DIR;
    Parameters.binType = Parameters.BIN_ORIENTATION;
  } else if (val==2) {   
    Parameters.nBins = 4;  //Flow separation based on edge orientation with a evasive blending function
    Parameters.histBlendMode = Parameters.BLEND_EVASIVE;
    Parameters.colorMode = Parameters.CM_DIR;
    Parameters.binType = Parameters.BIN_ORIENTATION;
  } else
    return;

  Parameters.colorSaturation = (Parameters.nBins == 1 ?  0.0f : 0.60f);

  reset();
  println("Criteria: " + val);
}

public void setSigma(float val) {
  Parameters.kernelSigma = val * Parameters.mapSize;  
  needsReset=true;
}

public void setHStep(float val) {
  Parameters.hStep = val;
  needsReset=true;
}

public void setEdgeSmoothing(float val) {
  Parameters.smoothingAlpha = val;
  needsReset=true;
}

public void setSeparationFactor(float val) {
  Parameters.separationFactor = val * Parameters.mapSize;
  reset();
}

public void setRepellence(float val) {
  Parameters.blendPenalty = val;
  needsReset=true;
}


public void setInvertColors(int val) {  
  Parameters.invertColors = (val != 0);
  redraw();
}

public void setUseWeight(int val) {
  Parameters.useStrokeWeight = (val != 0);
  redraw();
}

public void setEdgeAlpha(float val) {
  Parameters.edgeAlpha = val;
  redraw();
}

public void saveImage() {
  save("lr.png");
}

public void saveImageHR() {
  PGraphics pg = draw2pg(Parameters.screenFactor*4); 
  pg.save("hr.png");
}




