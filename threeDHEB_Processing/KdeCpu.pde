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
  KdeCpu: CPU-based implementation of density maps.
*/ 

public class KdeCpu {
  float [][]map;


  KdeCpu() {
    clearMap();
  }


  public float[][] getCopyOfMap() {
    float [][] mapCopy = new float[map.length][];
    for (int i = 0; i < map.length; i++)
      mapCopy[i] = map[i].clone();
    return mapCopy;
  }

  public void clearMap() {
    map =  new float[Parameters.mapHeight][Parameters.mapWidth];
  }

  //Bresenham's algorithm
  public void drawLine(int x0, int y0, int x1, int y1, float w) {
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);

    int sx = x0 < x1 ? 1 : -1; 
    int sy = y0 < y1 ? 1 : -1; 

    int err = dx-dy;    
    int currentX = x0;
    int currentY = y0;


    while (true) {      
      if ((currentX == x1 && currentY == y1) || currentX<0 || currentY<0 || currentX>=Parameters.mapWidth || currentY >=Parameters.mapHeight) {
        return;
      }

      map[currentY][currentX] += w;


      int e2 = err<<1;
      if (e2 > -dy) {
        err -= dy;
        currentX += sx;
      }

      if (e2 < dx) {
        err +=  dx;
        currentY += sy;
      }
    }
  }

  void applyFilter(ImageFilter aFilter) {
    aFilter.processImage(map);
  }


  void normalizeMapByMax() {
    divide(getMax());
  }

  float getMax() {
    float maxVal = max(map[0]);   
    for (int j=1; j<Parameters.mapHeight; j++) {    
      float tmpMax = max(map[j]);      
      if (tmpMax > maxVal) maxVal = tmpMax;
    }
    return maxVal;
  }

  void divide(float val) {
    for (int j=0; j<Parameters.mapHeight; j++)
      for (int i=0; i<Parameters.mapWidth; i++)
        map[j][i] /= val;
  }


  void normalizeMapByMaxMin() {
    float maxVal = max(map[0]);
    float minVal = min(0, min(map[0]));   
    for (int j=1; j<Parameters.mapHeight; j++) {    
      final float tmpMax = max(map[j]);      
      if (tmpMax > maxVal) maxVal = tmpMax;
      final float tmpMin = min(map[j]);      
      if (tmpMin < minVal) minVal = tmpMin;
    }

    // println("@ Max val: " + maxVal);

    for (int j=0; j<Parameters.mapHeight; j++)
      for (int i=0; i<Parameters.mapWidth; i++)
        map[j][i] = (map[j][i] - minVal) / (maxVal - minVal);
  }

  
  float correctX(float x) {
    return min(max(x,2),Parameters.mapWidth-3);     
  }    
  
  float correctY(float y) {
    return min(max(y,2),Parameters.mapHeight-3);     
  }      
  
  boolean outOfMap(int x, int y) {
    return (x<2||y<2||x>=Parameters.mapWidth-2||y>=Parameters.mapHeight-2);
  }    

  float getDensityAtPoint(int x, int y) {
    if (outOfMap(x,y))
      return -MAX_FLOAT;

    return map[y][x];
  }


  float getLocalDensity(int x, int y) {
    return map[y][x];
  }
  
  float[]  getGradientAtPoint(int x, int y) {
    try {  
      final float[] g = {
        map[y+1][x]-map[y-1][x], 
        map[y][x+1]-map[y][x-1]
      };    
      return g;
    }
    catch (Exception e) {    
      e.printStackTrace();
      Log.fatalError("EXCEPTION: " + x + " ; " + y);
      return null;
    }
  }


  float[] getNormalizedGradientAtPoint(int x, int y) {  
    float[] g = getGradientAtPoint( x, y);  
    final float gmag = mag(g[0], g[1]);
    //println("G: " + g[0] + "  " +  g[1] +  "  "  + gmag);
    if (gmag > 0) {
      g[0] /= gmag;
      g[1] /= gmag;
    }
    return g;
  }


  PImage makeImage(float hue, float sat) {
    int k = 0;
    PImage img = createImage(Parameters.mapWidth, Parameters.mapHeight, RGB);
    img.loadPixels();    
    for (int j = 0; j < Parameters.mapHeight; j++) {
      for (int i = 0; i < Parameters.mapWidth; i++) {      
        final float c = map[j][i];
        img.pixels[k++] = color(hue, sat, c, 0.15);
      }
    }
    img.updatePixels();       

    return img;
  }

  PImage makeImage(float hue) {
    return makeImage(hue, 1.0f);
  }

}

