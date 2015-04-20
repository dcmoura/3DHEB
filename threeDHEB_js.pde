// Author: Daniel C. Moura (daniel.c.moura@gmail.com)
// April 2015
// version 0.1


Graph graph=null;
int curIt = 0;
boolean needsReset = false;
boolean isJavaScript = true; //choose false if you are running a local Java app


public void loadGraph(String graphLocation) {   
  if (isJavaScript || graph == null) {  
    graph = new Graph();  
    graph.loadGraph(graphLocation);
  }

  size(round(Parameters.mapWidth*Parameters.screenFactor), 
  round(Parameters.mapHeight*Parameters.screenFactor), 
  (isJavaScript ? P2D : P3D));

  Log.log("Res: " + width + " x " + height);
  Log.log("Hist: " + Parameters.mapWidth + " x " + Parameters.mapHeight);

  graph.prepareBundlind(); 

  graph.makeDensityMap(); //remove this line?

  Parameters.dataFilename = graphLocation;
  curIt = 0;

  redraw();
}


public void reset() {
  needsReset=false;
  loadGraph(Parameters.dataFilename);
}

public void iteration(int n) {
  if (needsReset && curIt > 0) { 
    reset();
    needsReset = false;
  }

  for (int i=0; i<n; i++) {
    curIt++;
    graph.makeDensityMap(); //remove this line
    graph.updateEdges();
  }

  Log.log("H: " + graph.h + "  Hrel: " + graph.hrel  + "  Sigma: " + graph.sigma);

  if (isJavaScript)
    redraw();
}


public void setCriteria(int val) {
  if (val==0) {    
    Parameters.nBins = 1;
    Parameters.histBlendMode = Parameters.NO_BLEND;
    Parameters.colorMode = Parameters.CM_HUE;
    Parameters.binType = Parameters.BIN_TARGET_MATRIX;
  } else if (val==1) {    
    Parameters.nBins = 8;
    Parameters.histBlendMode = Parameters.BLEND_MODAL;
    Parameters.colorMode = Parameters.CM_DIR;
    Parameters.binType = Parameters.BIN_ORIENTATION;
  } else if (val==2) {   
    Parameters.nBins = 4;
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
  Parameters.hMax = Parameters.sigma2h * Parameters.kernelSigma;
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




void setup() {  
  Parameters.validate();

  if (isJavaScript)
    size(800, 800, P2D);
  frameRate(Parameters.vizFrameRate);
  colorMode(HSB, 1.0, 1.0, 1.0, 1.0);

  loadGraph(Parameters.dataFilename);

  if (isJavaScript)
    noLoop();
}



PGraphics draw2pg(int sf) { 
  PGraphics pg = createGraphics(width*sf, height*sf, JAVA2D);
  pg.colorMode(HSB, 1.0, 1.0, 1.0, 1.0);
  pg.beginDraw();
  pg.colorMode(HSB, 1.0, 1.0, 1.0, 1.0);
  pg.background(0);
  if (curIt==0)
    graph.drawGraph(pg, 0.25, sf, false, true);
  else 
    graph.drawGraph(pg, Parameters.edgeAlpha, sf, Parameters.useStrokeWeight, !true);

  if (Parameters.invertColors) 
    pg.filter(INVERT);  
  pg.endDraw();
  
  return pg;
  
}

void draw() { 
  if (!isJavaScript) 
    iteration(1);
    
  background(0);
  
  PGraphics pg = draw2pg(Parameters.screenFactor); 
  image(pg, 0,0);    
}

// Author: Daniel C. Moura (daniel.c.moura@gmail.com)
// April 2015
// version 0.1



class Edge {
  int sourceIdx, targetIdx;
  Node source, target;
  float weightOrig, weight;
  String id;
  Graph gparent;

  float[] x;
  float[] y;
  float distanceOrig, distanceNorm;

  int[] bins;
  boolean firstResample=true;
  boolean firstBining = true;

  int clusterIdx=0;

  Edge(Node source, Node target, float weightOrig, String id, Graph gparent) {
    this.source = source;
    this.target = target;
    sourceIdx = source.idx;
    targetIdx = target.idx;
    this.weightOrig = weightOrig;    
    this.id = id;    
    this.gparent = gparent;   
    this.distanceOrig = dist(source.xOrig, source.yOrig, target.xOrig, target.yOrig); //!!    
    if (!Parameters.pinNodes) {
      source.edges.add(this);
      target.edges.add(this);
    }
  }



  void normalize(float minW, float maxW, float minD, float maxD) {    
    //weight = (weightOrig-minW) / (maxW-minW);
    weight = (weightOrig) / (maxW);
    distanceNorm = (distanceOrig-minD) / (maxD-minD);
  }


  void initControlPoints() {
    //final int minSpace = (removeDist*9 + splitDist*1)/10; //tries to delay splitting
    //final float minSpace = (Parameters.removeDist*9 + Parameters.splitDist*1)/10; //tries to delay splitting
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

  // float version
  void applyGradientToControlPoints() {
    //Object o=null; o.getClass(); 

    final float halfBins = Parameters.nBins * 0.5f;

    //stepMag = new float[x.length]; //new!!

    final int idxMargin = (Parameters.pinNodes ? 1 : 0); 
    for (int i=idxMargin; i<x.length-idxMargin; i++) { 
      final KdeCpu map = gparent.kde[bins[min(i, bins.length-1)]];
      final int curX = (int) (x[i]); //gradient coulb be improvede by using the real coords
      final int curY = (int) (y[i]); 



      float[] g = map.getNormalizedGradientAtPoint(curX, curY);

      //if (i==1)
      //  println(" " + g[0] + " " + g[1]);


      float h = gparent.h; //step magnitude
      float stepY = g[0] * h;
      float stepX = g[1] * h;
      boolean update = true;

      if (!Parameters.blindStep) {
        if (map.getDensityAtPoint(curX, curY) <= -MAX_FLOAT)
          Log.log("Alert! Point @ " + curX + " , " + curY);
        
        final float originalDensity = map.getDensityAtPoint(curX, curY);
        //int newX = (int)(x[i]+stepX);
        //int newY = (int)(y[i]+stepY);
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

  void calculateOrientationBins() {     
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

/**
 * Copyright (c) 2011, The University of Southampton and the individual contributors.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *   *   Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 *
 *   *  Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 *   *  Neither the name of the University of Southampton nor the names of its
 *   contributors may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * Fast approximate Gaussian smoothing using repeated fast box filtering.
 * 
 * @author Jonathon Hare (jsh2@ecs.soton.ac.uk)
 * 
 
 @Reference(
 type = ReferenceType.Inproceedings,
 author = { "Kovesi, P." },
 title = "Fast Almost-Gaussian Filtering",
 year = "2010",
 booktitle = "Digital Image Computing: Techniques and Applications (DICTA), 2010 International Conference on",
 pages = { "121", "125" },
 month = "Dec",
 customData = {
 "keywords", "Gaussian processes;approximation theory;band-pass filters;image processing;Gaussian bandpass filters;fast almost-Gaussian filtering;image averaging;integral images;log-Gabor filters;separable moving average filters;summed area tables;symmetric transfer function;Approximation methods;Bandwidth;Computer vision;Frequency domain analysis;Laplace equations;Pixel;Transfer functions;Difference of Gaussian filtering;Gaussian smoothing",
 "doi", "10.1109/DICTA.2010.30"
 })
 */


//Code modified by Daniel C. Moura    
//daniel.c.moura@gmail.com


interface ImageFilter {
  public void processImage(final float[][] image);
  public float getSigma();
}  


public class FastGaussianConvolve implements ImageFilter {
  private  int n;
  private  int m;
  private AverageBoxFilter wlBox;
  private AverageBoxFilter wuBox;
  private float sigma;

  /**
   * Construct an {@link FFastGaussianConvolve} to approximate blurring with a
   * Gaussian of standard deviation sigma.
   * 
   * @param sigma
   *            Standard deviation of the approximated Gaussian
   * @param n
   *            Number of filtering iterations to perform (usually between 3
   *            and 6)
   */
  public FastGaussianConvolve(float sigma, int n) {
    this.sigma = sigma;
    //println("\n\nNEW SIGMA: " + sigma + "\n\n");
    if (sigma < 1.8) {
      // std.devs of less than 1.8 are not well approximated.
      n=0;
      m=0;
      Log.fatalError("Low sigma");
    } else {
      final float ss = sigma * sigma;
      final double wIdeal = sqrt((12.0 * ss / n) + 1.0);
      final int wl = (((int) wIdeal) % 2 == 0) ? (int) wIdeal - 1 : (int) wIdeal;
      final int wu = wl + 2;

      this.n = n;
      this.m = round((12 * ss - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4));

      this.wlBox = new AverageBoxFilter(wl);
      this.wuBox = new AverageBoxFilter(wu);
    }
  }


  public void processImage(final float[][] image) {
    for (int i = 0; i < m; i++)
      wlBox.processImage(image);
    for (int i = 0; i < n - m; i++)
      wuBox.processImage(image);
  }

  public float getSigma() {
    return sigma;
  }
}






/**
 * A rectangular averaging convolution operator (often known as a Box filter).
 * For efficiency, this is implemented using a {@link SummedAreaTable} rather
 * than through an actual convolution.
 * 
 * @author Jonathon Hare (jsh2@ecs.soton.ac.uk)
 */
public class AverageBoxFilter {
  private int bwidth;
  private int bheight;

  /**
   * Construct the averaging operator with a kernel of the given dimensions.
   * 
   * @param width
   *            width of the kernel
   * @param height
   *            height of the kernel
   */
  public AverageBoxFilter(int width, int height) {
    this.bwidth = width;
    this.bheight = height;
  }

  /**
   * Construct the averaging operator with a square kernel of the given
   * dimension.
   * 
   * @param dim
   *            The width and height of the box
   */
  public AverageBoxFilter(int dim) {
    this(dim, dim);
  }


  public void processImage(float[][] image) {
    // shortcut trivial case
    //if (this.bheight == 1 && this.bwidth == 1)
    //  return;

    //final SummedAreaTable sat = new SummedAreaTable(image, hw, hh);
    //sat.analyseImage(image);

    final int hw = int(bwidth / 2);
    final int hh = int(bheight / 2);    

    final SummedAreaTable sat = new SummedAreaTable(image);

    for (int y = 0; y < Parameters.mapHeight; y++) {
      final int sy = Math.max(0, y - hh);
      final int ey = Math.min(Parameters.mapHeight, y + hh + 1);
      final int dy = ey-sy;

      for (int x = 0; x < Parameters.mapWidth; x++) {
        final int sx = Math.max(0, x - hw);
        final int ex = Math.min(Parameters.mapWidth, x + hw + 1);
        final int area = (ex - sx) * dy;        
        image[y][x] = sat.calculateArea(sx, sy, ex, ey) / area;
        ;
      }
    }
  }
}





/**
 * Implementation of an Integral Image or Summed Area Table.
 * <p>
 * See http://en.wikipedia.org/wiki/Summed_area_table and
 * http://research.microsoft
 * .com/en-us/um/people/viola/Pubs/Detect/violaJones_IJCV.pdf
 * <p>
 * Basically, this provides an efficient way to find the sum of all pixels in a
 * rectangular area of an image.
 * 
 * @author Jonathon Hare (jsh2@ecs.soton.ac.uk)
 */
public class SummedAreaTable {
  /**
   * The SAT iData
   */
  public float[][] iData;


  /**
   * Construct an empty SAT
   */
  public SummedAreaTable() {
    iData = new float[Parameters.mapHeight+1][Parameters.mapWidth + 1];
  }

  /**
   * Construct a SAT from the provided image
   * 
   * @param image
   *            the image
   */
  public SummedAreaTable(float[][] image) {
    computeTable(image);
  }

  protected void computeTable(float[][] image) {
    iData = new float[Parameters.mapHeight+1][Parameters.mapWidth + 1];

    for (int y = 0; y < Parameters.mapHeight; y++) {
      for (int x = 0; x < Parameters.mapWidth; x++) {
        iData[y + 1][x + 1] = image[y][x] +
          iData[y + 1][x] +
          iData[y][x + 1] -
          iData[y][x];
      }
    }
  }

  /**
   * Calculate the sum of pixels in the image used for constructing this SAT
   * within the rectangle defined by (x1,y1) [top-left coordinate] and (x2,y2)
   * [bottom- right coordinate]
   * 
   * @param x1
   *            x1
   * @param y1
   *            y1
   * @param x2
   *            x2
   * @param y2
   *            y2
   * @return sum of pixels in given rectangle
   */
  public float calculateArea(final int xx1, final int yy1, final int xx2, final int yy2) {

    final float iA = iData[yy1][xx1];
    final float iB = iData[yy1][xx2];
    final float iC = iData[yy2][xx2];
    final float iD = iData[yy2][xx1];

    return iA + iC - iB - iD;
  }

  /**
   * Calculate the sum of pixels in the image used for constructing this SAT
   * within the given rectangle
   * 
   * @param r
   *            rectangle
   * @return sum of pixels in given rectangle
   
   public float calculateArea(Rectangle r) {
   return calculateArea(Math.round(r.x), Math.round(r.y), Math.round(r.x + r.width), Math.round(r.y + r.height));
   }
   */

  /*
   * (non-Javadoc)
   * 
   * @see
   * org.openimaj.image.analyser.ImageAnalyser#analyseImage(org.openimaj.image
   * .Image)
   */

  public void analyseImage(float[][] image) {
    computeTable(image);
  }
}

// Author: Daniel C. Moura (daniel.c.moura@gmail.com)
// April 2015
// version 0.1

import java.util.*;

class Graph {
  protected Node[] nodes;
  protected Edge[] edges;
  protected color[] colors; //color of each cluster
  protected int[] clusterOrder;
  final int nThreads = 1;
  protected float[][][] tmpVert;  
  protected ImageFilter ffilter = null;

  protected float h, hrel, sigma;

  KdeCpu[] kde;
  float[][] blendingWeights;

  float maxX, maxY, minX, minY, maxW, minW, minD, maxD, ar;
  

  Profiler[] prof = { 
    new Profiler("TOTAL B."), 
    new Profiler("RESAMPLE"), 
    new Profiler("DENSITTY"), 
    new Profiler("E UPDATE"), 
    new Profiler("N UPDATE"), 
    new Profiler("DRAWING ")
    };




    Graph() {    
      kde = null;
      h = Parameters.hMax;
      hrel = 1.0f;
      sigma = Parameters.kernelSigma;
    }


  void prepareBundlind() {
    normalizeVals();
        
    kde = new KdeCpu[Parameters.nBins];
    for (int i=0; i< kde.length; i++)
      kde[i] = new KdeCpu();


    setBlendMode(Parameters.histBlendMode);

    Log.log("Init: end");
  }




  void makeDensityMap() { 
    //upgrade to executorService?

    prof[0].start();//total

    final int edgesPerThread = ceil(edges.length*1.0f/nThreads);
    Thread[] t = new Thread[nThreads];

    prof[1].start();

    for (int i=0; i<edges.length; i++) {
      final Edge edge = edges[i];            
      edge.resampleControlPoints();
      if (Parameters.nBins > 1)                          
        edge.calculateOrientationBins();
    }

    prof[1].end();    



    prof[2].start();    
    setKernelSigma(sigma);

    for (int j=0; j< (kde.length); j++) {
      final int jj = j;

      kde[jj].clearMap();      
      for (int i=0; i<edges.length; i++)      
        drawEdgeOnKDEmap(edges[i], jj);
      kde[jj].applyFilter(ffilter);
    }

    prof[2].end();
  }




  protected void drawEdgeOnKDEmap(final Edge edge, int mapIdx) {
    if (Parameters.nBins == 1) {
      final KdeCpu map = kde[0];
      for (int j=mapIdx; j<edge.x.length-1; j+=nThreads)
        map.drawLine((int)(edge.x[j]), (int)(edge.y[j]), (int)(edge.x[j+1]), (int)(edge.y[j+1]), edge.weight);
    } else {      
      final KdeCpu map = kde[mapIdx];
      final float[] w = blendingWeights[mapIdx];

      if (Parameters.binType == Parameters.BIN_ORIENTATION) {
        for (int j=0; j<edge.x.length-1; j++) {
          final int edgeBin = edge.bins[j];                
          if (w[edgeBin] != 0)
            map.drawLine(round(edge.x[j]), round(edge.y[j]), round(edge.x[j+1]), round(edge.y[j+1]), edge.weight * w[edgeBin]);
        }
      } else {
        float linkw = w[edge.clusterIdx]; 
        if (linkw != 0.0f) {
          linkw *= edge.weight;
          for (int j=0; j<edge.x.length-1; j++)
            map.drawLine(round(edge.x[j]), round(edge.y[j]), round(edge.x[j+1]), round(edge.y[j+1]), linkw);
        }
      }
    }
  }


  //----

  void updateEdges() {            

    prof[3].start();

    for (int i=0; i<edges.length; i++) {
      final Edge edge = edges[i];
      edge.applyGradientToControlPoints();
      edge.smoothControlPoints(Parameters.smoothingIterations);
    }

    prof[3].end();


    h *= Parameters.hStep;
    hrel *= Parameters.hStep;
    sigma *= Parameters.hStep;

    prof[0].end();//total

    if (Profiler.logOn) {
      println("-----------*-----------");
      for (int i=0; i<prof.length; i++) {
        prof[i].inc();
        prof[i].print();
      }
      println("-----------=-----------");
    }
  }


  void drawDensityMap() {    
    //blendMode(WAS_BLEND);
    globalyNormalizeMaps();
    for (int i=0; i< kde.length; i++) {        
      drawDensityMap(i);
    }
  }

  void drawDensityMap(int i) {    
    //kde[i].normalizeMapByMax(); 
    image(kde[i].makeImage(i*1.0f/kde.length, Parameters.colorSaturation ), 0, 0, width, height);
  }


  void globalyNormalizeMaps() {
    float maxVal = -MAX_FLOAT;
    for (int i=0; i< kde.length; i++) {
      final float tmpMax = kde[i].getMax(); 
      if (tmpMax > maxVal) maxVal = tmpMax;
    }
    for (int i=0; i< kde.length; i++)
      kde[i].divide(maxVal);
  }

  void setKernelSigma(float newSigma) {    
    newSigma = max(1.8, newSigma);  
    if (ffilter == null || newSigma != ffilter.getSigma()) { 
      ffilter = new FastGaussianConvolve(newSigma, Parameters.convolutionSteps);
    }
  }

  void setBlendMode(int histBlendMode) {
    blendingWeights = (new MapBlender()).createBlender(histBlendMode).makeBlendingMatrix();
  }

  void drawGraph(PGraphics pg, float edgeAlpha, float s, boolean includeWeights, boolean drawNodes) {   
    prof[5].start(); 
    float hueVal;

    globalyNormalizeMaps();

    

    if (drawNodes) {
      pg.noStroke();    
      //fill(0xFF00FFFF, 0.5);
      pg.fill(0.5, 1-Parameters.colorSaturation, 1, 1);
      //blendMode(WAS_BLEND);
      for (int i=0; i<nodes.length; i++) {       
        Node node = nodes[i];
        pg.ellipse(node.x*s, node.y*s, 3*s, 3*s);
     
      }
    }


    pg.noFill();    
    pg.strokeWeight(1);
    pg.strokeCap(SQUARE);
    pg.smooth();

    for (int i=0; i<edges.length; i++) {      
      Edge edge = edges[i];

      //beginShape();        
      for (int j=0; j<edge.x.length; j++) {
        final int bin = (j==0?0:j-1);

        float dens=1.0;
        if (includeWeights) {
          dens = kde[edge.bins[bin]].getLocalDensity(round(edge.x[j]), round(edge.y[j]));
          if (dens <0) dens = 0;
          
          pg.strokeWeight((Parameters.strokeMaxWeight-Parameters.strokeMinWeight) * dens + Parameters.strokeMinWeight);
        }


        switch (Parameters.colorMode) {
        
        case Parameters.CM_HUE: //hue based  */

          hueVal = edge.bins[bin]*1.0f/Parameters.nBins;                              
          pg.stroke(hueVal, Parameters.colorSaturation, 1.0, edgeAlpha);

          break;


        case Parameters.CM_DIR: //direction

          float p = float(j)/(edge.x.length-1f);
          
          color r = min(255, round(max(0, p*2f) * 255f)); 
          color b = min(255, round(max(0, 2f*(1f-p)) * 255f));
          color g = int(min(r,b)*0.75);

          int c2 = 0xFF000000 | r << 16 | g << 8 |  b;

          pg.stroke(c2, edgeAlpha);
          break;

        default:
          pg.stroke(color(1, edgeAlpha));
        }



        //vertex(edge.x[j]*s, edge.y[j]*s);
        if (j>0)
          pg.line(edge.x[j-1]*s, edge.y[j-1]*s, edge.x[j]*s, edge.y[j]*s);
      }        

      //endShape();
    }
    prof[5].end();
  }




  float getAspectRatio() {
    return ar;
  }

  float calcAspectRatio() {
    ar = (maxY-minY)/(maxX-minX);

    if (ar <= 1.0f) 
    {
      Parameters.mapWidth = Parameters.mapSize;
      Parameters.mapHeight = round(Parameters.mapWidth * graph.getAspectRatio());
    } else {
      Parameters.mapHeight = Parameters.mapSize;
      Parameters.mapWidth = round(Parameters.mapHeight / graph.getAspectRatio());
    }

    //println("AR: " + ar);
    return ar;
  }

  void loadGEXF(String filename) {
    XML xml = loadXML(filename);      
    XML xmlGraph = xml.getChild("graph");
    XML[] xmlNodes = xmlGraph.getChild("nodes").getChildren("node");
    XML[] xmlEdges = xmlGraph.getChild("edges").getChildren("edge");

    Log.log("Loading " + xmlNodes.length + " nodes...");
    nodes = new Node[xmlNodes.length];
    Hashtable<String, Integer> nodeId2IdxTable = new Hashtable<String, Integer>();
    maxX = -MAX_FLOAT;
    maxY = -MAX_FLOAT;
    minX = MAX_FLOAT;
    minY = MAX_FLOAT;          
    for (int i=0; i<xmlNodes.length; i++) {
      String id = xmlNodes[i].getString("id");
      float x = xmlNodes[i].getChild("viz:position").getFloat("x");
      float y = xmlNodes[i].getChild("viz:position").getFloat("y");

      //println(id + ": " + x + " ; " + y);

      nodes[i] = new Node(x, y, id, i, this);
      nodeId2IdxTable.put(id, i);  

      if (x>maxX) maxX=x;
      if (x<minX) minX=x;
      if (y>maxY) maxY=y;
      if (y<minY) minY=y;
    }
    calcAspectRatio();

/*
    println("min_lon = " + minX/1e0);
    println("max_lon = " + maxX/1e0);
    println("min_lat = " + minY/1e0);
    println("max_lat = " + maxY/1e0);
*/

    Log.log("Loading " + xmlEdges.length + " edges...");
    edges = new Edge[xmlEdges.length];  
    maxW = maxD = -MAX_FLOAT;
    minW = minD =  MAX_FLOAT;
    for (int i=0; i<xmlEdges.length; i++) {
      String id = xmlEdges[i].getString("id", "");
      String sourceId = xmlEdges[i].getString("source");
      String targetId = xmlEdges[i].getString("target");
      float weight = xmlEdges[i].getFloat("weight", 1.0f);

      try {
        edges[i] = new Edge(nodes[nodeId2IdxTable.get(sourceId)], 
        nodes[nodeId2IdxTable.get(targetId)], 
        weight, id, this);
      }
      catch (Exception e) {
        println("Error loading");
        edges[i] = new Edge(nodes[0], 
        nodes[0], 
        0, id, this);
      }
      //nodeId2IdxTable.put(id, i);  
      if (weight>maxW) maxW=weight;
      if (weight<minW) minW=weight;
      if (edges[i].distanceOrig>maxD) maxD=edges[i].distanceOrig;
      if (edges[i].distanceOrig<minD) minD=edges[i].distanceOrig;
    }


    Log.log("Graph loaded.");
  }

  void loadGraphML(String filename) {
    XML xml = loadXML(filename);      
    XML xmlGraph = xml.getChild("graph");
    XML[] xmlNodes = xmlGraph.getChildren("node");
    XML[] xmlEdges = xmlGraph.getChildren("edge");

    Log.log("Loading " + xmlNodes.length + " nodes...");
    nodes = new Node[xmlNodes.length];
    //Hashtable<String, Integer> nodeId2IdxTable = new Hashtable<String, Integer>();
    HashMap nodeId2IdxTable = new HashMap();
    maxX = -MAX_FLOAT;
    maxY = -MAX_FLOAT;
    minX = MAX_FLOAT;
    minY = MAX_FLOAT;          
    for (int i=0; i<xmlNodes.length; i++) {
      String id = xmlNodes[i].getString("id");      
      float x = -MAX_FLOAT;
      float y = -MAX_FLOAT;

      XML[] xmlNodeData = xmlNodes[i].getChildren("data");
      for (int j=0; j<xmlNodeData.length; j++) {
        String s = xmlNodeData[j].getString("key");
        if (s.equals("x"))
          x = float(xmlNodeData[j].getContent());
        else if (s.equals("y"))
          y = (Parameters.invertY?-1f:1f) * float(xmlNodeData[j].getContent());
      }

      nodes[i] = new Node(x, y, id, i, this);
      nodeId2IdxTable.put(id, i);  

      if (x>maxX) maxX=x;
      if (x<minX) minX=x;
      if (y>maxY) maxY=y;
      if (y<minY) minY=y;
    }
    calcAspectRatio();


    Log.log("Loading " + xmlEdges.length + " edges...");
    edges = new Edge[xmlEdges.length];  
    maxW = maxD = -MAX_FLOAT;
    minW = minD =  MAX_FLOAT;
    for (int i=0; i<xmlEdges.length; i++) {
      String id = xmlEdges[i].getString("id");
      String sourceId = xmlEdges[i].getString("source");
      String targetId = xmlEdges[i].getString("target");
      float weight = 1.0f;

      XML[] xmlEdgeData = xmlEdges[i].getChildren("data");
      for (int j=0; j<xmlEdgeData.length; j++) {
        String s = xmlEdgeData[j].getString("key");
        if (s.equals("value") || s.equals("weight"))
          weight = float(xmlEdgeData[j].getContent());
          //weight = xmlEdgeData[j].getFloatContent();
      }


      try {
        edges[i] = new Edge(
          nodes[(Integer) (nodeId2IdxTable.get(sourceId))], 
          nodes[(Integer) (nodeId2IdxTable.get(targetId))], 
        weight, id, this);
      }
      catch (Exception e) {
        println("Error loading");
        edges[i] = new Edge(nodes[0], 
        nodes[0], 
        0, id, this);
      }
      //nodeId2IdxTable.put(id, i);  
      if (weight>maxW) maxW=weight;
      if (weight<minW) minW=weight;
      if (edges[i].distanceOrig>maxD) maxD=edges[i].distanceOrig;
      if (edges[i].distanceOrig<minD) minD=edges[i].distanceOrig;
    }


    Log.log("Graph loaded.");
    Log.log("Filename: " + filename);
  }

  void loadTxt(String filename) {        
    Log.log("Loading txt graph...");
    ArrayList<Node> tmpNodes = new ArrayList<Node>(50000);
    LinkedList<Edge> tmpEdges = new LinkedList<Edge>();
    Hashtable<String, Integer> nodeId2IdxTable = new Hashtable<String, Integer>();

    maxX = -MAX_FLOAT;
    maxY = -MAX_FLOAT;
    minX = MAX_FLOAT;
    minY = MAX_FLOAT;

    maxW = maxD = -MAX_FLOAT;
    minW = minD =  MAX_FLOAT;

    int i=0;
    Scanner scanner;
    try {
      scanner =  new Scanner(new File(filename));//, ENCODING.name());
    }
    catch (Exception e) {      
      Log.fatalError("Error reading file.");
    }
  }

  void makeRandomGraph(int nNodes, int nEdges) {
    Log.log("Ganerating graph with " + nNodes + " nodes and " + nEdges + " edges...");
    randomSeed(1);
    nodes = new Node[nNodes];    
    maxX = -MAX_FLOAT;
    maxY = -MAX_FLOAT;
    minX = MAX_FLOAT;
    minY = MAX_FLOAT;

    for (int i=0; i<nNodes; i++) {
      final String id = "";       
      final float x = random(1.0f);
      final float y = random(0.5);

      nodes[i] = new Node(x, y, id, i, this);

      if (x>maxX) maxX=x;
      if (x<minX) minX=x;
      if (y>maxY) maxY=y;
      if (y<minY) minY=y;
    }
    calcAspectRatio();



    edges = new Edge[nEdges];  
    maxW = maxD = -MAX_FLOAT;
    minW = minD =  MAX_FLOAT;
    for (int i=0; i<nEdges; i++) {
      int k=0;
      final String id = "";
      final float weight = random(1);
      final Node source = nodes[int(random(nNodes-1))];
      Node target = nodes[int(random(nNodes-1))];
      //while ( source.dist2Node(target) > 0.25 && k++ < 3) %accepts nodes > 0.25 with prob 25%
      if (source.dist2Node(target) > 0.25)
        target = nodes[int(random(nNodes-1))];      

      edges[i] = new Edge(source, target, weight, id, this);

      if (weight>maxW) maxW=weight;
      if (weight<minW) minW=weight;
      if (edges[i].distanceOrig>maxD) maxD=edges[i].distanceOrig;
      if (edges[i].distanceOrig<minD) minD=edges[i].distanceOrig;
    }

    Log.log("Graph loaded.");
  }

  void makeRandomGraph(int nEdges) {
    makeRandomGraph( round(sqrt((nEdges/0.1))), nEdges);
  }

  void loadGraph(String filename) {    
    //println(filename);
    final String[] parts = filename.split("\\.");
    final String ext = parts[parts.length-1];
    if (ext.equals("graphml"))
      loadGraphML(filename);
    else if (ext.equals("gexf"))
      loadGEXF(filename);
    else if (ext.equals("txt"))
      loadTxt(filename);
    else if (ext.equals("random"))
      switch (parts.length) {
      case 2: 
        makeRandomGraph(int(parts[0])); 
        break;
      case 3: 
        makeRandomGraph(int(parts[0]), int(parts[1])); 
        break;
      } 
    else    
      loadGraphML(filename);
  }



  private void normalizeVals() {
    float xl = graph.maxX-graph.minX;
    float yl = graph.maxY-graph.minY;



    minX-=Parameters.margin*xl;
    maxX+=Parameters.margin*xl;
    minY-=Parameters.margin*yl;
    maxY+=Parameters.margin*yl;        

/*
    println("min_lon = " + minX/1e0);
    println("max_lon = " + maxX/1e0);
    println("min_lat = " + minY/1e0);
    println("max_lat = " + maxY/1e0);
*/

    calcAspectRatio();

    for (int i=0; i<nodes.length; i++)
      nodes[i].calcCoords(minX, maxX, minY, maxY);
    for (int i=0; i<edges.length; i++)
      edges[i].normalize(minW, maxW, minD, maxD);
  }


  int getNumberOfCPs() {
    int n = 0;
    for (int i=0; i<edges.length; i++)
      n += edges[i].x.length;
    return n;
  }
}
// Author: Daniel C. Moura (daniel.c.moura@gmail.com)
// April 2015
// version 0.1


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

  public void drawLine(int x0, int y0, int x1, int y1, float w) {
    int dx = abs(x1 - x0);
    int dy = abs(y1 - y0);

    int sx = x0 < x1 ? 1 : -1; 
    int sy = y0 < y1 ? 1 : -1; 

    int err = dx-dy;    
    int currentX = x0;
    int currentY = y0;


    while (true) {
      //map[currentY][currentX] += w;

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

  /*
  float getLocalDensity(int x, int y) {
   float maxVal = -MAX_FLOAT;
   
   if (x<2||y<2||x>=Parameters.mapWidth-2||y>=Parameters.mapHeight-2)    
   return maxVal;
   
   final int s=0; //TODO: test s=0;
   for (int j=y-s;j<=y+s;j++)
   for (int i=x-s;i<=x+s;i++)
   if (map[j][i] > maxVal) maxVal = map[j][i];
   
   return maxVal;      
   }
   */


  float getLocalDensity(int x, int y) {
    return map[y][x];
  }
  /*
  float[]  getGradientAtPoint(int x, int y) {
   //println("GCalc (" + x+ " , " + y + "): " + map[y][x] + " , " + map[y-1][x] + " , " + map[y+1][x] + " , " + map[y][x-1] + " , " + map[y][x+1]);
   //Object o = null; o.getClass();
   float[] g=new float[2];
   try {      
   //central difference approximation of the derivative
   g[0] = map[y+1][x]-map[y-1][x];
   g[1] = map[y][x+1]-map[y][x-1];
   return g;
   }
   catch (Exception e) {
   println("EXCEPTION: " + x + " ; " + y);
   e.printStackTrace();
   //      return g;
   } 
   }
   */

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


  void normalizeMapDebug() {
    HashSet<Float> vals = new HashSet<Float>();

    float maxVal = -MAX_FLOAT;    
    for (int j=0; j<Parameters.mapHeight; j++) {
      for (int i=0; i<Parameters.mapWidth; i++) {
        if (map[j][i] > maxVal) maxVal = map[j][i];      
        vals.add(new Float(map[j][i]));
      }
    }    

    println("@ Max val: " + maxVal);
    println("@ N vals: " + vals.size());

    for (int j=0; j<Parameters.mapHeight; j++) {
      for (int i=0; i<Parameters.mapWidth; i++) {
        map[j][i] /= maxVal;
      }
    }
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

// Author: Daniel C. Moura (daniel.c.moura@gmail.com)
// April 2015
// version 0.1


static class Log {
  final static boolean logOn = !true;
  
  public static void log(String s) {
    if (logOn) println(s);  
  }
  
  public static void fatalError(String s) {
    log(s);    
  }
}


class Profiler {
  private long startTime = 0;
  private long endTime   = 0;
  private long accTime   = 0;
  private long n         = 0;
  private String descr   = "";
  public static final boolean logOn = Log.logOn;

  Profiler(String descr) {
    this.descr = descr;
  }

  public void start(){
    if (!logOn) return;
    Log.log(descr + " : begin");
    this.startTime = millis();
  }

  public void end() {
    if (!logOn) return;
    this.endTime   = millis();
    //n++;
    accTime += this.endTime - this.startTime;
    Log.log(descr + " : end");
  }
  
  public void inc() {
    n++;
  }

  public float getAvgTime() {
    if (!logOn) return 0.0f;
    return accTime*1.0f/n;
  }
  
  public long getTotal() {
    if (!logOn) return 0;
    return accTime;
  }
  
  public String toString() {
    if (!logOn) return "";
    return "Profiling \t" + descr + "  (nruns=" + n + "): " + getAvgTime() + "   Total: " + getTotal();
  }
  
  public void print() {
    println(toString());
  }
}
// Author: Daniel C. Moura (daniel.c.moura@gmail.com)
// April 2015
// version 0.1


class MapBlender {
  private final int nBins = Parameters.nBins;
  
  MapBlender() {      
  }
  
  float[][] makeBlendingMatrix() {
    float[][] w = new float[nBins][nBins];
    for (int i=0;i<nBins;i++) 
      for (int j=0;j<nBins;j++)
        w[i][j] = blendFunction(i,j);
    return w;
  } 
  
  float blendFunction(int i, int j) {
    return (i==j ? 1.0f : 0.0f); //default blending (does not blend) 
  }
  
  MapBlender createBlender(int histBlendMode) {    
    switch (histBlendMode) {    
      case Parameters.NO_BLEND:            return this;
      case Parameters.BLEND_EVASIVE:       return new EvasiveMapBlender();
      case Parameters.BLEND_MODAL:
                                           return new ModalMapBlender();                                           
    }
    return null;
  }

}




class EvasiveMapBlender extends MapBlender {
  final float penalty = Parameters.blendPenalty/(Parameters.nBins-1);
  final float attraction = 1.0f;
  
  float blendFunction(int i, int j) {
    return (i==j ? attraction : -penalty); 
  }
  
}


class ModalMapBlender extends MapBlender {
  final float penalty = Parameters.blendPenalty;
  final float attraction = 1.0f;
  
  float blendFunction(int i, int j) {    
    final float f = cos((j-i)*TWO_PI/Parameters.nBins) - Parameters.modalBlendOffset;
    return (f>0?f*attraction:f*penalty);  
  }
  
}

// Author: Daniel C. Moura (daniel.c.moura@gmail.com)
// April 2015
// version 0.1



class Node {
  float xOrig, yOrig;
  float x, y;
  float xNorm, yNorm;
  String id;
  int idx;
  Graph gparent;
  ArrayList<Edge> edges;
  
  
  
  Node(float xOrig, float yOrig, String id, int idx, Graph gparent) {
    this.xOrig = xOrig;
    this.yOrig = yOrig;
    this.id = id;
    this.idx = idx;
    this.gparent = gparent;
    this.x=-1;
    this.y=-1;    
    if (!Parameters.pinNodes)
      edges = new ArrayList<Edge>();
  }
  
  
  void calcCoords(float minX, float maxX, float minY, float maxY) {
    xNorm = (xOrig-minX) / (maxX-minX);
    yNorm = (yOrig-minY) / (maxY-minY);
    x =  round((Parameters.mapWidth-6) * xNorm)+3; //always ensure a 3px margin
    y =  round((Parameters.mapHeight-6) * yNorm)+3;
    
    //println(id + ": " + x + " ; " + y); 
  }
  
  void recomputeNodePositionFromEdges() {
    float xacc=0, yacc=0;
    for (Edge edge : edges) {
        final float[] pos = edge.getNodePosition(this);
        yacc += pos[0];
        xacc += pos[1];
    }
    x = min(max(xacc / edges.size(), 1), Parameters.mapWidth-3);
    y = min(max(yacc / edges.size(), 1), Parameters.mapHeight-3);
    for (Edge edge : edges) {
        edge.updateNodePosition(this);        
    }
    
  }
  
  float dist2Node(Node a) {
    return dist(x,y,a.x,a.y);
  }
  
  
}
// Author: Daniel C. Moura (daniel.c.moura@gmail.com)
// April 2015
// version 0.1


static class Parameters {
  //CONSTANTS
  static final int NO_BLEND                  =  0;
  static final int BLEND_EVASIVE             =  1;
  static final int BLEND_MODAL               = 13;
  
  static final int BIN_ORIENTATION       = 10;
  static final int BIN_SOURCE_MATRIX     = 20; 
  static final int BIN_TARGET_MATRIX     = 21;
  static final int BIN_CLUSTERING        = 25;  

  
  
  static final int CM_HUE     = 3;
  static final int CM_DIR     = 5;  

  //data source
  static String dataFilename = "" + 1000 + ".random";
  //static final String dataFilename = "" + 200000  + "." + 100000 + ".random";  
  //static String dataFilename = "airlines_original.graphml";
  //static String dataFilename = "migrations.graphml";


  //screen
  static final int vizFrameRate = 1;
  static final int mapSize = 800;
  static final float margin = 0.005*3;
  static int mapWidth;
  static int mapHeight;
  static final boolean invertY = !true;

  //edge control points resampling
  static final float removeDist = mapSize * 0.005f;// * 4;//  *  1.48 ; //adjust this value 1.48 for random KDE, 4.0 for random
  static final float splitDist = removeDist * 3.0f; //adjust this value  

  static float separationFactor = mapSize * 0.005 *3; //proposed

  //bundling
  static int nBins = 1;
  static int binType = BIN_TARGET_MATRIX;
  //static int binType = BIN_ORIENTATION;
  //static final int binType = BIN_CLUSTERING;
  

  static final boolean pinNodes = true;

  //density estimation
  static int histBlendMode = NO_BLEND; //BLEND_MODAL;//BLEND_EVASIVE;
  static float kernelSigma = 2*0.015 * mapSize; //0.04 * mapSize; //1*0.015 * mapSize;//10.0f;  
  static int convolutionSteps = 3; //0 means perfect conv; btw 3..6 for fast conv
  static float blendPenalty = 0.25;
  static float modalBlendOffset = 0.1;

  //optimization 
  static float sigma2h = 2.0f;
  static float hMax = sigma2h*kernelSigma; //2*0.0125 * mapSize;  //10*1.5;//15.0f*3;
  static float hStep = 0.9f; //0.8f; //0.7

  static final boolean blindStep = false; //0.5 .. 0.9


  //smoothing
  static int smoothingIterations = 5*2;
  static float smoothingAlpha = 0.25*1 ;

  //drawing
  static int screenFactor = 1;
  static float edgeAlpha = 3.0f/255f; //5.0f/255;//0.05f;
  static boolean useStrokeWeight = true;
  static float strokeMinWeight =  1.0f*screenFactor;
  static float strokeMaxWeight = 1*10.0f*screenFactor;
  static boolean drawVoronoi = false;
  static float colorSaturation = (nBins == 1 ?  0.0f : 0.60f);
  static int colorMode = CM_HUE;
  static boolean invertColors = false;





  static void validate() {
    if (histBlendMode/10 > 0 &&  //GENERIC BLENDERS
    histBlendMode/10 != binType/10) //SPECIFIC BLENDERS
      Log.fatalError("Invalid parameters: blend mode not compatible with bin type");

    switch (histBlendMode) {
    case BLEND_EVASIVE:
      if (blindStep)
        println("Warning: using blind step with repelent forces, may generate out of map errors");
      break;
    }  

    if (kernelSigma<2)
      Log.fatalError("Invalid parameters: numberof convolution steps should be >= 2.0");

    if (binType == BIN_ORIENTATION && separationFactor == 0.0f)
      println("Warning: 0 separation factor for orientation bins");
    
  }
}



