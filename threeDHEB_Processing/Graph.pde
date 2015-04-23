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
  Graph: Data structure that stores the graph.
 */

class Graph {
  protected Node[] nodes; 
  protected Edge[] edges;

  protected color[] colors; //color of each cluster
  protected int[] clusterOrder;
  final int nThreads = 1;
  protected ImageFilter ffilter = null;

  protected float h, hrel, sigma; 

  KdeCpu[] kde;               //density maps
  float[][] blendingWeights;  //defined by the blending function

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
      h = Parameters.hMax();
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


    prof[0].start();//total

    final int edgesPerThread = ceil(edges.length*1.0f/nThreads);
    Thread[] t = new Thread[nThreads];

    prof[1].start();

    for (int i=0; i<edges.length; i++) {
      final Edge edge = edges[i];            
      edge.resampleControlPoints();
      if (Parameters.nBins > 1)                          
        edge.calculateBins();
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

    globalyNormalizeMaps();    

    if (drawNodes) {
      pg.noStroke();          
      pg.fill(0.5, 1-Parameters.colorSaturation, 1, 1);

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

      if (!Parameters.isJavaScript) 
        pg.beginShape();
        
      for (int j=0; j<edge.x.length; j++) {
        final int bin = (j==0?0:j-1);

        float dens=1.0;
        if (includeWeights) {
          dens = kde[edge.bins[bin]].getLocalDensity(round(edge.x[j]), round(edge.y[j]));
          if (dens <0) 
            dens = 0;

          pg.strokeWeight((Parameters.strokeMaxWeight-Parameters.strokeMinWeight) * dens + Parameters.strokeMinWeight);
        }


        switch (Parameters.colorMode) {

          case Parameters.CM_HUE: //hue based  */
  
            float hueVal = edge.bins[bin]*1.0f/Parameters.nBins;                              
            pg.stroke(hueVal, Parameters.colorSaturation, 1.0, edgeAlpha);
  
            break;
  
  
          case Parameters.CM_DIR: //direction
  
            float p = float(j)/(edge.x.length-1f);
  
            color r = min(255, round(max(0, p*2f) * 255f)); 
            color b = min(255, round(max(0, 2f*(1f-p)) * 255f));
            color g = int(min(r, b)*0.75);
  
            int c2 = 0xFF000000 | r << 16 | g << 8 |  b;
  
            pg.stroke(c2, edgeAlpha);
            break;
  
          default:
            pg.stroke(color(1, edgeAlpha));
        }
        
        if (!Parameters.isJavaScript)
          pg.vertex(edge.x[j]*s, edge.y[j]*s);
        else if (j>0)
          //varying width and colour seems not to work well in javascript
          //so, edges are rendered as several lines, instead of a polyline
          pg.line(edge.x[j-1]*s, edge.y[j-1]*s, edge.x[j]*s, edge.y[j]*s);
      }        

      if (!Parameters.isJavaScript)
        pg.endShape();
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


  /* TO DO: make compatible with javascript 
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
   */

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

  /* TO DO: make compatible with javascript
   
   void loadTxt(String filename) {        
   println("Loading txt graph...");
   ArrayList<Node> tmpNodes = new ArrayList<Node>(50000);
   LinkedList<Edge> tmpEdges = new LinkedList<Edge>();
   Hashtable<String, Integer> nodeId2IdxTable = new Hashtable<String, Integer>();
   
   maxX = -Float.MAX_VALUE;
   maxY = -Float.MAX_VALUE;
   minX = Float.MAX_VALUE;
   minY = Float.MAX_VALUE;
   
   maxW = maxD = -Float.MAX_VALUE;
   minW = minD =  Float.MAX_VALUE;
   
   int i=0;
   Scanner scanner;
   try {
   scanner =  new Scanner(new File(filename));//, ENCODING.name());
   }
   catch (Exception e) {
   println("Error reading file.");  
   throw new RuntimeException(e);
   }  
   while (scanner.hasNextLine ()) {
   final String line = scanner.nextLine();
   
   //println(line);
   final String []toks = line.split(" ");
   //println("L: " + toks.length + "-> " + line);
   //if (1==1)return; 
   if (toks[0].equals("node")) {        
   final String id = toks[1];
   final float x = Float.parseFloat(toks[2]);
   final float y = (Parameters.invertY?-1f:1f) * Float.parseFloat(toks[3]);
   tmpNodes.add(new Node(x, y, id, i, this));
   nodeId2IdxTable.put(id, i);  
   if (x>maxX) maxX=x;
   if (x<minX) minX=x;
   if (y>maxY) maxY=y;
   if (y<minY) minY=y;
   i++;
   
   //println("L: " + toks.length + "-> " + line);
   //if (1==1)return;
   } else if (toks[0].equals("edge")) {
   final String id = "";//toks[1];
   ;
   final String sourceId = toks[2-1];
   final String targetId = toks[3-1];
   final float weight = 1.0f;
   
   try {
   final Edge edge = new Edge(tmpNodes.get(nodeId2IdxTable.get(sourceId)), 
   tmpNodes.get(nodeId2IdxTable.get(targetId)), 
   weight, id, this);
   tmpEdges.add( edge );
   if (weight>maxW) maxW=weight;
   if (weight<minW) minW=weight;
   if (edge.distanceOrig>maxD) maxD=edge.distanceOrig;
   if (edge.distanceOrig<minD) minD=edge.distanceOrig;
   }
   catch (Exception e) {
   println("Error loading edge " + id);
   }
   }
   }     
   scanner.close(); 
   
   nodes = tmpNodes.toArray(new Node[0]);
   edges = tmpEdges.toArray(new Edge[0]);
   
   calcAspectRatio();
   
   println("Loaded " + nodes.length + " nodes.");
   println("Loaded " + edges.length + " edges.");
   }
   */

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
    /*else if (ext.equals("gexf"))
     loadGEXF(filename);
     else if (ext.equals("txt"))
     loadTxt(filename);*/
    else if (ext.equals("random"))
      switch (parts.length) {
      case 2: 
        makeRandomGraph(int(parts[0])); 
        break;
      case 3: 
        makeRandomGraph(int(parts[0]), int(parts[1])); 
        break;
      } else    
      loadGraphML(filename);

    prepareBundlind();
  }



  private void normalizeVals() {
    float xl = graph.maxX-graph.minX;
    float yl = graph.maxY-graph.minY;



    minX-=Parameters.margin*xl;
    maxX+=Parameters.margin*xl;
    minY-=Parameters.margin*yl;
    maxY+=Parameters.margin*yl;        

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

