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
  Parameters: most of the parameters and global vars of the algorithm and application 
*/


static class Parameters {
  //CONSTANTS
  static final int NO_BLEND              =  0;
  static final int BLEND_EVASIVE         =  1;
  static final int BLEND_MODAL           = 13;
  
  static final int BIN_ORIENTATION       = 10;
  static final int BIN_SOURCE_MATRIX     = 20; 
  static final int BIN_TARGET_MATRIX     = 21;  

  static final int CM_HUE                = 3;
  static final int CM_DIR                = 5;  
  
  //choose false if you are running a local Java app
  static final boolean isJavaScript = false;

  //data source
  static String dataFilename = "" + 1000 + ".random"; //geerates a random graph 
  //static final String dataFilename = "" + 200000  + "." + 100000 + ".random";  
  //static String dataFilename = "airlines_original.graphml";
  //static String dataFilename = "migrations.graphml";  
  

  //screen
  static final int vizFrameRate = 30;
  static final int mapSize      = 800;
  static final float margin     = 0.01;
  static final boolean invertY  = false;
  static int mapWidth;
  static int mapHeight;
  
  //edge control points resampling
  static final float removeDist  = mapSize * 0.005f;
  static final float splitDist   = removeDist * 3.0f;  

  //intial separation factor for spliting flows
  static float separationFactor = mapSize * 0.015; //proposed

  //bundling criteria
  static int nBins = 1;  
  static int binType = BIN_ORIENTATION;  
  
  //density estimation
  static int histBlendMode       = NO_BLEND; //BLEND_MODAL;//BLEND_EVASIVE;
  static float kernelSigma       = 0.015 * mapSize;   
  static int convolutionSteps    = 3; //iterations of the gaussian approximation function (btw 3..6) 
  static float blendPenalty      = 0.25;
  static float modalBlendOffset  = 0.1;

  //optimization 
  static float sigma2h           = 2.0f; //tipically inf [0.5, 3]
  static float hMax()            { return sigma2h*kernelSigma; }
  static float hStep             = 0.9f; //in [0.5, 1.0[
  static final boolean blindStep = false;
  
  //edge smoothing
  static int smoothingIterations = 10;
  static float smoothingAlpha    = 0.25;

  //drawing
  static int screenFactor        = 1;
  static float edgeAlpha         = 3f/255f;
  static boolean useStrokeWeight = true;
  static float strokeMinWeight   =  1.0f*screenFactor;
  static float strokeMaxWeight   = 1*10.0f*screenFactor;
  static float colorSaturation   = (nBins == 1 ?  0.0f : 0.60f);
  static int colorMode           = CM_HUE;
  static boolean invertColors    = false;





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
