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
  MapBlender: defines interactions between edges taking into account their bins
*/ 

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
    return (i==j ? 1.0f : 0.0f); //default blending (edges of different bins are ignored) 
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



//repelles edges from different bins
class EvasiveMapBlender extends MapBlender {
  final float penalty = Parameters.blendPenalty/(Parameters.nBins-1);
  final float attraction = 1.0f;
  
  float blendFunction(int i, int j) {
    return (i==j ? attraction : -penalty); 
  }
  
}

//attracttion/repellence when bins are modal (e.g. angles)
//maximum attraction when angle between bins is 0, maximum reppelence when angle is 180 deg
class ModalMapBlender extends MapBlender {
  final float penalty = Parameters.blendPenalty;
  final float attraction = 1.0f;
  
  float blendFunction(int i, int j) {    
    final float f = cos((j-i)*TWO_PI/Parameters.nBins) - Parameters.modalBlendOffset;
    return (f>0?f*attraction:f*penalty);  
  }
  
}

