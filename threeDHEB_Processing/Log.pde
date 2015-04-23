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
  Utility classes for logging and profiling. 
*/ 

static class Log {
  final static boolean logOn = true;
  
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
    //Log.log(descr + " : begin");
    this.startTime = millis();
  }

  public void end() {
    if (!logOn) return;
    this.endTime   = millis();    
    accTime += this.endTime - this.startTime;
    //Log.log(descr + " : end");
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
