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


/*
Code adapted to Processing by Daniel C. Moura    
daniel.c.moura@gmail.com
*/


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

