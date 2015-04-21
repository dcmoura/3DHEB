#3DHEB
==========
3DHEB is a graph bundling algorithm that bundles edges of a graph to reveal connectivity patterns. 
3DHEB allows to tackle large graphs (tested on graphs with 5M edges) and has low computation times when compared to other bundling algorithms. 
In addition, it enables to specify the bundling criteria. Here we include a flow criteria based on edge orientation in addition to the default criteria based on edge proximity.

## Algorithm description
The algorithm is described here: [arxiv.org/abs/1504.02687](http://arxiv.org/abs/1504.02687).

## Implementation
Currently, a basic implementation of the algorithm is provided. The code is fully compatible with Processing (Java and JavaScript). 
This code is a simplified version of the code used in the [paper](http://arxiv.org/abs/1504.02687). Features requiring Java libraries ere not included 
(e.g. multi-threading, clustering and colour tables) as well as features that do not work well in processing.js 
(e.g. polylines with varying width and colour and hardware accelerated rendering). 
Thus, this implementation is slower than the original and has not got as many features, but still enables to performance very nice bundles out-of-the-box! 

## Online App
A [JavaScript app](http://dcmoura.github.io/3DHEB/) is provided for making your own graph bundles.
You may use one of the sample graphs or your own graph in GraphML format. 
The app includes an option to generate high resolution images that you may include in your publications. 
If you do so, please [cite me](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1504.02687). Thanks! 

[![Online app](http://dcmoura.github.io/3DHEB/img2/js_app.png)](http://dcmoura.github.io/3DHEB/)


## Gallery 
Some images produced with this code, and an animation based on an extension of this code.

* Animation based on 3DHEB:
[![Animation based on 3DHEB](http://dcmoura.github.io/3DHEB/img2/video_poster_3dheb_2.png)](https://vimeo.com/danielcmoura/3dheb)

* US Migrations Flow - 10K edges:
![US Migrations Flow](http://dcmoura.github.io/3DHEB/img/migrations_flow.png)

* US Airlines Flow - 1K edges:
![US Airlines Flow](http://dcmoura.github.io/3DHEB/img/airlines_flow.png)

* Pattern1 from the [University of Florida Sparse Matrix Collection](http://www.cise. ufl.edu/research/sparse/matrices/) - 5M edges:
<img src=http://dcmoura.github.io/3DHEB/img/pattern1_1.png width="67%" align="middle"/> 



