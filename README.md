Species centrality and uniqueness in ecological networks
======================================

## Introduction

In the context of network ecology, the structural importance of a species for an ecosystem is usually assessed via one or more measures of centrality: species involved in a larger number of energetic pathways can be key to food web stability.

Another important facet of the role of a species in an ecosystem is the uniqueness of its food-web role. The local extinction of distinctively unique species may result in a greater functional diversity loss in the community, with great cascading effects.

Some evidence seems to support a positive correlation between uniqueness and centrality, some other goes against this hypothesis: the conundrum may depend on the effect of network shape, dimension and other ecological variables.

In this data project we will test the hypothesis that central species are ecologically unique. In this presentation we focus on the analysis a publicly available dataset of 14 food webs. 

We consider node degree as a basic, common, centrality measure and we define the ecological uniqueness of a species in a food web as its averaged distance to the other species in the network-trait space determined by the food web [Dalla Riva and Stouffer, 2016, Oikos](http://gvdr.github.io/science/Rdpg.

## Vignettes

In this repository there are two vignettes in the `Vignettes` folder.

`read_pajek` explains the translation from a common network format (_Pajek projects_) to the more handy _igraph_s R objects.

`Uniqueness_degree` guide the fitting of a linear correlation test between ecological uniqueness and degree centrality.

## Limitations:
This are the vignettes, not the manuscript ;-) Most ecological and mathematical details are left to the manuscript itself, and the analytic flow presented here is a selection of the wider analysis supporting the manuscript. In fact, in the manuscript we consider a variety of classical central measures, define a new one, explore the _contribution_ of a species to the functional diversity of the community and look at evolutionary aspects. Moreover, we consider other datasets with both unipartite and bipartite webs, requiring ad hoc pre-processing. Finally, we valide the results with thorough null model simulations. However, this presentation is _complete_ in the sense that goes from the raw webs---as they appear in the wild on the web---to a satisfactory, yet partial, result.

Author:

* [Giulio V. Dalla Riva][gvdr]

[gvdr]: http://gvdr.github.io/