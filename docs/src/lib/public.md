# Public Documentation

Documentation for `PhyloNetworks`'s public (exported) interface.

See [Internal Documentation](@ref) for documentation on internal functions.

```@contents
Pages = ["public.md"]
```

## Index

```@index
Pages = ["public.md"]
```

## types

```@autodocs
Modules = [PhyloNetworks]
Private = false
Order   = [:type]
```

## utilities

```@docs
tipLabels
sorttaxa!
printEdges
printNodes
summarizeDataCF
directEdges!
checkroot!
preorder!
cladewiseorder!
rootatnode!
rootonedge!
hybridatnode!
setLength!
setGamma!
deleteleaf!
deleteHybridThreshold!
rotate!
getindex(::TraitSimulation, ::Symbol)
getNodeAges
pairwiseTaxonDistanceMatrix
biconnectedComponents
blobDecomposition
treeedgecomponents
getlabels
nparams
mapindividuals
nni!
```

## data and topology read/write

```@docs
readfastatodna
readTopology
readTopologyLevel1
readInputTrees
readMultiTopology
readNexusTrees
readSnaqNetwork
readTrees2CF
countquartetsintrees
readTableCF
readTableCF!
writeTableCF
readBootstrapTrees
writeSubTree!
writeTopology
writeMultiTopology
hybridlambdaformat
mapAllelesCFtable
```

## network inference

```@docs
snaq!
topologyMaxQPseudolik!
topologyQPseudolik!
fittedQuartetCF
bootsnaq
calibrateFromPairwiseDistances!
undirectedOtherNetworks
nj
```
## network Comparisons

```@docs
majorTree
minorTreeAt
displayedTrees
displayedNetworkAt!
hardwiredClusters
hardwiredCluster
hardwiredClusterDistance
treeEdgesBootstrap
hybridDetection
summarizeHFdf
hybridBootstrapSupport
```

## continuous trait evolution

```@docs
simulate
shiftHybrid
getShiftEdgeNumber
getShiftValue
phyloNetworklm
sigma2_estim
mu_estim
lambda_estim
ancestralStateReconstruction
expectations
predint
expectationsPlot
predintPlot
descendenceMatrix
regressorShift
regressorHybrid
sharedPathMatrix
vcv
```

## discrete trait evolution

```@docs
parsimonySoftwired
parsimonyGF
Q
randomTrait
randomTrait!
fitdiscrete
maxParsimonyNet
nstates
stationary
empiricalDNAfrequencies
phyLiNC
```
