```@setup phylinc
using PhyloNetworks
mkpath("../assets/figures")
fastafile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1_net.fasta")
startingtreefile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1_net.treefile")
```
This page provides examples and step-by-step instructions to guide analysts through estimating a network with phyLiNC (phylogenetic Likelihood Network from Concatenated data). We describe how to select an evolutionary model, group multiple samples or individuals, and run phyLiNC in parallel.
# Estimate a Network from Concatenated DNA Sequences   
Starting at a given starting tree or network, [`phyLiNC`] (@ref) uses a hill-climbing algorithm to search for the network that best fits the given DNA sequences. 
The search starts at or within a Geometric-distributed number of NNIs of the starting tree or network. Then, it proceeds with a hill-climbing search made up of nearest neighbor interchange moves, root change moves, add hybrid, and remove hybrid proposals. 
Throughout topology optimization, it optimizes evolutionary rates, rate variation across sites, branch lengths and inheritance inheritance weights. 
This search strategy runs ten times, by default, and returns the best network.

To estimate a phylogenetic network with [`phyLiNC`] (@ref), we need three pieces of information: concatenated DNA sequence data, a starting network, and a substitution model.

## Data
The DNA sequences should be in fasta format. We can find a small example called "h1.fasta" in the PhyloNetworks package in the examples folder. There is no need to read the data file into Julia ahead of time. To estimate a network, we simply need a path to the data file.
```julia
using PhyloNetworks
fastafile = joinpath(pkgdir(PhyloNetworks), "examples","h1_net.fasta")
```

## Starting Tree or Network
A starting network can be obtained using the analyst's knowledge of the field or with a tree-building tool like IQTree. It should be in Newick string format, as shown in the example "h1.tree" in the PhyloNetworks package in the examples folder. The starting tree or network cannot have more reticulations than the final estimated network: it should have at most `maxhybrid` reticulations. (If `maxhybrid` is not provided, then the starting network can have at most one reticulation.) The starting network can be of any level.

After selecting a starting tree or network, read it into the Julia session with PhyloNetwork's `ReadTopology()` function.
```julia
startingtreefile = joinpath(pkgdir(PhyloNetworks), "examples","h1_net.treefile")
startingtree = readTopology(startingtreefile)
```

If any branch lengths are missing in the input network, phyLiNC estimates the starting branch lengths using pairwise distances. Otherwise, it uses the input branch lengths as starting branch lengths, only unzipping all reticulations (as described in assumptions paragraph below).

## Estimate a Network
Now that we have data and a starting tree or network, we can estimate a phylogenetic network with the following [`phyLiNC`] (@ref) call, where the first argument is the starting tree, the second is the fasta file path, and the third is a choice of substitution model, either [`:HKY85`] (@ref) or [`JC69`] (@ref) (explained in detail in Modeling Evolutionary Rates, below).
```julia
fastafile = joinpath(pkgdir(PhyloNetworks), "examples","h1_net.fasta")
startingtreefile = joinpath(pkgdir(PhyloNetworks), "examples","h1_net.treefile")
startingtree = readTopology(startingtreefile)
netobject = PhyloNetworks.phyLiNC(startingtree, fastafile, :HKY85)
```
Note: Because Julia is a JIT (just-in-time) compiled programming language, the first use of a function takes longer than subsequent calls.

# Beyond Defaults: Customize Optimization with Optional Arguments
Only three arguments are required to estimate a network with PhyLiNC: a starting tree or network, a fasta file path, and a substitution model symbol. However, in many cases, analysts may want to move beyond this simple example for a more customized estimation. 
We use two types of optional arguments: positional and keyword arguments. Positional arguments appear before the semi-colon in the argument call and are identified by their position. Keyword arguments are placed after the semi-colon and must be appear with their name and an equal sign: maxhybrid = 1.

Positional Arguments (default value in parentheses):

- symbol for the model of rate variation across sites (`:noRV` for no rate variation). Use `:G` or `:Gamma` for Gamma-distributed rates, `:I` or `:Inv` for a proportion of invariable sites, or `:GI` or `:GammaInv` for a combination. See [`RateVariationAcrossSites`] (@ref).
- integer (4) for the number of categories to use in estimating
  evolutionary rates using a discretized gamma model. When allowing for rate
  variation, four categories is standard. With 1 category,
  no rate variation is assumed.

Optional keyword arguments for topology estimation and logging (default value in parentheses):

- `speciesfile` (""): path to a csv file with samples in rows and two columns:
  species (column 1), individual (column 2).
  Include this file to group individuals by species.
- `filename` ("phyLiNC"): root name for the output files (.out, .err). If empty (""), files are not created, progress log goes to the screen only (standard out).
- `maxhybrid` (1): the maximum number of hybridizations to include in the estimated network (allows one hybridization by default). 
Note: the starting network should have `maxhybrid` or fewer reticulations.
- `nruns` (10): number of independent starting points for the search
- `no3cycle` (true): prevents 3-cycles, which are (almost) not identifiable
- `nohybridladder` (true): prevents hybrid ladder in network. If true,
  the input network must not have hybrid ladders.
  Warning: Setting this to true will avoid most hybrid ladders, but some can still occur in some cases when deleting hybrid edges. (This might be replaced with an option to avoid all non-tree-child networks in the future.)


Optional keyword arguments to control the search and output (default value in parentheses):

- `seed` (0 for randomly-chosen seed): seed to replicate a given search
- `nreject` (75): Controls when to stop proposing new network topologies and, ultimately, the precision of the estimation.
This gives the maximum number of times that new topologies are proposed and rejected in a row. Lower values of nreject result in a less thorough but faster search. A higher value will lead to more precise estimation and higher computing times. When running test cases, use a small nreject such as 5.
- `probST` (0.5): probability to use net as the starting topology for each given run. If probST < 1, the starting topology is k NNI moves away from net, where k is drawn from a geometric distribution:
  <img src="https://render.githubusercontent.com/render/math?math=p(1-p)^{k}">
  with success probability p = probST.
- `maxmoves` (100): maximum number of topology moves before branch lengths, hybrid weight values, evolutionary rates, and rate variation parameters are reestimated.
- `verbose` (true): set to false to turn off screen output
- `alphamin` (0.02): minimum value for shape parameter alpha in rate variation across sites model.
- `alphamax` (50.0): maximum value for shape parameter alpha in rate variation across sites model.


Optional keyword arguments to control the optimization of branch lengths and inheritance weights (default value in parentheses):

- `ftolRel` (1e-6) and `ftolAbs` (1e-6): relative and absolute differences of the network score between the current and proposed parameters
- `xtolRel` (1e-5) and `xtolAbs` (1e-5): relative and absolute differences between the current and proposed parameters.
Greater values will result in a less thorough but faster search. These parameters are used when evaluating candidate networks only.
Regardless of these arguments, once a final topology is chosen, branch lengths are optimized using stricter tolerances (1e-10, 1e-12, 1e-10, 1e-10) for better estimates.


The following code block shows an example using both positional and keyword optional arguments.
```julia
obj_custom = PhyloNetworks.phyLiNC(startingtree, fastafile, :JC69, :G, 2; maxhybrid=3, nreject=5, verbose=true, seed=123)
```

# Modeling Evolutionary Rates
PhyLiNC allows analysts to choose from two Markov models for evolution and four ways to model rate variation. Together, these can be combined for eight ways to model evolutionary rates.

## Markov Model Substitution Models
There are two nucleic acid substitution models to choose from: [`JC69`] (@ref) and [`HKY85`] (@ref) . For each model, the transition matrix [`Q`] (@ref) is normalized to have an average of 1 mutation per unit of time.

- `:JC69` The Jukes & Cantor (1969) nucleic acid substitution model assumes a single evolutionary rate for all transitions between the characters A, C, G, and T.

- `:HKY85` The HKY85 nucleic acid substitution model lets transitions occur at the different rate than transversions (Hasegawa et al., 1985).
The model is made up of two components: a transition/transversion ratio, which we call κ or rate, and a vector of probabilities, called π. The rate presents the transition/transversion ratio: κ = α/β. The vector π shows the base frequencies of each character state (A, C, G, T) in the data. These are estimated from the data and must sum to 1.0.


PhyLiNC optimizes rate parameters throughout network estimation. After estimation is complete, it present the final transition rate matrix [`Q`] (@ref). The model is stored in the object that the function phyLiNC returns, which we have called netobj here. To explore the model further, use the following commands.
```julia
model = netobject.model
nstates(model)
getlabels(netobject.model)
```

## Rate Variation and Invariable Sites
By default, the model does not allow rates to vary across sites.
However, we can allow for more complexity by including an extra symbol argument in the function call. 
This turns either of the substitution models above into the "+G" or "+I" version: e.g. HKY85+G or JC69+I. 
To allow rates to vary across sites, we can add a discrete Gamma model (+G), developed by Yang et al (1994). To account for invariable sites, we can use the +I model developed by Hasegawa et al (1985). We can combine the two types of rate variation (+G+I, Gu1995) but this is discouraged (Jia2014). 

Using rate variation increases the number of parameters by one (+G or +I) or by two (+G+I).
Because the mean of the desired distribution or rates is 1, we use a Gamma distribution with shape parameter α and scale parameter θ=1/α (so β=α). 

For example, to allow for variable rates across sites using the discrete Gamma model, we add the `:G` symbol to our function call.
We can specify the number of categories for our discrete Gamma model by including the number after the rate variation model symbol. Four categories is standard. (If one category is specified, then rates are assumed not to vary.)
```julia
obj_G = PhyloNetworks.phyLiNC(startingtree, fastafile, :JC69, :G, 4)
```

Our model is well-suited to data with a high percentages of invariable sites. 
To model invariant sites, add +I to the Markov model by running phyLiNC with the `:I` symbol.
```julia
obj_I = PhyloNetworks.phyLiNC(startingtree, fastafile, :JC69, :I)
```

After estimating the network, phyLiNC returns the rate variation model in addition to the Markov model.
We can explore the details of our network's model in this way.
```julia
ratemodel = obj_G.ratemodel
nparams(ratemodel)
```

# Group Multiple Alleles or Individuals within a Species
We allow analysts to fix alleles or individuals within a species group. 
This constrains the output network to connect this set of individuals to one polytomy node representing the base of the species.

To group individuals within a species in this way, include a species mapping file path as a keyword argument. This file groups individuals by species and should have the following format: a csv file with samples in rows and two columns named "species" and "individual".
An example can be found in the examples folder ("mappingIndividuals.csv").

With a population-level starting network or tree, a fasta file path with a line for each individual, and a species mapping file path, we can estimate a network with grouped individuals with this call:
```julia
populationtreefile = joinpath(pkgdir(PhyloNetworks), "examples", "population_level.treefile")
populationtree = readTopology(populationtreefile)
individualfasta = joinpath(pkgdir(PhyloNetworks), "examples", "individuals.fasta")
mappingfile = joinpath(pkgdir(PhyloNetworks), "examples", "mappingIndividuals.csv")
obj_multipleindiv = PhyloNetworks.phyLiNC(populationtree, individualfasta, :JC69; maxhybrid=2, speciesfile=mappingfile)
```
# Estimate Networks in Parallel

By default, phyLiNC estimates ten networks, then returns the network with the highest likelihood. To speed this process, we can compute these ten runs in parallel.
For example, if our machine has multiple cores, we can tell Julia to add ten worker processes by starting Julia with the command 
```julia
julia -p 10
```
This starts Julia with eleven processes: one main and ten workers. Alternatively, we can start Julia as usual, then add processes. If we are unsure how many we have added, we can check using nworkers().
```julia
using Distributed
addprocs(10)
nworkers()
```
This allows PhyLiNC to estimate a network more quickly, completing the ten runs in parallel using these ten worker processes. 

# Assumptions
When estimating a network with phyLiNC, we make five major assumptions.
- Our model does not account for incomplete lineage sorting. 
- The sites within the DNA data are independent. 
- The final network will not include parallel edges, also known as two-cycles. These are not or nearly not-identifiable under our model, so we do not attempt to infer them. 
- The length of the edge below a reticulation is not identifiable, so we infer a network with reticulations unzipped: edges below reticulations are set to zero and hybrid edges (parental lineages) lengths are increased accordingly. 
- A hybrid node can have only two parent edges, meaning that each hybridization event can only occur between two species. However, we do allow for multiple hybridization events in succession (with keyword argument nohybridladder=false) for more complex hybridization events.
