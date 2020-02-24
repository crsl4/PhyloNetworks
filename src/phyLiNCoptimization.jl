# any change to these constants must be documented in phyLiNC!
const moveweights = Distributions.aweights([0.4, 0.2, 0.2, 0.2])
const likAbsAddHybLiNC = 0.5
const likAbsDelHybLiNC = -0.1
const alphaRASmin = 0.05
const alphaRASmax = 50.0

"""
    phyLiNC!(net::HybridNetwork, fastafile::String, modSymbol::Symbol;
             maxhybrid=1::Int64, no3cycle=true::Bool,
             nohybridladder=true::Bool, maxmoves=100::Int64,
             nreject=75::Int64, nruns=10::Int64,
             filename="phyLiNC"::AbstractString, verbose=false::Bool,
             seed=0::Int64, probST=0.5::Float64, NLoptMethod=:LD_MMA::Symbol,
             ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
             xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64,
             constraints=TopologyConstraint[]::Vector{TopologyConstraint},
             alphamin=alphaRASmin::Float64, alphamax=alphaRASmax::Float64)

Estimate a phylogenetic network from concatenated DNA data using
maximum likelihood, ignoring incomplete lineage sorting
(phyLiNC: phylogenetic Likelihood Network from Concatenated data).
The network is constrained to have `maxhybrid` reticulations at most,
but can be of any level.
The search starts at (or near) the network `net`,
using a local hill-climbing search to optimize the topology
(nearest-neighbor interchange moves, add hybridizations,
and remove hybridizations). Also optimized are evolutionary rates,
amount of rate variation across sites, branch lengths and inheritance γs.
This search strategy is run `nruns` times, and the best of the `nruns`
networks is returned.

Return a [`StatisticalSubstitutionModel`](@ref) object, say `obj`, which
contains the estimated network in `obj.net`.

The length of the edge below a reticulation is not identifiable.
Therefore, phyLiNC estimates the canonical version of the network: with
reticulations **unzipped**: edges below reticulations are set to 0, and
hybrid edges (parental lineages) have estimated lengths that are
increased accordingly.

Optional arguments include (default value in parenthesis):
- `nruns` (10): number of independent starting points for the search
- `filename` ("phyLiNC"): root name for the output files (`.out`, `.err`).
  If empty (""), files are *not* created, progress log goes to the screen only
  (standard out).
- `maxhybrid` (1): maximum number of hybridizations allowed
- `no3cycle` (true): prevents 3-cycles, which are (almost) not
  identifiable
- `nohybridladder` (true): prevents hybrid ladder in network. If true,
  the input network must not have hybrid ladders.
- `maxmoves` (100): maximum number of topology moves before branch lengths,
  hybrid γ values, evolutionary rates, and rate variation parameters are
  reestimated.
- `verbose` (false): if true, print information about the numerical optimization
- `seed` (default 0 to get it from the clock): seed to replicate a given search
- `probST` (0.5): probability to use `net` as the starting topology
  for each given run. If probST < 1, the starting topology is k NNI moves
  away from `net`, where k is drawn from a geometric distribution: p (1-p)ᵏ,
  with success probability p = `probST`.
- `constraints` (none): topology constraints to meet during the search,
  such as constrained clades or species groups.
  Created using [`TopologyConstraint`] (@ref)

The following optional arguments control when to stop the optimization of branch
lengths and gamma values on each individual candidate network. Defaults in
parentheses.
- `ftolRel` (1e-6) and `ftolAbs` (1e-6): relative and absolute differences of the
  network score between the current and proposed parameters
- `xtolRel` (1e-2) and `xtolAbs` (1e-3): relative and absolute differences between the
  current and proposed parameters.
Greater values will result in a less thorough but faster search. These parameters
are used when evaluating candidate networks only.

The following optional arguments control when to stop proposing new
network topologies:

- `nreject` (75): maximum number of times that new topologies are
  proposed and rejected in a row.
- `liktolAbs` (1e-6): the proposed network is accepted if its score is better
  than the current score by at least `liktolAbs`.
Lower values of `nreject` and greater values of `liktolAbs` would
result in a less thorough but faster search.

fixit: delete the description of liktolAbs: unused option
fixit: add as options and describe `likAbsAddHybLiNC` and `likAbsDelHybLiNC`
fixit: remove all defaults from phyLiNCone!, to avoid different defaults and
to avoid an incorrect documentation about defaults.
"""
function phyLiNC!(net::HybridNetwork, fastafile::String, modSymbol::Symbol;
                  maxhybrid=1::Int64, no3cycle=true::Bool,
                  nohybridladder=true::Bool,
                  constraints=TopologyConstraint[]::Vector{TopologyConstraint},
                  kwargs...)

    # create starting object for all runs
    obj = StatisticalSubstitutionModel(net, fastafile, modSymbol, maxhybrid)
    # check_matchtaxonnames (inside constructor) renumbers nodes so we recalc constraint numbers here
    # fixit: let the user provide a data frame with 2 columns: species and individuals,
    # perhaps another data frame for clade constraints,
    # then build the constraints inside phyLiNC here.
    updateconstraints!(constraints, obj.net)
    checknetwork_LiNC!(obj.net, maxhybrid, no3cycle, nohybridladder, constraints)
    startingBL!(obj.net, true, obj.trait, obj.siteweight) # true: to unzip

    phyLiNC!(obj; maxhybrid=maxhybrid, no3cycle=no3cycle, nohybridladder=nohybridladder,
                  kwargs...)
end
function phyLiNC!(obj::SSM;
                  maxhybrid=1::Int64, no3cycle=true::Bool,
                  nohybridladder=true::Bool, maxmoves=100::Int64, nreject=75::Int64,
                  nruns=10::Int64, filename="phyLiNC"::AbstractString, verbose=false::Bool,
                  seed=0::Int64, probST=0.5::Float64, NLoptMethod=:LD_MMA::Symbol,
                  ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
                  xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64,
                  constraints=TopologyConstraint[]::Vector{TopologyConstraint},
                  alphamin=alphaRASmin::Float64, alphamax=alphaRASmax::Float64)
    writelog = true
    writelog_1proc = true
    if filename != ""
        julialog = string(filename,".log")
        logfile = open(julialog,"w")
        juliaout = string(filename,".out")
        if Distributed.nprocs() == 1
            writelog_1proc = true
            juliaerr = string(filename,".err")
            errfile = open(juliaerr,"w")
        end
    else
      writelog = false
      logfile = stdout
    end
    str = """optimization of topology, BL and inheritance probabilities using:
              maxhybrid = $(maxhybrid),
              tolerance parameters: ftolRel=$(ftolRel), ftolAbs=$(ftolAbs),
                                    xtolAbs=$(xtolAbs), xtolRel=$(xtolRel).
              max number of consecutive failed proposals = $(nreject)"
             """
    str *= (writelog ? "filename for files: $(filename)\n" : "no output files\n")
    str *= "BEGIN: $(nruns) runs starting near $(writeTopology(obj.net))\n"
    if Distributed.nprocs()>1
        str *= "       using $(Distributed.nprocs()) processors\n"
    end
    if writelog
      write(logfile,str)
      flush(logfile)
    end
    print(stdout,str)
    print(stdout, Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s") * "\n")
    # if 1 proc: time printed to logfile at start of every run, not here.

    if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) # seed based on clock
    end
    if writelog
      write(logfile,"\nmain seed $(seed)\n")
      flush(logfile)
    else print(stdout,"\nmain seed $(seed)\n"); end
    Random.seed!(seed)
    seeds = [seed; round.(Integer,floor.(rand(nruns-1)*100000))]

    if writelog && !writelog_1proc
        for i in 1:nruns # workers won't write to logfile
            write(logfile, "seed: $(seeds[i]) for run $(i)\n")
        end
        flush(logfile)
    end

    tstart = time_ns()


    bestnet = Distributed.pmap(1:nruns) do i # for i in 1:nruns
        logstr = "seed: $(seeds[i]) for run $(i), $(Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s"))\n"
        print(stdout, logstr)
        msg = "\nBEGIN PhyLiNC for run $(i), seed $(seeds[i]) and maxhybrid $(maxhybrid)"
        if writelog_1proc # workers can't write on streams opened by master
            write(logfile, logstr * msg)
            flush(logfile)
        end
        verbose && print(stdout, msg)
        GC.gc()
        try
            objcopy = deepcopy(obj)
            # todo in future, only deepcopy net and other fields of obj that could be modified by phyLiNCone!
            # refactor to use SharedArrays for data?
            best = phyLiNCone!(objcopy, maxhybrid, no3cycle,
                            nohybridladder, maxmoves, nreject, verbose,
                            seeds[i], probST, constraints, NLoptMethod, ftolRel, ftolAbs,
                            xtolRel, xtolAbs, alphamin, alphamax)
            logstr *= "\nFINISHED PhyLiNC for run $(i), -loglik of best $(best.loglik)\n"
            verbose && print(stdout, logstr)
            if writelog_1proc
                logstr = writeTopology(best.net)
                logstr *= "\n---------------------\n"
                write(logfile, logstr)
                flush(logfile)
            end
            return best
        catch err
            msg = "\nERROR found on PhyLiNC for run $(i) seed $(seeds[i]): $(err)\n"
            logstr = msg * "\n---------------------\n"
            if writelog_1proc
                write(logfile, logstr)
                flush(logfile)
                write(errfile, msg)
                flush(errfile)
            end
            @warn msg # returns: nothing
        end
    end
    tend = time_ns() # in nanoseconds
    telapsed = round(convert(Int64, tend-tstart) * 1e-9, digits=2) # in seconds
    writelog_1proc && close(errfile)
    msg = "\n" * Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s")
    if writelog
        write(logfile, msg)
    elseif verbose
        print(stdout, msg)
    end
    filter!(n -> n !== nothing, bestnet) # remove "nothing", failed runs
    if length(bestnet) > 0
        ind = sortperm([n.loglik for n in bestnet])
        bestnet = bestnet[ind]
        maxNet = bestnet[1]::StatisticalSubstitutionModel # tell type to compiler
    else
        error("all runs failed")
    end
    return maxNet
end

"""
    phyLiNCone!(obj::SSM, maxhybrid=1::Int64, no3cycle=true::Bool,
                nohybridladder=true::Bool, maxmoves=100::Int64,
                nreject=75::Int64, verbose=false::Bool, seed=0::Int64,
                probST=0.5::Float64,
                constraints=TopologyConstraint[]::Vector{TopologyConstraint},
                NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
                ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64,
                xtolAbs=xAbsBL::Float64, alphamin=alphaRASmin::Float64,
                alphamax=alphaRASmax::Float64)

Estimate one phylogenetic network (or tree) from concatenated DNA data,
like [`phyLiNC!`](@ref), but doing one run only, and taking as input an
StatisticalSubstitutionModel object `obj`. The starting network is `obj.net`
and is assumed to meet all the requirements.

See [`phyLiNC!`](@ref) for optional arguments, which are positional arguments
here, and keyword arguments in `phyLiNC!`.
"""
function phyLiNCone!(obj::SSM, maxhybrid=1::Int64, no3cycle=true::Bool,
                    nohybridladder=true::Bool,
                    maxmoves=100::Int64, nreject=75::Int64, verbose=false::Bool,
                    seed=0::Int64, probST=0.5::Float64,
                    constraints=TopologyConstraint[]::Vector{TopologyConstraint},
                    NLoptMethod=:LD_MMA::Symbol,
                    ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
                    xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64,
                    alphamin=alphaRASmin::Float64, alphamax=alphaRASmax::Float64)

    Random.seed!(seed)
    if probST < 1.0 # modify starting tree by k nni moves (if possible), k=0 or more
        numNNI = rand(Geometric(probST)) # number of NNIs follows a geometric distribution: p (1 - p)^k
        for i in 1:numNNI
            nni_LiNC!(obj, no3cycle, nohybridladder, verbose, constraints)
        end
        # todo write to logfile here, pass logfile argument from phyLiNC! function
        # writelog && suc && write(logfile," changed starting topology by $numNNI attempted NNI move(s)\n")
    end
    # rough optimization of rates and alpha:
    fit!(obj; optimizeQ=true, optimizeRVAS=true, maxeval=20)

    done = false
    while !done # break out of this loop only is if nnmoves < mmaxmoves and rejections < nreject.
        done = optimizestructure!(obj, maxmoves, maxhybrid, no3cycle, nohybridladder,
                                  nreject, verbose, constraints)
        fit!(obj; optimizeQ=true, optimizeRVAS=true, ftolRel=1e-2, ftolAbs=1e-2,
             xtolRel=1e-1, xtolAbs=1e-2)
        optimizeBL_LiNC!(obj, obj.net.edge, verbose, 100, NLoptMethod,
                    ftolRel, ftolAbs, xtolRel, xtolAbs)
        optimizeallgammas_LiNC!(obj, verbose, 100, NLoptMethod,
                           ftolRel, ftolAbs, xtolRel, xtolAbs)
    end
    return obj
end

"""
    checknetwork_LiNC!(net::HybridNetwork, maxhybrid::Int64, no3cycle::Bool,
        nohybridladder::Bool,
        constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Check that `net` is an adequate starting network before phyLiNC:
Remove nodes of degree 2 (possibly including the root) and
unzip the network: set all edges below hybrid nodes to length zero.
According to user-given options, also check for 3-cycles, hybrid ladders,
and max number of hybrids.

```jldoctest
julia> maxhybrid = 3;

julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> preorder!(net) # for correct unzipping in checknetwork_LiNC!

julia> PhyloNetworks.checknetwork_LiNC!(net, maxhybrid, true, true)
HybridNetwork, Rooted Network
8 edges
8 nodes: 4 tips, 1 hybrid nodes, 3 internal tree nodes.
tip labels: A, B, C, D
((A:2.0,(B:0.0)#H1:1.1::0.9):1.5,(C:0.6,#H1:2.0::0.1):1.0,D:2.5);

```
"""
function checknetwork_LiNC!(net::HybridNetwork, maxhybrid::Int64, no3cycle::Bool,
    nohybridladder::Bool,
    constraints=TopologyConstraint[]::Vector{TopologyConstraint})

    if maxhybrid > 0
        net.numTaxa >= 3 ||
            error("cannot estimate hybridizations in topologies with 3 or fewer tips: $(net.numTaxa) tips here.")
    end
    # checks for polytomies, constraint violations, nodes of degree 2
    checkspeciesnetwork!(net, constraints) ||
        error("The species or clade constraints are not satisfied in the starting network.")
    if no3cycle
        !contain3cycles(net, no3cycle) || @warn("Options indicate there should
        be no 3-cycles in the returned network, but the input network contains
        one or more 3-cycles (after removing any nodes of degree 2 including the
        root). These 3-cycles have been removed.")
    end
    if nohybridladder
        !hashybridladder(net) || error("Options indicate there should be no
        hybrid ladders in the returned network, but the input network contains
        one or more hybrid ladders.")
    end
    if length(net.hybrid) > maxhybrid
        error("Options indicate a maximum of $(maxhybrid) reticulations, but
        the input network contains $(length(net.hybrid)) hybrid nodes. Please
        increase maxhybrid to $(length(net.hybrid)) or provide an input network
        with $(maxhybrid) or fewer reticulations.")
    end
    unzip_canonical!(net)
    return net
end

"""
    optimizestructure!(obj::SSM, maxmoves::Int64, maxhybrid::Int64,
        no3cycle::Bool, nohybridladder::Bool, nreject=75::Int64,
        verbose=false::Bool,
        constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Alternate nni moves, hybrid moves, and root changes. Optimizes local branch
lengths and hybrid gammas after each move, then decides whether or not to accept
the move by comparing likelihoods. After adding or removing a hybrid, updates
`SSM` object's displayed `trees` and their attributes.
Return `done` boolean indicating if number of rejections was met.

The percent of nni moves, hybrid moves, and root changes to be performed is
0.5, 0.3, and 0.2 respectively.

For a description of optional arguments, see [`phyLiNC`](@ref).

Assumptions:
- `checknetworkbeforeLiNC` and `discrete_corelikelihood!` have been called on
  `obj.net`.
- starting with a network without 2- and 3- cycles
  (checked by `checknetworkbeforeLiNC`)

Note: When removing a hybrid edge, always removes the minor edge.
"""
function optimizestructure!(obj::SSM, maxmoves::Int64, maxhybrid::Int64,
    no3cycle::Bool, nohybridladder::Bool, nreject=75::Int64,
    verbose=false::Bool,
    constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    nmoves = 0
    rejections = 0
    while nmoves < maxmoves && rejections < nreject # both should be true to continue
        currLik = obj.loglik
        movechoice = sample(["nni", "addhybrid", "deletehybrid", "root"], moveweights)
        if movechoice == "nni"
            nmoves += 1
            result = nni_LiNC!(obj, no3cycle,  nohybridladder, verbose, constraints)
            if isnothing(result) # no nni moves possible
                verbose && println("There are no nni moves possible in this network.")
            elseif result # move successful and accepted
                rejections = 0
            else # move made, rejected, and undone
                rejections += 1 # reset
            end
        elseif movechoice in ["addhybrid", "deletehybrid"]  # perform hybrid move
            if maxhybrid == 0
                @debug("The maximum number of hybrids allowed is $maxhybrid,
                so hybrid moves are not legal on this network.")
                #TODO in future update moveprobability here as in SNaQ
            elseif length(obj.net.hybrid) == 0
                movechoice = "addhybrid"
            elseif length(obj.net.hybrid) == maxhybrid
                movechoice = "deletehybrid"
            elseif length(obj.net.hybrid) > maxhybrid # this should never happen
                error("""The network has more hybrids than allowed. maxhybrid =
                 $maxhybrid, but network has $(obj.net.hybrid) hybrids.""")
            end # either move is possible
            if movechoice == "addhybrid"
                added = addhybridedgeLiNC!(obj, currLik, maxhybrid, no3cycle,
                            nohybridladder, verbose, constraints)
                nmoves += 1
                if isnothing(added)
                    verbose && println("Cannot add a hybrid to the network.")
                    movechoice = "add hybrid (unsuccessful attempt)"
                elseif added
                    rejections = 0 # reset
                else
                    movechoice = "add hybrid (but deleted afterward)"
                    rejections += 1
                end
            else # delete hybrid
                deleted = deletehybridedgeLiNC!(obj, currLik, maxhybrid,
                        no3cycle, nohybridladder, verbose, constraints)
                nmoves += 1
                if isnothing(deleted)
                    verbose && println("""Cannot delete a hybrid to the network
                     without violating a topology constraint.""")
                    movechoice = "delete hybrid (unsuccessful attempt)"
                elseif deleted
                    rejections = 0 # reset
                else
                    movechoice = "delete hybrid (but added back)"
                    rejections += 1
                end
            end
        else # change root (doesn't affect likelihood)
            originalroot = obj.net.root
            changednet = moveroot!(obj.net, constraints)
            nmoves += 1
            if !changednet
                @debug("Cannot perform a root change move on current network.")
                #? reduce likelihood of root move?
            end
        end
        verbose && println("""loglik = $(loglikelihood(obj)) after move of type
        $movechoice, $nmoves total moves, and $rejections rejected moves""")
    end
    return rejections >= nreject
end

"""
    nni_LiNC!(obj::SSM, no3cycle::Bool, nohybridladder::Bool,
             verbose::Bool,
             constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Loop over possible edges for a nearest-neighbor interchange move until one is
found. Performs move and compares the original and modified likelihoods to
decide whether to accept the move or not.
Return true if move accepted, false if move rejected. Return nothing if there
are no nni moves possible in the network.

Assumptions:
- called by [`optimizestructure!`](@ref) or [`phyLiNC!`](@ref)
"""
function nni_LiNC!(obj::SSM, no3cycle::Bool, nohybridladder::Bool,
                  verbose::Bool,
                  constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    currLik = obj.loglik
    edgefound = false
    blacklist = Edge[]
    while !edgefound # randomly select interior edge
        if length(blacklist) == length(obj.net.edge)
            return nothing
        end
        remainingedges = setdiff(obj.net.edge, blacklist) # remove already tried edges
        eindex = Random.rand(1:length(remainingedges))
        e1 = remainingedges[eindex]
        if !(e1 in blacklist) # else go back to top of nni while loop
            undoinfo = nni!(obj.net,e1,nohybridladder,no3cycle,constraints)
            if !isnothing(undoinfo)
                edgefound = true
                discrete_corelikelihood!(obj)
                optimizelocalBL_LiNC!(obj, e1, verbose)
                optimizelocalgammas_LiNC!(obj, e1, verbose)
                if obj.loglik - currLik < likAbs
                    nni!(undoinfo...) # undo move
                    return false # result = false
                else
                    return true # rejections = 0 # resets to zero
                end
            else # if move unsuccessful, search for edge until successful
                push!(blacklist, e1)
            end
        end
    end
end

"""
    addhybridedgeLiNC!(obj::SSM, currLik::Float64, maxhybrid::Int64,
        no3cycle::Bool, nohybridladder::Bool, verbose::Bool,
        constraints::Vector{TopologyConstraint})

Completes checks, adds hybrid in a random location, updates SSM object, and
optimizes branch lengths and gammas locally as part of PhyLiNC optimization.

Return true if accepted add hybrid move. If move not accepted, return false.
If cannot add a hybrid, return nothing.

Assumptions:
- called by [`optimizestructure!`](@ref)
"""
function addhybridedgeLiNC!(obj::SSM, currLik::Float64, maxhybrid::Int64,
    no3cycle::Bool, nohybridladder::Bool, verbose::Bool,
    constraints::Vector{TopologyConstraint})
    orignet = deepcopy(obj.net) # hold old network in case we remove new hybrid
        #? in future use the same memory space for this every time?
    result = addhybridedge!(obj.net, nohybridladder, no3cycle, constraints)
    if !isnothing(result)
        newhybridnode, newhybridedge = result
        updateSSM!(obj)
        discrete_corelikelihood!(obj)
        optimizelocalBL_LiNC!(obj, newhybridedge, verbose)
        optimizelocalgammas_LiNC!(obj, newhybridedge, verbose)
        if obj.loglik - currLik < likAbsAddHybLiNC # improvement too small or negative: undo
            obj.net = orignet
            updateSSM!(obj)
            discrete_corelikelihood!(obj)
            return false
        else
            return true
        end
    else
        return nothing
    end
end

"""
    deletehybridedgeLiNC!(obj::SSM, currLik::Float64, maxhybrid::Int64,
        no3cycle::Bool, nohybridladder::Bool, verbose::Bool,
        constraints::Vector{TopologyConstraint})

Deletes a random hybrid edge, completes checks, and updates SSM object as part of
PhyLiNC optimization.

Return true if accepted delete hybrid move. If move not accepted, return false.

Assumptions:
- called by [`optimizestructure!`](@ref) which does some checks.
"""
function deletehybridedgeLiNC!(obj::SSM, currLik::Float64, maxhybrid::Int64,
    no3cycle::Bool, nohybridladder::Bool, verbose::Bool,
    constraints::Vector{TopologyConstraint})
    hybridnode = obj.net.hybrid[Random.rand(1:length(obj.net.hybrid))]
    minorhybridedge = getMinorParentEdge(hybridnode)
    if length(constraints) > 0 # check constraints
        edgefound = false
        blacklist = Node[]
        while !edgefound
            if length(blacklist) < obj.net.numHybrids
                verbose && println("There are no delete hybrid moves possible in this network.")
            end
            hybridnode = obj.net.hybrid[Random.rand(1:length(obj.net.hybrid))]
            if !(hybridnode in blacklist)
                edgefound = true
                for c in constraints
                    if minorhybridedge == c.edge # edge to remove is stem edge
                        push!(blacklist, hybridnode)
                        edgefound = false
                        break # out of constraint loop
                    end
                end
            end
        end
        if !edgefound # tried all hybrids
            return nothing
        end
    end
    savededges, savedhybridedges = savelocalBLgamma(obj.net, minorhybridedge)
    setGamma!(minorhybridedge, 0.0)
    discrete_corelikelihood!(obj)
    optimizelocalBL_LiNC!(obj, minorhybridedge, verbose)
    optimizelocalgammas_LiNC!(obj, minorhybridedge, verbose)
    if obj.loglik - currLik > likAbsDelHybLiNC # -0.1: loglik can decrease for parsimony
        deletehybridedge!(obj.net, minorhybridedge, false, true) # don't keep nodes; unroot
        updateSSM!(obj, true; constraints=constraints)
        discrete_corelikelihood!(obj)
        return true
    else # keep hybrid
        resetlocalBLgamma!(obj.net, savededges, savedhybridedges)
        obj.loglik = currLik
        return false
    end
end
"""
    updateSSM!(obj::SSM, renumber=false::Bool;
               constraints=TopologyConstraint[]::Vector{TopologyConstraint})

After adding or removing a hybrid, displayed trees will change. Updates
the displayed tree list. Return SSM object.

if `renumber`, reorder edge and internal node numbers. Only need
to renumber after deleting a hybrid (which could remove edges and nodes
from the middle of the edge and node lists).

Assumptions:
- The SSM object has cache arrays of size large enough, that is,
  the constructor [`StatisticalSubstitutionModel`](@ref) was previously
  called with maxhybrid equal or greater than in `obj.net`.
  `obj.priorltw` is not part of the "cache" arrays.

Warning:
Does not update the likelihood.

```jldoctest
julia> maxhybrid = 3;

julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69, maxhybrid);

julia> PhyloNetworks.checknetwork_LiNC!(obj.net, maxhybrid, true, true);

julia> using Random; Random.seed!(432);

julia> PhyloNetworks.addhybridedge!(obj.net, obj.net.edge[8], obj.net.edge[1], true, 0.0, 0.4);

julia> PhyloNetworks.updateSSM!(obj);

julia> writeTopology(obj.net)
"(((B:0.0)#H1:1.1::0.9,(A:1.0)#H2:1.0::0.6):1.5,(C:0.6,#H1:2.0::0.1):1.0,(D:1.25,#H2:0.0::0.4):1.25);"
```
"""
function updateSSM!(obj::SSM, renumber=false::Bool;
                   constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    if renumber # traits are in leaf.number order, so leaf nodes not reordered
        resetNodeNumbers!(obj.net; checkPreorder=false, type=:internalonly)
        resetEdgeNumbers!(obj.net, false) # verbose=false
        updateconstraints!(constraints, obj.net)
    end
    # extract displayed trees
    obj.displayedtree = displayedTrees(obj.net, 0.0; nofuse=true)
    nnodes = length(obj.net.node)
    for tree in obj.displayedtree
        preorder!(tree) # no need to call directEdges! before: already done on net
        #core likelihood uses nodes_changed to traverse tree in post-order
        # length(tree.nodes_changed) == nnodes ||
        #     error("displayed tree with too few nodes: $(writeTopology(tree))")
        # length(tree.edge) == length(obj.net.edge)-obj.net.numHybrids ||
        #     error("displayed tree with too few edges: $(writeTopology(tree))")
        # allow this because, in some cases, we remove a hybrid node during displayedtrees because it has no children.
    end
    # log tree weights: sum log(γ) over edges, for each displayed tree
    obj.priorltw = inheritanceWeight.(obj.displayedtree)
    @debug begin
        all(!ismissing, obj.priorltw) ? "" :
        "one or more inheritance γ's are missing or negative. fix using setGamma!(network, edge)"
    end
    return obj
end

## Optimize Branch Lengths and Gammas ##
"""
    startingBL!(net::HybridNetwork, unzip::Bool,
                trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
                siteweight=ones(length(trait[1]))::AbstractVector{Float64})

Calibrate branch lengths in `net` by minimizing the mean squared error
between the JC-adjusted pairwise distance between taxa, and network-predicted
pairwise distances, using [`calibrateFromPairwiseDistances!`](@ref).
`siteweight[k]` gives the weight of site (or site pattern) `k` (default: all 1s).
`unzip` = true sets all edges below a hybrid node to length zero.

Assumptions:

- all species have the same number of traits (sites): `length(trait[i])` constant
- `trait[i]` is for leaf with `node.number = i` in `net`, and
  `trait[i][j] = k` means that leaf number `i` has state index `k` for trait `j`.
  These indices are those used in a substitution model:
  kth value of `getlabels(model)`.
- Hamming distances are < 0.75 with four states, or < (n-1)/n for n states.
  If not, all pairwise hamming distances are scaled by `.75/(m*1.01)` where `m`
  is the maximum observed hamming distance, to make them all < 0.75.
"""
function startingBL!(net::HybridNetwork, unzip::Bool,
        trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
        siteweight=ones(length(trait[1]))::AbstractVector{Float64})
    nspecies = net.numTaxa
    M = zeros(Float64, nspecies, nspecies) # pairwise distances initialized to 0
    # count pairwise differences, then multiply by pattern weight
    ncols = length(trait[1]) # assumption: all species have same # columns
    length(siteweight) == ncols ||
      error("$(length(siteweight)) site weights but $ncols columns in the data")
    for i in 2:nspecies
        species1 = trait[i]
        for j in 1:(i-1)
            species2 = trait[j]
            for col in 1:ncols
                if !(ismissing(species1[col]) || ismissing(species2[col])) &&
                    (species1[col] != species2[col])
                    M[i, j] += siteweight[col]
                end
            end
            M[j,i] = M[i,j]
        end
    end
    Mp = M ./ sum(siteweight) # to get proportion of sites, for each pair

    # estimate pairwise evolutionary distances using extended Jukes Cantor model
    nstates = mapreduce(x -> maximum(skipmissing(x)), max, trait)
    maxdist = (nstates-1)/nstates
    Mp[:] = Mp ./ max(maxdist, maximum(Mp*1.01)) # values in [0,0.9901]: log(1-Mp) well defined
    dhat = - maxdist .* log.( 1.0 .- Mp)

    taxonnames = [net.leaf[i].name for i in sortperm([n.number for n in net.leaf])]
    # taxon names: to tell the calibration that row i of dhat if for taxonnames[i]
    # ASSUMPTION: trait[i][j] = trait j for taxon at node number i: 'node.number' = i
    calibrateFromPairwiseDistances!(net, dhat, taxonnames,
        forceMinorLength0=false, ultrametric=false)

    if unzip
        unzip_canonical!(net)
    end
    return net
end


"""
    optimizelocalBL_LiNC!(obj::SSM, edge::Edge, verbose)

Optimize branch lengths in `net` locally around `edge`. Update all edges that
share a node with `edge` (including itself).
Constrains branch lengths to zero below hybrid nodes.
Return vector of updated `edges`.

Used after `nni!` or `addhybridedge!` moves to update local branch lengths.

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
HybridNetwork, Rooted Network
9 edges
9 nodes: 4 tips, 1 hybrid nodes, 4 internal tree nodes.
tip labels: A, B, C, D
(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);


julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);

julia> obj.net.edge[4]
PhyloNetworks.Edge:
 number:4
 length:1.5
 attached to 2 node(s) (parent first): 6 8


julia> PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[4], true);

julia> obj.net.edge[4]
PhyloNetworks.Edge:
 number:4
 length:0.0
 attached to 2 node(s) (parent first): 6 8
```
"""
function optimizelocalBL_LiNC!(obj::SSM, edge::Edge, verbose=false::Bool)
    edges = Edge[]
    for n in edge.node
        for e in n.edge # all edges sharing a node with `edge` (including self)
            if !(e in edges)
                push!(edges, e)
            end
        end
    end
    optimizeBL_LiNC!(obj, edges, verbose, 10, :LD_MMA, fRelBL, fAbsBL,
                    xRelBL, xAbsBL) # maxeval = 10
    return edges
end

"""
    optimizeBL_LiNC!(obj::SSM, edges::Vector{Edge},
                verbose=false::Bool, maxeval=1000::Int64,
                NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
                ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64,
                xtolAbs=xAbsBL::Float64)

Optimize branch lengths for edges in vector `edges`.
Constrains branch lengths to zero below hybrid nodes.
Return vector of updated `edges`.

Assumption: None of the branch length are negative.
"""
function optimizeBL_LiNC!(obj::SSM, edges::Vector{Edge},
    verbose=false::Bool, maxeval=1000::Int64, NLoptMethod=:LD_MMA::Symbol,
    ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64,
    xtolAbs=xAbsBL::Float64)

    if !isempty(obj.net.hybrid)
        constrainededges = unzip_canonical!(obj.net)
        edges = setdiff(edges, constrainededges) # edges - constrainededges
    end
    counter = [0]
    function loglikfunBL(lengths::Vector{Float64}, grad::Vector{Float64})
        counter[1] += 1
        setlengths!(edges, lengths) # set lengths in order of vector `edges`
        res = discrete_corelikelihood!(obj)
        verbose && println("loglik: $res, branch lengths: $(lengths)")
        isempty(grad) || error("gradient not implemented")
        return res
    end
    # set-up optimization object for BL parameter
    NLoptMethod=:LN_COBYLA # no gradient
    # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
    nparBL = length(edges)
    optBL = NLopt.Opt(NLoptMethod, nparBL)
    NLopt.ftol_rel!(optBL,ftolRel) # relative criterion
    NLopt.ftol_abs!(optBL,ftolAbs) # absolute criterion
    NLopt.xtol_rel!(optBL,xtolRel)
    NLopt.xtol_abs!(optBL,xtolAbs)
    NLopt.maxeval!(optBL, maxeval) # max number of iterations
    # NLopt.maxtime!(optBL, t::Real)
    NLopt.lower_bounds!(optBL, zeros(length(edges)))
    NLopt.max_objective!(optBL, loglikfunBL)
    fmax, xmax, ret = NLopt.optimize(optBL, getlengths(edges)) # get lengths in order of edges vector
    verbose && println("""BL: got $(round(fmax, digits=5)) at
    BL = $(round.(xmax, digits=5)) after $(counter[1]) iterations
    (return code $(ret))""")
    return edges
end

"""
    optimizeallgammas_LiNC!(obj::SSM,
                        verbose::Bool, NLoptMethod::Symbol,
                        ftolRel::Float64, ftolAbs::Float64,
                        xtolRel::Float64, xtolAbs::Float64)

Optimize all gammas in a network. Creates a list containing one parent edge
per hybrid then calls optimizegammas_LiNC! on that list.
"""
function optimizeallgammas_LiNC!(obj::SSM,
    verbose::Bool, maxeval::Int64, NLoptMethod::Symbol, ftolRel::Float64,
    ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64)

    edges = [getMajorParentEdge(h) for h in obj.net.hybrid]
    if isempty(edges)
        @debug "no gammas to optimize, $(length(obj.net.hybrid)) hybrids"
    else
        optimizegammas_LiNC!(obj, edges, verbose, maxeval, NLoptMethod,
                ftolRel, ftolAbs, xtolRel, xtolAbs)
    end
end

"""
    optimizelocalgammas_LiNC!(obj::SSM, edge::Edge, verbose=false::Bool)

Optimize gammas in `net` locally around `edge`. Update all edges that share a
node with `edge` (including itself). Does not include edges' partners because
the `setGamma!` updates partners automatically.

Return modified edges.

Used after `nni!` or `addhybridedge!` moves to update local gammas.

Assumptions:
- correct `isChild1` field for `edge` and for hybrid edges
- no in-coming polytomy: a node has 0, 1 or 2 parents, no more

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);

julia> obj.net.hybrid[1].edge
3-element Array{PhyloNetworks.Edge,1}:
 PhyloNetworks.Edge:
 number:2
 length:1.0
 attached to 2 node(s) (parent first): 9 2

 PhyloNetworks.Edge:
 number:3
 length:0.1
 major hybrid edge with gamma=0.9
 attached to 2 node(s) (parent first): 8 9

 PhyloNetworks.Edge:
 number:6
 length:1.0
 minor hybrid edge with gamma=0.1
 attached to 2 node(s) (parent first): 7 9

julia> PhyloNetworks.optimizelocalgammas_LiNC!(obj, obj.net.hybrid[1].edge[2], true);

julia> obj.net.hybrid[1].edge
3-element Array{PhyloNetworks.Edge,1}:
 PhyloNetworks.Edge:
 number:2
 length:1.0
 attached to 2 node(s) (parent first): 9 2

 PhyloNetworks.Edge:
 number:3
 length:0.1
 major hybrid edge with gamma=0.9000000000259154
 attached to 2 node(s) (parent first): 8 9

 PhyloNetworks.Edge:
 number:6
 length:1.0
 minor hybrid edge with gamma=0.0999999999740846
 attached to 2 node(s) (parent first): 7 9

````
"""
function optimizelocalgammas_LiNC!(obj::SSM, edge::Edge, verbose::Bool)
    edges = Edge[]
    for n in edge.node
        for e in n.edge # edges that share a node with `edge` (including self)
            if e.hybrid && !(e in edges) && !(getPartner(e) in edges)
                push!(edges, e)
            end
        end
    end
    if isempty(edges)
        @debug "no local gammas to optimize around edge $edge"
    else
        optimizegammas_LiNC!(obj, edges, verbose, 10, :LD_MMA, fRelBL,
            fAbsBL, xRelBL, xAbsBL) # maxeval = 20
    end
end

"""
    optimizegammas_LiNC!(obj::SSM, edges::Vector{Edge},
                    verbose=false::Bool, maxeval=1000::Int64,
                    NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
                    ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64,
                    xtolAbs=xAbsBL::Float64)

Optimize gammas for hybrid edges in vector `edges`.
Return vector of updated `edges`.

Assumption: `edges` vector does not contain hybrid partners.

Warning: Do not call directly.
Instead use [`optimizelocalgammas_LiNC!`](@ref) or [`optimizeallgamma!`](@ref).
"""
function optimizegammas_LiNC!(obj::SSM, edges::Vector{Edge},
    verbose=false::Bool, maxeval=1000::Int64,
    NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
    ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64)

    counter = [0]
    function loglikfungamma(gammas::Vector{Float64}, grad::Vector{Float64})
        counter[1] += 1
        setmultiplegammas!(edges, gammas)
        res = discrete_corelikelihood!(obj)
        # verbose && println("loglik: $res, gammas: $(gammas)")
        isempty(grad) || error("gradient not implemented")
        return res
    end
    # set-up optimization object for gamma parameter
    NLoptMethod=:LN_COBYLA # no gradient
    # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
    npargamma = length(edges)
    optgamma = NLopt.Opt(NLoptMethod, npargamma)
    NLopt.ftol_rel!(optgamma,ftolRel) # relative criterion
    NLopt.ftol_abs!(optgamma,ftolAbs) # absolute criterion
    NLopt.xtol_rel!(optgamma,xtolRel)
    NLopt.xtol_abs!(optgamma,xtolAbs)
    NLopt.maxeval!(optgamma, maxeval) # max number of iterations
    # NLopt.initial_step!(optgamma, 0.05) # step size
    # NLopt.maxtime!(optgamma, t::Real)
    NLopt.lower_bounds!(optgamma, zeros(Float64, npargamma))
    NLopt.upper_bounds!(optgamma, ones(Float64, npargamma))
    counter[1] = 0
    NLopt.max_objective!(optgamma, loglikfungamma)
    fmax, xmax, ret = NLopt.optimize(optgamma, [e.gamma for e in edges])
    verbose && println("gamma: got $(round(fmax, digits=5)) at
    $(round.(xmax, digits=5)) after $(counter[1]) iterations
    (return code $(ret))")
    return edges
end

## Prep Functions ##

"""
    savelocalBLgamma(net::HybridNetwork, edge::Edge)

Saves local branch lengths and gammas before they're optimized so they can be
reset.
Return a tuple of two dictionaries holding edges and branch lengths
and edges and gammas, respectively. Each dictionary is of this format:
Dict{Edge, Float64}
"""
function savelocalBLgamma(net::HybridNetwork, edge::Edge)
    localedges = Dict{Edge, Float64}()
    for n in edge.node
        for e in n.edge # all edges sharing a node with `edge` (including self)
            if !(haskey(localedges, e))
                localedges[e] = e.length
            end
        end
    end
    localhybridedges = Dict{Edge, Float64}()
    for n in edge.node
        for e in n.edge # edges that share a node with `edge` (including self)
            if e.hybrid && !haskey(localhybridedges, e) && !haskey(localhybridedges, getPartner(e))
                localhybridedges[e] = e.gamma
            end
        end
    end
    return localedges, localhybridedges
end

"""
    resetlocalBLgamma!(net::HybridNetwork, localedges::Dict{Edge, Float64},
                        localhybridedges::Dict{Edge, Float64})

Reset local branch lengths and gammas to undo a local optimization.
Return net.
"""
function resetlocalBLgamma!(net::HybridNetwork, localedges::Dict{Edge, Float64},
    localhybridedges::Dict{Edge, Float64})
    for edge in keys(localedges) # set lengths
        edge.length = localedges[edge]
    end
    for hybridedge in keys(localhybridedges) # set gammas
        hybridedge.gamma = localhybridedges[hybridedge]
    end
    return net
end
