#=
constraints:
- species groups e.g. as polytomies (type 1), or
- clades in the major tree (type 2), or
- clades in *some* displayed tree (type 3, not implemented).
outgroups "grades" can be coded as ingroup clades

move types, all under topological constraints:
- NNI = nearest neighbor interchange for undirected network,
  similar to NNIs for rooted networks defined in
  Gambette, van Iersel, Jones, Lafond, Pardi and Scornavacca (2017),
  Rearrangement moves on rooted phylogenetic networks.
  PLOS Computational Biology 13(8):e1005611.
- move the root: needed for optimization of a semi-directed network,
  because the root position affects the feasibility of the
  NNIs starting from a BR configuration (bifurcation -> reticulation).

- remove a hybrid edge: `deletehybridedge!` in file deleteHybrid.jl,
  but does not check for clade constraints
- add a hybrid edge: in file addHybrid.jl
- change the direction of a hybrid edge

in `moves_snaq.jl`, functions are tailored to level-1 networks;
here, functions apply to semidirected networks of all levels.
=#

"""
Type for various topological constraints, such as:

1. a set of taxa forming a clade in the major tree
2. a set of individuals belonging to the same species
3. a set of taxa forming a clade in any one of the displayed trees
"""
mutable struct TopologyConstraint
    "type of constraint: 1=species, 2=clade in major tree. 3=clade in some displayed tree: not yet implemented."
    type::UInt8
    "names of taxa in the constraint (clade / species members)"
    taxonnames::Vector{String}
    """
    node numbers of taxa in the constraint (clade / species members).
    warning: interpretation dependent on the internal network representation
    """
    taxonnums::Set{Int}
    "stem edge of the constrained group"
    edge::Edge
    "crown node, that is, child node of stem edge"
    node::Node
end

"""
    TopologyConstraint(type::UInt8, taxonnames::Vector{String}, net::HybridNetwork)

Create a topology constraint from user-given type, taxon names, and network.
There are 3 types of constraints:

- type 1: A set of tips that forms a species.
- type 2: A set of tips that forms a clade in the major tree.
  Note that the root matters. Constraining a set of species to be an outgroup
  is equivalent to constraining the ingroup to form a clade.
- type 3: A set of tips that forms a clade in any one of the displayed trees.

Note: currently, with type-1 constraints, hybridizations are prevented
from coming into or going out of a species group.

# examples

```jldoctest
julia> net = readTopology("(((2a,2b),(((((1a,1b,1c),4),(5)#H1),(#H1,(6,7))))#H2),(#H2,10));");

julia> c_species1 = PhyloNetworks.TopologyConstraint(0x01, ["1a","1b","1c"], net)
Species constraint, on tips: 1a, 1b, 1c
 stem edge number 7
 crown node number -9

julia> c_species2 = PhyloNetworks.TopologyConstraint(0x01, ["2a","2b"], net)
Species constraint, on tips: 2a, 2b
 stem edge number 3
 crown node number -4

julia> nni!(net , net.edge[3], true, true, [c_species2]) === nothing # we get nothing: would break species 2
true

julia> # the following gives an error, because 4,5 do not form a clade in "net"
       # PhyloNetworks.TopologyConstraint(0x02, ["4","5"], net)

julia> c_clade145 = PhyloNetworks.TopologyConstraint(0x02, ["1a","1b","1c","4","5"], net)
Clade constraint, on tips: 1a, 1b, 1c, 4, 5
 stem edge number 12
 crown node number -7

julia> nni!(net , net.edge[12], true, true, [c_species1, c_clade145]) # we get nothing (failed NNI): would break the clade
```
"""
function TopologyConstraint(type::UInt8, taxonnames::Vector{String}, net::HybridNetwork)

    ntax_clade = length(taxonnames)
    ntax_clade >= 2 ||
        error("there must be 2 or more taxon name in a clade or species constraint.")
    taxonnums = Set{Int}()
    indices = findall(in(taxonnames), [n.name for n in net.leaf])
    length(indices) == ntax_clade || error("taxa cannot be matched to nodes. Check for typos.")
    outsideclade = setdiff([n.name for n in net.leaf], taxonnames)
    # find the set of node numbers that correspond to taxonnames
    for (i,taxon) in enumerate(taxonnames)
        index = findfirst(n -> n.name == taxon, net.leaf)
        !isnothing(index) || error("taxon $taxon is not found in the network. Check for typos.")
        push!(taxonnums, net.leaf[index].number) # note: not ordered as in taxonnames
    end
    # get interior ancestor edges in major tree (no hybrids)
    matrix = hardwiredClusters(majorTree(net), vcat(taxonnames, outsideclade))
    # look for row with ones in relevant columns, zeros everywhere else (or zeros there and ones everywhere else)
    edgenum = 0 # 0 until we find the stem edge
    comparator = zeros(Int8, size(matrix)[2]-2)
    comparator[1:ntax_clade] = ones(Int8, ntax_clade)
    # go through matrix row by row to find match with comparator. Return number
    for i in 1:size(matrix, 1)
        if matrix[i,2:size(matrix)[2]-1] == comparator # found the mrca!
            edgenum = matrix[i, 1]
            break
        end
    end
    if edgenum == 0
        comparator .= 1 .- comparator
        for i in 1:size(matrix)[1]
            if matrix[i,2:size(matrix)[2]-1] == comparator
                error("""The clade given is not rooted correctly, making it a grade instead of a clade.
                You can re-try after modifying the root of your network.""")
            end
        end
        error("The taxa given do not form a clade in the network")
    end
    edgei = findfirst(e -> e.number == edgenum, net.edge)
    edgei !== nothing || error("hmm. hardwiredClusters on the major tree got an edge number not in the network")
    stemedge = net.edge[edgei]
    mrcanode = getChild(stemedge)
    TopologyConstraint(type, taxonnames, taxonnums, stemedge, mrcanode)
end
function Base.show(io::IO, obj::TopologyConstraint)
    str = (obj.type == 0x01 ? "Species constraint, on tips: " : "Clade constraint, on tips: " )
    str *= join(obj.taxonnames, ", ")
    str *= "\n stem edge number $(obj.edge.number)\n crown node number $(obj.node.number)"
    print(io, str)
end

"""
    constraintviolated(network::HybridNetwork,
                       constraints::Vector{TopologyConstraint})

True if `network` violates one (or more) of the constraints of type 1
(individuals in a species group) or type 2 (must be clades in the major tree).
Warning: constraints of type 3 are not implemented.
"""
function constraintviolated(net::HybridNetwork, constraints::Vector{TopologyConstraint})
    # fixit next PR: add option to give a vector of constraint types to check,
    #                then only check these constraint types
    if isempty(constraints)
        return false
    end # avoids extracting major tree when no constraint
    tree = majorTree(net)
    for con in constraints # checks directionality of stem edge hasn't changed
        if con.type in [0x01, 0x02]
            getChild(con.edge) === con.node || return true
            tei = findfirst(e -> e.number == con.edge.number, tree.edge)
            tei !== nothing ||
                error("hmm. edge number $(con.edge.number) was not found in the network's major tree")
            treeedge = tree.edge[tei]
            des = descendants(treeedge) # vector of node numbers, descendant tips only by default
            Set(des) == con.taxonnums || return true
        end
    end
    return false
end

"""
    updateconstraints!(constraints::Vector{TopologyConstraint}, net::HybridNetwork)

Update the set `taxonnum` in each constraint, assuming that the stem edge and
the crown node are still correct, and that their descendants are still correct.
May be needed if the node and edge numbers were modified by
[`resetNodeNumbers!`](@ref) or [`resetEdgeNumbers!`](@ref).

Warning: does *not* check that the names of leaves with numbers in `taxonnum`
are `taxonnames`.

```jldoctest
julia> net = readTopology("(((2a,2b),(((((1a,1b,1c),4),(5)#H1),(#H1,(6,7))))#H2),(#H2,10));");

julia> c_species1 = PhyloNetworks.TopologyConstraint(0x01, ["1a","1b","1c"], net)
Species constraint, on tips: 1a, 1b, 1c
 stem edge number 7
 crown node number -9

julia> c_species1.taxonnums
Set{Int64} with 3 elements:
  5
  4
  3

julia> c_clade145 = PhyloNetworks.TopologyConstraint(0x02, ["1a","1b","1c","4","5"], net)
Clade constraint, on tips: 1a, 1b, 1c, 4, 5
 stem edge number 12
 crown node number -7

julia> PhyloNetworks.resetNodeNumbers!(net)

julia> net.node[4].number = 111;

julia> PhyloNetworks.updateconstraints!([c_species1, c_clade145], net)

julia> c_species1
Species constraint, on tips: 1a, 1b, 1c
 stem edge number 7
 crown node number 21

julia> c_species1.taxonnums
Set{Int64} with 3 elements:
  5
  4
  111
```
"""
function updateconstraints!(constraints::Vector{TopologyConstraint}, net::HybridNetwork)
    if isempty(constraints)
        return nothing
    end # avoids extracting major tree when no constraint
    tree = majorTree(net)
    for con in constraints
      if con.type in [0x01, 0x02]
        getChild(con.edge) === con.node ||
            error("the stem edge and crown node have been disconnected")
        tei = findfirst(e -> e.number == con.edge.number, tree.edge)
        tei !== nothing ||
            error("stem edge number $(con.edge.number) was not found in the network's major tree")
        des = descendants(tree.edge[tei]) # vector of node numbers, descendant tips only by default
        d1 = setdiff(con.taxonnums, des)  # set
        d2 = setdiff(des, con.taxonnums)  # array
        length(d1) == length(d2) || error("missing or extra taxa in the clade / species.")
        replace!(con.taxonnums, Dict(tn => d2[i] for (i,tn) in enumerate(d1))...)
      end
    end
    return nothing
end

"""
    updateconstraintfields!(constraints::Vector{TopologyConstraint}, net::HybridNetwork)

Update fields stem `edge` and crown `node` to match the given `net`.

Assumes that the constraints are still met in `net`,
and that nodes & edges are numbered identically in `net` as in the network
used to create all `constraints`.

fixit: remove the assumption that constraints are still met, since
an NNI near the crown of a constrained clade might change the crown node
and / or the stem edge (u and v exchange).
"""
function updateconstraintfields!(constraints::Vector{TopologyConstraint}, net::HybridNetwork)
    for con in constraints
        num = con.edge.number
        con.edge = net.edge[findfirst([e.number == num for e in net.edge])]
        num = con.node.number
        con.node = net.node[findfirst([n.number == num for n in net.node])]
    end
end


#=
nice to add: option to modify branch lengths during NNI
=#

"""
    nni!(net::HybridNetwork, e::Edge, nohybridladder::Bool=true, no3cycle::Bool=true,
         constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Attempt to perform a nearest neighbor interchange (NNI) around edge `e`,
randomly chosen among all possible NNIs (e.g 3, sometimes more depending on `e`)
satisfying the constraints, and such that the new network is a DAG. The number of
possible NNI moves around an edge depends on whether the edge's parent/child nodes
are tree or hybrid nodes. This is calculated by [`nnimax`](@ref).

The option `no3cycle` forbids moves that would create a 3-cycle in the network.
When `no3cycle` = false, 2-cycle and 3-cycles may be generated.

Note that the defaults values are for positional (not keyword) arguments, so
two or more arguments can be used, but in a specific order: nni!(net, e) or
nni!(net, e, nohybridladder),
nni!(net, e, nohybridladder, no3cycle),
nni!(net, e, nohybridladder, no3cycle, contraints).

Assumptions:
- The starting network does not have 3-cycles, if `no3cycle=true`.
  No check for the presence of 2- and 3-cycles in the input network.
- The edges' field `isChild1` is correct in the input network. (This field
  will be correct in the output network.)

Output:
information indicating how to undo the move or `nothing` if all NNIs failed.

# examples

```jldoctest
julia> str_network = "(((S8,S9),(((((S1,S2,S3),S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));";

julia> net = readTopology(str_network);

julia> using Random; Random.seed!(321);

julia> undoinfo = nni!(net, net.edge[3], true, true); # true's to avoid hybrid ladders and 3-cycles

julia> writeTopology(net)
"(((S8,(((((S1,S2,S3),S4),(S5)#H1),(#H1,(S6,S7))))#H2),S9),(#H2,S10));"

julia> nni!(undoinfo...);

julia> writeTopology(net) == str_network # net back to original topology: the NNI was "undone"
true
```
"""
function nni!(net::HybridNetwork, e::Edge, nohybridladder::Bool=true, no3cycle::Bool=true,
              constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    for con in constraints
        # reject NNI if the focus edge is the stem of a constraint.
        # sufficient to maintain star species constraints
        # but not sufficient for clades.
        if con.edge === e
            return nothing
        end
    end
    nnirange = 0x01:nnimax(e) # 0x01 = 1 but UInt8 instead of Int
    nnis = Random.shuffle(nnirange)
    for nummove in nnis # iterate through all possible NNIs, but in random order
        moveinfo = nni!(net, e, nummove, nohybridladder, no3cycle)
        !isnothing(moveinfo) || continue # to next possible NNI
        # if constraintviolated(net, constraints)
        #     # fixit for next PR: not needed for species constraints,
        #     # not sufficient for clades if γ / isMajor modified later, or if
        #     # constraints met after NNI but stem / crows need to be moved
        #     nni!(moveinfo...) # undo the previous NNI
        #     continue          # try again
        # end
        return moveinfo
    end
    return nothing # if we get to this point, all NNIs failed
end

"""
    nnimax(e::Edge)

Return the number of NNI moves around edge `e`,
assuming that the network is semi-directed.
Return 0 if `e` is not internal, and more generally if either node
attached to `e` does not have 3 edges.

Output: UInt8 (e.g. `0x02`)
"""
function nnimax(e::Edge)
    length(e.node[1].edge) == 3 || return 0x00
    length(e.node[2].edge) == 3 || return 0x00
    # e.hybrid and hybrid parent: RR case, 2 possible NNIs only
    # e.hybrid and tree parent:   BR case, 3 or 6 NNIs if e is may contain the root
    # e not hybrid, hyb parent:   RB case, 4 NNIs
    # e not hybrid, tree parent:  BB case, 2 NNIs if directed, 8 if undirected
    hybparent = getParent(e).hybrid
    n = (e.hybrid ? (hybparent ? 0x02 : (e.containRoot ? 0x06 : 0x03)) : # RR & BR
                    (hybparent ? 0x04 : (e.containRoot ? 0x08 : 0x02)))  # RB & BB
    return n
end

"""
    nni!(net::HybridNetwork, uv::Edge, nummove::UInt8,
         nohybridladder::Bool, no3cycle::Bool)

Modify `net` with a semi-directed nearest neighbor interchange (sNNI) around
edge `uv`. Return the information necessary to undo the sNNI, or `nothing` if
the move was not successful (such as if the resulting graph was not acyclic
(not a DAG) or if the focus edge is adjacent to a polytomy). If the move fails,
the network is not modified.

`nummove` specifies which of the available sNNIs is performed.

Expanded from the rooted-NNI (rNNI) options in Gambette et al. (2017), fig. 8:
* BB: 2 moves, both to BB, if directed edges. 8 moves if undirected.
* RR: 2 moves, both to RR.
* BR: 3 moves, 1 RB & 2 BRs, if directed. 6 moves if e is undirected.
* RB: 4 moves, all 4 BRs.
The extra sNNI options come from assuming a semi-directed network, in
constrast to Gambette et al (2017)'s rooted networks. On a semi-directed network,
there may be a choice in how to direct the edges that may contain the root:
e.g. choice of e=uv versus vu and choice of labelling adjacent nodes as
α/β (BB) or as α/γ (BR).

`nohybridladder` = true prevents moves that would create a hybrid ladder in
the network, that is, 2 consecutive hybrid nodes (one parent of the other).
`no3cycle` = true prevents NNI moves that would make a 3-cycle, and assumes
that the input network does not have any 2- or 3-cycles.
If `no3cycle` is false, 3-cycles can be generated, but NNIs generating
2-cycles are prevented.

The edge field `isChild1` is assumed to be correct according to the `net.root`.
"""
function nni!(net::HybridNetwork, uv::Edge, nummove::UInt8,
              nohybridladder::Bool, no3cycle::Bool)
    nmovemax = nnimax(uv)
    if nmovemax == 0x00 # uv not internal, or polytomy e.g. species represented by polytomy with >2 indiv
        return nothing
    end
    nummove <= nmovemax || error("nummove $(nummove) must be <= $nmovemax for edge number $(uv.number)")
    u = getParent(uv)
    v = getChild(uv)
    ## TASK 1: grab the edges adjacent to uv: αu, βu, vγ, vδ and adjacent nodes
    # get edges αu & βu connected to u
    # α = u's major parent if it exists, β = u's last child or minor parent
    if u.hybrid
        αu = getMajorParentEdge(u)
        βu = getMinorParentEdge(u)
        α = getParent(αu)
        β = getParent(βu)
    else # u may not have any parent, e.g. if root node
        # pick αu = parent edge if possible, first edge of u (other than uv) otherwise
        labs = [edgerelation(e, u, uv) for e in u.edge]
        pti = findfirst(isequal(:parent), labs) # parent index: there should be 1 or 0
        ci = findall(  isequal(:child), labs)  # vector. there should be 1 or 2
        if pti === nothing
            length(ci) == 2 || error("node $(u.number) should have 2 children other than node number $(v.number)")
            pti = popfirst!(ci)
            αu = u.edge[pti]
            α = getChild(αu)
        else
            αu = u.edge[pti]
            α = getParent(αu)
        end
        βu = u.edge[ci[1]]
        β = getChild(βu)
    end
    # get edges vδ & vγ connected to v
    # δ = v's last child, γ = other child or v's parent other than u
    labs = [edgerelation(e, v, uv) for e in v.edge]
    if v.hybrid # then v must have another parent edge, other than uv
        vγ = v.edge[findfirst(isequal(:parent), labs)] # γ = getParent(vδ) ?this should be vγ right?
            #on net_hybridladder edge 1 get ERROR: ArgumentError: invalid index: nothing of type Nothing
        vδ = v.edge[findfirst(isequal(:child), labs)]
        γ = getParent(vγ)
        δ = getChild(vδ)
    else
        ci = findall(isequal(:child), labs)
        length(ci) == 2 || error("node $(v.number) should have 2 children")
        vγ = v.edge[ci[1]] # γ = getChild(vγ)
        vδ = v.edge[ci[2]]
        γ = getChild(vγ)
        δ = getChild(vδ)
    end
    ## TASK 2: semi-directed network swaps:
    ## swap u <-> v, α <-> β if undirected and according to move number
    if !u.hybrid && uv.containRoot # case BB or BR, with uv that may contain the root
        nmoves = 0x03 # number of moves in the rooted (directed) case
        if !v.hybrid # BB, uv and all adjacent edges are undirected: 8 moves
            if nummove > 0x04  # switch u & v with prob 1/2
                nummove -= 0x04  # now in 1,2,3,4
                (u,v) = (v,u)
                (α,αu,β,βu,γ,vγ,δ,vδ) = (γ,vγ,δ,vδ,α,αu,β,βu)
            end
            nmoves = 0x02
        end
        # now nmoves should be in 1, ..., 2*nmoves: 1:4 (BB) or 1:6 (BR) #? should this be 2*nummoves?
        # next: switch α & β with probability 1/2
        if nummove > nmoves
            nummove -= nmoves # now nummoves in 1,2 (BB) or 1,2,3 (BR)
            (α,αu,β,βu) = (β,βu,α,αu)
        end
    end
    ## TASK 3: rooted network swaps, then do the NNI
    if u.hybrid # RR or RB cases
        # swap α & β: reduce NNIs from 2 to 1 (RR) or from 4 to 2 (RB)
        if nummove == 0x02 || nummove == 0x04
            (α,αu,β,βu) = (β,βu,α,αu)
        end
        if !v.hybrid && nummove > 0x02 # case RB, moves 3 or 4
            (γ,vγ,vδ,δ) = (δ,vδ,vγ,γ)
        end
        # detach β and graft onto vδ
        if no3cycle # don't create a 3-cycle
            if problem4cycle(α,γ, β,δ) return nothing; end
        else # don't create a 2-cycle
            if α === γ || β === δ return nothing; end
        end
        res = nni!(αu, u, uv, v, vδ)
    else # u = bifurcation: BB or BR cases
        if !v.hybrid # case BB: 2 options
            if nummove == 0x02
                (γ,vγ,vδ,δ) = (δ,vδ,vγ,γ)
            end
            # detach β and graft onto vδ
            if no3cycle
                if problem4cycle(α,γ, β,δ) return nothing; end
            elseif α === γ || β === δ return nothing
            end
            res = nni!(αu, u, uv, v, vδ)
        else # case BR: 3 options
            # DAG check
            # if α parent of u, check if β -> γ (directed or undirected)
            # if undirected with β parent of u, check α -> γ
            # if undirected and u = root:
            #   moves 1 and 3 will fail if α -> γ
            #   moves 2 and 3 will fail if β -> γ
            # nummove 3 always creates a nonDAG when u is the root
            αparentu = getChild(αu)===u
            βparentu = getChild(βu)===u
            if αparentu
                if γ === β || isdescendant(γ, β) return nothing; end
            elseif βparentu
                if γ === α || isdescendant(γ, α) return nothing; end
            else # cases when u is root
                if nummove == 0x01 && (γ === α || isdescendant(γ, α))
                    return nothing
                elseif nummove == 0x02 && (γ === β || isdescendant(γ, β))
                    return nothing
                elseif nummove == 0x03 # fail: not DAG
                    return nothing
                end
            end
            if nummove == 0x01 # graft γ onto α
                if no3cycle
                    if problem4cycle(α,γ, β,δ) return nothing; end
                elseif α === γ || β === δ return nothing
                end
                if nohybridladder && α.hybrid return nothing; end
                res = nni!(vδ, v, uv, u, αu)
            elseif nummove == 0x02 # graft γ onto β
                if nohybridladder && β.hybrid return nothing; end
                if no3cycle
                    if problem4cycle(β,γ, α,δ) return nothing; end
                elseif β === γ || α === δ return nothing
                end
                res = nni!(vδ, v, uv, u, βu)
            else # nummove == 0x03
                if αparentu # if α->u: graft δ onto α
                    if no3cycle
                        if problem4cycle(α,δ, β,γ) return nothing; end
                    elseif β === γ || α === δ return nothing
                    end
                    if nohybridladder && α.hybrid return nothing; end
                    res = nni!(vγ, v, uv, u, αu)
                else # if β->u: graft δ onto β
                    if no3cycle
                        if problem4cycle(α,γ, β,δ) return nothing; end
                    elseif α === γ || β === δ return nothing
                    end
                    # we could check nohybridladder && β.hybrid but no need because
                    # we only get here if uv.containRoot so β cannot be hybrid
                    res = nni!(vγ, v, uv, u, βu)
                # case when u is root already excluded (not DAG)
                end
            end
        end
    end
    return res
end

"""
    nni!(αu, u, uv::Edge, v, vδ)
    nni!(αu,u,uv,v,vδ, flip::Bool, inner::Bool, indices)

Transform a network locally around the focus edge `uv` with
the following NNI, that detaches u-β and grafts it onto vδ:

    α - u -- v ------ δ
        |    |
        β    γ

    α ------ v -- u - δ
             |    |
             γ    β

`flip` boolean indicates if the uv edge was flipped
`inner` boolean indicates if edges αu and uv both point toward node u,
i.e. α->u<-v<-δ. If this is true, we flip the hybrid status of αu and vδ.

`indices` give indices for nodes and edges u_in_αu, αu_in_u, vδ_in_v, and v_in_vδ.
These are interpreted as:

    u_in_αu: the index for u in the edge αu
    αu_in_u: the index for αu in node u
    vδ_in_v: the index for vδ in node v
    v_in_vδ: the index for v in edge vδ

**Warnings**:

- *No* check of assumed adjacencies
- Not implemented for cases that are not necessary thanks to symmetry,
  such as cases covered by `nni!(vδ, v, uv, u, αu)` or `nni!(βu, u, v, vγ)`.
  More specifically, these cases are not implemented (and not checked):
  * u not hybrid & v hybrid
  * u hybrid, v not hybrid, α -> u <- v -> δ
- Because of this, `nni(αu,u,uv,v,vδ, ...)` should not be used directly;
  use instead `nni!(net, uv, move_number)`.
- nni!(undoinfo...) restores the topology, but edges below hybrid nodes
  will have length 0.0 even if they didn't before.

Node numbers and edge numbers are not modified.
Edge `uv` keeps its direction unchanged *unless* the directions were
`α -> u -> v -> δ` or `α <- u <- v <- δ`,
in which case the direction of `uv` is flipped.

The second version's input has the same signature as the output, but
will undo the NNI more easily. This means that if `output = nni!(input)`,
then `nni!(output...)` is valid and undoes the first operation.

Right now, branch lengths are not modified except when below a hybrid node.
Future versions might implement options to modify branch lengths.
"""
function nni!(αu::Edge, u::Node, uv::Edge, v::Node, vδ::Edge)
    # find indices to detach αu from u, and to detach v from vδ
    u_in_αu = findfirst(n->n===u, αu.node)
    αu_in_u = findfirst(e->e===αu, u.edge)
    vδ_in_v = findfirst(e->e===vδ, v.edge)
    v_in_vδ = findfirst(n->n===v, vδ.node)
    # none of them should be 'nothing' --not checked
    αu_child = getChild(αu)
    uv_child = getChild(uv)
    vδ_child = getChild(vδ)
    # flip the direction of uv if α->u->v->δ or α<-u<-v<-δ
    flip = (αu_child === u && uv_child === v && vδ_child !== v) ||
             (αu_child !== u && uv_child === u && vδ_child === v)
    inner = (flip ? false : αu_child === u && vδ_child === v )
    # fixit: calculate new edge lengths here; give them as extra arguments to nni! below
    res = nni!(αu,u,uv,v,vδ, flip,inner, u_in_αu,αu_in_u,vδ_in_v,v_in_vδ)
    return res
end
function nni!(αu::Edge, u::Node, uv::Edge, v::Node, vδ::Edge,
              flip::Bool, inner::Bool,
              u_in_αu::Int, αu_in_u::Int, vδ_in_v::Int, v_in_vδ::Int)
              # extra arguments: new edge lengths
    # fixit: extract old edge lengths; update to new lengths; return old lengths as extra items
    ## TASK 1: do the NNI swap: to detach & re-attach at once
    (αu.node[u_in_αu], u.edge[αu_in_u], v.edge[vδ_in_v], vδ.node[v_in_vδ]) =
        (v, vδ, αu, u)
    ## TASK 2: update edges' directions
    if flip
        uv.isChild1 = !uv.isChild1
    end
    ## TASK 3: update hybrid status of nodes (#? is it ever nodes?)
             # & edges, and edges containRoot field
    if u.hybrid && !v.hybrid
        if flip # RB, BR 1 or BR 2': flip hybrid status of uv and (αu or vδ)
            if uv.hybrid # BR 1 or BR 2'
                vδ.hybrid = true
                vδ.gamma = uv.gamma
                vδ.isMajor = uv.isMajor
                # length switches: for keeping the network unzipped. could be option later
                αu.length = uv.length
                uv.length = 0.0
                uv.hybrid = false
                uv.gamma = 1.0
                uv.isMajor = true
                norootbelow!(uv)
            else # RB
                uv.hybrid = true
                uv.gamma = αu.gamma
                uv.isMajor = αu.isMajor
                # length switches: for keeping the network unzipped. could be option later
                uv.length = vδ.length
                vδ.length = 0.0
                αu.hybrid = false
                αu.gamma = 1.0
                αu.isMajor = true
                if αu.containRoot
                    allowrootbelow!(v, αu)
                end
            end
        else # assumes α<-u or v<-δ --- but not checked
            # not implemented: α->u<-v->δ: do uv.hybrid=false and γv.hybrid=true
            if inner # BR 3 or 4': α->u<-v<-δ: flip hybrid status of αu and vδ
                # no lengths to switch in these cases
                vδ.hybrid = true
                vδ.gamma = αu.gamma
                vδ.isMajor = αu.isMajor
                αu.hybrid = false
                αu.gamma = 1.0
                αu.isMajor = true
                # containRoot: nothing to update
            else # BR 1' or 2: α<-u<-v->δ : hybrid status just fine
                # length switches: for keeping the network unzipped. could be option later
                αu.length = vδ.length
                vδ.length = 0.0
                norootbelow!(vδ)
                if uv.containRoot
                    allowrootbelow!(αu)
                end
            end
        end
    elseif u.hybrid && v.hybrid # RR: hybrid edges remain hybrids, but switch γs
        # we could have -αu-> -uv-> -vδ-> or <-αu- <-uv- <-vδ-
        hyb = ( getChild(αu) == v ? αu : vδ ) # uv's hybrid parent edge: either αu or vδ
        (uv.gamma,   hyb.gamma)   = (hyb.gamma,   uv.gamma)
        (uv.isMajor, hyb.isMajor) = (hyb.isMajor, uv.isMajor)
    end
    # throw error if u not hybrid & v hybrid? (does not need to be implemented)
    return vδ, u, uv, v, αu, flip, inner, v_in_vδ, αu_in_u, vδ_in_v, u_in_αu
end

"""
    problem4cycle(β::Node, δ::Node, α::Node, γ::Node)

Check if the focus edge uv has a 4 cycle that could lead to a 3 cycle
after an chosen NNI. Return true if there is a problem 4 cycle, false if none.
"""
function problem4cycle(β::Node, δ::Node, α::Node, γ::Node)
    isconnected(β, δ) || isconnected(α, γ)
end

"""
    checkspeciesnetwork!(network::HybridNetwork,
                         constraints::Vector{TopologyConstraint})

Check that the network satisfies a number of requirements:
- no polytomies, other than at species constraints: throws an error otherwise
- no unnecessary nodes: fuse edges at nodes with degree two, including at the root
  (hence the bang: `net` may be modified)
- topology constraints are met.

Output: true if all is good, false if one or more clades are violated.
"""
function checkspeciesnetwork!(net::HybridNetwork, constraints::Vector{TopologyConstraint})
    for n in net.node # check for no polytomies: all nodes should have up to 3 edges
        if length(n.edge) > 3 # polytomies allowed at species constraints only
            coni = findfirst(c -> c.type == 1 && c.node === n, constraints)
            coni !== nothing ||
                error("The network has a polytomy at node number $(n.number). Please resolve.")
        end
    end
    removedegree2nodes!(net)
    return !constraintviolated(net, constraints)
end

"""
    mapindividuals(net::HybridNetwork, mappingFile::String)

Return a network expanded from `net`, where species listed in the mapping file
are replaced by individuals mapped to that species.
If a species has only 1 individual, the name of the leaf for that species is
replaced by the name of its one individual representative.
If a species has 2 or more individuals, the leaf for that species is
expanded into a "star" (polytomy if 3 or more individuals) with a tip for each
individual. If a species is in the network but not listed in the mapping file,
the tip for that species is left as is. Species listed in the mapping file
but not present in the network are ignored.

The mapping file should be readable by `CSV.File` and contain two columns:
one for the species names and one for the individual (or allele) names.
fixit: make this function more flexible by accepting column names

Output: individual-level network and vector of species constraint(s).

# examples

```jldoctest
julia> species_net = readTopology("(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));");

julia> filename = joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "mappingIndividuals.csv");

julia> filename |> read |> String |> print # to see what the mapping file contains
species,individual
S1,S1A
S1,S1B
S1,S1C

julia> individual_net, species_constraints = mapindividuals(species_net, filename);

julia> writeTopology(individual_net, internallabel=true)
"(((S8,S9),(((((S1A,S1B,S1C)S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"

julia> species_constraints
1-element Vector{PhyloNetworks.TopologyConstraint}:
 Species constraint, on tips: S1A, S1B, S1C
 stem edge number 4
 crown node number 3
```
"""
function mapindividuals(net::HybridNetwork, mappingFile::String)
    mappingDF = DataFrame(CSV.File(mappingFile); copycols=false)
    specieslist = unique(mappingDF[:, 1])
    individualnet = deepcopy(net)
    constraints = TopologyConstraint[]
    for species in specieslist
        individuals = mappingDF[mappingDF.species .== species, 2]
        res = addindividuals!(individualnet, species, individuals)
        res === nothing || # nothing if species not found in the network, or 0,1 individuals
            push!(constraints, res)
    end
    return individualnet, constraints
end

"""
    addindividuals!(net::HybridNetwork, species::AbstractString, individuals::Vector)

Add individuals to their species leaf as a star, and return the
corresponding species constraint. `nothing` is returned in 2 cases:

- if `individuals` contains only 1 individual, in which case the
  name of the leaf for that species is replaced by the name of
  its one individual representative
- if `species` is not found in the network

Spaces in individuals' names are eliminated, see [`cleantaxonname`](@ref).

Called by [`mapindividuals`](@ref).

# examples

```jldoctest
julia> net = readTopology("(S8,(((S1,S4),(S5)#H1),(#H1,S6)));");

julia> PhyloNetworks.addindividuals!(net, "S1", ["S1A", "S1B", "S1C"])
Species constraint, on tips: S1A, S1B, S1C
 stem edge number 2
 crown node number 2

julia> writeTopology(net) # 3 new nodes, S1 now internal: not a tip
"(S8,((((S1A,S1B,S1C)S1,S4),(S5)#H1),(#H1,S6)));"
```
"""
function addindividuals!(net::HybridNetwork, species::AbstractString, individuals::Vector{String})
    length(individuals) > 0 || return nothing
    species = strip(species)
    speciesnodeindex = findfirst(l -> l.name == species, net.node)
    if speciesnodeindex === nothing
        @warn "species $species not found in network (typo in mapping file?). will be ignored."
        return nothing
    end
    if length(individuals) == 1 # change species name into name of individual representative
        net.node[speciesnodeindex].name = cleantaxonname(individuals[1])
        return nothing
    end
    # there are 2+ individuals: replace the species tip by a star, will need a species constraint
    taxnames = map(cleantaxonname, individuals)
    for tax in taxnames
        addleaf!(net, net.node[speciesnodeindex], tax);
    end
    return TopologyConstraint(0x01, taxnames, net)
end

"""
    cleantaxonname(taxonname::AbstractString)

Return a String with leading and trailing spaces removed, and interior spaces
replaced by underscores: good for using as tip names in a network without
causing future error when reading the newick description of the network.
"""
function cleantaxonname(taxonname::AbstractString)
    tname = string(strip(taxonname)) # SubString if we don't do string()
    m = match(r"\s", tname)
    if m !== nothing
        @warn """Spaces in "$tname" may cause errors in future network readability: replaced by _"""
        tname = replace(tname, r"\s+" => "_")
    end
    return tname
end


"""
    moveroot!(net::HybridNetwork, constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Move the root to a randomly chosen non-leaf node that is different from
the current root, and not within a constraint clade or species.
Output: `true` if successul, `nothing` otherwise.
"""
function moveroot!(net::HybridNetwork, constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    newrootrandomorder = Random.shuffle(1:length(net.node))
    oldroot = net.root
    for newrooti in newrootrandomorder
        newrooti != oldroot || continue
        newrootnode = net.node[newrooti]
        newrootnode.leaf && continue # try next potential new root if current is a leaf
        # Check the new root is NOT at the top or inside contraint clade or within pecies
        newrootfound = true
        for con in constraints
            if con.node === newrootnode || isdescendant(newrootnode, con.node) # assumes the constraint is met
                newrootfound = false
                break # out of constraint loop
            end
        end
        newrootfound || continue # to next potential new root
        # Check for no directional conflict
        try
            net.root = newrooti
            directEdges!(net)
            # warning: assumes that the previous root has degree 3
            # otherwise, need to delete previous root by fusing its 2 adjacent edges
        catch e # RootMismatch error if root placement is below a hybrid
            isa(e, RootMismatch) || rethrow(e)
            net.root = oldroot
            directEdges!(net) # revert edges' directions to match original rooting
            continue # to next potential new root
        end
        # if we get here: newrooti passed all the checks
        return true
    end
    # if we get here: none of the root positions worked
    return nothing
end

"""
    fliphybrid!(net::HybridNetwork, minor=true::Bool, nohybridladder=false::Bool,
                constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Cycle through hybrid nodes in random order until an admissible flip is found.
At this hybrid node, flip the indicated hybrid parent edge (minor or major).

If an admissible flip is found, return the tuple: newhybridnode, flippededge, oldchildedge.
Otherwise, return nothing.

The flip can be undone with
`fliphybrid!(net, newhybridnode, minor, constraints)`, or
`fliphybrid!(net, newhybridnode, !flippededge.isMajor, constraints)`
more generally, such as if the flipped edge had its γ modified after the
original flip.

*Warnings*
- if the root needed to be reset and if the original root was of degree 2,
  then a node of degree 2 remains in the modified network.
- undoing the flip may not recover the original root in case
  the root position was modified during the original flip.
"""
function fliphybrid!(net::HybridNetwork, minor=true::Bool,
                     nohybridladder=false::Bool,
                     constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    hybridindex = Random.shuffle(1:length(net.hybrid)) # indices in net.hybrid
    while !isempty(hybridindex) # all minor edges
        i = pop!(hybridindex)
        undoinfo = fliphybrid!(net, net.hybrid[i], minor, nohybridladder, constraints)
        if !isnothing(undoinfo) # if it failed, the network was not changed
            return undoinfo
        end # else, continue until explored all i in hybridindex
    end
    return nothing
end

"""
    fliphybrid!(net::HybridNetwork, hybridnode::Node, minor=true::Bool,
                nohybridladder=false::Bool,
                constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Flip the direction of a single hybrid edge:
the minor parent edge of `hybridnode` by default,
or the major parent edge if `minor` is false.
The parent node of the hybrid edge becomes the new hybrid node.
The former hybrid edge partner is converted to a tree edge (with γ=1),
and `hybridnode` becomes a tree node.

For the flip to be admissible, the new network must be a semi-directed
phylogenetic network: with a root such that the rooted version is a DAG.
If `nohybridladder` is false (default), the flip may create a hybrid ladder
If `nohybridladder` is true and if the flip would create a hybrid ladder,
then the flip is not admissible.
A hybrid ladder is when a hybrid child of another hybrid.

The new hybrid partner is an edge adjacent to the new hybrid node,
such that the flip is admissible (so it must be a tree edge).
The flipped edge retains its original γ.
The new hybrid edge is assigned inheritance 1-γ.

Output:
`(newhybridnode, flippededge, oldchildedge)` if the flip is admissible,
`nothing` otherwise.

The network is unchanged if the flip is not admissible.
If the flip is admissible, the root position may be modified, and
the direction of tree edges (via `isChild1`) is modified accordingly. If the
root needs to be modified, then the new root is set to the old hybrid node.

The index of the new hybrid node in `net.hybrid` is equal to that of the
old `hybridnode`.

Warning: Undoing this move may not recover the original root if
the root position was modified.
"""
function fliphybrid!(net::HybridNetwork, hybridnode::Node, minor=true::Bool,
                     nohybridladder=false::Bool,
                     constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    #= for species constraints, there is nothing to check, because hybrids cannot point into or come out of the group
    for con in constraints
        if con.type == # 0x02/0x03 for types 2 and 3, need to consider more cases
            return nothing
        end
    end
    =#
    runDirectEdges = false
    edgetoflip, edgetokeep = minor ?
        (getMinorParentEdge(hybridnode), getMajorParentEdge(hybridnode)) :
        (getMajorParentEdge(hybridnode), getMinorParentEdge(hybridnode))
    oldchildedge = getChildEdge(hybridnode)
    newhybridnode = getParent(edgetoflip)
    if newhybridnode.hybrid # already has 2 parents: cannot had a third.
        return nothing
    end
    ## choose newhybridedge ##
    p2 = getParent(edgetokeep) # parent node
    isdesc = Bool[]
    for e in newhybridnode.edge # is p2 undirected descendant of nhn via this edge?
        isp2desc = false
        #= isp2desc = true means:
        there is a semi-directed path  newhybridnode -e-> neighbor --...-> p2
        so: flipping hybridedge without making e the new partner would:
        - make e a child of newhybridnode, and
        - create a directed cycle: p2 -hybridedge-> newhybridnode -e-> neibr --...-> p2
        so we HAVE to make e the new partner hybrid edge if its "isp2desc" is true
        =#
        if e !== edgetoflip # e cannot be hybrid --e-> nhn because earlier check
            neibr = getOtherNode(e, newhybridnode) # neighbor of new hybrid node via e
            if !neibr.leaf && (p2 === neibr ||
                isdescendant_undirected(p2, neibr, e))
                isp2desc = true
            end
        end
        push!(isdesc, isp2desc)
    end
    sum_isdesc = sum(isdesc)
    if sum_isdesc > 1 # if 2+ edges have a semi-directed path to p2,
        # flipping either would create a semi-directed cycle via the other(s)
        return nothing
    end
    if sum_isdesc == 0 # no risk to create directed cycle.
        # choose the current parent edge of nhn to: keep its current direction
        # and keep reaching all nodes from a single root
        newhybridedge = getMajorParentEdge(newhybridnode)
    else
        newhybridedge = newhybridnode.edge[findfirst(isdesc)]
        if newhybridedge.hybrid # cannot flip another hybrid edge, but we needed
            return nothing      # to calculate its isp2desc for total # of true's
        end
    end
    if nohybridladder
        # case when edgetoflip would be bottom rung of ladder
        if getOtherNode(newhybridedge, newhybridnode).hybrid
            return nothing
        end
        # case when edgetoflip would be the top rung: happens when
        # newhybridnode = center point of a W structure initially
        for e in newhybridnode.edge
            if e !== edgetoflip && e.hybrid # newhybridedge is NOT hybrid, so far
                return nothing
            end
        end
    end
    # check if newhybridedge's direction needs to be flipped. If so, reroot.
        # includes case when newhybridnode was the original root
    if getChild(newhybridedge) !== newhybridnode
        # choose a new root: after the flip, all edges connected to the root
            # must be able to containRoot
        if edgetoflip.containRoot && edgetokeep.containRoot
            newroot = hybridnode
        else # set parent node of newhybridedge as new root
            p3 = getOtherNode(newhybridedge, newhybridnode)
            all(e.containRoot for e in p3.edge) || return nothing
            newroot = p3
        end
        net.root = findfirst(n -> n === newroot, net.node)
        # then flip newhybridedge, make it point towards newhybridnode
        newhybridedge.isChild1 = !newhybridedge.isChild1
        runDirectEdges = true # in most cases, many edges will need to be flipped
    end
    # change hybrid status and major status for nodes and edges
    hybridnode.hybrid = false
    edgetokeep.hybrid = false
    newhybridnode.hybrid = true
    newhybridedge.hybrid = true
    newhybridedge.isMajor = minor
    edgetokeep.isMajor = true
    # update node order to keep isChild1 attribute of hybrid edges true
    edgetoflip.isChild1 = !edgetoflip.isChild1 # just switch
    # give new hybridnode a name
    newhybridnode.name = hybridnode.name
    hybridnode.name = "" # remove name from former hybrid node
    # update gammas
    newhybridedge.gamma = edgetokeep.gamma
    edgetokeep.gamma = 1.0
    # update hybrids in network (before directEdges!)
    hybridindex = findfirst(n -> n === hybridnode, net.hybrid)
    net.hybrid[hybridindex] = newhybridnode
    if runDirectEdges # direct edges only if root is moved
        directEdges!(net)
    else # if not, only update containRoot attributes
        # norootbelow in child edges of newhybridedge
        for ce in newhybridnode.edge
            ce !== newhybridedge || continue # skip e
            getParent(ce) === newhybridnode || continue # skip edges that aren't children of cn
            norootbelow!(ce)
        end
        # allow root for edges below old hybrid node
        if edgetokeep.containRoot #else, already forbidden below due to hybrid above
            allowrootbelow!(hybridnode, edgetokeep) # old hybrid node
        end
    end
    return newhybridnode, edgetoflip, oldchildedge
end

"""
    nnineighbors(net::HybridNetwork, nohybridladder::Bool, no3cycle::Bool,
                 constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Return the immediate sNNI neighborhood of `net` as an array
of the neighbor networks' Newick strings.

Warning:
- Because different NNIs can lead to the same topology, the list of
neighbors may contain duplicates. Therefore, the length of the list should not
be taken as the number of neighbors. For the unique neighbors, see
`uniqueneighbornets` function below.
- creates one deep copy of net
- clade constraints are not fully implemented yet --species constraints only.

"""
function nnineighbors(net::HybridNetwork, nohybridladder::Bool, no3cycle::Bool,
                      constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    neighbors = String[]
    neibr = deepcopy(net) # do NOT modify the input network
    for e in neibr.edge
        if any(con.edge === e for con in constraints) # if edge is a stem
            continue # skip to next edge
        end
        for nummove in 0x01:nnimax(e) # empty if nnimax(e) <= 0x00
            moveinfo = nni!(neibr, e, nummove, nohybridladder, no3cycle)
            !isnothing(moveinfo) || continue # move didn't work, skip to next NNI
            neibr.root == net.root || error("root changed by NNI")
            #= cecile:
            did you check that NNIs don't never move the root?

            why have keep the multiple occurences of the same neighbor?
                why not remove the duplicates here, and
                delete the function uniqueneighbornets? example:
            duplicate = findlast(n -> hardwiredClusterDistance(n,neibr,true)==0, neighbors)
            if isnothing(duplicate) then push!(neighbors, deepcopy(neibr)).
            else don't push the duplicate.
            findlast: because the last neighbors pushed were from NNIs along
            the same edge, so most likely to match the duplicate.
            =#
            push!(neighbors, writeTopology(neibr))
            #= cecile: why push the string version? Every downstream use
                       transforms it back to a HybridNetwork object.
            =#
            nni!(moveinfo...) # undo this nni
            # not done: confirm that neibr is unchanged after undoing nni
            # hardwiredClusterDistance(neibr, net, true) == 0 || error("original net changed")
        end
    end
    return neighbors
end

"""
    uniqueneighbornets(net::HybridNetwork, nohybridladder::Bool, no3cycle::Bool,
                 constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Return all networks obtained by applying to `net` one
nearest neighbor interchange (NNI). If the same topology can be obtained
by several NNIs (possibly with different sets of branch lengths),
it is listed once only.

Assumes that rooting does not change during an NNI.
"""
function uniqueneighbornets(net::HybridNetwork, nohybridladder::Bool,
                            no3cycle::Bool,
                            constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    neighbors = nnineighbors(net, nohybridladder, no3cycle, constraints)
    # cecile: why not do this below?
    # function isequal(n1, n2)
    #     hardwiredClusterDistance(n1,n2,true) == 0
    # end
    # return unique(neighbors)
    ncount = length(neighbors)
    hwcds = zeros(Int, ncount, ncount) # rooted distances between neighbors
    for i in 1:ncount
        for j in 1:ncount
            hwcds[i,j] = hardwiredClusterDistance(readTopology(neighbors[i]),
                                            readTopology(neighbors[j]), true)
        end
    end
    # find rows with duplicate HWCD patterns
    duplicates = nonunique(DataFrame(hwcds)) # true entries indicate duplicate rows
    duplicateindices = findall(x->x==1, duplicates)
    for di in duplicateindices # confirm duplicate net has hwcd == 0 with >= 1 other net
        any(i-> i == 0, hwcds[di, 1:end .!= di]) ||
            error("Duplicate net does not have a HWCD of zero with another net")
    end
    deleteat!(neighbors, duplicateindices)
    return neighbors
end
