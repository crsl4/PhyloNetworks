# functions to read/write networks topologies
# originally in functions.jl
# Claudia March 2015

# ------------------------------- Read topology ------------------------------------------

function peekchar(s::IOStream)
    Base.peekchar(s)
end

function peekchar(s::IOBuffer)
    mark(s)
    c = read(s,Char)
    reset(s)
    return c
end

# aux function to advance stream in readSubtree
# input: s IOStream/IOBuffer
function advance!(s::IO, c::Char, numLeft::Array{Int64,1})
    c = peekchar(s)
    if(Base.eof(s))
        error("Tree ends prematurely while reading subtree after left parenthesis $(numLeft[1]).")
    end
    return read(s,Char)
end


# aux function to read all digits of taxon name
# it allows names with letters and numbers
# it also reads # as part of the name and returns pound=true
# it returns the node name as string as well to check if it exists already (as hybrid)
# warning: treats digit taxon numbers as strings to avoid repeated node numbers
function readNum(s::IO, c::Char, net::HybridNetwork, numLeft::Array{Int64,1})
    pound = 0;
    if(isalnum(c) || isValidSymbol(c) || c == '#')
        pound += (c == '#') ? 1 : 0
        num = read(s,Char)
        c = peekchar(s)
        while(isalnum(c) || isValidSymbol(c) || c == '#')
            d = read(s,Char)
            num = string(num,d)
            if(d == '#')
                pound += 1;
               c = peekchar(s);
               if(isalnum(c))
                   if(c != 'L' && c != 'H' && c != 'R')
                       warn("Expected H, R or LGT after # but received $(c) in left parenthesis $(numLeft[1]).")
                   end
               else
                   a = readall(s);
                   error("Expected name after # but received $(c) in left parenthesis $(numLeft[1]). Remaining is $(a).")
               end
            end
            c = peekchar(s);
        end
        if(pound == 0)
            return size(net.names,1)+1, num, false
        elseif(pound == 1)
            return size(net.names,1)+1, num, true
        else
            a = readall(s);
            error("strange node name with $(pound) # signs. remaining is $(a).")
        end
    else
        a = readall(s);
        error("Expected int digit, alphanum or # but received $(c). remaining is $(a).");
    end
end


# aux function to read floats like length
function readFloat(s::IO, c::Char)
    founde = false
    if(isdigit(c))
        num = read(s,Char);
        c = peekchar(s);
        while(isdigit(c) || c == '.' || c == 'e' || c == '-')
            d = read(s,Char);
            num = string(num,d);
            c = peekchar(s);
        end
        try
            float(num)
        catch
            error("problem with number read $(num), not a float number")
        end
        return float(num)
    else
        a = readall(s);
        error("Expected float digit after : but received $(c). remaining is $(a).");
    end
end

# aux function to identify if a symbol in taxon name is valid
# symbol cannot be: () [] : ; ' , . space \t \r \n
# according to richnewick.pdf
function isValidSymbol(c::Char)
    return !isspace(c) && c != '(' && c != ')' && c != '[' && c != ']' && c != ':' && c != ';' && c != ',' #&& c != '.'
end

# aux function to read subtree
# input: s IOStream/IOBuffer
# warning: reads additional info :length:bootstrap:gamma
# warning: allows for name of internal nodes without # after: (1,2)A,...
# warning: warning if hybrid edge without gamma value, warning if gamma value (ignored) without hybrid edge
# modified from original Cecile c++ code to allow polytomies
function readSubtree!(s::IO, parent::Node, numLeft::Array{Int64,1}, net::HybridNetwork, hybrids::Array{ASCIIString,1}, index::Array{Int64,1})
    c = peekchar(s)
    e = nothing;
    hasname = false; # to know if the current node has name
    pound = false;
    if(c =='(')
       numLeft[1] += 1
       DEBUGC && println(numLeft)
       n = Node(-1*numLeft[1],false);
       cind = 1;
       keepon = true;
       c = read(s,Char)
       while (keepon)
           bl = readSubtree!(s,n,numLeft,net,hybrids,index)
           c = advance!(s,c,numLeft)
           if (c == ')')
               keepon = false
           elseif (c != ',')
               a = readall(s);
               error("Expected right parenthesis after left parenthesis $(numLeft[1]) but read $(c). The remainder of line is $(a).")
           end
       end
        c = peekchar(s);
        if(isalnum(c) || isValidSymbol(c) || c == '#') # internal node has name
            hasname = true;
            num,name,pound = readNum(s,c,net,numLeft);
            n.number = num;
            c = peekchar(s);
            if(!pound)
                warn("internal node with name without it being a hybrid node. node name might be meaningless after tree modifications.")
            end
        end
    elseif(isalnum(c) || isValidSymbol(c) || c == '#')
        hasname = true;
        bl = true;
        num,name,pound = readNum(s,c,net,numLeft)
        n = Node(num,true);
    else
        a = readall(s);
        error("Expected beginning of subtree but read $(c) after left parenthesis $(numLeft[1]), remaining is $(a).");
    end
    if(pound) # found pound sign in name
        n.hybrid = true;
        DEBUGC && println("encontro un hybrid $(name).")
        DEBUGC && println("hybrids list tiene size $(size(hybrids,1))")
        if(in(name,hybrids))
            DEBUGC && println("dice que $(name) esta en hybrids")
            ind = getIndex(name,hybrids);
            other = net.node[index[ind]];
            DEBUGC && println("other is leaf? $(other.leaf), n is leaf? $(n.leaf)")
            if(!n.leaf && !other.leaf)
                error("both hybrid nodes are internal nodes: successors of the hybrid node must only be included in the node list of a single occurrence of the hybrid node.")
            elseif(n.leaf)
                e = Edge(net.numEdges+1);
                e.hybrid = true
                e.isMajor = false;
                pushEdge!(net,e);
                setNode!(e,[other,parent]);
                setEdge!(other,e);
                setEdge!(parent,e);
            else # !n.leaf
                if(size(other.edge,1) == 1) #other should be a leaf
                    DEBUGC && println("other is $(other.number), n is $(n.number), edge of other is $(other.edge[1].number)")
                    otheredge = other.edge[1];
                    otherparent = getOtherNode(otheredge,other);
                    DEBUGC && println("parent of other is $(otherparent.number)")
                    removeNode!(other,otheredge);
                    deleteNode!(net,other);
                    setNode!(otheredge,n);
                    setEdge!(n,otheredge);
                    e = Edge(net.numEdges+1);
                    e.hybrid = true
                    setNode!(e,[n,parent]);
                    setEdge!(n,e);
                    setEdge!(parent,e);
                    pushNode!(net,n);
                    pushEdge!(net,e);
                    n.number = other.number;
                    n.name = other.name;
                else
                    error("strange: node $(other.number) is a leaf hybrid node so it should have only one edge and it has $(size(other.edge,1))")
                end
            end
        else
            DEBUGC && println("dice que $(name) no esta en hybrids")
            n.hybrid = true;
            push!(net.names,string(name));
            n.name = string(name);
            DEBUGC && println("aqui vamos a meter a $(name) en hybrids")
            push!(hybrids,string(name));
            pushNode!(net,n);
            push!(index,size(net.node,1));
            e = Edge(net.numEdges+1);
            e.hybrid = true
            n.leaf ? e.isMajor = false : e.isMajor = true
            pushEdge!(net,e);
            setNode!(e,[n,parent]);
            setEdge!(n,e);
            setEdge!(parent,e);
        end
        e.containRoot = !e.hybrid
    else
        if(hasname)
            push!(net.names,string(name));
            n.name = string(name)
        end
        pushNode!(net,n);
        e = Edge(net.numEdges+1);
        pushEdge!(net,e);
        setNode!(e,[n,parent]);
        setEdge!(n,e);
        setEdge!(parent,e);
    end
    c = peekchar(s);
    if(c == ':')
        c = read(s,Char);
        c = peekchar(s);
        if(isdigit(c))
            length = readFloat(s,c);
            setLength!(e,length);
            c = peekchar(s);
            if(c == ':')
                c = read(s,Char);
                c = peekchar(s);
                if(isdigit(c))
                    length = readFloat(s,c); #bootstrap value
                    c = peekchar(s);
                    if(c == ':')
                        c = read(s, Char);
                        c = peekchar(s);
                        if(isdigit(c))
                            length = readFloat(s,c); #gamma
                            if(!e.hybrid)
                                warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(length) ignored")
                            else
                                setGamma!(e,length, false, true);
                            end
                        else
                            warn("third colon : without gamma value after in $(numLeft[1]) left parenthesis, ignored")
                        end
                    else
                        e.hybrid ? error("hybrid edge $(e.number) read but without gamma value in left parenthesis $(numLeft[1])") : nothing
                    end
                elseif(c == ':')
                    c = read(s, Char);
                    c = peekchar(s);
                    if(isdigit(c))
                        length = readFloat(s,c); #gamma
                        if(!e.hybrid)
                            warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(length) ignored")
                        else
                            setGamma!(e,length, false, true);
                        end
                    else
                        warn("third colon : without gamma value after in $(numLeft[1]) left parenthesis, ignored.")
                    end
                else
                    warn("second colon : read without any double in left parenthesis $(numLeft[1]), ignored.")
                end
            end
        elseif(c == ':')
            e.length = -1.0 # do not use setLength because it does not allow BL too negative
            c = read(s,Char);
            c = peekchar(s);
            if(isdigit(c))
                length = readFloat(s,c); #bootstrap value
                c = peekchar(s);
                if(c == ':')
                    c = read(s, Char);
                    c = peekchar(s);
                    if(isdigit(c))
                        length = readFloat(s,c); #gamma
                        if(!e.hybrid)
                            warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(length) ignored")
                        else
                            setGamma!(e,length, false, true);
                        end
                    else
                        warn("third colon : without gamma value after in $(numLeft[1]) left parenthesis, ignored")
                    end
                else
                    e.hybrid ? warn("hybrid edge $(e.number) read but without gamma value in left parenthesis $(numLeft[1])") : nothing
                end
            elseif(c == ':')
                c = read(s, Char);
                c = peekchar(s);
                if(isdigit(c))
                    length = readFloat(s,c); #gamma
                    if(!e.hybrid)
                        warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(length) ignored")
                    else
                        setGamma!(e,length, false, true);
                    end
                else
                    warn("third colon : without gamma value after in left parenthesis number $(numLeft[1]), ignored")
                end
            else
                warn("second colon : read without any double in left parenthesis $(numLeft[1]), ignored.")
            end
        else
            warn("one colon read without double in left parenthesis $(numLeft[1]), ignored.")
        end
    else
        e.length = -1.0 # do not use setLength because it does not allow BL too negative
    end
    return true
end


# function to read topology from parenthetical format
# input: file name or tree in parenthetical format
# calls readTopology(s::IO)
# warning: crashes if file name starts with (
function readTopology(input::AbstractString,verbose::Bool)
    if(input[1] == '(') #it is a tree
       s = IOBuffer(input)
    else
        try
            s = open(input)
        catch
            error("Could not find or open $(input) file");
        end
       s = open(input)
    end
    net = readTopology(s,verbose)
    return net
end

"""
`readTopology(file name); readTopology(tree in parenthetical format)`

function to read tree or network topology from parenthetical format.
Input: text file or parenthetical format directly
"""
readTopology(input::AbstractString) = readTopology(input,true)

function readTopology(s::IO,verbose::Bool)
    net = HybridNetwork()
    line = readuntil(s,";");
    if(line[end] != ';')
        error("file does not end in ;")
    end
    seekstart(s)
    c = peekchar(s)
    numLeft = [0]; # made Array to make it mutable
    hybrids = ASCIIString[];
    index = Int64[];
    if(c == '(')
       numLeft[1] += 1;
       #println(numLeft)
       n = Node(-1*numLeft[1],false);
       c = read(s,Char)
       b = false;
       while(c != ';')
           b |= readSubtree!(s,n,numLeft,net,hybrids,index)
           c = read(s,Char);
           if(eof(s))
               error("Tree ended while reading in subtree beginning with left parenthesis number $(numLeft[1]).")
           elseif(c == ',')
               continue;
           elseif(c == ')')
               c = peekchar(s);
               if(c == ':')
                   while(c != ';')
                       c = read(s,Char)
                   end
               end
           end
       end
       if(size(n.edge,1) == 1) # root has only one child
           edge = n.edge[1]; # assume it has only one edge
           child = getOtherNode(edge,n);
           removeEdge!(child,edge);
           net.root = getIndex(child,net);
           deleteEdge!(net,edge);
       else
           pushNode!(net,n);
           net.root = getIndex(n,net);
       end
    else
       error("Expected beginning of tree with ( but received $(c) instead")
    end
    storeHybrids!(net)
    checkNumHybEdges!(net)
    ## if(verbose)
    ##     any([(e.length == -1.0 && e.istIdentifiable) for e in net.edge]) && println("edges lengths missing, so assigned default value of -1.0") #fixit: not best approach, better to add a flag inside readSubTree, careful bool not modified inside function, need bool array
    ## end
    net.isRooted = true
    return net
end

readTopology(s::IO) = readTopology(s,true)

# aux function to send an error if the number of hybrid attached to every
# hybrid node are 0,1 or >2
# to be used in readTopology
# need to be run after storeHybrids
function checkNumHybEdges!(net::HybridNetwork)
    if(!isTree(net))
        isempty(net.hybrid) && error("net.hybrid should not be empty for this network")
        for(n in net.hybrid)
            hyb = sum([e.hybrid?1:0 for e in n.edge]);
            if(hyb > 2)
                error("hybrid node $(n.number) has more than two hybrid edges attached to it: polytomy that cannot be resolved without intersecting cycles.")
            elseif(hyb == 1)
                if(net.numHybrids == 1)
                    error("only one hybrid node number $(n.number) with name $(net.names[n.number]) found with one hybrid edge attached")
                else
                    error("current hybrid node $(n.number) with name S(net.names[n.number]) has only one hybrid edge attached. there are other $(net.numHybrids-1) hybrids out there but this one remained unmatched")
                end
            elseif(hyb == 0)
                if(length(n.edge) == 0)
                    error("strange hybrid node $(n.number) attached to 0 edges")
                elseif(length(n.edge) == 1)
                    n.leaf || error("hybrid node $(n.number) has only one tree edge attached and it is not a leaf")
                elseif(length(n.edge) >= 2)
                    warn("hybrid node $(n.number) is not connected to any hybrid edges, it was transformed to tree node")
                    n.hybrid = false;
                end
            end
        end
    end
end


# aux function to solve a polytomy
# warning: chooses one resolution at random
function solvePolytomyRecursive!(net::HybridNetwork, n::Node)
    if(size(n.edge,1) == 4)
        edge1 = n.edge[1];
        edge2 = n.edge[2];
        edge3 = n.edge[3];
        edge4 = n.edge[4];
        removeEdge!(n,edge3);
        removeEdge!(n,edge4);
        removeNode!(n,edge3);
        removeNode!(n,edge4);
        ednew = Edge(net.numEdges+1,0.0);
        max_node = maximum([e.number for e in net.node]);
        n1 = Node(max_node+1,false,false,[edge3,edge4,ednew]);
        setEdge!(n,ednew);
        setNode!(edge3,n1);
        setNode!(edge4,n1);
        setNode!(ednew,[n,n1]);
        pushNode!(net,n1);
        pushEdge!(net,ednew);
    else
        edge1 = n.edge[1];
        removeEdge!(n,edge1);
        solvePolytomyRecursive!(net,n);
        setEdge!(n,edge1);
    end
end

# function to solve a polytomy among tree edges recursively
function solvePolytomy!(net::HybridNetwork, n::Node)
    !n.hybrid || error("cannot solve polytomy in a hybrid node $(n.number).")
    while(size(n.edge,1) > 3)
        solvePolytomyRecursive!(net,n);
    end
end

# aux function to add a child to a leaf hybrid
function addChild!(net::HybridNetwork, n::Node)
    n.hybrid || error("cannot add child to tree node $(n.number).")
    ed1 = Edge(net.numEdges+1,0.0);
    n1 = Node(size(net.names,1)+1,true,false,[ed1]);
    setEdge!(n,ed1);
    setNode!(ed1,[n,n1]);
    pushNode!(net,n1);
    pushEdge!(net,ed1);
end
# aux function to expand the children of a hybrid node
function expandChild!(net::HybridNetwork, n::Node)
    if(n.hybrid)
        suma = sum([!e.hybrid?1:0 for e in n.edge]);
        #println("create edge $(net.numEdges+1)")
        ed1 = Edge(net.numEdges+1,0.0);
        n1 = Node(size(net.names,1)+1,false,false,[ed1]);
        #println("create node $(n1.number)")
        hyb = Edge[];
        for(i in 1:size(n.edge,1))
            !n.edge[i].hybrid ? push!(hyb,n.edge[i]) : nothing
        end
        #println("hyb tiene $([e.number for e in hyb])")
        for(e in hyb)
            #println("se va a borrar a $(e.number)")
            removeEdge!(n,e);
            removeNode!(n,e);
            setEdge!(n1,e);
            setNode!(e,n1);
        end
        #println("now node $(n1.number) has the edges $([e.number for e in n1.edge])")
        setEdge!(n,ed1);
        setNode!(ed1,[n,n1]);
        pushNode!(net,n1);
        pushEdge!(net,ed1);
        if(size(n1.edge,1) > 3)
            solvePolytomy!(net,n1);
        end
    else
        error("cannot expand children of a tree node.")
    end
end

# function to clean topology after readTopology
# looks for:
# TREE:
# - all tree edges must have gamma=1. fixit: cannot point out which doesn't,
#   only shows error.
# - internal nodes with only 2 edges and solves this case.
# - polytomies and choose one resolution at random, issuing a warning
# NETWORK:
# - number of hybrid edges per hybrid node:
#   if 0,1: error (with warning in old functions)
#   if >2: error of hybrid polytomy
#   if 2: check number of tree edges
# - number of tree edges per hybrid node:
#   if 0: leaf hybrid, add child
#   if >1: expand child
#   if 1: check values of gamma:
# - gammas: need to sum to one and be present.
#   error if they do not sum up to one
#   default values of 0.1,0.9 if not present
# leaveRoot=true: leaves the root even if it has only 2 edges (for plotting), default=false
function cleanAfterRead!(net::HybridNetwork, leaveRoot::Bool)
    mod(sum([!e.hybrid?e.gamma:0 for e in net.edge]),1) == 0 ? nothing : error("tree (not network) read and some tree edge has gamma different than 1")
    nodes = copy(net.node)
    for(n in nodes)
        if(isNodeNumIn(n,net.node))
            if(size(n.edge,1) == 2)
                if(!n.hybrid)
                    if(!leaveRoot || !isEqual(net.node[net.root],n)) #if n is the root
                        deleteIntNode!(net,n);
                    end
                else
                    hyb = sum([e.hybrid?1:0 for e in n.edge]);
                    if(hyb == 1)
                        deleteIntNode!(net,n);
                    end
                end
            end
            if(!n.hybrid)
                if(size(n.edge,1) > 3)
                    DEBUG && warn("polytomy found in node $(n.number), random resolution chosen")
                    solvePolytomy!(net,n);
                end
                hyb = sum([e.hybrid?1:0 for e in n.edge]);
                if(hyb == 1)
                    n.hasHybEdge == true;
                elseif(hyb > 1)
                    warn("strange tree node $(n.number) with more than one hybrid edge, intersecting cycles maybe")
                end
            else
                hyb = sum([e.hybrid?1:0 for e in n.edge]);
                tre = sum([!e.hybrid?1:0 for e in n.edge]);
                if(hyb > 2)
                    error("hybrid node $(n.number) has more than two hybrid edges attached to it: polytomy that cannot be resolved without intersecting cycles.")
                elseif(hyb == 1)
                    hybnodes = sum([n.hybrid?1:0 for n in net.node]);
                    if(hybnodes == 1)
                        error("only one hybrid node number $(n.number) with name $(net.names[n.number]) found with one hybrid edge attached")
                    else
                        error("current hybrid node $(n.number) with name S(net.names[n.number]) has only one hybrid edge attached. there are other $(hybnodes-1) hybrids out there but this one remained unmatched")
                    end
                elseif(hyb == 0)
                    warn("hybrid node $(n.number) is not connected to any hybrid edges, it was transformed to tree node")
                    n.hybrid = false;
                else # 2 hybrid edges
                    if(tre == 0) #hybrid leaf
                        warn("hybrid node $(n.number) is a leaf, so we add an extra child")
                        addChild!(net,n);
                    elseif(tre > 1)
                        warn("hybrid node $(n.number) has more than one child so we need to expand with another node")
                        expandChild!(net,n);
                    end
                    suma = sum([e.hybrid?e.gamma:0 for e in n.edge]);
                    if(suma == 2)
                        warn("hybrid edges in read network without gammas")
                        println("hybrid edges for hybrid node $(n.number) do not contain gamma value, set default: 0.9,0.1")
                        for(e in n.edge)
                            if(e.hybrid)
                                (!e.isMajor) ? setGamma!(e,0.1, false, true) : setGamma!(e,0.9, false, true)
                            end
                        end
                    elseif(suma != 1)
                        ed1 = nothing
                        ed2 = nothing
                        for(e in n.edge)
                            if(e.hybrid)
                                isa(ed1,Void) ? ed1=e : ed2=e
                            end
                        end
                        if(ed1.gamma < 1 && ed2.gamma < 1) #both gammas were set, but contradictory
                            error("hybrid edges for hybrid node $(n.number) have gammas that do not sum up to one: $(ed1.gamma),$(ed2.gamma)")
                        elseif(ed1.gamma < 1)
                            warn("only one hybrid edge of hybrid node $(n.number) has gamma value $(ed1.gamma) set, the other edge will be assigned $(1-ed1.gamma).")
                            setGamma!(ed2,1-ed1.gamma, false);
                        else
                            warn("only one hybrid edge of hybrid node $(n.number) has gamma value $(ed2.gamma) set, the other edge will be assigned $(1-ed2.gamma).")
                            setGamma!(ed1,1-ed2.gamma, false);
                        end
                    elseif(suma == 1)
                        for(e in n.edge)
                            if(e.hybrid)
                                if(approxEq(e.gamma,0.5))
                                    e.isMajor = false
                                    break
                                else
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

cleanAfterRead!(net::HybridNetwork) = cleanAfterRead!(net,false)

# function to search for the hybrid nodes in a read network after cleaning it
# and store this information as a network's attribute
function storeHybrids!(net::HybridNetwork)
    flag = true;
    try
        hybrid = searchHybridNode(net)
    catch
        #warn("topology read is a tree as it has no hybrid nodes")
        flag = false;
    end
    if(flag)
        hybrid = searchHybridNode(net);
        net.hybrid = hybrid;
        net.numHybrids = size(hybrid,1);
    end
end

# function to update the read topology after reading
# it will go over the net.hybrid array and check each
# of the hybridization events defined to update:
# - in cycle
# - contain root
# - gammaz
# it uses updateAllNewHybrid! function that
# returns: success (bool), hybrid, flag, nocycle, flag2, flag3
# if tree read, also check that contain root is true for all, ishybrid and hashybedge is false
# warning: needs to have run storeHybrids! before
# warning: it will stop when finding one conflicting hybrid
function updateAllReadTopology!(net::HybridNetwork)
    if(isTree(net))
        #warn("not a network read, but a tree as it does not have hybrid nodes")
        all((e->e.containRoot), net.edge) ? nothing : error("some tree edge has contain root as false")
        all((e->!e.hybrid), net.edge) ? nothing : error("some edge is hybrid and should be all tree edges in a tree")
        all((n->!n.hasHybEdge), net.node) ? nothing : error("some tree node has hybrid edge true, but it is a tree, there are no hybrid edges")
    else
        if(!net.cleaned)
            for(n in net.hybrid)
                success,hyb,flag,nocycle,flag2,flag3 = updateAllNewHybrid!(n,net,false,true,false)
                if(!success)
                    warn("current hybrid $(n.number) conflicts with previous hybrid by intersecting cycles: $(!flag), nonidentifiable topology: $(!flag2), empty space for contain root: $(!flag3), or does not create a cycle (probably problem with the root placement): $(nocycle).")
                    #net.cleaned = false
                end
            end
            DEBUG && println("before update partition")
            DEBUG && printPartitions(net)
            for(n in net.hybrid) #need to updatePartition after all inCycle
                nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
                updatePartition!(net,nodesInCycle)
                DEBUG && println("after updating partition for hybrid node $(n.number)")
                DEBUG && printPartitions(net)
            end
        end
    end
end

# cleanAfterReadAll includes all the step to clean a network after read
function cleanAfterReadAll!(net::HybridNetwork, leaveRoot::Bool)
    DEBUG && println("cleanBL -----")
    cleanBL!(net)
    DEBUG && println("cleanAfterRead -----")
    cleanAfterRead!(net,leaveRoot)
    DEBUG && println("updateAllReadTopology -----")
    updateAllReadTopology!(net) #fixit: it could break if leaveRoot = true (have not checked it), but we need to updateContainRoot
    if(!leaveRoot)
        DEBUG && println("parameters -----")
        parameters!(net)
    end
    DEBUG && println("check root placement -----")
    checkRootPlace!(net)
    net.node[net.root].leaf && warn("root node $(net.node[net.root].number) is a leaf, so when plotting net, it can look weird")
    net.cleaned = true #fixit: set to false inside updateAllReadTopology if problem encountered
    net.isRooted = false
end

cleanAfterReadAll!(net::HybridNetwork) = cleanAfterReadAll!(net,false)

# function to read a topology from file name/tree directly and update it
# by calling updateAllReadTopology after
# leaveRoot=true if the root will not be deleted even if it has only 2 edges
# used for plotting (default=false)
# warning: if leaveRoot=true, net should not be used outside plotting, things will crash
function readTopologyUpdate(file::AbstractString, leaveRoot::Bool,verbose::Bool)
    DEBUG && println("readTopology -----")
    net = readTopology(file,verbose)
    cleanAfterReadAll!(net,leaveRoot)
    return net
end

readTopologyUpdate(file::AbstractString) = readTopologyUpdate(file, false, true)
readTopologyUpdate(file::AbstractString,verbose::Bool) = readTopologyUpdate(file, false, verbose)

"""
`readTopologyLevel1(filename); readTopologyLevel1(parenthetical format)`

same as readTopology, reads a tree or network from parenthetical
format, but this function enforces the necessary conditions for any
starting topology in SNaQ: non-intersecting cycles, no polytomies,
unrooted. It sets any missing branch length to 1.0.
"""
readTopologyLevel1(file::AbstractString) = readTopologyUpdate(file, false, true)


# aux function to check if the root is placed correctly, and re root if not
# warning: it needs updateCR set
function checkRootPlace!(net::HybridNetwork)
    if(!canBeRoot(net.node[net.root]))
        #warn("root node $(net.node[net.root].number) placement is not ok, we will change it to the first found node that agrees with the direction of the hybrid edges")
        for(i in 1:length(net.node))
            if(canBeRoot(net.node[i]))
                net.root = i
                break
            end
        end
    end
    canBeRoot(net.node[net.root]) || error("tried to place root, but couldn't. root is node $(net.node[net.root])")
end

    # --------------------------- write topology -------------------------------------
# function to write a node and its descendants in
#parenthetical format
# di=true in densdroscope format, names=true, prints names
function writeSubTree!(s::IOBuffer, n::Node, parent::Edge,di::Bool,names::Bool, printID::Bool)
    if((parent.hybrid && !parent.isMajor) || n.leaf)
        if(names)
            if(n.name != "")
                print(s,n.name) #instead of number
            else
                if(parent.hybrid)
                    print(s,"#H")
                    print(s,n.number)
                end
            end
        else
            if(parent.hybrid)
                print(s,"#H")
            end
            print(s,n.number)
        end
    else
        print(s,"(")
        numChildren = length(n.edge) - 1
        for(e in n.edge)
            if(!isEqual(e,parent) && !(e.hybrid && parent.hybrid))
                child = getOtherNode(e,n)
                writeSubTree!(s,child,e, di,names, printID)
                if(!parent.hybrid) #cecile handles this step differently: if(parent.hybrid) numChildren-=1 end
                    numChildren -= 1
                    if(numChildren > 0)
                        print(s,",")
                    end
                end
            end
        end
        print(s,")")
        if(parent.hybrid)
            if(!names)
                print(s,string("#H",n.number))
            else
                if(n.name != "")
                    print(s,n.name)
                else
                    print(s,string("#H",n.number))
                end
            end
        end
    end
    printBL = false
    if(printID) #which BL to print
        if(!n.leaf)
            if(parent.istIdentifiable && parent.length >= 0.0) #we do not want to print BL of non-id edges
                print(s,string(":",round(parent.length,3)))
                printBL = true
            end
        end
    else
        if(parent.length >= 0.0) #print all BL != -1
            print(s,string(":",round(parent.length,3)))
            printBL = true
        end
    end
    if(parent.hybrid && !di && !n.isBadDiamondI)
        if(!printBL)
            if(parent.gamma != 1.0)
                print(s,string(":::",round(parent.gamma,3)))
            end
        else
            if(parent.gamma != 1.0)
                print(s,string("::",round(parent.gamma,3)))
            end
        end
    end
end

# function to writeTopology
# if string=true, returns a string with network in parenthetical format
#                 ow returns the IOBuffer object
# need as input HybridNetwork, since QuartetNetwork does not have root
# input di=true if written for Dendroscope (without gammas)
# names=true, writes the names instead of node numbers, default true
# outgroup: place the root in the external edge of this taxon if possible,
# if none given, placed the root wherever possible
# printID=true, only print identifiable BL, default false (only true inside snaq)
function writeTopology(net0::HybridNetwork, di::Bool, string::Bool, names::Bool,outgroup::AbstractString, printID::Bool)
    s = IOBuffer()
    net = deepcopy(net0) #writeTopology needs containRoot, but should not alter net0
    if(net.numBad > 0)
        println("net has $(net.numBad) bad diamond I, gammas and some branch lengths are not identifiable, and therefore, meaningless")
    end
    if(net.numNodes == 1)
        print(s,string(net.node[net.root].number,";"))
    else
        if(!isTree(net) && !net.cleaned)
            DEBUG && println("net not cleaned inside writeTopology, need to run updateContainRoot")
            for(n in net.hybrid)
                flag,edges = updateContainRoot!(net,n)
                flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
            end
        end
        print(s,"(")
        updateRoot!(net,outgroup)
        #DEBUG && printEverything(net)
        CHECKNET && canBeRoot(net.node[net.root])
        degree = length(net.node[net.root].edge)
        for(e in net.node[net.root].edge)
            writeSubTree!(s,getOtherNode(e,net.node[net.root]),e,di,names, printID)
            degree -= 1
            if(degree > 0)
                print(s,",")
            end
        end
        print(s,");")
    end
    outgroup != "none" && undoRoot!(net) #to delete the node with only two edges
    if(string)
        return bytestring(s)
    else
        return s
    end
end

#writeTopology(net::HybridNetwork) = writeTopology(net,false, true,true,"none") #not needed because of last function definition
writeTopology(net::HybridNetwork,printID::Bool) = writeTopology(net,false, true,true,"none",printID)
writeTopology(net::HybridNetwork,outgroup::AbstractString) = writeTopology(net,false, true,true,outgroup,true)
writeTopology(net::HybridNetwork,di::Bool,outgroup::AbstractString) = writeTopology(net,di, true,true,outgroup,true)

"""
`writeTopology(net::HybridNetwork)`

writes the parenthetical format of a HybridNetwork object with many optional arguments:

- di=true: write in format for Dendroscope (default false)
- names=false: write the leaf nodes numbers instead of taxon names (default true)
- outgroup (string): name of outgroup to root the tree/network
- printID=true, only print branch lengths for identifiable egdes according to the snaq estimation procedure (default false)
"""
writeTopology(net::HybridNetwork; di=false::Bool, string=true::Bool, names=true::Bool,outgroup="none"::AbstractString, printID=false::Bool) = writeTopology(net, di, string, names,outgroup,printID)

# function to check if root is well-placed
# and look for a better place if not
# searches on net.node because net.root is the index in net.node
# if we search in net.edge, we then need to search in net.node
function updateRoot!(net::HybridNetwork, outgroup::AbstractString)
    checkroot = false
    if(outgroup == "none")
        DEBUG && println("no outgroup defined")
        checkroot = true
    else
        println("outgroup defined $(outgroup)")
        try
            index = getIndex(true,[isequal(outgroup,n.name) for n in net.node])
        catch
            error("outgroup $(outgroup) not in net.names $(net.names)")
        end
        index = getIndex(true,[isequal(outgroup,n.name) for n in net.node])
        node = net.node[index]
        node.leaf || error("outgroup $(outgroup) is not a leaf in net")
        length(net.node[index].edge) == 1 || error("strange leaf $(outgroup), node number $(net.node[index].number) with $(length(net.node[index].edge)) edges instead of 1")
        edge = net.node[index].edge[1]
        if(edge.containRoot)
            DEBUG && println("creating new node in the middle of the external edge $(edge.number) leading to outgroup $(node.number)")
            othernode = getOtherNode(edge,node)
            removeEdge!(othernode,edge)
            removeNode!(othernode,edge)
            max_edge = maximum([e.number for e in net.edge]);
            max_node = maximum([e.number for e in net.node]);
            newedge = Edge(max_edge+1) #fixit: maybe this edge not identifiable, need to add that check
            newnode = Node(max_node+1,false,false,[edge,newedge])
            if(net.cleaned && !isTree(net))
                part = whichPartition(net,edge)
                push!(net.partition[part].edges,newedge)
            end
            setNode!(edge,newnode)
            setNode!(newedge,newnode)
            setEdge!(othernode,newedge)
            setNode!(newedge,othernode)
            pushEdge!(net,newedge)
            pushNode!(net,newnode)
            t = edge.length
            setLength!(edge,t/2)
            setLength!(newedge,t/2)
            net.root = length(net.node) #last node is root
        else
            warn("external edge $(net.node[index].edge[1].number) leading to outgroup $(outgroup) cannot contain root, root placed wherever")
            checkroot = true
        end
    end
    if(checkroot && !isTree(net))
        checkRootPlace!(net)
    end
end

# function to check if a node could be root
# by the containRoot attribute of edges around it
function canBeRoot(n::Node)
    !n.hybrid || return false
    #!n.hasHybEdge || return false #need to allow for some reason, check ipad notes
    !n.leaf || return false
    return any([e.containRoot for e in n.edge])
end

# function to delete the extra node created in updateRoot
# this extra node is needed to be able to compare networks with the distance function
# but if left in the network, everything crashes (as everything assumes three edges per node)
# fromUpdateRoot=true if called after updateRoot (in which case leaf has to be a leaf), ow it is used in readTopUpdate
function undoRoot!(net::HybridNetwork, fromUpdateRoot::Bool)
    if(length(net.node[net.root].edge) == 2)
        root = net.node[net.root]
        leaf = getOtherNode(root.edge[1],root).leaf ? getOtherNode(root.edge[1],root) : getOtherNode(root.edge[2],root)
        (fromUpdateRoot && leaf.leaf) || error("root should have one neighbor leaf which has to be the outgroup defined")
        deleteIntLeafWhile!(net,root,leaf);
    end
end

undoRoot!(net::HybridNetwork) = undoRoot!(net, true)


# function to read the .out file from snaq (optTopRuns) function
function readOutfile(file::AbstractString)
    try
        s = open(file)
    catch
        error("Could not find or open $(file) file");
    end
    s = open(file)
    line = readline(s)
    DEBUG && println("line read $(line)")
    c = line[1]
    if(c == '(')
       println("Estimated network from file $(file): $(line)")
       net = readTopologyUpdate(line)
    else
        error("output file $(filename).out does not contain a tree in the first line, instead it has $(line)")
    end
    return net
end

"""
`readSnaqNetwork(output file)`

function to read the estimated network from an .out file generated by the snaq function
"""
readSnaqNetwork(file::AbstractString) = readOutfile(file)


# function to change branch lengths of -1.0 to 1.0 for starting topology
function cleanBL!(net::HybridNetwork)
    ##println("missing branch lengths will be set to 1.0")
    for(e in net.edge)
        if(e.length == -1.0)
            setLength!(e,1.0)
        end
    end
end
