__precompile__()

module PhyloNetworks

    # stdlib (standard libraries)
    using Dates
    using Distributed
    using LinearAlgebra: diag, I, logdet, norm, LowerTriangular, mul!, lmul!, Diagonal, cholesky, BLAS
    # alternative: drop support for julia v1.4, because LinearAlgebra.rotate! requires julia v1.5
    # using LinearAlgebra # bring all of LinearAlgebra into scope
    # import LinearAlgebra.rotate! # allow re-definition of rotate!
    using Printf: @printf, @sprintf
    using Random
    using Statistics: mean, quantile, median

    # other libraries, indicate compatible version in Project.toml
    using BioSequences
    using BioSymbols
    using Combinatorics: combinations
    using CSV
    using DataFrames # innerjoin new in v0.21
    using DataStructures # for updateInCycle with priority queue
    using Distributions #for RateVariationAcrossSites
    using FASTX
    using GLM # for the lm function
    using NLopt # for branch lengths optimization
    using StaticArrays
    using StatsBase # sample, coef etc.
    using StatsFuns # logsumexp, logaddexp, various cdf
    using StatsModels # re-exported by GLM. for ModelFrame ModelMatrix Formula etc

    import Base: show
    import GLM: ftest

    const DEBUGC = false # even more debug messages
    global CHECKNET = false # for debugging only

    export ftest
    export
        ## Network Definition
        HybridNetwork,
        DataCF,
        Quartet,
        readTopology,
        readTopologyLevel1,
        tipLabels,
        writeTopology,
        writeSubTree!,
        hybridlambdaformat,
        deleteleaf!,
        printEdges,
        printNodes,
        sorttaxa!,
        ## SNAQ
        readTrees2CF,
        countquartetsintrees,
        readTableCF,
        readTableCF!,
        writeTableCF,
        mapAllelesCFtable,
        readInputTrees,
        readNexusTrees,
        summarizeDataCF,
        snaq!,
        readSnaqNetwork,
        topologyMaxQPseudolik!,
        topologyQPseudolik!,
        ## Network Manipulation
        # getParent, getParents, getMajorParentEdge, getMinorParentEdge, getChildren,
        # functions above: first rename them throughout to be consistent with other packages, like:
        # parent child parents children parentmajor parentminor ancestor sibling offspring
        rootatnode!,
        rootonedge!,
        directEdges!,
        preorder!,
        cladewiseorder!,
        fittedQuartetCF,
        rotate!,
        setLength!,
        setGamma!,
        deleteHybridThreshold!,
        displayedTrees,
        majorTree,
        minorTreeAt,
        displayedNetworkAt!,
        hardwiredClusters,
        hardwiredCluster,
        hardwiredCluster!,
        hardwiredClusterDistance,
        biconnectedComponents,
        blobDecomposition!,
        blobDecomposition,
        mapindividuals,
        nni!,
        checkroot!,
        treeedgecomponents,
        ## Network Bootstrap
        treeEdgesBootstrap,
        hybridDetection,
        summarizeHFdf,
        hybridBootstrapSupport,
        bootsnaq,
        readBootstrapTrees,
        writeMultiTopology,
        readMultiTopologyLevel1,
        readMultiTopology,
        hybridatnode!,
        undirectedOtherNetworks,
        ## Network Calibration
        getNodeAges,
        pairwiseTaxonDistanceMatrix,
        calibrateFromPairwiseDistances!,
        ## Network PCM
        phyloNetworklm,
        PhyloNetworkLinearModel,
        simulate,
        TraitSimulation,
        ParamsBM,
        ParamsMultiBM,
        ShiftNet,
        shiftHybrid,
        getShiftEdgeNumber,
        getShiftValue,
        sharedPathMatrix,
        descendenceMatrix,
        regressorShift,
        regressorHybrid,
        ancestralStateReconstruction,
        ReconstructedStates,
        sigma2_estim,
        mu_estim,
        lambda_estim,
        expectations,
        expectationsPlot,
        predint,
        predintPlot,
        vcv,
        ## Discrete Trait PCM
        parsimonySoftwired,
        parsimonyGF,
        maxParsimonyNet,
        TraitSubstitutionModel,
        EqualRatesSubstitutionModel,
        BinaryTraitSubstitutionModel,
        TwoBinaryTraitSubstitutionModel,
        JC69, HKY85,
        nstates,
        Q,
        getlabels,
        nparams,
        RateVariationAcrossSites,
        randomTrait,
        randomTrait!,
        fitdiscrete,
        readfastatodna,
        stationary,
        empiricalDNAfrequencies,
        phyLiNC,
        # neighbor joining
        nj

    include("types.jl")
    include("auxiliary.jl")
    include("update.jl")
    include("undo.jl")
    include("addHybrid_snaq.jl")
    include("addHybrid.jl")
    include("deleteHybrid.jl")
    include("moves_snaq.jl")
    include("moves_semidirected.jl")
    include("readwrite.jl")
    include("readData.jl")
    include("snaq_optimization.jl")
    include("pseudolik.jl")
    include("descriptive.jl")
    include("manipulateNet.jl")
    include("bootstrap.jl")
    include("multipleAlleles.jl")
    include("compareNetworks.jl")
    include("traits.jl")
    include("parsimony.jl")
    include("pairwiseDistanceLS.jl")
    include("interop.jl")
    include("substitutionModels.jl")
    include("graph_components.jl")
    include("traitsLikDiscrete.jl")
    include("deprecated.jl")
    include("nj.jl")
    include("phyLiNCoptimization.jl")

end #module
