include("bad_general_pairwise_alignment.jl")
include("custom_2D_stats.jl")
export kmerChainMatching, KmerMatch, Move, general_pairwise_aligner, Gap, make_match_score_matrix
using BioSequences
using Kmers

A = LongDNA{2}("ACGGTTAGCGCGCAAGGTCGATGTGTGTGTGTG")
B = LongDNA{2}("TCGGTTACGCGCAAGGTCGATGAGTGTGTGTG")

function getCorrelationScore(s::SampleMetrics2D)
    return sqrt(s.n) * s.correlation
end

#Needs refinement
function getTransitionScore(x, y, k)
    return max(y.posA - x.posA, y.posB - x.posB) - k
end

struct KmerMatch
    posA::Int64
    posB::Int64
end

struct Endpoint
    x::Int64
    y::Int64
    id::Int64
    isBeginning::Bool
end

function kmerChainMatching(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, moves::Array{Move}, affine_gap::Float64, kmerLength::Int64 = 12) 
    return kmerChainMatching(A, B, make_match_score_matrix(match_score, mismatch_score), moves, affine_gap, kmerLength)
end
function kmerChainMatching(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, moves::Array{Move}, affine_gap::Float64, kmerLength::Int64 = 12)
    k = kmerLength
    kmerDict = Dict{UInt128, Array{Int64}}()
    repetitionTolerance = 5

    m = length(A)
    n = length(B)
    #Initialize List of Kmer Matches
    kmerMatches = KmerMatch[]
    #Initialize Kmer Dictionary from A
    #k = min(m-2, n-2, floor(Int64, log2(m*n)), 30)
    for i in 1:m-k+1
        kmer = hash(A[i:i+k-1])
        if kmer in keys(kmerDict)
            push!(kmerDict[kmer], i)
        else
            kmerDict[kmer] = [i]
        end
    end

    #Compare B with dictionary
    diagonals = fill(typemin(Int64), m+n) #diagonals[i] is the kmer in diagonal i furthest to the right. (used to find overlaps in matches)
    for iB in 1 : n-k+1
        kmer = B[iB : iB+k-1]
        hashedKmer = hash(kmer)
        if !haskey(kmerDict, hashedKmer)
           continue
        end 
        if length(kmerDict) < repetitionTolerance
            continue
        end
        for iA in kmerDict[hashedKmer]
            if diagonals[iA - iB + n + 1] + k <= iA
                push!(kmerMatches, KmerMatch(iA, iB))
            end
            diagonals[iA - iB + n + 1] = iA
        end
    end
    matchCount = length(kmerMatches)
    #Prune set of kmers  
    endpoints = Array{Endpoint}(undef,matchCount * 2)
    for i in 1 : matchCount
        match = kmerMatches[i]
        endpoints[2 * i - 1] = Endpoint(match.posA, match.posB, i, true)
        endpoints[2 * i] = Endpoint(match.posA + k, match.posB + k, i, false)
    end
    sort!(endpoints, by = e -> (e.x, e.y, e.isBeginning))
    
    chains = Vector{Vector{Int64}}()
    chainIds = Array{Int64}(undef, matchCount)
    bestTransitions = Array{Int64}(undef, matchCount)
    bestScores = Array{Int64}(undef, matchCount)
    occupiedChains = Vector{Bool}()
    chainNum = 0
    transitions = Vector{Vector{Int64}}()
    for e in endpoints
        if e.isBeginning
            
            #Find the closest chain it can connect with
            bestChain = -1
            bestChainY = -1
            for j in 1 : chainNum
                if !occupiedChains[j]
                    chainY = kmerMatches[last(chains[j])].posB + k
                    if chainY < e.y && chainY > bestChainY
                        bestChain = j
                        bestChainY = chainY
                    end
                end
            end

            #Add it to the chain or create a new one
            if bestChain == -1
                chainNum += 1
                bestChain = chainNum
                push!(transitions, fill(0, chainNum)) 
                for j in 1 : chainNum-1
                    push!(transitions[j], 0)
                end
                push!(chains, [])
                push!(occupiedChains, true)
            else
                occupiedChains[bestChain] = true
            end
            
            #Find best transition using the chains
            chainIds[e.id] = bestChain
            bestScores[e.id] = getTransitionScore(KmerMatch(-k,-k), KmerMatch(e.x, e.y), k)
            bestTransitions[e.id] = e.id
            for j in 1 : chainNum
                t = transitions[bestChain][j]
                while t + 1 <= length(chains[j]) && kmerMatches[chains[j][t + 1]].posB + k <= e.y
                    t += 1
                end
                transitions[bestChain][j] = t
                if t != 0
                    candidate = chains[j][t]
                    score = getTransitionScore(kmerMatches[candidate], KmerMatch(e.x, e.y), k) + bestScores[candidate]
                    if score < bestScores[e.id]
                        bestScores[e.id] = score
                        bestTransitions[e.id] = candidate
                    end
                end
            end
        else
            occupiedChains[chainIds[e.id]] = false
            push!(chains[chainIds[e.id]], e.id)
        end
    end
    for i in 1 : chainNum
        @show i
        println(map(x -> [kmerMatches[x].posA, kmerMatches[x].posB], chains[i]))
    end
    for i in 1 : matchCount
        println(kmerMatches[i], " -> ", kmerMatches[bestTransitions[i]], "  ", bestScores[i])
    end

    #Back propagation
    kmerPath = Vector{KmerMatch}()
    if matchCount > 0
        kmer = argmin([bestScores[i] + getTransitionScore(kmerMatches[i], KmerMatch(m, n), k) for i in 1 : matchCount])
        while true
            push!(kmerPath, kmerMatches[kmer])
            if bestTransitions[kmer] == kmer
                break
            else
                kmer = bestTransitions[kmer]
            end
        end
        reverse!(kmerPath)
    end
    @show kmerPath
    

    prevA = -k+1
    prevB = -k+1
    #Join kmers using needleman Wunsch
    #To first Kmer
    result = [LongDNA{4}(""), LongDNA{4}("")]
    for kmer in kmerPath
        if !(kmer.posA == prevA + k && kmer.posB == prevB + k)
            result .*= general_pairwise_aligner(A[prevA + k : kmer.posA - 1], B[prevB + k : kmer.posB - 1], match_score_matrix, moves, affine_gap)
        end
        result .*= [A[kmer.posA : kmer.posA + k - 1], B[kmer.posB : kmer.posB + k - 1]]
        prevA = kmer.posA
        prevB = kmer.posB
    end
    result .*= general_pairwise_aligner(A[prevA + k : m], B[prevB + k : n], match_score_matrix, moves, affine_gap)

    return (result[1], result[2])
end

# A = LongDNA{2}("ACGGTTAGCGCGCAAGGTCGATGTGTGTGTGTGTG")
# B = LongDNA{2}("TCGGTTACGCGCAAGGTCGATGAGTGTGTGTG")
# mismatch_score = 0.5
# match_score = 0.0
# kmerLength = 5
# affine_score = 0.5
# moveset = [Move((1, 1), 0), Move((1, 0), 1), Move((0, 1), 1), Move((3, 3), 0), Move((3, 0), 2), Move((0, 3), 2)]
# alignment = kmerChainMatching(A, B, match_score, mismatch_score, moveset, affine_score, kmerLength)
# println(alignment[1])
# println(alignment[2])