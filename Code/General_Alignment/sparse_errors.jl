include("bad_general_pairwise_alignment.jl")
include("custom_2D_stats.jl")
export kmerMatching, KmerMatch, Move, general_pairwise_aligner, Gap, make_match_score_matrix
using BioSequences
using Kmers

A = LongDNA{2}("ACGGTTAGCGCGCAAGGTCGATGTGTGTGTGTG")
B = LongDNA{2}("TCGGTTACGCGCAAGGTCGATGAGTGTGTGTG")

function CorrelationScore(s::SampleMetrics2D)
    return sqrt(s.n) * s.correlation
end

struct KmerMatch
    posA::Int
    posB::Int
end

function kmerMatching(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, moves::Array{Move}, affine_gap::Float64, kmerLength::Int64 = 12) 
    return kmerMatching(A, B, make_match_score_matrix(match_score, mismatch_score), moves, affine_gap, kmerLength)
end
function kmerMatching(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, moves::Array{Move}, affine_gap::Float64, kmerLength::Int64 = 12)
    k = kmerLength
    kmerDict = Dict{UInt128, Array{Int64}}()

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
        kmer = hash(B[iB : iB+k-1])
        if haskey(kmerDict, kmer)
            for iA in kmerDict[kmer]
                if diagonals[iA - iB + n + 1] + k <= iA
                    push!(kmerMatches, KmerMatch(iA, iB))
                else
                    println("deleted kmer $iA, $iB")
                end
                diagonals[iA - iB + n + 1] = iA
            end
        end
    end
    @show kmerMatches

    #Prune set of kmers
    matchCount = length(kmerMatches)
    deletionFlags = BitArray([0 for i in 1:matchCount])
    println(deletionFlags)
    kmerMetrics = GetSampleMetrics2D(map(x -> x.posA, kmerMatches), map(x -> x.posB, kmerMatches))
    corScore = CorrelationScore(kmerMetrics)
    @show corScore
    for i in 1 : matchCount
        bestDeletion = -1
        if kmerMetrics.n <= 2
            break
        end
        for j in 1 : matchCount
            if !deletionFlags[j]
                kmer = kmerMatches[j]
                newCorScore = CorrelationScore(RemoveVector(kmerMetrics, kmer.posA, kmer.posB))
                if newCorScore > corScore
                    corScore = newCorScore
                    bestDeletion = j
                end
            end
        end
        if bestDeletion == -1
            break
        else
            deletionFlags[bestDeletion] = true
            kmer = kmerMatches[bestDeletion]
            kmerMetrics = RemoveVector(kmerMetrics, kmer.posA, kmer.posB)
            corScore = CorrelationScore(kmerMetrics)
            println("Deleted kmer $(A[kmer.posA:kmer.posA+k-1]) at position $(kmer.posA), $(kmer.posB)")
            println("New correlation score = $corScore")
        end
    end

    #Ensure that kmers are compatible
    filteredKmerMatches = KmerMatch[]
    prevA = -k
    prevB = -k
    for i in 1 : matchCount
        match = kmerMatches[i]
        if !deletionFlags[i] && prevA + k <= match.posA && prevB + k <= match.posB
            prevA = match.posA
            prevB = match.posB
            push!(filteredKmerMatches, match)
        end
    end

    @show kmerMatches
    prevA = -k+1
    prevB = -k+1
    #Join kmers using needleman Wunsch
    #To first Kmer
    result = [LongDNA{4}(""), LongDNA{4}("")]
    for kmer in filteredKmerMatches
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
# kmerMatching(A, B, match_score, mismatch_score, moveset, affine_score, kmerLength)