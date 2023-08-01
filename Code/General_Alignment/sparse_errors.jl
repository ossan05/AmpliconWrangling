include("bad_general_pairwise_alignment.jl")
include("custom_2D_stats.jl")
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

function initiate_kmerMatching(A::LongDNA{2}, B::LongDNA{2})
    mismatch_score = 0.5
    match_score = 0.0
    kmerLength = 30
    affine_score = 0.5
    matches = [Move(1, 0), Move(3, 0)]
    gaps = [Move(1, 1), Move(3, 2)]
    kmerMatching(A, B, match_score, mismatch_score, matches, gaps, gaps, affine_score, kmerLength)
end

function kmerMatching(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, match_moves::Vector{Move}, hgap_moves::Vector{Move}, vgap_moves::Vector{Move}, affine_gap::Float64, kmerLength::Int64 = 12) 
    return kmerMatching(A, B, make_match_score_matrix(match_score, mismatch_score), match_moves, hgap_moves, vgap_moves, affine_gap, kmerLength)
end

function kmerMatching(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, match_moves::Vector{Move}, hgap_moves::Vector{Move}, vgap_moves::Vector{Move}, affine_gap::Float64, kmerLength::Int64 = 12)
    k = kmerLength
    kmerDict = Dict{UInt128, Array{Int64}}()
    repetitionTolerance = 5

    m, n = length(A), length(B)
    #Initialize List of Kmer Matches
    kmerMatches = KmerMatch[]
    #Initialize Kmer Dictionary from A
    for i in 1:m-k+1
        kmer = hash(A[i : i + k - 1])
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
        if length(kmerDict[hashedKmer]) > repetitionTolerance
            continue
        end
        for iA in kmerDict[hashedKmer]
            if diagonals[iA - iB + n + 1] + k <= iA
                push!(kmerMatches, KmerMatch(iA, iB))
                diagonals[iA - iB + n + 1] = iA
            end
        end
    end
    @show kmerMatches
    #Prune set of kmers
    matchCount = length(kmerMatches)
    remaining_matches = matchCount
    deletionFlags = BitArray([0 for i in 1:matchCount])
    kmerMetrics = GetSampleMetrics2D(map(x -> x.posA, kmerMatches), map(x -> x.posB, kmerMatches))
    corScore = CorrelationScore(kmerMetrics)
    deletion_candidates = Array{Tuple{Float64, Int64}}(undef, matchCount)
    pruning_steps = 10
    max_pruning_fraction = 0.15
    for i in 1 : pruning_steps
        if kmerMetrics.n <= 2
            break
        end
        for j in 1 : matchCount
            if !deletionFlags[j]
                kmer = kmerMatches[j]
                deletion_candidates[j] = (CorrelationScore(RemoveVector(kmerMetrics, kmer.posA, kmer.posB)), j)
            end
        end
        max_pruning_num = floor(Int64,remaining_matches*max_pruning_fraction + 1.0)
        println(corScore)
        sort!(deletion_candidates, rev = true)
        cand = map(x -> (x[1], kmerMatches[x[2]]), deletion_candidates)
        @show cand
        prev_remaining_matches = remaining_matches
        for j in 1 : max_pruning_num
            if deletion_candidates[j][1] > corScore && remaining_matches > 3
                kmer = deletion_candidates[j][2]
                println("deleted ", kmerMatches[kmer].posA, " ", kmerMatches[kmer].posB)
                deletionFlags[kmer] = true
                kmerMetrics = RemoveVector(kmerMetrics, kmerMatches[kmer].posA, kmerMatches[kmer].posB)
                corScore = CorrelationScore(kmerMetrics)
                remaining_matches -= 1
            else
                break
            end
        end
        if remaining_matches == prev_remaining_matches
            break
        end
    end
    @show kmerMatches[filter(i -> !deletionFlags[i], 1 : matchCount)]
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
    @show filteredKmerMatches

    prevA = -k+1
    prevB = -k+1
    #Join kmers using needleman Wunsch
    #To first Kmer
    result = [LongDNA{4}(""), LongDNA{4}("")]
    for kmer in filteredKmerMatches
        if !(kmer.posA == prevA + k && kmer.posB == prevB + k)
            result .*= general_pairwise_aligner(A[prevA + k : kmer.posA - 1], B[prevB + k : kmer.posB - 1], match_score_matrix, match_moves, hgap_moves, vgap_moves, affine_gap)
        end
        result .*= [A[kmer.posA : kmer.posA + k - 1], B[kmer.posB : kmer.posB + k - 1]]
        prevA = kmer.posA
        prevB = kmer.posB
    end
    result .*= general_pairwise_aligner(A[prevA + k : m], B[prevB + k : n], match_score_matrix, match_moves, hgap_moves, vgap_moves, affine_gap)
    return (result[1], result[2])
end

# #A = LongDNA{2}("ACGGTTAGCGCGCAAGGTCGATGTGTGTGTGTGTG")
# #B = LongDNA{2}("TCGGTTACGCGCAAGGTCGATGAGTGTGTGTG")
# A, B = generate_seq(3000)

# mismatch_score = 1.0
# match_score = 0.0
# kmerLength = 5
# affine_score = 0.5
# match_moves = [Move(1, 0.0), Move(3, 0.0)]
# gap_moves = [Move(3, 2.0), Move(1, 1.0)]
# kmerMatching(A, B, match_score, mismatch_score, match_moves, gap_moves, gap_moves, affine_score, kmerLength)