#include("sparse_errors.jl")
include("bad_general_pairwise_alignment.jl")
include("sequence_generator.jl")
include("sparse_errors_2.jl")
using BioSequences
using Kmers
using NextGenSeqUtils
using BenchmarkTools

#A = LongDNA{2}("ACGGTTAGCGCGCAAGGTCGATGTGTGTGTGTG")
#B = LongDNA{2}("TCGGTTACGCGCAAGGTCGATGAGTGTGTGTG")
mismatch_score = 0.5
match_score = 0.0
kmerLength = 6
affine_score = 0.5

# good = 0
# bad = 0
# for i in 1 : 100
#     A, B = generate_seq(50)
    
#     alignment1 = kmerMatching(A, B, match_score, mismatch_score, moveset, affine_score, kmerLength)
#     alignment2 = general_pairwise_aligner(A, B, match_score, mismatch_score, moveset, affine_score)
#     println("kmerAlign:")
#     println(alignment1[1])
#     println(alignment1[2])
#     println("generalAlign:")
#     println(alignment2[1])
#     println(alignment2[2])

#     if reduce((a, b) -> a && b, alignment1 .== alignment2)
#         global good += 1
#     else
#         global bad += 1
#     end
# end
# @show good
# @show bad

# A, B = generate_seq(20)
# @benchmark triplet_nw_align(String(A), String(B))
# @benchmark initiate_general_pairwise_aligner(A, B)
# @benchmark affine_nw_align(String(A), String(B))

# @benchmark (kmer_seeded_align($(String(A)), $(String(B)); wordlength = 15, skip = 5))
# @benchmark (initiate_kmerMatching($A, $B))

function calculate_score(A, B, match_score, mismatch_score, moveset, affine_score)
    if length(A) != length(B)
        error("Aligned sequences must be of equal lengths.")
    end
    n = length(A)
    score = 0
    gap = 1
    streak = 0
    frame = 1

    for i in axes(A)
        if A[i] == DNA_Gap && B[i] == DNA_Gap
            score += gap
            streak += 1
            if streak % 3 == 0
                score -= frame
            end
            if streak > 1
                score -= affine
            end
        elseif A[i] != B[i]
            score += mis
            streak = 0
        else
            streak = 0
        end

        if A[i] == DNA_Gap
            score += gap
            streak += 1
            if streak % 3 == 0
                score -= frame
            end
            if streak > 1
                score -= affine
            end
        elseif A[i] == DNA_Gap && A[i] != B[i]
            score += mis
            streak = 0
        else
            streak = 0
        end
    end
    return score
end

A, B = generate_seq(20)


align1 = kmer_seeded_align(String(A), String(B))
align2 = kmerChainMatching(A, B, match_score, mismatch_score, moveset, affine_score, kmerLength)


println(align1[1])
println(align1[2])
println(calculate_score(align1...))
println(align2[1])
println(align2[2])
println(calculate_score(align2...))
#println(calculate_score(initiate_kmerMatching(A, B)...))