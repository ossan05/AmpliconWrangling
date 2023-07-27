include("sparse_errors.jl")
include("bad_general_pairwise_alignment.jl")
include("sequence_generator.jl")
using BioSequences
using Kmers
using BenchmarkTools

#A = LongDNA{2}("ACGGTTAGCGCGCAAGGTCGATGTGTGTGTGTG")
#B = LongDNA{2}("TCGGTTACGCGCAAGGTCGATGAGTGTGTGTG")
mismatch_score = 0.5
match_score = 0.0
kmerLength = 5
affine_score = 0.5
moveset = [Move((1, 1), 0), Move((1, 0), 1), Move((0, 1), 1), Move((3, 3), 0), Move((3, 0), 2), Move((0, 3), 2)]





good = 0
bad = 0
for i in 1 : 100
    A, B = generate_seq(50)
    
    alignment1 = kmerMatching(A, B, match_score, mismatch_score, moveset, affine_score, kmerLength)
    alignment2 = general_pairwise_aligner(A, B, match_score, mismatch_score, moveset, affine_score)
    println("kmerAlign:")
    println(alignment1[1])
    println(alignment1[2])
    println("generalAlign:")
    println(alignment2[1])
    println(alignment2[2])

    if reduce((a, b) -> a && b, alignment1 .== alignment2)
        global good += 1
    else
        global bad += 1
    end
end
@show good
@show bad