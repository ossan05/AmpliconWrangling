include("sparse_errors.jl")
include("sparse_errors_2.jl")
include("bad_general_pairwise_alignment.jl")
include("new_general_pairwise_alignment.jl")

using CodecXz
using BioSequences, FASTX, NextGenSeqUtils

#Because the old one is broken by a FASTX update
function read_fasta1(filename)
    stream = XzDecompressorStream(open(filename))
    decompressed = FASTA.Reader(stream)
    seqs = [String(sequence(record)) for record in decompressed]
    close(stream)
    return seqs
end

# println(read_fasta1("sequences.fasta.xz"))
seqs = read_fasta1("sequences.fasta.xz")
modSeqs = LongDNA{2}.(filter.(x -> contains("ACGT", x), seqs))


mismatch_score = 2.0
match_score = 0.0
kmerLength = 30
affine_score = 0.5
matches = [Move(1, 0)]
gaps = [Move(1, 2.0)]

chain_rel_scores = Vector{Float64}()
correlation_rel_scores = Vector{Float64}()
old_rel_scores = Vector{Float64}()

chain_misses = 0
correlation_misses = 0
old_misses = 0

for i in 1 : 2
    println("iteration ", i)
    A = rand(modSeqs)
    typeof(A)
    B = rand(modSeqs)
    corr = kmerMatching(A, B, match_score, mismatch_score, matches, gaps, gaps, affine_score, kmerLength)
    chain = kmerChainMatching(A, B, match_score, mismatch_score, matches, gaps, gaps, kmerLength)
    #old = kmer_seeded_align(SCring(A), String(B);wordlength = kmerLength, skip = Int64(kmerLength/3))
    a = [corr, chain]#, old]
    b = [0.0, 0.0]#, 0.0]
    for j in axes(a)
        x = a[j]
        @show typeof(x)
        @show typeof(x[1])
        b[j] = alignment_score(x[1], x[2], make_match_score_matrix(match_score, mismatch_score), matches, gaps, gaps)
    end
    println(b)
end
#A, B = kmerMatching(modSeqs[20], modSeqs[19], match_score, mismatch_score, matches, gaps, gaps, affine_score, kmerLength)
#A, B = kmerChainMatching(modSeqs[20], modSeqs[19], match_score, mismatch_score, matches, gaps, gaps, kmerLength)
#A, B = kmer_seeded_align(String(modSeqs[20]), String(modSeqs[19]);wordlength = kmerLength, skip = Int64(kmerLength/3))
println(kmer_matches, alignment_score(A, B, make_match_score_matrix(match_score, mismatch_score), matches, gaps, gaps))
#println(alignment_score(LongDNA{4}.(A), LongDNA{4}.(B), make_match_score_matrix(match_score, mismatch_score), matches, gaps, gaps))