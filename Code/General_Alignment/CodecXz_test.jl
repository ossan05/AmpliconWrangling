include("sparse_errors.jl")
include("sparse_errors_2.jl")
include("bad_general_pairwise_alignment.jl")
include("new_general_pairwise_alignment.jl")

using CodecXz
using BioSequences, FASTX

#Because the old one is broken by a FASTX update
function read_fasta(filename)
    stream = XzDecompressorStream(open(filename))
    decompressed = FASTA.Reader(stream)
    seqs = [String(sequence(record)) for record in decompressed]
    close(stream)
    return seqs
end

# println(read_fasta("sequences.fasta.xz"))
seqs = read_fasta("sequences.fasta.xz")
modSeqs = LongDNA{2}.(filter.(x -> contains("ACGT", x), seqs))

println(length(modSeqs))

mismatch_score = 0.5
match_score = 0.0
kmerLength = 30
affine_score = 0.5
matches = [Move(1, 0), Move(3, 0)]
gaps = [Move(1, 1), Move(3, 2)]

#A, B = kmerMatching(modSeqs[1], modSeqs[2], match_score, mismatch_score, matches, gaps, gaps, affine_score, kmerLength)
#kmerChainMatching(modSeqs[1], modSeqs[2], match_score, mismatch_score, matches, gaps, gaps, affine_score, kmerLength)
A, B = general_pairwise_aligner(modSeqs[1], modSeqs[2], make_match_score_matrix(match_score, mismatch_score), matches, gaps, gaps, affine_score) 
println(alignment_score(A, B, make_match_score_matrix(match_score, mismatch_score), matches, gaps, gaps, affine_score))
