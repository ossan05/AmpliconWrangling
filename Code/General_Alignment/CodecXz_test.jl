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

bests = [0,0,0]
sums = [0,0,0]
n = 50
times = [0.0, 0.0, 0.0]
for i in 1 : n
    println("iteration ", i)
    A = rand(modSeqs)
    B = rand(modSeqs)
    time1 = time()
    corr = kmerMatching(A, B, match_score, mismatch_score, matches, gaps, gaps, affine_score, kmerLength)
    time2 = time()
    chain = kmerChainMatching(A, B, match_score, mismatch_score, matches, gaps, gaps, affine_score, kmerLength)
    time3 = time()
    old = kmer_seeded_align(String(A), String(B);wordlength = kmerLength, skip = Int64(kmerLength÷ 3))
    time4 = time()
    times .+= [time2-time1, time3-time2, time4-time3]
    @show times ./ i

    a = [corr, chain, old]
    b = [0.0, 0.0, 0.0]
    for j in eachindex(a)
        x = a[j]
        b[j] = alignment_score(LongDNA{4}(x[1]), LongDNA{4}(x[2]), make_match_score_matrix(match_score, mismatch_score), matches, gaps, gaps, affine_score)
    end
    best = minimum(b)
    for j in 1 : 3
        if b[j] == best
            bests[j] += 1
        end
    end
    sums .+= b
    @show bests
    @show sums./ i
end
println("avrage correlation ", sums[1] ./ n)
println("avrage chain       ", sums[2] ./ n)
println("avrage old         ", sums[3] ./ n)
println("bests correlation ", bests[1])
println("bests chain       ", bests[2])
println("bests old         ", bests[3])
#A, B = kmerMatching(modSeqs[20], modSeqs[19], match_score, mismatch_score, matches, gaps, gaps, affine_score, kmerLength)
#A, B = kmerChainMatching(modSeqs[20], modSeqs[19], match_score, mismatch_score, matches, gaps, gaps, kmerLength)
#A, B = kmer_seeded_align(String(modSeqs[20]), String(modSeqs[19]);wordlength = kmerLength, skip = Int64(kmerLength/3))
println(kmer_matches, alignment_score(A, B, make_match_score_matrix(match_score, mismatch_score), matches, gaps, gaps))
#println(alignment_score(LongDNA{4}.(A), LongDNA{4}.(B), make_match_score_matrix(match_score, mismatch_score), matches, gaps, gaps))