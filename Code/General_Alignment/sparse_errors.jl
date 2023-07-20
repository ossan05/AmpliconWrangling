using BioSequences

A = dna"ACGGTTAGCGCGC"
B = dna"ACCCGTGA"
m = length(A)
using BioSequences

sequence = DNASequence("ATCGATCGATCG")
k = 3

kmers = kmer(sequence, k)
println("All 3-mers: ", kmers)

kmer_counts = count(kmers)
println("K-mer counts: ", kmer_counts)