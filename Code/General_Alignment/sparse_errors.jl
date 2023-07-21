using BioSequences
using Kmers

A = LongDNA{2}("ACGGTTAGCGCGCAAGGTCGATGTGTGTGTGTG")
B = LongDNA{2}("TCGGTTACGCGCAAGGTCGATGAGTGTGTGTG")
@show A[2,7].data[0]
kmerDict = Dict{Int64, Array{Int64}}()

m = length(A)
n = length(B)

k = min(m-2, n-2, log2(m*n), 30)
if k <= 1
    k = 2
end

for i in 1:m-k+1
    kmer=A[i:i+k-1].data[0]
    if haskey(kmerDict, kmer)
        push(kmerDict[kmer], i)
    else
        kmerDict[kmer] = [i]
    end
end
@show kmerDict
