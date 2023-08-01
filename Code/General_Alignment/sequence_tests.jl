include("sequence_generator.jl")
include("sparse_errors.jl")
using BenchmarkTools, NextGenSeqUtils
len = 3000
upper = 10

seq = Vector{Any}()
for i in 1:upper
    A, B = generate_seq(3000, 0.001i, 0.002i, 0.00005i, 0.0001i, 4)
    println(length(A))
    b = @benchmark initiate_kmerMatching($A, $B)
    push!(seq, median(b))
end

println(seq)