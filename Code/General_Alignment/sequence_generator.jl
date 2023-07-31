using BioSequences
using Distributions

# function check_divergence(A::LongSequence{2}, B::LongSequence{2})
#     divergent_bases = max(length(A), length(B)) - min(length(A), length(B))

#     for i in min(length(A), length(B))
#         if A[i] != B[i]
#             divergent_bases += 1
#         end
#     end

#     return divergent_bases
# end


function generate_seq(seq_len::Int64, long_indel_frequency::Float64, indel_frequency::Float64, substitution_frequency::Float64, frameshift_frequency::Float64, long_length::Int64)
    seq = randseq(DNAAlphabet{2}(), SamplerUniform(dna"ACGT"), seq_len)
    dna = LongDNA{2}("ACGT")
   
    # mutations

    mutated_seq = LongDNA{2}()

    for i in seq
        push!(mutated_seq, i)
    end

    n = seq_len

    indels = Poisson(n * indel_frequency)
    frameshifts = indels * frameshift_frequency
    long_indels  = Poisson(indel_count * long_indel_frequency)
    substitutions  = Poisson(n * substitution_frequency)

    Geometric(1/long_length)

    for i in 1:substitutions
        mutated_seq[rand(1:n)] = dna[rand(findall(x->x != mutated_seq[i], dna))]
    end

    for i in 1:indels
        mutated_seq[rand]

    while x > 0
        indel_pos = rand(1:n)
        if rand(Float64) <long_indel_frequencyg && x ÷ 7 > 1
            indel_len = rand(1:x ÷ 7)
        else
            indel_len = 1
        end

        if rand(1:2) == 1
            if rand(1:100) == 21
                other = mutated_seq[indel_pos + 1:end]
                mutated_seq = append!(mutated_seq[1:indel_pos], randseq(DNAAlphabet{2}(), SamplerUniform(dna"ACGT"), Int64(indel_len)))
                append!(mutated_seq, other)
                x -= indel_len
            else
                other = mutated_seq[indel_pos + 1:end]
                indel = randseq(DNAAlphabet{2}(), SamplerUniform(dna"ACGT"), Int64(indel_len * 3))
                mutated_seq = append!(mutated_seq[1:indel_pos], indel)
                append!(mutated_seq, other)
                x -= indel_len * 3
            end
        elseif n > indel_len * 3
            if rand(1:100) == 21
                if indel_len + indel_pos - 1 >= n
                    deleteat!(mutated_seq, Int64(n - indel_len):n)
                elseif indel_len % 2 != 0
                    deleteat!(mutated_seq, Int64(indel_pos - indel_len ÷ 2):Int64(indel_pos + indel_len ÷ 2))   
                else
                    deleteat!(mutated_seq, Int64(indel_pos - indel_len / 2):Int64(indel_pos + indel_len / 2 - 1))
                end
                x -= indel_len
            else
                if indel_pos - 1 + indel_len * 3 >= n
                    deleteat!(mutated_seq, Int64(n - indel_len * 3):n)
                elseif indel_pos - (indel_len * 3) ÷ 2 == 0
                    deleteat!(mutated_seq, 1:indel_len * 3)
                elseif indel_len * 3 % 2 != 0
                    deleteat!(mutated_seq, Int64(indel_pos - (indel_len * 3) ÷ 2):Int64(indel_pos + (indel_len * 3) ÷ 2))
                else
                    deleteat!(mutated_seq, Int64(indel_pos - (indel_len * 3) / 2):Int64(indel_pos + (indel_len * 3) / 2 - 1))
                end
                x -= indel_len * 3
            end
        end
        n = length(mutated_seq)
    end

    # copy
    copy_chance = rand(1:3) == 3
    if copy_chance
        copy_len = n ÷ rand(10:20)
        copy_ind = rand(1:n - copy_len)

        other = mutated_seq[copy_ind:end]
        mutated_seq = append!(mutated_seq[1:copy_ind - 1], mutated_seq[copy_ind:copy_ind + copy_len])
    end

    # sequencing errors

    for i in 1:n
        if rand(1:100) == 7
            mutated_seq[i] = dna[rand(findall(x->x != mutated_seq[i], dna))]
        end
    end

    return seq, mutated_seq
end