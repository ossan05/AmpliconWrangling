using BioSequences

function generate_seq(seq_len::Int)
    seq = randseq(DNAAlphabet{2}(), SamplerUniform(dna"ACGT"), seq_len)
    dna = LongDNA{2}("ACGT")
   
    # mutations

    mutated_seq = LongDNA{2}()

    for i in seq
        push!(mutated_seq, i)
    end

    n = seq_len

    x = ceil(seq_len / 8)
    long_indel_frequency  = 10

    while x > 0
        indel_pos = rand(1:n)
        if rand(1:long_indel_frequency) == 1 && x ÷ 7 > 1
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
    copy_chance = rand(1:1) == 3
    if copy_chance
        copy_len = n ÷ rand(10:20)
        copy_ind = rand(1:n - copy_len)
        println(copy_ind)
        @show mutated_seq

        other = mutated_seq[copy_ind:end]
        mutated_seq = append!(mutated_seq[1:copy_ind - 1], mutated_seq[copy_ind:copy_ind + copy_len])
        append!(mutated_seq, other)
        @show mutated_seq
    end

    # sequencing errors

    for i in 1:n
        if rand(1:100) == 7
            mutated_seq[i] = dna[rand(findall(x->x != mutated_seq[i], dna))]
        end
    end

    return seq, mutated_seq
end

A, B = generate_seq(50)
println(A, "\n", B)
