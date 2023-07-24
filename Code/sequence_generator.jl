using BioSequences

function generate_seq(seq_len::Int)
    seq = dna""
    dna = dna"ACGT"

    for i in 1:seq_len
        push!(seq, dna[rand(1:4)])
    end

    # mutations

    mutated_seq = seq

    n = seq_len

    while n < seq_len * 1.2
        dein = rand(1:n)
        
        if dein <= n - 2
            if rand(1:2) == 1
                if rand(1:100) == 21
                    other = mutated_seq[dein + 1:end]
                    mutated_seq = push!(mutated_seq[1:dein], dna[rand(1:4)])
                    append!(mutated_seq, other)
                else
                    other = mutated_seq[dein + 1:end]
                    mutated_seq = push!(mutated_seq[1:dein], dna[rand(1:4)], dna[rand(1:4)], dna[rand(1:4)])
                    append!(mutated_seq, other)
                end
            elseif n > seq_len
                if rand(1:100) == 21
                    mutated_seq = deleteat!(mutated_seq, dein)
                else
                    for i in 1:3
                        mutated_seq = deleteat!(mutated_seq, dein)
                    end
                end
            end
        elseif dein == n - 1
            if rand(1:2) == 1
                if rand(1:100) == 21
                    mutated_seq = push!(utated_seq[1:dein], dna[rand(1:4)], mutated_seq[end])
                else
                    mutated_seq = push!(mutated_seq[1:dein], dna[rand(1:4)], dna[rand(1:4)], dna[rand(1:4)], mutated_seq[end])
                end
            elseif n > seq_len
                if rand(1:100) == 21
                    deleteat!(mutated_seq, dein)
                else
                    for i in 1:3
                        deleteat!(mutated_seq, dein - 1)
                    end
                end
            end
        else
            if rand(1:2) == 1
                if rand(1:100) == 21
                    mutated_seq = push!(mutated_seq[1:dein], dna[rand(1:4)])
                else
                    mutated_seq = push!(mutated_seq[1:dein], dna[rand(1:4)], dna[rand(1:4)], dna[rand(1:4)])
                end
            elseif n > seq_len
                if rand(1:100) == 21
                    deleteat!(mutated_seq, dein)
                else
                    for i in 1:3
                        deleteat!(mutated_seq, dein - 2)
                    end
                end
            end
        end
        n = length(mutated_seq)
    end

    # sequencing errors

    for i in 1:n
        if rand(1:10) == 7
            mutated_seq[i] = dna[rand(findall(x->x != i, dna))]
        end
    end

    return seq, mutated_seq
end

A, B = generate_seq(7)
println(A, "\n", B)