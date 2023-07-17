function generate_seq(seq_len::Int)
    dna = "ACGT"
    seq = ""

    for i in 1:seq_len
        seq *= rand(dna)
    end

    # mutations

    seq_copy1 = ""
    dein = seq_len
    while length(seq_copy1) + seq_len - dein < seq_len * 1.2
        dein = rand(1:seq_len)
        if rand(1:2) == 1
            if rand(1:100) == 21
                seq_copy1 *= seq[1:dein] * rand(dna)
            else
                seq_copy1 *= seq[1:dein] * rand(dna) * rand(dna) * rand(dna) 
            end
        else
            if rand(1:100) == 21
                seq_copy1 *= seq[1:dein] * rand(dna) 
            else
                seq_copy1 *= seq[1:dein] * rand(dna) * rand(dna) * rand(dna)
            end
        end
    end

    seq_copy1 *= seq[dein + 1:end]

    # sequencing errors

    seq_copy2 = collect(seq_copy1)
    for i in 1:length(seq_copy2)
        if rand(1:10) == 7
            seq_copy2[i] = dna[rand(findall(x->x != seq_copy2[i], dna))]
        end
    end

    return seq, join(seq_copy2, "")
end