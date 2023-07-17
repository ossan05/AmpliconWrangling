# using FASTX

ma_penalty = 2

tri_penalty_a = 1
tri_penalty_b = 0.5

gap_penalty_a = 1
gap_penalty_b = 0.5


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

function insert_blanks(pos::Vector{Int}, original::String)
    spos = sort(pos)
    push!(spos, typemax(Int64))
    l = length(original)
    l2 = length(pos)
    j = 1
    result=""
    
    for i in 1:l
        while i == spos[j]
            j+=1
            result *= '_'
        end
        result *= original[i]
    end
    while j <= l2 && l+1 == spos[j]
        result *= '_'
        j+=1
    end
    return result
end

function Wunsch(A::String,B::String)
    n, m = (length(A), length(B))
    score_matrix = zeros(Float64, n+1, m+1)
    move_matrix = zeros(Int, n+1, m+1)

    a, b = (1, 0.5)
    gap_penalty = 1
    tri_bonus = 0.5


    for i in 2:n + 1
        move_matrix[i, 1] = i - 1
        score_matrix[i, 1] = i - 1
    end
    
    for i in 2:m + 1
        move_matrix[1, i] = i - 1
        score_matrix[1, i] = i - 1
    end

    # for i in 2:n + 1
    #     if (i - 1) % 3 == 0
    #         score_matrix[i, 1] = gap_penalty + b * (i - 2) - (i - 1) / 3 * tri_bonus
    #     else
    #         score_matrix[i, 1] = gap_penalty + b * (i - 2)
    #     end
    # end
    
    # for i in 2:m + 1
    #     if (i - 1) % 3 == 0
    #         score_matrix[1, i] = gap_penalty + b * (i - 2) - (i - 1) / 3 * tri_bonus
    #     else
    #         score_matrix[1, i] = gap_penalty + b * (i - 2)
    #     end
    # end

    for i in 2:n+1
        for j in 2:m+1
            if A[i-1] == B[j-1] 
                score_matrix[i,j] = score_matrix[i-1, j-1]
                move_matrix[i, j] = move_matrix[i-1, j-1] + 1
            elseif score_matrix[i-1, j-1] + ma_penalty <= min(score_matrix[i-1, j], score_matrix[i, j-1]) + gap_penalty
                score_matrix[i,j] = score_matrix[i-1, j-1] + ma_penalty
                move_matrix[i, j] = move_matrix[i-1, j-1] + 1
            elseif score_matrix[i-1, j] < score_matrix[i, j-1]
                score_matrix[i, j] = score_matrix[i-1, j] + gap_penalty
                move_matrix[i, j] = move_matrix[i-1, j] + 1
            else
                score_matrix[i, j] = score_matrix[i, j-1] + gap_penalty
                move_matrix[i, j] = move_matrix[i, j-1] + 1
            end

            if i > 3 && move_matrix[i-3, j] % 3 == 0 && score_matrix[i-3, j] + tri_penalty_a < score_matrix[i, j]
                score_matrix[i, j] = score_matrix[i-3, j] + tri_penalty_a
                move_matrix[i, j] = move_matrix[i-3, j] + 3
            end
            if j > 3 && move_matrix[i, j-3] % 3 == 0 && score_matrix[i, j-3] + tri_penalty_a < score_matrix[i, j]
                score_matrix[i, j] = score_matrix[i, j-3] + tri_penalty_a
                move_matrix[i, j] = move_matrix[i, j-3] + 3
            end
        end
    end

    # display(move_matrix)

    # display(score_matrix[5:end, 3:end])
    show(stdout, "text/plain", score_matrix)
    show(stdout, "text/plain", move_matrix)

    i, j = (n + 1, m + 1)

    gap1 = Vector{Int}()
    gap2 = Vector{Int}()

    while i > 1 || j > 1
        if i == 1
            j -= 1
            push!(gap1, i)
        elseif j == 1
            i -= 1
            push!(gap2, j)
        elseif score_matrix[i, j] == score_matrix[i-1, j] + gap_penalty_a 
            i -= 1
            push!(gap2, j)
        elseif score_matrix[i,j] == score_matrix[i,j - 1] + gap_penalty_a 
            j -= 1
            push!(gap1, i)
        elseif i > 3 && score_matrix[i,j] == score_matrix[i-3,j] + tri_penalty_a
            i-=3
            push!(gap2, j, j, j)
        elseif j > 3 && score_matrix[i,j] == score_matrix[i,j-3] + tri_penalty_a
            j-=3
            push!(gap1, i, i, i)
        else
            i -= 1
            j -= 1
        end
    end

    println(insert_blanks(gap1, A), "\n" ,insert_blanks(gap2, B))
    # return score_matrix[n + 1, m + 1], sort!(gap1), sort!(gap2)
    # return insert_blanks(gap1, A) ,insert_blanks(gap2, B)
end

A, B = generate_seq(10)

Wunsch("CAGCTGCGAG", "CAGCGAG")