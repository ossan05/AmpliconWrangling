# using FASTX

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

ma_penalty = 1
tri_penalty = 2
gap_penalty = 2
two_penalty = 2

function Wunsch(A::String,B::String)
    n, m = (length(A), length(B))
    matrix = zeros(Float64, n+1, m+1)

    for i in 2:n + 1
        matrix[i, 1] = matrix[i-1, 1] + gap_penalty
        if i > 2 && matrix[i-2, 1] + two_penalty < matrix[i, 1]
            matrix[i, 1] = matrix[i-2, 1] + two_penalty
        end
        if i > 3 && matrix[i-3, 1] + tri_penalty < matrix[i, 1]
            matrix[i, 1] = matrix[i-3, 1] + tri_penalty
        end
    end

    for i in 2:m + 1
        matrix[1, i] = matrix[1, i-1] + gap_penalty
        if i > 2 && matrix[1, i-2] + two_penalty < matrix[1, i]
            matrix[1, i] = matrix[1, i-2] + two_penalty
        end
        if i > 3 && matrix[1, i-3] + tri_penalty < matrix[1, i]
            matrix[1, i] = matrix[1, i-3] + tri_penalty
        end
    end

    for i in 2:n+1
        for j in 2:m+1
            if A[i-1] == B[j-1]
                matrix[i,j] = matrix[i-1, j-1]
            else
                matrix[i, j] = min(
                    matrix[i-1, j-1] + ma_penalty,
                    matrix[i, j-1] + gap_penalty,
                    matrix[i-1, j] + gap_penalty
                    )
            end
            # if j > 3
            #     if i > 3
            #         if A[i-3:i-1] == B[j-3:j-1]

            if i > 2 && matrix[i-2, j] + two_penalty < matrix[i, j]
                matrix[i, j] = matrix[i-2, j] + two_penalty
            end
            if j > 2 && matrix[i, j-2] + two_penalty < matrix[i, j]
                matrix[i, j] = matrix[i, j-2] + two_penalty
            end

            if i > 3 && matrix[i-3, j] + tri_penalty < matrix[i, j]
                matrix[i, j] = matrix[i-3, j] + tri_penalty
            end
            if j > 3 && matrix[i, j-3] + tri_penalty < matrix[i, j]
                matrix[i, j] = matrix[i, j-3] + tri_penalty
            end
            if j > 3 && i > 3 && A[i-3:i-1] == B[j-3:j-1]
                matrix[i, j] = matrix[i-3, j-3]
            end
        end
    end

    display(matrix)

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
        # elseif matrix[i, j] == matrix[i-1, j-1] || matrix[i, j] == matrix[i-1, j-1] + ma_penalty
        #     i -= 1
        #     j -= 1
        elseif matrix[i, j] == matrix[i-1, j] + gap_penalty
            i -= 1
            push!(gap2, j)
        elseif matrix[i,j] == matrix[i,j - 1] + gap_penalty
            j -= 1
            push!(gap1, i)
        elseif i > 2 && matrix[i, j] == matrix[i-2, j] + two_penalty
            i -= 2
            push!(gap2, j, j)
        elseif j > 2 && matrix[i, j] == matrix[i, j-2] + two_penalty
            j -= 2
            push!(gap1, i, i)
        elseif i > 3 && matrix[i,j] == matrix[i-3,j] + tri_penalty
            i-=3
            push!(gap2, j, j, j)
        elseif j > 3 && matrix[i,j] == matrix[i,j-3] + tri_penalty
            j-=3
            push!(gap1, i, i, i)
        elseif j > 3 && i > 3 && any([matrix[i, j] == matrix[i-3, j-3], matrix[i, j] == matrix[i-3, j-3] + ma_penalty, matrix[i, j] == matrix[i-3, j-3] + ma_penalty * 2])
            i -= 3
            j -= 3
        else
            i -= 1
            j -= 1
        end
    end

    # show(stdout, "text/plain", matrix[end-36:end-34, end-36:end-34])
    # println()

    #println(insert_blanks(gap1, A), "\n" ,insert_blanks(gap2, B))
    # return matrix[n + 1, m + 1], sort!(gap1), sort!(gap2)
    return insert_blanks(gap1, A), insert_blanks(gap2, B)
end

A, B = generate_seq(30)

println(A, "\n", B, "\n")

A2, B2 = Wunsch(A, B)

println(A2, "\n", B2)