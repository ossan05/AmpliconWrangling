# using FASTX

ma_penalty = 2
tri_penalty = 1
gap_penalty = 1
fs_penalty = 1

# penalty - aff_penalty = b
aff_penalty = 0.5

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
    matrix = zeros(Float64, n+1, m+1)


    # for i in 2:n + 1
    #     matrix[i, 1] = gap_penalty * (i-1)
    # end
    
    # for i in 2:m + 1
    #     matrix[1, i] = gap_penalty * (i-1)
    # end

    for i in 2:n + 1
        matrix[i, 1] = gap_penalty + aff_penalty * (i-1)
    end
    
    for i in 2:m + 1
        matrix[1, i] = gap_penalty + aff_penalty * (i-1)
    end

    for i in 2:n+1
        for j in 2:m+1
            if A[i-1] == B[j-1]
                matrix[i,j] = matrix[i-1, j-1] + fs_penalty
            else
                matrix[i, j] = min(
                    matrix[i-1, j-1] + ma_penalty,
                    matrix[i, j-1] + gap_penalty + fs_penalty,
                    matrix[i-1, j] + gap_penalty + fs_penalty
                    )
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
            # det kan ändras till att kolla tillbaka som trigap
            if any([matrix[i, j-1] == matrix[i, j-2] + gap_penalty, matrix[i, j-1] == matrix[i, j-2] + aff_penalty]) && matrix[i, j-1] + aff_penalty < matrix[i, j]
                matrix[i, j] = matrix[i, j-1] + aff_penalty 
            elseif any([matrix[i-1, j] == matrix[i-2, j] + gap_penalty, matrix[i-1, j] == matrix[i-2, j] + aff_penalty]) && matrix[i-1, j] + aff_penalty < matrix[i, j]
                matrix[i, j] = matrix[i, j-1] + aff_penalty 
            end
        end
    end

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
        elseif matrix[i, j] == matrix[i-1, j-1] + fs_penalty || matrix[i, j] == matrix[i-1, j-1] + ma_penalty
            i -= 1
            j -= 1
        elseif matrix[i, j] == matrix[i-1, j] + gap_penalty + fs_penalty
            i -= 1
            push!(gap2, j)
        elseif matrix[i,j] == matrix[i,j - 1] + gap_penalty + fs_penalty
            j -= 1
            push!(gap1, i)
        elseif i > 3 && matrix[i,j] == matrix[i-3,j] + tri_penalty
            i-=3
            push!(gap2, j, j, j)
        elseif j > 3 && matrix[i,j] == matrix[i,j-3] + tri_penalty
            j-=3
            push!(gap1, i, i, i)
        else
            i -= 3
            j -= 3
        end
    end

    println(insert_blanks(gap1, A), "\n" ,insert_blanks(gap2, B))
    # return matrix[n + 1, m + 1], sort!(gap1), sort!(gap2)
    # return insert_blanks(gap1, A) ,insert_blanks(gap2, B)
end

# A, B = generate_seq(80)

Wunsch("CGAGTGTCCCATAAAGAGATCTAGGCTGCTATGATGGATATGACTGCTGGCGCGGACTGTTAATCTTATTCAAGTGAATT", "CGAGTGTCCCATAAAGAGTTTAGGCTGCTATGATGGATATGACTGATGGAGCGGTCTGTTTATCTTATCGTCGAGTGTTCTATAAAGAAATCTACGCTGCTATTCTGGCTATGACCGCTGCAAGCGCGGACTGTTAATCTTATTCAAGTGAATT")