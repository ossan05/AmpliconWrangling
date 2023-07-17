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
    gap_penalty = 1
    ma_penalty = 2
    tri_penalty = 1

    for i in 2:n+1
        matrix[i, 1] = i - 1
    end
    
    for i in 2:m+1
        matrix[1, i] = i - 1
    end

    for i in 2:n+1
        for j in 2:m+1
            if A[i-1] == B[j-1]
                matrix[i,j] = matrix[i-1, j-1]
            else
                matrix[i, j] = min(
                    matrix[i-1, j-1] + ma_penalty,
                    matrix[i, j-1]+gap_penalty,
                    matrix[i-1, j]+gap_penalty
                    )
                if i > 3 && matrix[i-3, j]+tri_penalty < matrix[i,j]
                    matrix[i,j] = matrix[i-3, j]+tri_penalty
                end
                if j > 3 && matrix[i, j-3]+tri_penalty < matrix[i,j]
                    matrix[i,j] = matrix[i, j-3]+tri_penalty
                end
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
        elseif matrix[i, j] == matrix[i-1, j] + gap_penalty
            i -= 1
            push!(gap2, j)
        elseif matrix[i,j] == matrix[i,j - 1] + gap_penalty
            j -= 1
            push!(gap1, i)
        elseif i > 3 && matrix[i,j] == matrix[i-3,j] + tri_penalty
            i-=3
            push!(gap2, j, j, j)
        elseif j > 3 && matrix[i,j] == matrix[i,j-3] + tri_penalty
            j-=3
            push!(gap1, i, i, i)
        else
            i -= 1
            j -= 1
        end
    end
    return (insert_blanks(gap1, A), insert_blanks(gap2, B))
end

function Hirschberg(seq1, seq2)
    m = length(seq1)
    n = length(seq2)
    if m < 5 || n < 5
        return Wunsch(seq1, seq2)
    end
    toMid_sb = lastrow_scoreboard(seq1[1 : m ÷ 2], seq2)
    fromMid_sb = reverse(lastrow_scoreboard(reverse(seq1[m ÷ 2 : m]), seq2))
    display(toMid_sb .+ fromMid_sb)
    mid = findmin(toMid_sb .+ fromMid_sb)[2]
    println(mid)
    return Hirschberg(seq1[1 : m ÷ 2], seq2[1 : mid]) .* Hirschberg(seq1[m ÷ 2 : m], seq2[mid : n])
end

function lastrow_scoreboard(seq1, seq2, )
    m = length(seq1)
    n = length(seq2)
    rows = fill(typemax(Int) ÷ 2, (3, n + 1))
    cyc = [3, 1, 2] #[2, 3, 1] [3, 1, 2]
    rows[1,1] = 0
    for i in 2:n+1
        rows[1, i] = rows[1, i-1] + gap_penalty
    end
    for i in 2:m+1
        rows[cyc[1], 1] = rows[cyc[2], 1] + gap_penalty 
        for j in 2:n+1
            rows[cyc[1], j] = min(
                rows[cyc[2], j-1] + ma_penalty,
                rows[cyc[1], j-1] + gap_penalty,
                rows[cyc[2], j] + gap_penalty,
                rows[cyc[1], j] + tri_penalty
            )
            if j > 3
                rows[cyc[1], j] = min(rows[cyc[1], j], rows[cyc[1], j-3] + tri_penalty)
            end
            if seq1[i-1] == seq2[j-1]
                rows[cyc[1], j] = min(rows[cyc[1], j], rows[cyc[2], j-1])
            end
        end
        tmp = cyc[3]
        cyc[3]=cyc[2]
        cyc[2]=cyc[1]
        cyc[1]= tmp
    end
    return rows[cyc[2], : ]
end
gap_penalty = 1
ma_penalty = 2
tri_penalty = 1

println(Hirschberg("AATCGGT", "CGTTATT"))