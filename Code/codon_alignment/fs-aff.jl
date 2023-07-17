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

fs_penalty = 0.5  # 1.4 + 1.2 = 2.6 || 0.5 + 2 * 1.4 = 3.3

ma_penalty = 1
tri_penalty = 2
gap_penalty = 1 + fs_penalty
aff_penalty = 0.4 + fs_penalty
tri_aff = (aff_penalty - fs_penalty) * 3 # affine for tri_penalty (1.2)

function Wunsch(A::String,B::String)
    n, m = (length(A), length(B))
    matrix = zeros(Float64, n+1, m+1)

    matrix[2, 1] = gap_penalty
    matrix[1, 2] = gap_penalty

    for i in 3:n + 1
        matrix[i, 1] = matrix[i-1, 1] + aff_penalty
        if i == 4
            matrix[i, 1] = tri_penalty
        elseif i > 4
            matrix[i, 1] = matrix[i-3, 1] + tri_aff
        end
    end
    
    for i in 3:m + 1
        matrix[1, i] = matrix[1, i-1] + aff_penalty
        if i > 4
            matrix[1, i] = matrix[1, i-3] + tri_aff
        elseif i == 4
            matrix[1, i] = tri_penalty
        end
    end

    for i in 2:n+1
        for j in 2:m+1
            if A[i-1] == B[j-1]
                matrix[i,j] = matrix[i-1, j-1] + fs_penalty
            else
                matrix[i, j] = min(
                    matrix[i-1, j-1] + ma_penalty + fs_penalty,
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

            # affine gap check for know if tri_aff possible
            if i > 6 && any([matrix[i-3, j] == matrix[i-4, j] + gap_penalty,
                matrix[i-3, j] == matrix[i-4, j] + aff_penalty,
                matrix[i-3, j] == matrix[i-6, j] + tri_penalty,
                matrix[i-3, j] == matrix[i-6, j] + tri_aff]) && matrix[i-3, j] + tri_aff < matrix[i, j]

                matrix[i, j] = matrix[i-3, j] + tri_aff
            elseif i > 4 && any([matrix[i-3, j] == matrix[i-4, j] + gap_penalty,
                matrix[i-3, j] == matrix[i-4, j] + aff_penalty]) && matrix[i-3, j] + tri_aff < matrix[i, j]

                matrix[i, j] = matrix[i-3, j] + tri_aff
            end

            if j > 6 && any([matrix[i, j-3] == matrix[i, j-4] + gap_penalty,
                matrix[i, j-3] == matrix[i, j-4] + aff_penalty,
                matrix[i, j-3] == matrix[i, j-6] + tri_penalty,
                matrix[i, j-3] == matrix[i, j-6] + tri_aff]) && matrix[i, j-3] + tri_aff < matrix[i, j]

                matrix[i, j] = matrix[i, j-3] + tri_aff
            elseif j > 4 && any([matrix[i, j-3] == matrix[i, j-4] + gap_penalty,
                matrix[i, j-3] == matrix[i, j-4] + aff_penalty]) && matrix[i, j-3] + tri_aff < matrix[i, j]

                matrix[i, j] = matrix[i, j-3] + tri_aff
            end
            # amount of triplet mismatches going back three steps
            if j > 3 && i > 3
                mismatches = 0
                for l in 1:3
                    if A[i-l] != B[j-l]
                        mismatches += 1
                    end
                end
            end
            
            # above could be stored in a matrix instead
            if j > 3 && i > 3 && matrix[i-3, j-3] + mismatches * ma_penalty <= matrix[i, j]
                matrix[i, j] = matrix[i-3, j-3] + mismatches * ma_penalty
            end

            # kolla ifall det är affine gap för ett 1-gap
            if i > 4
                if any([matrix[i-1, j] == matrix[i-2, j] + gap_penalty,
                    matrix[i-1, j] == matrix[i-2, j] + aff_penalty,
                    matrix[i-1, j] == matrix[i-4, j] + tri_penalty,
                    matrix[i-1, j] == matrix[i-4, j] + tri_aff]) && matrix[i-1, j] + aff_penalty < matrix[i, j]

                    matrix[i, j] = matrix[i-1, j] + aff_penalty
                end
            elseif i > 2 && any([matrix[i-1, j] == matrix[i-2, j] + gap_penalty,
                matrix[i-1, j] == matrix[i-2, j] + aff_penalty]) && matrix[i-1, j] + aff_penalty < matrix[i, j]

                matrix[i, j] = matrix[i-1, j] + aff_penalty
            end

            if j > 4
                if any([matrix[i, j-1] == matrix[i, j-2] + gap_penalty,
                    matrix[i, j-1] == matrix[i, j-2] + aff_penalty,
                    matrix[i, j-1] == matrix[i, j-4] + tri_penalty,
                    matrix[i, j-1] == matrix[i, j-4] + tri_aff]) && matrix[i, j-1] + aff_penalty < matrix[i, j]
                    matrix[i, j] = matrix[i, j-1] + aff_penalty
                end
            elseif j > 2 && any([matrix[i, j-1] == matrix[i, j-2] + gap_penalty,
                matrix[i, j-1] == matrix[i, j-2] + aff_penalty]) && matrix[i, j-1] + aff_penalty < matrix[i, j]

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
        # elseif matrix[i, j] == matrix[i-1, j-1] || matrix[i, j] == matrix[i-1, j-1] + ma_penalty
        #     i -= 1
        #     j -= 1
        elseif i > 6 && any([matrix[i-3, j] == matrix[i-4, j] + gap_penalty,
            matrix[i-3, j] == matrix[i-4, j] + aff_penalty,
            matrix[i-3, j] == matrix[i-6, j] + tri_penalty,
            matrix[i-3, j] == matrix[i-6, j] + tri_aff]) && matrix[i, j] == matrix[i-3, j] + tri_aff 
        
            i-=3
            push!(gap2, j, j, j)

        elseif i > 4 && any([matrix[i-3, j] == matrix[i-4, j] + gap_penalty,
            matrix[i-3, j] == matrix[i-4, j] + aff_penalty]) && matrix[i, j] == matrix[i-1, j] + tri_aff

            i -= 3
            push!(gap2, j, j, j)
        elseif j > 6 && any([matrix[i, j-3] == matrix[i, j-4] + gap_penalty,
            matrix[i, j-3] == matrix[i, j-4] + aff_penalty,
            matrix[i, j-3] == matrix[i, j-6] + tri_penalty,
            matrix[i, j-3] == matrix[i, j-6] + tri_aff]) && matrix[i, j] == matrix[i, j-3] + tri_aff

            j-= 3
            push!(gap1, i, i, i)
        elseif j > 4 && any([matrix[i, j-3] == matrix[i, j-4] + gap_penalty,
            matrix[i, j-3] == matrix[i, j-4] + aff_penalty]) && matrix[i, j] == matrix[i, j-3] + tri_aff

            j -= 3
            push!(gap1, i, i, i)

        elseif i > 4 && any([matrix[i-1, j] == matrix[i-2, j] + gap_penalty,
            matrix[i-1, j] == matrix[i-2, j] + aff_penalty,
            matrix[i-1, j] == matrix[i-4, j] + tri_penalty,
            matrix[i-1, j] == matrix[i-4, j] + tri_aff]) && matrix[i, j] == matrix[i-1, j] + aff_penalty

            i -= 1
            push!(gap2, j)
        elseif i > 2 && any([matrix[i-1, j] == matrix[i-2, j] + gap_penalty,
            matrix[i-1, j] == matrix[i-2, j] + aff_penalty]) && matrix[i, j] == matrix[i-1, j] + aff_penalty

            i -= 1
            push!(gap2, j)
        elseif j > 4 && any([matrix[i, j-1] == matrix[i, j-2] + gap_penalty,
            matrix[i, j-1] == matrix[i, j-2] + aff_penalty,
            matrix[i, j-1] == matrix[i, j-4] + tri_penalty,
            matrix[i, j-1] == matrix[i, j-4] + tri_aff]) && matrix[i, j] == matrix[i, j-1] + aff_penalty

            j -= 1
            push!(gap1, i)
        elseif j > 2 && any([matrix[i, j-1] == matrix[i, j-2] + gap_penalty,
            matrix[i, j-1] == matrix[i, j-2] + aff_penalty]) && matrix[i, j] == matrix[i, j-1] + aff_penalty

            j -= 1
            push!(gap1, i)
        elseif matrix[i, j] == matrix[i-1, j] + gap_penalty 
            i -= 1
            push!(gap2, j)
        elseif matrix[i, j] == matrix[i, j-1] + gap_penalty 
            j -= 1
            push!(gap1, i)
        elseif i > 3 && matrix[i, j] == matrix[i-3, j] + tri_penalty 
            i-=3
            push!(gap2, j, j, j)
        elseif j > 3 && matrix[i,j] == matrix[i,j-3] + tri_penalty
            j-=3
            push!(gap1, i, i, i)
        elseif j > 3 && i > 3 && any([matrix[i, j] == matrix[i-3, j-3],
            matrix[i, j] == matrix[i-3, j-3] + ma_penalty,
            matrix[i, j] == matrix[i-3, j-3] + ma_penalty * 2])
            i -= 3
            j -= 3

            while j > 3 && i > 3 && any([matrix[i, j] == matrix[i-3, j-3],
                matrix[i, j] == matrix[i-3, j-3] + ma_penalty,
                matrix[i, j] == matrix[i-3, j-3] + ma_penalty * 2])
                i -= 3
                j -= 3
            end

            while i > 1 && j > 1 && (matrix[i, j] == matrix[i-1, j-1] + fs_penalty || matrix[i, j] == matrix[i-1, j-1] + ma_penalty + fs_penalty) # && (fs || ma) maybe
                i -= 1
                j -= 1
            end
        else
            i -= 1
            j -= 1

            while j > 3 && i > 3 && any([matrix[i, j] == matrix[i-3, j-3],
                matrix[i, j] == matrix[i-3, j-3] + ma_penalty,
                matrix[i, j] == matrix[i-3, j-3] + ma_penalty * 2])
                i -= 3
                j -= 3
            end

            while i > 1 && j > 1 && (matrix[i, j] == matrix[i-1, j-1] + fs_penalty || matrix[i, j] == matrix[i-1, j-1] + ma_penalty + fs_penalty) # && (fs || ma) maybe
                i -= 1
                j -= 1
            end
        end
    end

    # show(stdout, "text/plain", matrix[end-36:end-34, end-36:end-34])
    # println()

    return insert_blanks(gap1, A), insert_blanks(gap2, B)
end

A, B = generate_seq(16)

# println(A, "\n", B, "\n")

A2, B2 = Wunsch("CCTCATTAAGGCGCCCCATCTTAAGGCTTCAGGC", "CCATATTAAGGTTGGC")

println(A2, "\n", B2)