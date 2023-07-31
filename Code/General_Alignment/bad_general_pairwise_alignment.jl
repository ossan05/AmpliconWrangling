using BioSequences

struct Move
    step::Int64
    score::Float64
end

function initiate_general_pairwise_aligner(seq1::LongDNA{2}, seq2::LongDNA{2})
    general_pairwise_aligner(seq1, seq2, .0, 0.5, [((1, 1), .0), ((1, 0), 1.0), ((0, 1), 1.0), ((3, 3), 0.0), ((3, 0), 2.0), ((0, 3), 2.0)], 0.5)
end
function toInt(x::NucleicAcid)
    trailing_zeros(reinterpret(UInt8,x))+1
end
toInt(DNA_A)

function make_match_score_matrix(match_score, mismatch_score)
    m = zeros(4, 4)
    for i in 1:4, j in 1:4
        if i == j
            m[i, j] = match_score
        else
            m[i, j] = mismatch_score
        end
    end
    return m
end

# match and mismatch matrix
function general_pairwise_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, moves::Vector{Tuple{Tuple{Int64, Int64}, Float64}}) 
    general_pairwise_aligner(A, B, make_match_score_matrix(match_score, mismatch_score), moves) 
end

# general_pairwise_aligner makes a match/mismatch matrix and takes the tuple with the moves and scores, 
# making them objects of type Move or Move and puting them into lists 
function general_pairwise_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, moves::Vector{Tuple{Tuple{Int64, Int64}, Float64}})
    # sort the moves into gaps
    # moves are all the moves
    # v - vertical, h - horizontal
    matches = Vector{Move}()
    vgaps = Vector{Move}()
    hgaps = Vector{Move}()

    # fill the vectors with gaps and moves
    biggest_match = -Inf
    biggest_vgap = -Inf
    biggest_hgap = -Inf
    for i ∈ moves
        if i[1][2] == i[1][1]
            push!(matches, Move(i[1][1], i[2]))
            if i[1][2] > biggest_match
                biggest_match = i[1][2]
            end
        elseif i[1][2] == 0
            push!(vgaps, Move(i[1][1], i[2]))
            if i[1][1] > biggest_vgap
                biggest_vgap = i[1][1]
            end
        else
            push!(hgaps, Move(i[1][2], i[2]))
            if i[1][2] > biggest_hgap
                biggest_hgap = i[1][2]
            end
        end
    end
    
    # make the alignment

    # Needleman-Wunsch alignment

    n, m = length(A), length(B)
    dp_matrix = fill(typemax(Float64),
                    n + 1 + max(biggest_match, biggest_vgap),
                    m + 1 + max(biggest_match, biggest_hgap))
    dp_matrix[1, 1] = .0

    # first row
    for i ∈ 2:m + 1, Move ∈ hgaps 
        if dp_matrix[1, i - Move.step] + Move.score < dp_matrix[1, i]
            dp_matrix[1, i] = dp_matrix[1, i - Move.step] + Move.score 
        end
    end

    # first column
    for i ∈ 2:n + 1, Move ∈ vgaps 
        if dp_matrix[i - Move.step, 1] + Move.score < dp_matrix[i, 1]
            dp_matrix[i, 1] = dp_matrix[i - Move.step, 1] + Move.score
        end
    end

    # whole dp matrix
    for i ∈ 2:n + 1
        for j ∈ 2:m + 1

            # finding the lowest score move

            # the lowest score of vertical gaps is asigned to the current position in the matrix
            for k ∈ vgaps
                if dp_matrix[i - k.step, j] + k.score < dp_matrix[i, j]
                    dp_matrix[i, j] = dp_matrix[i - k.step, j] + k.score
                end
            end

            # the lowest score of horizontal gaps is asigned to the current position in the matrix if the score is lower than the current one
            for k ∈ hgaps
                if dp_matrix[i, j - k.step] + k.score < dp_matrix[i, j]
                    dp_matrix[i, j] = dp_matrix[i, j - k.step] + k.score
                end
            end

            for k ∈ matches
                # calculate mismatch score
                s = 0
                for l ∈ i - k.step:i - 1
                    s += match_score_matrix[toInt(A[l]), toInt(B[l + j - i])]
                end

                # asign a new value to the matrix if the move score is lower
                if dp_matrix[i - k.step, j - k.step] + k.score + s < dp_matrix[i, j]
                    dp_matrix[i, j] = dp_matrix[i - k.step, j - k.step] + k.score + s 
                end
            end
        end
    end

    y = m + 1
    x = n + 1
    res_A = dna""
    res_B = dna""

    while x > 1 || y > 1
        if x == 1
            push!(res_A, DNA_Gap)
            push!(res_B, B[y - 1])
            y -= 1
        elseif y == 1
            push!(res_A, A[x - 1])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # iterate through vertical Move moves
            for k ∈ vgaps
                # check that the move isn't out of bounds
                if k.step < x
                    # opening Move
                    # check if the move lead to the current cell
                    if dp_matrix[x, y] == dp_matrix[x - k.step, y] + k.score
                        for i ∈ 1:k.step
                            push!(res_A, A[x - i])
                            push!(res_B, DNA_Gap)
                        end
                        x -= k.step
                        # break to stop iterating through moves
                        break
                    end
                end
            end

            # iterate through horizontal Move moves
            for k ∈ hgaps
                # check that the move isn't out of bounds
                if k.step < y
                    # opening Move
                    # check if the move lead to the current cell
                    if dp_matrix[x, y] == dp_matrix[x, y - k.step] + k.score
                        for i ∈ 1:k.step
                            push!(res_A, DNA_Gap)
                            push!(res_B, B[y - i])
                        end
                        y -= k.step
                        # break to stop iterating through moves
                        break
                    end
                end
            end

            # iterate through digonal match moves
            for k ∈ matches
                
                # check that the move isn't out of bounds
                if k.step < x && k.step < y

                    # calculate score of mismatches
                    s = 0
                    for i ∈ x - k.step:x - 1
                        s += match_score_matrix[toInt(A[i]), toInt(B[i + y - x])]
                    end

                    # check if the move lead to the current cell
                    if dp_matrix[x, y] == dp_matrix[x - k.step, y - k.step] + k.score + s
                        
                        # write the resulting sequences
                        # k.step and k.step is the same
                        for i ∈ 1:k.step
                            push!(res_A, A[x - i])
                            push!(res_B, B[y - i])
                        end
                        x -= k.step
                        y -= k.step

                        # break to stop iterating through moves
                        break
                    end
                end
            end
        end
    end

    return reverse(res_A), reverse(res_B)
end

# Another method for affine Move penalties
function general_pairwise_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, moves::Vector{Tuple{Tuple{Int64, Int64}, Float64}}, affine_gap::Float64) 
    
    match_score_matrix = zeros(4, 4)

    # this is when every mismatch between nucletides have the same penalty 
    for i ∈ 1:4, j ∈ 1:4
        if i == j
            match_score_matrix[i, j] = match_score
        else
            match_score_matrix[i, j] = mismatch_score
        end
    end
    general_pairwise_aligner(A, B, match_score_matrix, moves, affine_gap) 
end
function general_pairwise_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, moves::Vector{Tuple{Tuple{Int64, Int64}, Float64}}, affine_gap::Float64) 

    # sort the moves into gaps 
    matches = Vector{Move}()
    vgaps = Vector{Move}()
    hgaps = Vector{Move}()

    biggest_match = -Inf
    biggest_vgap = -Inf
    biggest_hgap = -Inf
    for i ∈ moves
        if i[1][2] == i[1][1]
            push!(matches, Move(i[1][1], i[2]))
            if i[1][2] > biggest_match
                biggest_match = i[1][2]
            end
        elseif i[1][2] == 0
            push!(vgaps, Move(i[1][1], i[2]))
            if i[1][1] > biggest_vgap
                biggest_vgap = i[1][1]
            end
        else
            push!(hgaps, Move(i[1][2], i[2]))
            if i[1][2] > biggest_hgap
                biggest_hgap = i[1][2]
            end
        end
    end

    # Needleman-Wunsch alignment

    n, m = length(A), length(B)

    # Elongated sequence and matrices for less comparisons
    A = LongDNA{2}("A")^(biggest_match - 1) * A
    B = LongDNA{2}("A")^(biggest_match - 1) * B

    vertical_offset = max(biggest_match, biggest_vgap) + 1
    horizontal_offset = max(biggest_match, biggest_hgap) + 1

    dp_matrix = fill(typemax(Float64),
                     n + vertical_offset,
                     m + horizontal_offset)

    dp_matrix[vertical_offset, horizontal_offset] = .0


    vaffine_matrix = fill(typemax(Float64),
                         n + vertical_offset,
                         m + horizontal_offset)

    haffine_matrix = fill(typemax(Float64),
                         n + vertical_offset,
                         m + horizontal_offset)

    vboundary = n + vertical_offset
    hboundary = m + horizontal_offset

    # first row
    for i ∈ 1 + horizontal_offset:hboundary
        for Move ∈ hgaps
            if dp_matrix[vertical_offset, i - Move.step] + affine_gap * Move.step < dp_matrix[vertical_offset, i]
                dp_matrix[vertical_offset, i] = dp_matrix[vertical_offset, i - Move.step] + affine_gap * Move.step
                haffine_matrix[vertical_offset, i] = dp_matrix[vertical_offset, i]

            elseif dp_matrix[vertical_offset, i - Move.step] + Move.score < dp_matrix[vertical_offset, i]
                dp_matrix[vertical_offset, i] = dp_matrix[vertical_offset, i - Move.step] + Move.score
                haffine_matrix[vertical_offset, i] = dp_matrix[vertical_offset, i]
            end
        end
    end

    # first column
    for i ∈ 1 + vertical_offset:vboundary
        for Move ∈ vgaps
            if dp_matrix[i - Move.step, horizontal_offset] + affine_gap * Move.step < dp_matrix[i, horizontal_offset]
                dp_matrix[i, horizontal_offset] = dp_matrix[i - Move.step, horizontal_offset] + affine_gap * Move.step
                vaffine_matrix[i, horizontal_offset] = dp_matrix[i, horizontal_offset]

            elseif dp_matrix[i - Move.step, horizontal_offset] + Move.score < dp_matrix[i, horizontal_offset]
                dp_matrix[i, horizontal_offset] = dp_matrix[i - Move.step, horizontal_offset] + Move.score
            end
        end
    end

    for i ∈ 1 + vertical_offset:vboundary
        for j ∈ 1 + horizontal_offset:hboundary

            # finds the best vertical move

            for k ∈ vgaps
                # affine Move
                if vaffine_matrix[i - k.step, j] + affine_gap * k.step < vaffine_matrix[i, j]
                    vaffine_matrix[i, j] = vaffine_matrix[i - k.step, j] + affine_gap * k.step
                end
                # opening Move
                if dp_matrix[i - k.step, j] + k.score < vaffine_matrix[i, j]
                    vaffine_matrix[i, j] = dp_matrix[i - k.step, j] + k.score
                end
            end

            # finds the best horizontal move

           
            for k ∈ hgaps
                # affine_gap
                if haffine_matrix[i, j - k.step] + affine_gap * k.step < haffine_matrix[i, j]
                    haffine_matrix[i, j] = haffine_matrix[i, j - k.step] + affine_gap * k.step
                end
                # opening Move
                if dp_matrix[i, j - k.step] + k.score < haffine_matrix[i, j]
                    haffine_matrix[i, j] = dp_matrix[i, j - k.step] + k.score
                end
            end
            # the fastest diagonal move

            for k ∈ matches
                s = 0
                iB = j - horizontal_offset + biggest_match - k.step

                for iA ∈ i - vertical_offset + biggest_match - k.step:i - vertical_offset - 1 + biggest_match
                    s += match_score_matrix[toInt(A[iA]), toInt(B[iB])]
                    iB += 1
                end 

                if dp_matrix[i - k.step, j - k.step] + k.score + s < dp_matrix[i, j]
                    dp_matrix[i, j] = dp_matrix[i - k.step, j - k.step] + k.score + s 
                end
            end

            # fastest move to dp_matrix
            dp_matrix[i, j] = min(dp_matrix[i, j], haffine_matrix[i, j], vaffine_matrix[i, j])
        end
    end

    # backtracking
    y = m + horizontal_offset
    x = n + vertical_offset
    res_A = LongDNA{4}("")
    res_B = LongDNA{4}("")


    while x > vertical_offset || y > horizontal_offset
        if x == vertical_offset
            push!(res_A, DNA_Gap)
            push!(res_B, B[y - horizontal_offset - 1 + biggest_match])
            y -= 1
        elseif y == horizontal_offset
            push!(res_A, A[x - vertical_offset - 1 + biggest_match])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # iterate through vertical Move moves
            for k ∈ vgaps
                # affine Move
                # check if the affine move lead to the current cell
                if dp_matrix[x, y] == vaffine_matrix[x - k.step, y] + affine_gap * k.step
                    
                    # write the result
                    for i ∈ 1:k.step
                        push!(res_A, A[x - i - vertical_offset + biggest_match])
                        push!(res_B, DNA_Gap)
                    end
                    x -= k.step
                    # how long the affine Move is
                    while x > k.step && vaffine_matrix[x, y] == vaffine_matrix[x - k.step, y] + affine_gap * k.step
                        for i ∈ 1:k.step
                            push!(res_A, A[x - i - vertical_offset + biggest_match])
                            push!(res_B, DNA_Gap)
                        end
                        x -= k.step
                    end

                    # include the opening Move
                    for i ∈ 1:k.step
                        push!(res_A, A[x - i - vertical_offset + biggest_match])
                        push!(res_B, DNA_Gap)
                    end
                    x -= k.step
                    # break to stop iterating through moves
                    break

                # opening Move
                # check if the move lead to the current cell
                elseif dp_matrix[x, y] == dp_matrix[x - k.step, y] + k.score
                    for i ∈ 1:k.step
                        push!(res_A, A[x - i - vertical_offset + biggest_match])
                        push!(res_B, DNA_Gap)
                    end
                    x -= k.step

                    # break to stop iterating through moves
                    break
                end
            end

            # iterate through horizontal Move moves
            for k ∈ hgaps
                # affine Move
                # check if the affine move lead to the current cell
                if dp_matrix[x, y] == haffine_matrix[x, y - k.step] + affine_gap * k.step
                    
                    # write the resulting sequences
                    for i ∈ 1:k.step
                        push!(res_A, DNA_Gap)
                        push!(res_B, B[y - horizontal_offset - i + biggest_match])
                    end
                    y -= k.step

                    # how long the affine Move is
                    while y > k.step && haffine_matrix[x, y] == haffine_matrix[x, y - k.step] + affine_gap * k.step
                        for i ∈ 1:k.step
                            push!(res_A, DNA_Gap)
                            push!(res_B, B[y - horizontal_offset - i + biggest_match])
                        end
                        y -= k.step
                    end
                    for i ∈ 1:k.step
                        push!(res_A, DNA_Gap)
                        push!(res_B, B[y - horizontal_offset - i + biggest_match])
                    end
                    y -= k.step
                    break

                # opening Move
                # check if the move lead to the current cell
                elseif dp_matrix[x, y] == dp_matrix[x, y - k.step] + k.score
                    for i ∈ 1:k.step
                        push!(res_A, DNA_Gap)
                        push!(res_B, B[y - horizontal_offset - i + biggest_match])
                    end
                    y -= k.step

                    # break to stop iterating through moves
                    break
                end
            end

            # iterate through digonal match moves
            for k ∈ matches
                
                s = 0
                iB = y - horizontal_offset + biggest_match - k.step

                for iA ∈ x - vertical_offset + biggest_match - k.step:x - vertical_offset - 1 + biggest_match
                    s += match_score_matrix[toInt(A[iA]), toInt(B[iB])]
                    iB += 1
                end

                if dp_matrix[x - k.step, y - k.step] + k.score + s < dp_matrix[x, y]
                    dp_matrix[x, y] = dp_matrix[x - k.step, y - k.step] + k.score + s 
                end

                # check if the move lead to the current cell
                if dp_matrix[x, y] == dp_matrix[x - k.step, y - k.step] + k.score + s
                    
                    # write the resulting sequences
                    # k.step and k.step is the same
                    for i ∈ 1:k.step
                        push!(res_A, A[x - vertical_offset - i + biggest_match])
                        push!(res_B, B[y - horizontal_offset - i + biggest_match])
                    end
                    x -= k.step
                    y -= k.step

                    # break to stop iterating through moves
                    break
                end
            end
        end
    end

    deleteat!(A, 1:biggest_match - 1)
    deleteat!(B, 1:biggest_match - 1)

    return reverse(res_A), reverse(res_B)
end

print(general_pairwise_aligner(LongDNA{2}("TTCGACTG"), LongDNA{2}("TACGACGACTG"), .0, 0.5, [((1, 1), .0), ((1, 0), 1.0), ((0, 1), 1.0), ((3, 3), 0.0), ((3, 0), 2.0), ((0, 3), 2.0)], 0.5))
