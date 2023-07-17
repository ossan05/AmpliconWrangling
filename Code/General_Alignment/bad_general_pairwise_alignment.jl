module BadGeneralAlignments
export general_pairwise_aligner, dp_alignment, Move


struct Move
    move::Tuple{Int16, Int16}
    score::Number
end

struct Gap
    move::Int16
    score::Number
end

# dp_alignment is the needleman-wunsch algorithm
function dp_alignment(A::String, B::String, vgaps::Vector{Gap}, hgaps::Vector{Gap}, moves::Vector{Move}, frameshift_score::Number, mismatch_score::Number)
    n, m = length(A), length(B)

    # initialize matrices
    dp_matrix = fill(typemax(Float64), n + 1, m + 1)
    dp_matrix[1, 1] = .0
    bt_matrix = fill((0, 0), n + 1, m + 1)
    ma_matrix = zeros(Float16, n, m)

    # Fill the first row in the dp_matrix
    for i ∈ 2:m + 1, gap ∈ hgaps 
        if gap.move[2] % 3 == 0
            if gap.move[2] < i && dp_matrix[1, i - gap.move[2]] + gap.score < dp_matrix[1, i]
                dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + gap.score
                bt_matrix[1, i] = (0, gap.move[2])
            end
        else
            if gap.move[2] < i && dp_matrix[1, i - gap.move[2]] + gap.score + frameshift_score < dp_matrix[1, i]
                dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + gap.score + frameshift_score
                bt_matrix[1, i] = (0, gap.move[2])
            end
        end
    end

    # Fill the first column in the dp_matrix
    for i ∈ 2:m + 1, gap ∈ vgaps 
        if gap.move[1] % 3 == 0
            if gap.move[1] < i && dp_matrix[i - gap.move[1], 1] + gap.score < dp_matrix[i, 1]
                dp_matrix[i, 1] = dp_matrix[i - gap.move[1], 1] + gap.score
                bt_matrix[i, 1] = (gap.move[1], 0)
            end
        else
            if gap.move[1] < i && dp_matrix[i - gap.move[1], 1] + gap.score + frameshift_score < dp_matrix[i, 1]
                dp_matrix[i, 1] = dp_matrix[i - gap.move[1], 1] + gap.score + frameshift_score
                bt_matrix[i, 1] = (gap.move[1], 0)
            end
        end
    end 

    for i ∈ 2:n + 1
        for j ∈ 2:m + 1

            # Store mismatch score in dp sized mismatch matrix
            if A[i - 1] != B[j - 1] 
                ma_matrix[1, j - 1] = mismatch_score # some other index for the row
            end

            # the minimum move is stored in an object and changed when a better move is found
            # constructed with the struct "Move"
            min_move = Move((0, 0), Inf)

            # calculate the best move by iterating over the moves 
            for k ∈ moves

                # k.move is a tuple e.g. (0, 3)

                # Checks a move is a multiple of three
                if (k.move[1] > 2 && k.move[1] % 3 == 0) || (k.move[2] > 2 && k.move[2] % 3 == 0 )

                    # if the move is a gap:
                    if k.move[1] != k.move[2] && k.move[1] < i && k.move[2] < j && dp_matrix[i - k.move[1], j - k.move[2]] + k.score < min_move.score 
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score)

                    # if the move is diagonal
                    elseif k.move[1] < i && k.move[2] < j

                        # the sum of the mismatch scores on the diagonal
                        s = 0
                        for l ∈ i - k.move[1]:i - 1
                            s += ma_matrix[l, l]
                        end                       # sum([ma_matrix[x, x] for x ∈ i - k.move:i - 1])

                        if dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s < min_move.score 
                            min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s)
                        end
                    end
                else # when the move is not a multiple of three

                    # check if it's a gap
                    if k.move[1] != k.move[2] && k.move[1] < i && k.move[2] < j && dp_matrix[i - k.move[1], j - k.move[2]] + k.score + frameshift_score < min_move.score 
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + frameshift_score)
                    # if it's a match (diagonal move)
                    elseif k.move[1] < i && k.move[2] < j
                        s = 0
                        for l ∈ i - k.move[1]:i - 1
                            s += ma_matrix[l, l]
                        end

                        # check if the diagonal move is less than the current best move
                        if dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s + frameshift_score < min_move.score 
                            min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s + frameshift_score)
                        end
                    end
                end
            end

            # put the new minimum score in the dp_matrix
            dp_matrix[i, j] = min_move.score

            # and the move in the backtracking matrix
            bt_matrix[i, j] = min_move.move
        end
    end

    # backtracking
    y = m + 1
    x = n + 1
    res_A = ""
    res_B = ""

    while x > 1 || y > 1
        if bt_matrix[x, y][1] == 0
            for i ∈ 1:bt_matrix[x, y][2]
                res_B *= '_'
                res_A *= A[x - i]
            end
        elseif bt_matrix[x, y][2] == 0
            for i ∈ 1:bt_matrix[x, y][1]
                res_A *= '_'
                res_B *= B[y - i]
            end
        else
            for i ∈ 1:bt_matrix[x, y][1]
                res_A *= A[x - i]
                res_B *= B[y - i]
            end
        end

        x_prev = x
        x -= bt_matrix[x, y][1]
        y -= bt_matrix[x_prev, y][2]
    end

    # return the resulting strings
    return reverse(res_A), reverse(res_B)
end

# same thing but with affine gap
function dp_alignment(A::String, B::String, vgaps::Vector{Move}, hgaps::Vector{Move}, moves::Vector{Move}, mismatch_score::Number, affine_gap::Number, frameshift_score::Number)
    n, m = length(A), length(B)

    dp_matrix = fill(typemax(Float64), n + 1, m + 1)
    dp_matrix[1, 1] = .0
    bt_matrix = fill((0, 0), n + 1, m + 1)
    ma_matrix = zeros(Float16, n, m)


    for i ∈ 2:m + 1, gap ∈ hgaps 
        if gap.move[2] % 3 == 0
            if gap.move[2] < i + 1 && dp_matrix[1, i - gap.move[2]] + affine_gap * gap.move[2] < dp_matrix[1, i]
                dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + affine_gap * gap.move[2]
                bt_matrix[1, i] = (0, gap.move[2])
            elseif gap.move[2] < i && dp_matrix[1, i - gap.move[2]] + gap.score / gap.move[2] + affine_gap * (gap.move[2] - 1) < dp_matrix[1, i]
                dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + gap.score / gap.move[2] + affine_gap * (gap.move[2] - 1)
                bt_matrix[1, i] = (0, gap.move[2])
            end
        else
            if gap.move[2] < i + 1 && dp_matrix[1, i - gap.move[2]] + affine_gap * gap.move[2] + frameshift_score < dp_matrix[1, i]
                dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + affine_gap * gap.move[2] + frameshift_score
                bt_matrix[1, i] = (0, gap.move[2])
            elseif gap.move[2] < i && dp_matrix[1, i - gap.move[2]] + gap.score / gap.move[2] + affine_gap * (gap.move[2] - 1) + frameshift_score < dp_matrix[1, i]
                dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + gap.score / gap.move[2] + affine_gap * (gap.move[2] - 1) + frameshift_score
                bt_matrix[1, i] = (0, gap.move[2])
            end
        end
    end

    for i ∈ 2:m + 1, gap ∈ vgaps 
        if gap.move[1] < i + 1
        if gap.move[1] < i && dp_matrix[i - gap.move[1], 1] + gap.score < dp_matrix[i, 1]
            dp_matrix[i, 1] = dp_matrix[i - gap.move[1], 1] + gap.score
            bt_matrix[i, 1] = (gap.move[1], 0)
        end
    end 


    for i ∈ 2:n + 1
        for j ∈ 2:m + 1

            if A[i - 1] != B[j - 1] 
                ma_matrix[1, j - 1] = mismatch_score # some other index for the row
            end

            min_move = Move((0, 0), Inf)

            for k ∈ vgaps
                if k.move[1] % 3 == 0
                    if k.move[1] < i + 1 && bt_matrix[i - k.move[1], j][2] == 0 && dp_matrix[i - k.move[1], j] + affine_gap * k.move[1] < min_move.score
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j] + affine_gap * k.move[1])
                    elseif k.move[1] < i && dp_matrix[i - k.move[1], j] + k.score < min_move.score # k.score / k.move[1] + affine_gap * (k.move[1] - 1)
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j] + affine_gap * k.move[1])
                    end
                else
                    if k.move[1] < i + 1 && bt_matrix[i - k.move[1], j][2] == 0 && dp_matrix[i - k.move[1], j] + affine_gap * k.move[1] + frameshift_score < min_move.score
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j] + affine_gap * k.move[1])
                    elseif k.move[1] < i && bt_matrix[i - k.move[1], j][2] == 0 && dp_matrix[i - k.move[1], j] + k.score / k.move[1] + affine_gap * (k.move[1] - 1) < min_move.score # k.score / k.move[1] + affine_gap * (k.move[1] - 1)
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j] + affine_gap * k.move[1])
                    end
                end
            end

            for k ∈ hgaps

            # iterate over gaps and matches seperately instead of moves to remove if statements
            for k ∈ moves
                if k.move[1] % 3 == 0 || k.move[2] % 3 == 0

                    if k.move[1] == 0 && k.move[2] < j + 1 && bt_matrix[i, j - k.move[1]][1] == 0 && dp_matrix[i, j - k.move[2]] + affine_gap * k.move[2] < min_move.score 
                        min_move = Move(k.move, dp_matrix[i, j - k.move[2]] + affine_gap * k.move[2])
                    elseif k.move[2] == 0 && k.move[1] < i + 1 && bt_matrix[i - k.move[1], j][2] == 0 && dp_matrix[i - k.move[1], j] + affine_gap * k.move[1] < min_move.score
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j] + affine_gap * k.move[1])
                    elseif k.move[1] == 0 && k.move[2] < j && bt_matrix[i, j - k.move[1]][1] == 0 && dp_matrix[i, j - k.move[2]] + affine_gap * k.move[2] < min_move.score 
                        min_move = Move(k.move, dp_matrix[i, j - k.move[2]] + affine_gap * k.move[2])
                    elseif k.move[1] < i && k.move[2] < j
                        s = 0
                        for l ∈ i - k.move[1]:i - 1
                            s += ma_matrix[l, l]
                        end                       # sum([ma_matrix[x, x] for x ∈ i - k.move:i - 1])

                        if dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s < min_move.score 
                            min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s)
                        end
                    end
            end

            dp_matrix[i, j] = min_move.score
            bt_matrix[i, j] = min_move.move
        end
    end

    y = m + 1
    x = n + 1
    res_A = ""
    res_B = ""

    while x > 1 || y > 1
        if bt_matrix[x, y][1] == 0
            for i ∈ 1:bt_matrix[x, y][2]
                res_B *= '_'
                res_A *= A[x - i]
            end
        elseif bt_matrix[x, y][2] == 0
            for i ∈ 1:bt_matrix[x, y][1]
                res_A *= '_'
                res_B *= B[y - i]
            end
        else
            for i ∈ 1:bt_matrix[x, y][1]
                res_A *= A[x - i]
                res_B *= B[y - i]
            end
        end

        x_prev = x
        x -= bt_matrix[x, y][1]
        y -= bt_matrix[x_prev, y][2]
    end

    return reverse(res_A), reverse(res_B)
end

# general_pairwise_aligner makes a match/mismatch matrix and takes the tuple with the moves and scores, 
# making them objects of type Move or Gap and puting them into lists 
function general_pairwise_aligner(A::String, B::String, match_score::Number, mismatch_score::Number, frameshift_score::Number, moves_with_penalties::Tuple) 
    # match/mismatch scores between ACGT and ACGT
    match_matrix = zeros(4, 4)

    # this is when every mismatch between nucletides have the same penalty 
    for i in 1:4, j in 1:4
        if i == j
            match_matrix[i, j] = match_score
        else
            match_matrix[i, j] = mismatch_score
        end
    end

    # sort the moves into gaps 
    # moves are all the moves
    # v - vertical, h - horizontal
    moves = Vector{Move}()
    vgaps = Vector{Move}()
    hgaps = Vector{Move}()

    # fill the vectors with gaps and moves
    for i ∈ 1:2:length(moves_with_penalties)

        # vertical gap
        if moves_with_penalties[i][2] == 0
            push!(vgaps, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))

        # horizontal gap
        elseif moves_with_penalties[i][1] == 0
            push!(hgaps, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
        end

        # all moves
        push!(moves, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
    end
    
    # make the alignment
    return dp_alignment(A, B, vgaps, hgaps, moves, frameshift_score, mismatch_score) 
end

# Same thing with affine gap
function general_pairwise_aligner(A::String, B::String, match_score::Number, mismatch_score::Number, moves_with_penalties::Tuple, affine_gap::Number) 
    # match/mismatch scores between ACGT and ACGT
    match_matrix = zeros(4, 4)

    for i in 1:4, j in 1:4
        if i == j
            match_matrix[i, j] = match_score
        else
            match_matrix[i, j] = mismatch_score
        end
    end

    # sort the moves into gaps 
    moves = Vector{Move}()
    vgaps = Vector{Move}()
    hgaps = Vector{Move}()

    for i ∈ 1:2:length(moves_with_penalties)
        if moves_with_penalties[i][2] == 0
            push!(vgaps, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
        elseif moves_with_penalties[i][1] == 0
            push!(hgaps, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
        end

        push!(moves, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
    end
    
    return dp_alignment(A, B, vgaps, hgaps, moves, mismatch_score, affine_gap)

end

end

# general_pairwise_aligner("ACGT", "ACGT", 0, 0.5, ((1, 1), 1, (1, 0), 2, (0, 1), 2, (3, 3), 0, (3, 0), 2, (0, 3), 2))