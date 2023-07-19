module BadGeneralAlignments
export general_pairwise_aligner, dp_alignment, Move


struct Move
    move::Tuple{Int64, Int64}
    score::Number
end

struct Gap
    move::Int64
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

    # Fill the first row ∈ the dp_matrix
    for i ∈ 2:m + 1, gap ∈ hgaps 
        if gap.move[2] % 3 == 0
            if gap.move[2] < i && dp_matrix[1, i - gap.move[2]] + gap.score < dp_matrix[1, i]
                dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + gap.score
                bt_matrix[1, i] = (0, gap.move[2])
            end
        else
            if gap.move[2] < i && dp_matrix[1, i - gap.move[2]] + gap.score < dp_matrix[1, i]
                dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + gap.score
                bt_matrix[1, i] = (0, gap.move[2])
            end
        end
    end

    # Fill the first column ∈ the dp_matrix
    for i ∈ 2:m + 1, gap ∈ vgaps 
        if gap.move[1] % 3 == 0
            if gap.move[1] < i && dp_matrix[i - gap.move[1], 1] + gap.score < dp_matrix[i, 1]
                dp_matrix[i, 1] = dp_matrix[i - gap.move[1], 1] + gap.score
                bt_matrix[i, 1] = (gap.move[1], 0)
            end
        else
            if gap.move[1] < i && dp_matrix[i - gap.move[1], 1] + gap.score < dp_matrix[i, 1]
                dp_matrix[i, 1] = dp_matrix[i - gap.move[1], 1] + gap.score
                bt_matrix[i, 1] = (gap.move[1], 0)
            end
        end
    end 

    for i ∈ 2:n + 1
        for j ∈ 2:m + 1

            # Store mismatch score ∈ dp sized mismatch matrix
            if A[i - 1] != B[j - 1] 
                ma_matrix[1, j - 1] = mismatch_score # some other index for the row
            end

            # the minimum move is stored ∈ an object and changed when a better move is found
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
                    if k.move[1] != k.move[2] && k.move[1] < i && k.move[2] < j && dp_matrix[i - k.move[1], j - k.move[2]] + k.score < min_move.score 
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score)
                    # if it's a match (diagonal move)
                    elseif k.move[1] < i && k.move[2] < j
                        s = 0
                        for l ∈ i - k.move[1]:i - 1
                            s += ma_matrix[l, l]
                        end

                        # check if the diagonal move is less than the current best move
                        if dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s < min_move.score 
                            min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s)
                        end
                    end
                end
            end

            # put the new minimum score ∈ the dp_matrix
            dp_matrix[i, j] = min_move.score

            # and the move ∈ the backtracking matrix
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

        x_copy = x
        x -= bt_matrix[x, y][1]
        y -= bt_matrix[x_copy, y][2]
    end

    # return the resulting strings
    return reverse(res_A), reverse(res_B)
end

# same thing but with affine gap

# general_pairwise_aligner makes a match/mismatch matrix and takes the tuple with the moves and scores, 
# making them objects of type Move or Gap and puting them into lists 
function general_pairwise_aligner(A::String, B::String, match_score::Number, mismatch_score::Number, frameshift_score::Number, moves_with_penalties::Tuple) 
    # match/mismatch scores between ACGT and ACGT
    match_matrix = zeros(4, 4)

    # this is when every mismatch between nucletides have the same penalty 
    for i ∈ 1:4, j ∈ 1:4
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
function general_pairwise_aligner(A::String, B::String, match_score::Number, mismatch_score::Number, frameshift_score::Number, moves_with_penalties::Tuple, affine_gap::Number) 
    # match/mismatch scores between ACGT and ACGT
    match_matrix = zeros(4, 4)

    for i ∈ 1:4, j ∈ 1:4
        if i == j
            match_matrix[i, j] = match_score
        else
            match_matrix[i, j] = mismatch_score
        end
    end

    # sort the moves into gaps 
    matches = Vector{Move}()
    vgaps = Vector{Gap}()
    hgaps = Vector{Gap}()

    for i ∈ 1:2:length(moves_with_penalties)
        if moves_with_penalties[i][2] == 0
            if moves_with_penalties[i][1] % 3 == 0
                push!(vgaps, Gap(moves_with_penalties[i][1], moves_with_penalties[i + 1]))
            else
                # move is not divisable by three, hence the added frameshift_score
                push!(vgaps, Gap(moves_with_penalties[i][1], moves_with_penalties[i + 1] + frameshift_score))
            end

        elseif moves_with_penalties[i][1] == 0
            if moves_with_penalties[i][2] % 3 == 0
                push!(hgaps, Gap(moves_with_penalties[i][2], moves_with_penalties[i + 1]))
            else
                # move is not divisable by three, hence the added frameshift_score                
                push!(hgaps, Gap(moves_with_penalties[i][2], moves_with_penalties[i + 1] + frameshift_score))
            end
        elseif moves_with_penalties[i][1] % 3 == 0
            push!(matches, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
        else
            # move is not divisable by three, hence the added frameshift_score
            push!(matches, Move(moves_with_penalties[i], moves_with_penalties[i + 1] + frameshift_score))
        end
    end


    # Needleman-Wunsch alignment

    n, m = length(A), length(B)

    dp_matrix = fill(typemax(Float64), n + 1, m + 1)
    dp_matrix[1, 1] = .0
    bt_matrix = fill((0, 0), n + 1, m + 1)
    vaffine_matrix = zeros(Float64, n + 1, m + 1)
    haffine_matrix = zeros(Float64, n + 1, m + 1)


    # first row
    for i ∈ 2:m + 1, gap ∈ hgaps 
        if gap.move + 1 < i && dp_matrix[1, i - gap.move] + affine_gap * gap.move < dp_matrix[1, i]
            dp_matrix[1, i] = dp_matrix[1, i - gap.move] + affine_gap * gap.move
            haffine_matrix[1, i] = dp_matrix[1, i]
            vaffine_matrix[1, i] = Inf
            bt_matrix[1, i] = (0, gap.move)
        elseif gap.move < i && dp_matrix[1, i - gap.move] + gap.score < dp_matrix[1, i]
            dp_matrix[1, i] = dp_matrix[1, i - gap.move] + gap.score 
            haffine_matrix[1, i] = dp_matrix[1, i]
            vaffine_matrix[1, i] = Inf
            bt_matrix[1, i] = (0, gap.move)
        end
    end

    # first column
    for i ∈ 2:n + 1, gap ∈ vgaps 
        if gap.move + 1 < i && dp_matrix[i - gap.move, 1] + affine_gap * gap.move < dp_matrix[i, 1]
            dp_matrix[i, 1] = dp_matrix[i - gap.move, 1] + affine_gap * gap.move
            vaffine_matrix[i, 1] = dp_matrix[i, 1]
            haffine_matrix[i, 1] = Inf
            bt_matrix[i, 1] = (gap.move, 0)
        elseif gap.move < i && dp_matrix[i - gap.move, 1] + gap.score / gap.move + affine_gap * (gap.move - 1) < dp_matrix[i, 1]
            dp_matrix[i, 1] = dp_matrix[i - gap.move, 1] + gap.score / gap.move + affine_gap * (gap.move - 1)
            vaffine_matrix[i, 1] = dp_matrix[i, 1]
            haffine_matrix[i, 1] = Inf
            bt_matrix[i, 1] = (gap.move, 0)
        end
    end

    for i ∈ 2:n + 1
        for j ∈ 2:m + 1

            # finds the best vertical move
            vaffine_matrix[i, j] = Inf

            for k ∈ vgaps
                if k.move < i
                    # affine gap
                    if vaffine_matrix[i - k.move, j] + affine_gap * k.move < vaffine_matrix[i, j]
                        vaffine_matrix[i, j] = vaffine_matrix[i - k.move, j] + affine_gap * k.move
                    end
                    # opening gap
                    if dp_matrix[i - k.move, j] + k.score < vaffine_matrix[i, j]
                        vaffine_matrix[i, j] = dp_matrix[i - k.move, j] + k.score
                    end
                end
            end

            # finds the best horizontal move
            haffine_matrix[i, j] = Inf

            for k ∈ hgaps
                if k.move < j
                    # affine_gap
                    if haffine_matrix[i, j - k.move] + affine_gap * k.move < haffine_matrix[i, j]
                        haffine_matrix[i, j] = haffine_matrix[i, j - k.move] + affine_gap * k.move
                    end
                    # opening gap
                    if dp_matrix[i, j - k.move] + k.score < haffine_matrix[i, j]
                        haffine_matrix[i, j] = dp_matrix[i, j - k.move] + k.score
                    end
                end
            end

            # the fastest diagonal move

            dp_matrix[i, j] = Inf

            for k ∈ matches
                if k.move[1] < i && k.move[2] < j
                    s = 0
                    for l ∈ i - k.move[1]:i - 1
                        if A[l] != B[l]
                            s += mismatch_score
                        end
                    end 
                    if dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s < dp_matrix[i, j]
                        dp_matrix[i, j] = dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s 
                        # asign the score to the matrix instead for all of these
                    end
                end
            end

            # fastest move to dp_matrix
            dp_matrix[i, j] = min(dp_matrix[i, j], haffine_matrix[i, j], vaffine_matrix[i, j])
        end
    end

    println()
    show(stdout, "text/plain", dp_matrix)
    
    println()
    show(stdout, "text/plain", haffine_matrix)
    
    println()
    show(stdout, "text/plain", vaffine_matrix)

    y = m + 1
    x = n + 1
    res_A = ""
    res_B = ""

    while x > 1 || y > 1
        if x == 1
            res_A *= '_'
            res_B *= B[y - 1]
            y -= 1
        elseif y == 1
            res_A *= A[x - 1]
            res_B *= '_'
            x -= 1
        else
            # iterate through vertical gap moves
            for k ∈ vgaps
                # check that the move isn't out of bounds
                if k.move < x
                    # affine gap
                    # check if the affine move lead to the current cell
                    if dp_matrix[x, y] == vaffine_matrix[x - k.move, y] + affine_gap * k.move
                        
                        # write the result
                        for i ∈ 1:k.move
                            res_A *= A[x - i]
                            res_B *= '_' 
                        end
                        x -= k.move
                        # how long the affine gap is
                        while x > k.move && vaffine_matrix[x, y] == vaffine_matrix[x - k.move, y] + affine_gap * k.move
                            for i ∈ 1:k.move
                                res_A *= A[x - i]
                                res_B *= '_' 
                            end
                            x -= k.move
                        end

                        # include the opening gap
                        for i ∈ 1:k.move
                            res_A *= A[x - i]
                            res_B *= '_' 
                        end
                        x -= k.move
                        @show x
                        # break to stop iterating through moves
                        break

                    # opening gap
                    # check if the move lead to the current cell
                    elseif dp_matrix[x, y] == dp_matrix[x - k.move, y] + k.score
                        for i ∈ 1:k.move
                            res_A *= A[x - i]
                            res_B *= '_' 
                        end
                        x -= k.move

                        @show x
                        # break to stop iterating through moves
                        break
                    end
                end
            end

            # iterate through horizontal gap moves
            for k ∈ hgaps
                # check that the move isn't out of bounds
                if k.move < y

                    # affine gap
                    # check if the affine move lead to the current cell
                    if dp_matrix[x, y] == haffine_matrix[x, y - k.move] + affine_gap * k.move
                        
                        # write the resulting sequences
                        for i ∈ 1:k.move
                            res_A *= '_'
                            res_B *= B[y - i]
                        end
                        y -= k.move

                        # how long the affine gap is
                        while y > k.move && haffine_matrix[x, y] == haffine_matrix[x, y - k.move] + affine_gap * k.move
                            for i ∈ 1:k.move
                                res_A *= '_'
                                res_B *= B[y - i]
                            end
                            y -= k.move
                        end
                        for i ∈ 1:k.move
                            res_A *= '_'
                            res_B *= B[y - i]
                        end
                        y -= k.move
                        @show y
                        break

                    # opening gap
                    # check if the move lead to the current cell
                    elseif dp_matrix[x, y] == dp_matrix[x, y - k.move] + k.score
                        for i ∈ 1:k.move
                            res_A *= '_'
                            res_B *= B[y - i]
                        end
                        y -= k.move

                        @show y
                        # break to stop iterating through moves
                        break
                    end
                end
            end

            # iterate through digonal match moves
            for k ∈ matches
                
                # check that the move isn't out of bounds
                if k.move[1] < x && k.move[2] < y

                    # calculate score of mismatches
                    s = 0
                    for i ∈ x - k.move[1]:x - 1
                        if A[i] != B[i]
                            s += mismatch_score
                        end
                    end

                    # check if the move lead to the current cell
                    if dp_matrix[x, y] == dp_matrix[x - k.move[1], y - k.move[2]] + k.score + s
                        
                        # write the resulting sequences
                        # k.move[1] and k.move[2] is the same
                        for i ∈ 1:k.move[1]
                            res_A *= A[x - i]
                            res_B *= B[y - i]
                        end
                        x -= k.move[1]
                        y -= k.move[2]

                        @show x, y
                        # break to stop iterating through moves
                        break
                    end
                end
            end
        end
    end



    
    # y = m + 1
    # x = n + 1
    # res_A = ""
    # res_B = ""

    # while x > 1 || y > 1
    #     if bt_matrix[x, y][1] == 0
    #         for i ∈ 1:bt_matrix[x, y][2]
    #             res_B *= '_'
    #             res_A *= A[x - i]
    #         end
    #     elseif bt_matrix[x, y][2] == 0
    #         for i ∈ 1:bt_matrix[x, y][1]
    #             res_A *= '_'
    #             res_B *= B[y - i]
    #         end
    #     else
    #         for i ∈ 1:bt_matrix[x, y][1]
    #             res_A *= A[x - i]
    #             res_B *= B[y - i]
    #         end
    #     end

    #     x_copy = x
    #     x -= bt_matrix[x, y][1]
    #     y -= bt_matrix[x_copy, y][2]
    # end

    return reverse(res_A), reverse(res_B)
end

end

# general_pairwise_aligner("TTCGACTG", "TACGACGACTG", 0, 0.5, 1, ((1, 1), 0, (1, 0), 1, (0, 1), 1, (3, 3), 0, (3, 0), 2, (0, 3), 2), 0.5)