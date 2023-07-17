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

# w/o affine gap
function dp_alignment(A::String, B::String, vgaps::Vector{Gap}, hgaps::Vector{Gap}, moves::Vector{Move}, frameshift_score::Number, mismatch_score::Number)
    n, m = length(A), length(B)

    dp_matrix = fill(typemax(Float64), n + 1, m + 1)
    dp_matrix[1, 1] = .0
    bt_matrix = fill((0, 0), n + 1, m + 1)
    ma_matrix = zeros(Float16, n, m)


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

            if A[i - 1] != B[j - 1] 
                ma_matrix[1, j - 1] = mismatch_score # some other index for the row
            end

            min_move = Move((0, 0), Inf)

            for k ∈ moves
                if k.move[1] % 3 == 0 || k.move[2] % 3 == 0
                    if k.move[1] != k.move[2] && k.move[1] < i && k.move[2] < j && dp_matrix[i - k.move[1], j - k.move[2]] + k.score < min_move.score 
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score)
                    elseif k.move[1] < i && k.move[2] < j
                        s = 0
                        for l ∈ i - k.move[1]:i - 1
                            s += ma_matrix[l, l]
                        end                       # sum([ma_matrix[x, x] for x ∈ i - k.move:i - 1])

                        if dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s < min_move.score 
                            min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s)
                        end
                    end
                else
                    if k.move[1] != k.move[2] && k.move[1] < i && k.move[2] < j && dp_matrix[i - k.move[1], j - k.move[2]] + k.score + frameshift_score < min_move.score 
                        min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + frameshift_score)
                    elseif k.move[1] < i && k.move[2] < j
                        s = 0
                        for l ∈ i - k.move[1]:i - 1
                            s += ma_matrix[l, l]
                        end

                        if dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s + frameshift_score < min_move.score 
                            min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + s + frameshift_score)
                        end
                    end
                end
            end

            dp_matrix[i, j] = min_move.score
            bt_matrix[i, j] = min_move.move


            # # Move left
            # if vgap.move[1] < i
            #     if vgap.move[1] > cyc[2] && dp_matrix[cyc[2] - vgap.move[1] + , j] + vgap.score < dp_matrix[cyc[2], j]
            #         dp_matrix[cyc[2], j] = dp_matrix[cyc[2] - vgap.move[1] + , j] + vgap.score
            #         bt_matrix[i, j] = (i - vgap.move[1], j)
            #     elseif vgap.move[1] < cyc[2] && dp_matrix[cyc[2] - vgap.move[1], j] + vgap.score < dp_matrix[cyc[2], j]
            #         dp_matrix[cyc[2], j] = dp_matrix[cyc[2] - vgap.move[1], j] + vgap.score
            #         bt_matrix[i, j] = (i - vgap.move[1], j)
            #     end
            # end

            # # Move right
            # for hgap ∈ hgaps
            #     if hgap.move[2] < j && dp_matrix[cyc[2], j  - hgap.move[2]] + hgap.score < dp_matrix[cyc[2], j]
            #         dp_matrix[cyc[2], j] = dp_matrix[cyc[2], j - hgap.move[2]] + hgap.score
            #         bt_matrix[i, j] = (i, j - hgap.move[2])
            #     end
            # end

            # # Move diagonaly
            # for match ∈ matches
            #     if match.move[1] < i && match.move[2] < j 
            #         if match.move[1] > cyc[2]
            #             if match.move[1] > cyc[2] && dp_matrix[cyc[2] - match.move[1] + , j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1] < dp_matrix[cyc[2], j]
            #                 dp_matrix[cyc[2], j] = dp_matrix[cyc[2] - match.move[1] + , j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1]
            #                 bt_matrix[i, j] = (i - match.move[1], j - match.move[2])
            #             elseif match.move[1] > cyc[2] && dp_matrix[cyc[2] - match.move[1] + , j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1] < dp_matrix[cyc[2], j]
            #                 dp_matrix[cyc[2], j] = dp_matrix[cyc[2] - match.move[1], j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1]
            #                 bt_matrix[i, j] = (i - match.move[1], j - match.move[2])
            #             end
            #         else
            #             if match.move[1] > cyc[2] && dp_matrix[cyc[2] - match.move[1], j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1] < dp_matrix[cyc[2], j]
            #                 dp_matrix[cyc[2], j] = dp_matrix[cyc[2] - match.move[1], j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1]
            #                 bt_matrix[i, j] = (i - match.move[1], j - match.move[2])
            #             elseif match.move[1] > cyc[2] && dp_matrix[cyc[2] - match.move[1], j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1] < dp_matrix[cyc[2], j]
            #                 dp_matrix[cyc[2], j] = dp_matrix[cyc[2] - match.move[1], j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1]
            #                 bt_matrix[i, j] = (i - match.move[1], j - match.move[2])
            #             end
            #         end
            #     end
            # end
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

# w/o affine gap
function general_pairwise_aligner(A::String, B::String, match_score::Number, mismatch_score::Number, frameshift_score::Number, moves_with_penalties::Tuple) 
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
    
    return dp_alignment(A, B, vgaps, hgaps, moves, mismatch_score) 
end


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