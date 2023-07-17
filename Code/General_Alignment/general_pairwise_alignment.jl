# add fs_penalty

module GeneralAlignments
export general_pairwise_aligner, dp_alignment, Move


struct Move
    move::Tuple{Int16, Int16}
    score::Number
end

function init_matrices(n::Int64, m::Int64, vgaps::Vector{Move}, hgaps::Vector{Move})
    dp_matrix = fill(Inf, n + 1, m + 1)
    dp_matrix[1, 1] = .0
    bt_matrix = fill((0, 0), n + 1, m + 1)

    for i ∈ 2:n + 1, gap ∈ vgaps
        if gap.move[1] < i && dp_matrix[i - gap.move[1], 1] + gap.score < dp_matrix[i, 1]
            dp_matrix[i, 1] = dp_matrix[i - gap.move[1], 1] + gap.score
            bt_matrix[i, 1] = (i - gap.move[1], 1)
        end
    end

    for i ∈ 2:m + 1, gap ∈ hgaps
        if gap.move[2] < i && dp_matrix[1, i - gap.move[2]] + gap.score < dp_matrix[1, i]
            dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + gap.score
            bt_matrix[1, i] = (1, i - gap.move[2])
        end
    end

    return dp_matrix, bt_matrix
end

function dp_alignment(A::String, B::String, vgaps::Vector{Move}, hgaps::Vector{Move}, moves::Vector{Move}, mismatch_score::Float64, long::Int64)
    n, m = length(A), length(B)

    dp_matrix = fill(typemax(Float64), long, m + 1)
    dp_matrix[1, 1] = .0
    bt_matrix = fill((0, 0), n + 1, m + 1)
    ma_matrix = zeros(Float16, long, m)


    for i ∈ 2:m + 1, gap ∈ hgaps
        if gap.move[2] < i && dp_matrix[1, i - gap.move[2]] + gap.score < dp_matrix[1, i]
            dp_matrix[1, i] = dp_matrix[1, i - gap.move[2]] + gap.score
            bt_matrix[1, i] = (1, gap.move[2])
        end
    end

    cyc = [i for i ∈ axes(dp_matrix, 1)] # 1, 2, 3...

    for i ∈ 2:n + 1


        # Asign value to first column in new row
        if vgap.move[1] < i 
            if vgap.move[1] > cyc[2] && dp_matrix[cyc[2] - vgap.move[1] + long, 1] + vgap.score < dp_matrix[cyc[2], 1]
                dp_matrix[cyc[2], 1] = dp_matrix[cyc[2] - vgap.move[1] + long, 1] + vgap.score
                bt_matrix[i, 1] = (i - vgap.move[1], 1)
            elseif vgap.move[1] < cyc[2] && dp_matrix[cyc[2] - vgap.move[1], 1] + vgap.score < dp_matrix[cyc[2], 1]
                dp_matrix[cyc[2], 1] = dp_matrix[cyc[2] - vgap.move[1], 1] + vgap.score
                bt_matrix[i, 1] = (i - vgap.move[1], 1)
            end
        end

        for j ∈ 2:m + 1
            if A[i - 1] != B[j - 1] 
                ma_matrix[1, j - 1] = mismatch_score # some other index for the row
            end

            min_move = Move((0, 0), Inf)

            for k ∈ moves
                if k.move[1] != k.move[2] && k.move[1] < i && k.move[2] < j && dp_matrix[i - k.move[1], j - k.move[2]] + k.score < min_move.score 
                    min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score)
                elseif k.move[1] < i && k.move[2] < j && dp_matrix[i - k.move[1], j - k.move[2]] + k.score + sum([ma_matrix[i, i] for i in 1:3]) < min_move.score 
                    min_move = Move(k.move, dp_matrix[i - k.move[1], j - k.move[2]] + k.score + sum([ma_matrix[i, i] for i in 1:3]))
                end
            end

            dp_matrix[i, j] = min_move.score
            bt_matrix[i, j] = min_move.move


            # # Move left
            # if vgap.move[1] < i
            #     if vgap.move[1] > cyc[2] && dp_matrix[cyc[2] - vgap.move[1] + long, j] + vgap.score < dp_matrix[cyc[2], j]
            #         dp_matrix[cyc[2], j] = dp_matrix[cyc[2] - vgap.move[1] + long, j] + vgap.score
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
            #             if match.move[1] > cyc[2] && dp_matrix[cyc[2] - match.move[1] + long, j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1] < dp_matrix[cyc[2], j]
            #                 dp_matrix[cyc[2], j] = dp_matrix[cyc[2] - match.move[1] + long, j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1]
            #                 bt_matrix[i, j] = (i - match.move[1], j - match.move[2])
            #             elseif match.move[1] > cyc[2] && dp_matrix[cyc[2] - match.move[1] + long, j - match.move[2]] + match.score + ma_matrix[cyc[1], j - 1] < dp_matrix[cyc[2], j]
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
        tmp = cyc[1]
        cyc[1] = cyc[2]
        cyc[2] = cyc[3]
        cyc[3] = tmp
    end

    return bt_matrix
end

function general_pairwise_aligner(A::String, B::String, match_score::Number, mismatch_score::Number, moves_with_penalties::Tuple) 
    # match/mismatch scores between ACGT and ACGT
    match_matrix = zeros(4, 4)

    for i in 1:4, j in 1:4
        if i == j
            match_matrix[i, j] = match_score
        else
            match_matrix[i, j] = mismatch_score
        end
    end

    # sort the moves into gaps aswell as 
    moves = Vector{Move}()
    vgaps = Vector{Move}()
    hgaps = Vector{Move}()
    matches = Vector{Move}()
    
    long = 0

    for i ∈ 1:2:length(moves_with_penalties)
        if moves_with_penalties[i][1] > long
            long = moves_with_penalties[i][1]
        end
        if moves_with_penalties[i][2] > long
            long = moves_with_penalties[i][2]
        end

        # gap and match
        if moves_with_penalties[i][1] == moves_with_penalties[i][2]
            push!(matches, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
        elseif moves_with_penalties[i][2] == 0
            push!(vgaps, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
        else
            push!(hgaps, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
        end

        push!(moves, Move(moves_with_penalties[i], moves_with_penalties[i + 1]))
    end

    # initiate dp_matrix, bt_matrix and ma_matrix
    # bt_matrix = init_matrices(n, m, vgaps, hgaps)
    bt_matrix = dp_alignment(A, B, vgaps, hgaps, moves, mismatch_score, long)

    display(bt_matrix)
    println(A, "\n", B)
end




function general_pairwise_aligner(A::String, B::String, match_score::Number, missmatch_score::Number, moves::Tuple, affine_gap::Number) 
    match_matrix = zeros(4, 4)

    for i ∈ 1:4, j ∈ 1:4
        if i == j
            match_matrix[i, j] = match_score
        else
            match_matrix[i, j] = missmatch_score
        end
    end

    n, m = length(A), length(B)

    dp_matrix = zeros(Float64, n + 1, m + 1)
    bt_matrix = zeros(Int64, n + 1, m + 1)


end



# println(A, "\n", B)


# general_pairwise_aligner("ACGT", "ACGT", 0, 0.5, ((1, 1), 1, (1, 0), 2, (0, 1), 2, (3, 3), 0, (3, 0), 2, (0, 3), 2))
end