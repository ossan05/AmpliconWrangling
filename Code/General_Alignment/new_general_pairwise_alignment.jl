using BioSequences

struct Move
    step::Int64
    score::Float64
end

function initiate_general_pairwise_aligner(seq1::LongDNA{2}, seq2::LongDNA{2})
    general_pairwise_aligner(seq1, seq2, .0, 0.5, [Move((1, 1), 0), Move((1, 0), 1), Move((0, 1), 1), Move((3, 3), 0), Move((3, 0), 2), Move((0, 3), 2)])
end
function toInt(x::NucleicAcid)
    trailing_zeros(reinterpret(UInt8,x))+1
end

function simple_match_score_matrix(mismatch_score)
    match_score_matrix = zeros(4, 4)
    for i in 1:4, j in 1:4
        if i != j
            match_score_matrix[i, j] = mismatch_score
        end
    end
    match_score_matrix
end

# match and mismatch matrix
function nw_aligner(A::LongDNA{2}, B::LongDNA{2}, mismatch_score, reg_moves::Array{Move}, gap_moves_a::Array{Move}, gap_moves_b::Array{Move})
    nw_aligner(A, B, simple_match_score_matrix(mismatch_score), reg_moves, gap_moves_a, gap_moves_b) 
end

function find_best_move(A::LongDNA{4}, B::LongDNA{4}, match_score_matrix::Array{Float64, 2}, reg_moves::Array{Move}, gap_moves_a::Array{Move}, gap_moves_b::Array{Move}, dp_matrix::Array{Float64, 2}, x::Int64, y::Int64)
    
    best_move_a = findmin(map(k -> dp_matrix[x - k.step, y] + k.score, gap_moves_a))
    best_move_b = findmin(map(k -> dp_matrix[x, y - k.step] + k.score, gap_moves_b))
    best_reg_move = findmin(map(k -> dp_matrix[x, y - k.step] + k.score + sum([match_score_matrix[toInt(A[x - t]), toInt(B[y - t])] for t in 1 : k.step]), reg_moves))
    
    scores = [best_move_a[1], best_move_b[1], best_reg_move[1]]
    moves = [gap_moves_a[best_move_a[2]], gap_moves_b[best_move_b[2]], reg_moves[best_reg_move[2]]]

    move_type = argmin(scores)
    return moves[move_type], scores[move_type], move_type 
end

# general_pairwise_aligner makes a match/mismatch matrix and takes the tuple with the moves and scores, 
# making them objects of type Move or Gap and puting them into lists 
function nw_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, reg_moves::Array{Move}, gap_moves_a::Array{Move}, gap_moves_b::Array{Move}) 
    # Needleman-Wunsch alignment
    max_move_a = max(maximum([k.step for k in gap_moves_a]), maximum([k.step for k in reg_moves]))
    max_move_b = max(maximum([k.step for k in gap_moves_a]), maximum([k.step for k in reg_moves]))
    
    m, n = length(A), length(B)

    A2 = repeat(LongDNA{4}([DNA_N]), max_move_a + 1) * A
    B2 = repeat(LongDNA{4}([DNA_N]), max_move_b + 1) * B
    
    dp_matrix = fill(typemax(Float64), m + max_move_a + 1, n + max_move_b + 1)
    dp_matrix[max_move_a + 1, max_move_b + 1] = .0

    # whole dp matrix
    for i ∈ max_move_a + 1 : m + max_move_a + 1
        for j ∈ max_move_b + 1 : n + max_move_b + 1
            dp_matrix[i, j] = find_best_move(A2, B2, match_score_matrix, reg_moves, gap_moves_a, gap_moves_b, dp_matrix, i, j)[2]
        end
    end

    x = m + max_move_a + 1 
    y = n + max_move_b + 1
    res_A = dna""
    res_B = dna""

    while x > max_move_a + 1 && y > max_move_b + 1
        best_move, best_score, move_type = find_best_move(A2, B2, match_score_matrix, reg_moves, gap_moves_a, gap_moves_b, dp_matrix, x, y)
        
        if move_type == 1 #gap move a
            for i ∈ 1 : best_move.step
                push!(res_A, A2[x - i])
                push!(res_B, DNA_Gap)
            end
            x -= best_move.step
        elseif move_type == 2
            for i ∈ 1 : best_move.step
                push!(res_A, DNA_Gap)
                push!(res_B, B2[y - i])
            end
            y -= best_move.step
        elseif move_type == 3
            for i ∈ 1 : best_move.step[1]
                push!(res_A, A2[x - i])
                push!(res_B, B2[y - i])
            end
            x -= best_move.step
            y -= best_move.step
        end
    end

    return reverse(res_A), reverse(res_B)
end

A = LongDNA{2}("TTCGACTG")
B = LongDNA{2}("TACGACGACTG") 
reg_moves = [Move(1, 0)]
gap_moves = [Move(3, 2), Move(1, 1)]
print(nw_aligner(A, B, 0.5, reg_moves, gap_moves, gap_moves))
#print(general_pairwise_aligner(LongDNA{2}("TTCGACTG"), LongDNA{2}("TACGACGACTG"), .0, 0.5, [Move((1, 1), 0), Move((1, 0), 1), Move((0, 1), 1), Move((3, 3), 0), Move((3, 0), 2), Move((0, 3), 2)], 0.5))
