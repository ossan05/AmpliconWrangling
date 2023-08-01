using BioSequences
include("bad_general_pairwise_alignment.jl")

# function initiate_general_pairwise_aligner(seq1::LongDNA{2}, seq2::LongDNA{2})
#     general_pairwise_aligner(seq1, seq2, .0, 0.5, [Move((1, 1), 0), Move((1, 0), 1), Move((0, 1), 1), Move((3, 3), 0), Move((3, 0), 2), Move((0, 3), 2)])
# end
# function toInt(x::NucleicAcid)
#     trailing_zeros(reinterpret(UInt8,x))+1
# end

# function simple_match_score_matrix(mismatch_score)
#     match_score_matrix = zeros(4, 4)
#     for i in 1:4, j in 1:4
#         if i != j
#             match_score_matrix[i, j] = mismatch_score
#         end
#     end
#     match_score_matrix
# end

# # match and mismatch matrix
# function nw_aligner(A::LongDNA{2}, B::LongDNA{2}, mismatch_score, reg_moves::Array{Move}, gap_moves_a::Array{Move}, gap_moves_b::Array{Move})
#     nw_aligner(A, B, simple_match_score_matrix(mismatch_score), reg_moves, gap_moves_a, gap_moves_b) 
# end

# function find_best_move(A::LongDNA{4}, B::LongDNA{4}, match_score_matrix::Array{Float64, 2}, reg_moves::Array{Move}, gap_moves_a::Array{Move}, gap_moves_b::Array{Move}, dp_matrix::Array{Float64, 2}, x::Int64, y::Int64)
    
#     best_move_a = findmin(map(k -> dp_matrix[x - k.step, y] + k.score, gap_moves_a))
#     best_move_b = findmin(map(k -> dp_matrix[x, y - k.step] + k.score, gap_moves_b))
#     best_reg_move = findmin(map(k -> dp_matrix[x - k.step, y - k.step] + k.score + sum([match_score_matrix[toInt(A[x - t]), toInt(B[y - t])] for t in 1 : k.step]), reg_moves))
    
#     scores = [best_move_a[1], best_move_b[1], best_reg_move[1]]
#     moves = [gap_moves_a[best_move_a[2]], gap_moves_b[best_move_b[2]], reg_moves[best_reg_move[2]]]

#     move_type = argmin(scores)
#     return moves[move_type], scores[move_type], move_type 
# end

# # general_pairwise_aligner makes a match/mismatch matrix and takes the tuple with the moves and scores, 
# # making them objects of type Move or Gap and puting them into lists 
# function nw_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, reg_moves::Array{Move}, gap_moves_a::Array{Move}, gap_moves_b::Array{Move}) 
#     # Needleman-Wunsch alignment
#     max_move_a = max([k.step for k in gap_moves_a]..., [k.step for k in reg_moves]...)
#     max_move_b = max([k.step for k in gap_moves_b]..., [k.step for k in reg_moves]...)
    
#     A2 = repeat(LongDNA{4}([DNA_N]), max_move_a + 1) * A
#     B2 = repeat(LongDNA{4}([DNA_N]), max_move_b + 1) * B
    
#     m, n = length(A2), length(B2)

#     dp_matrix = fill(Inf64, m, n)
#     dp_matrix[max_move_a + 1, max_move_b + 1] = .0

#     i = max_move_a + 1
#     j = max_move_b + 2
#     while i <= m
#         dp_matrix[i, j] = find_best_move(A2, B2, match_score_matrix, reg_moves, gap_moves_a, gap_moves_b, dp_matrix, i, j)[2]
#         j += 1
#         if j > n
#             j = max_move_b + 1
#             i += 1
#         end
#     end

#     x = m; y = n
#     res_A = dna""; res_B = dna""

#     while x > max_move_a + 1 || y > max_move_b + 1
#         best_move, best_score, move_type = find_best_move(A2, B2, match_score_matrix, reg_moves, gap_moves_a, gap_moves_b, dp_matrix, x, y)
#         if move_type == 1 #gap move a
#             for i ∈ 1 : best_move.step
#                 push!(res_A, A2[x - i])
#                 push!(res_B, DNA_Gap)
#             end
#             x -= best_move.step
#         elseif move_type == 2 #gap move b
#             for i ∈ 1 : best_move.step
#                 push!(res_A, DNA_Gap)
#                 push!(res_B, B2[y - i])
#             end
#             y -= best_move.step
#         elseif move_type == 3 #regular move
#             for i ∈ 1 : best_move.step
#                 push!(res_A, A2[x - i])
#                 push!(res_B, B2[y - i])
#             end
#             x -= best_move.step
#             y -= best_move.step
#         end
#     end

#     return reverse(res_A), reverse(res_B)
# end

function alignment_score(A::LongDNA{4}, B::LongDNA{4}, match_score, mismatch_score, match_moves::Array{Move}, vgap_moves::Array{Move}, hgap_moves::Array{Move}, affine_score = -1)
    return alignment_score(A, B, make_match_score_matrix(match_score, mismatch_score), match_moves, vgap_moves, hgap_moves, affine_score)
end

function alignment_score(A::LongDNA{4}, B::LongDNA{4}, match_score_matrix::Array{Float64, 2}, match_moves::Array{Move}, vgap_moves::Array{Move}, hgap_moves::Array{Move}, affine_score = -1)
    if length(A) != length(B)
        error("Aligned sequences must be of equal lengths.")
    end
    
    max_move = maximum([k.step for k in vcat(vgap_moves, hgap_moves, match_moves)])
    
    n = length(A)
    
    no_gaps = filter(i -> A[i] != DNA_Gap && B[i] != DNA_Gap, 1 : n)
    match_score_total = sum(i -> match_score_matrix[toInt(A[i]), toInt(B[i])], no_gaps)

    A2 = repeat(LongDNA{4}([DNA_N]), max_move) * A
    B2 = repeat(LongDNA{4}([DNA_N]), max_move) * B

    dp = fill(Inf64, n+max_move)
    dp[max_move] = 0
    for i in max_move + 1 : n + max_move
        possible_match_moves = filter(k -> !any(t -> A2[i - t] == DNA_Gap || B2[i - t] == DNA_Gap, 0 : k.step - 1), match_moves)
        possible_vgap_moves = filter(k -> !any(t -> A2[i - t] != DNA_Gap, 0 : k.step - 1), vgap_moves)
        possible_hgap_moves = filter(k -> !any(t -> B2[i - t] != DNA_Gap, 0 : k.step - 1), hgap_moves)
        possible_affine_moves = Vector{Move}()
        if affine_score >= 0 && (A2[i] == A2[i - 1] == DNA_Gap || B2[i] == B2[i - 1] == DNA_Gap)
            push!(possible_affine_moves, Move(1, affine_score))
        end
        dp[i] = minimum(k -> dp[i - k.step] + k.score, vcat(possible_match_moves, possible_vgap_moves, possible_hgap_moves, possible_affine_moves))
    end
    return dp[n + max_move] + match_score_total
end


# A = LongDNA{4}("CGG-G---")
# B = LongDNA{4}("-ACCGCTG")
# reg_moves = [Move(1, 0)]
# gap_moves = [Move(3, 2), Move(1, 1)]
# alignment = (A, B)#general_pairwise_aligner(LongDNA{2}("TTCGACTG"), LongDNA{2}("TACGACGACTG"), .0, 0.5, [Move((1, 1), 0), Move((1, 0), 1), Move((0, 1), 1), Move((3, 3), 0), Move((3, 0), 2), Move((0, 3), 2)], 0.5)
# affine_score = 0.5
# # println(alignment[1])
# # println(alignment[2])
# # @show alignment_score(alignment[1], alignment[2], 0.0, 0.5, reg_moves, gap_moves, gap_moves, affine_score)