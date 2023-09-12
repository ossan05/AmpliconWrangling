using BioSequences

struct Move
    step::Int64
    score::Float64
end

# Convert NucleicAcid to integer A -> 1, C -> 2, G -> 3, T -> 4
function toInt(x::NucleicAcid)
    trailing_zeros(reinterpret(UInt8,x)) + 1
end

function simple_match_penalty_matrix(match_score, mismatch_score)
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
function nw_align(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move}) 
    nw_align(A, B, simple_match_penalty_matrix(match_score, mismatch_score), match_moves, vgap_moves, hgap_moves) 
end

#Needleman Wunsch alignment without affine scoring
function nw_align(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move})
    n, m = length(A), length(B)
    
    #Margins in the dp_matrix streamlines the code by avoiding boundschecking
    vertical_offset = maximum(k -> k.step, vcat(match_moves, vgap_moves)) + 1
    horizontal_offset = maximum(k -> k.step, vcat(match_moves, hgap_moves)) + 1

    #The matrix and the sequences are expanded according to the offsets
    vboundary = n + vertical_offset
    hboundary = m + horizontal_offset
    A2 = LongDNA{2}("A")^(vertical_offset - 1) * A
    B2 = LongDNA{2}("A")^(horizontal_offset - 1) * B

    #Initialize DP matrix
    #The cell at [x + vertical_offset, y + vertical_offset] is the score of the best alignment of A[1 : x] with B[1 : y]
    dp_matrix = fill(Inf64, vboundary, hboundary) 

    # Assign score 0 to the empty alignment
    dp_matrix[vertical_offset, horizontal_offset] = .0

    # first row
    for i ∈ 1 + horizontal_offset:hboundary, Move ∈ hgap_moves 
        if dp_matrix[vertical_offset, i - Move.step] + Move.score < dp_matrix[vertical_offset, i]
            dp_matrix[vertical_offset, i] = dp_matrix[vertical_offset, i - Move.step] + Move.score 
        end
    end

    # first column
    for i ∈ 1 + vertical_offset:vboundary, Move ∈ vgap_moves 
        if dp_matrix[i - Move.step, horizontal_offset] + Move.score < dp_matrix[i, horizontal_offset]
            dp_matrix[i, horizontal_offset] = dp_matrix[i - Move.step, horizontal_offset] + Move.score
        end
    end

    # whole dp matrix
    for j ∈ 1 + horizontal_offset:hboundary
        for i ∈ 1 + vertical_offset:vboundary
            # finding the lowest score move

            # the lowest score of vertical gaps is asigned to the current position in the matrix
            for k ∈ vgap_moves
                if dp_matrix[i - k.step, j] + k.score < dp_matrix[i, j]
                    dp_matrix[i, j] = dp_matrix[i - k.step, j] + k.score
                end
            end

            # the lowest score of horizontal gaps is asigned to the current position in the matrix if the score is lower than the current one
            for k ∈ hgap_moves
                if dp_matrix[i, j - k.step] + k.score < dp_matrix[i, j]
                    dp_matrix[i, j] = dp_matrix[i, j - k.step] + k.score
                end
            end

            for k ∈ match_moves
                # calculate mismatch score
                s = 0
                iB = j - k.step

                for iA ∈ i - k.step:i - 1
                    s += match_score_matrix[toInt(A2[iA]), toInt(B2[iB])]
                    iB += 1
                end 

                # asign a new value to the matrix if the move score is lower
                if dp_matrix[i - k.step, j - k.step] + k.score + s < dp_matrix[i, j]
                    dp_matrix[i, j] = dp_matrix[i - k.step, j - k.step] + k.score + s 
                end
            end
        end
    end

    # Backtracking
    x = vboundary
    y = hboundary
    res_A = LongDNA{4}("")
    res_B = LongDNA{4}("")


    while x > vertical_offset || y > horizontal_offset
        if x == vertical_offset # first row
            push!(res_A, DNA_Gap)
            push!(res_B, B2[y - 1])
            y -= 1
        elseif y == horizontal_offset # first column
            push!(res_A, A2[x - 1])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # iterate through vertical moves
            for k ∈ vgap_moves
                
                # check if the move leads to the current cell
                if dp_matrix[x, y] == dp_matrix[x - k.step, y] + k.score
                    for i ∈ 1:k.step
                        push!(res_A, A2[x - i])
                        push!(res_B, DNA_Gap)
                    end
                    x -= k.step
                    break
                end
            end

            # iterate through horizontal moves
            for k ∈ hgap_moves
                
                if dp_matrix[x, y] == dp_matrix[x, y - k.step] + k.score
                    for i ∈ 1:k.step
                        push!(res_A, DNA_Gap)
                        push!(res_B, B2[y - i])
                    end
                    y -= k.step
                    break
                end
            end

            # iterate through diagonal (match) moves
            for k ∈ match_moves
                
                #calculate total (mis-)match score from move
                s = sum(t -> match_score_matrix[toInt(A2[x-t]), toInt(B2[y-t])], 1 : k.step)
                
                # check if the move leads to the current cell
                if dp_matrix[x, y] == dp_matrix[x - k.step, y - k.step] + k.score + s
                    
                    # write the resulting sequences
                    for i ∈ 1:k.step
                        push!(res_A, A2[x - i])
                        push!(res_B, B2[y - i])
                    end
                    x -= k.step
                    y -= k.step
                    break
                end
            end
        end
    end
    return reverse(res_A), reverse(res_B)
end

#Needleman Wunsch alignment with affine scoring
function nw_align(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move}, extension_score::Float64) 
    nw_align(A, B, simple_match_penalty_matrix(match_score, mismatch_score), match_moves, vgap_moves, hgap_moves, extension_score) 
end

function nw_align(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move}, extension_score::Float64)

    n, m = length(A), length(B)

    if (extension_score < 0)
        # Do non-affine NW
        return nw_align(A, B, match_score_matrix, match_moves, vgap_moves, hgap_moves)
    end

    # Offset indeces to avoid bounds-checking
    vertical_offset = maximum(k -> k.step, vcat(match_moves, vgap_moves)) + 1
    horizontal_offset = maximum(k -> k.step, vcat(match_moves, hgap_moves)) + 1
    vboundary = n + vertical_offset
    hboundary = m + horizontal_offset

    # Length of sequences and matrices are increased according to offset
    A2 = LongDNA{2}("A")^(vertical_offset - 1) * A
    B2 = LongDNA{2}("A")^(horizontal_offset - 1) * B

    # Initialize DP matrix
    # The cell at [x + vertical_offset, y + vertical_offset] is the score of the best alignment of A[1 : x] with B[1 : y]
    dp_matrix = fill(Inf64, vboundary, hboundary)

    # Assign score 0 to the empty alignment
    dp_matrix[vertical_offset, horizontal_offset] = .0

    # Affine moves requires two extra DP matrices
    vaffine_matrix = fill(Inf64, vboundary, hboundary)
    haffine_matrix = fill(Inf64, vboundary, hboundary)

    # first row (optional speedup)
    for i ∈ 1 + horizontal_offset:hboundary
        for Move ∈ hgap_moves
            dp_matrix[vertical_offset, i] = min(
                dp_matrix[vertical_offset, i],
                haffine_matrix[vertical_offset, i - Move.step] + extension_score * Move.step,
                dp_matrix[vertical_offset, i - Move.step] + Move.score
            )
        end
    end

    # first column (optional speedup)
    for i ∈ 1 + vertical_offset:vboundary
        for Move ∈ vgap_moves
            dp_matrix[i, horizontal_offset] = min(
                dp_matrix[i, horizontal_offset],
                dp_matrix[i - Move.step, horizontal_offset] + extension_score * Move.step,
                dp_matrix[i - Move.step, horizontal_offset] + Move.score
            )
        end
    end

    #Main DP -step
    for j ∈ 1 + horizontal_offset : hboundary
        for i ∈ 1 + vertical_offset : vboundary

            # find the best diagonal move
            for k ∈ match_moves
                mismatch_sum = sum(t -> match_score_matrix[toInt(A2[i - t]), toInt(B2[j - t])], 1 : k.step)
                dp_matrix[i, j] = min(dp_matrix[i, j], dp_matrix[i - k.step, j - k.step] + k.score + mismatch_sum)
            end

            # finds the best vertical move
            for k ∈ vgap_moves
                # affine move
                vaffine_matrix[i, j] = min(
                    vaffine_matrix[i, j],
                    vaffine_matrix[i - k.step, j] + extension_score * k.step,
                    dp_matrix[i - k.step, j] + k.score
                )
            end

            # finds the best horizontal move
            for k ∈ hgap_moves
                # affine gap
                haffine_matrix[i, j] = min(
                    haffine_matrix[i, j],
                    haffine_matrix[i, j - k.step] + extension_score * k.step,
                    dp_matrix[i, j - k.step] + k.score
                )
            end

            #find best move overall
            dp_matrix[i, j] = min(dp_matrix[i, j], haffine_matrix[i, j], vaffine_matrix[i, j])
        end
    end

    # Backtracking
    res_A = LongDNA{4}("")
    res_B = LongDNA{4}("")
    
    # Start at the final cell
    x = vboundary
    y = hboundary
    
    # Flags for affine moves
    must_move_ver = false
    must_move_hor = false

    while x > vertical_offset || y > horizontal_offset
        if x == vertical_offset # first row
            push!(res_A, DNA_Gap)
            push!(res_B, B2[y - 1])
            y -= 1
        elseif y == horizontal_offset # first col
            push!(res_A, A2[x - 1])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # record previous position
            px = x
            py = y

            if !must_move_hor
    
                # iterate through vertical moves
                for k ∈ vgap_moves

                    # check if the move leads to the current cell
                    current_score = must_move_ver ? vaffine_matrix[x, y] : dp_matrix[x, y]
                    can_move_affine::Bool = (current_score == vaffine_matrix[x - k.step, y] + extension_score * k.step)
                    can_move_regular::Bool = (current_score == dp_matrix[x - k.step, y] + k.score)

                    if can_move_affine || can_move_regular
                        
                        # record the path
                        for i ∈ 1 : 2*k.step   
                            push!(res_A, A2[x - i])
                            push!(res_B, DNA_Gap)
                        end
                        x -= k.step

                        # constrain next move
                        must_move_ver = !can_move_regular
                        
                        break
                    end
                end
            end

            if !must_move_ver

                # iterate through horizontal Move moves
                for k ∈ hgap_moves
                    
                    # check if the move leads to the current cell
                    current_score = must_move_hor ? haffine_matrix[x, y] : dp_matrix[x, y]
                    can_move_affine::Bool = (current_score == haffine_matrix[x, y - k.step] + extension_score * k.step)
                    can_move_regular::Bool = (current_score == dp_matrix[x, y - k.step] + k.score)

                    if can_move_affine || can_move_regular
                        
                        # record the path
                        for i ∈ 1:k.step
                            push!(res_A, DNA_Gap)
                            push!(res_B, B2[y - i])
                        end
                        y -= k.step

                        # constrain next move
                        must_move_hor = !can_move_regular

                        break
                    end
                end
            end

            if !must_move_hor && !must_move_ver

                # iterate through digonal match moves
                for k ∈ match_moves

                    #calculate total (mis-)match score 
                    s = sum(t -> match_score_matrix[toInt(A2[x-t]), toInt(B2[y-t])], 1 : k.step)
                    
                    # check if the move leads to the current cell
                    if dp_matrix[x, y] == dp_matrix[x - k.step, y - k.step] + k.score + s
                        
                        # record the path
                        for i ∈ 1:k.step
                            push!(res_A, A2[x - i])
                            push!(res_B, B2[y - i])
                        end
                        x -= k.step
                        y -= k.step

                        break
                    end
                end
            end

            # if no move was found

            if px == x && py == y
                error("Backtracking failed")
            end
        end
    end
    return reverse(res_A), reverse(res_B)
end

# A = LongDNA{2}("CCGACCCCGATTCCCGTTA")
# B = LongDNA{2}("GACCTTCACCGTTA")
# match_moves = [Move(1, 0.0), Move(3, 0.0)]
# gap_moves = [Move(3, 2.0), Move(1, 1.0)]
# alignment = nw_align(A, B, .0, 0.5, match_moves, gap_moves, gap_moves)
# println(alignment[1])
# println(alignment[2])