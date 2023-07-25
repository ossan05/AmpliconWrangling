using BioSequences

struct Move
    step::Tuple{Int64, Int64}
    score::Float64
end

struct Gap
    step::Int64
    score::Float64
end

function initiate(seq1::LongDNA{2}, seq2::LongDNA{2})
    general_pairwise_aligner(seq1, seq2, .0, 0.5, [Move((1, 1), 0), Move((1, 0), 1), Move((0, 1), 1), Move((3, 3), 0), Move((3, 0), 2), Move((0, 3), 2)], 0.5)
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
function general_pairwise_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, moves::Array{Move}) 
    general_pairwise_aligner(A, B, make_match_score_matrix(match_score, mismatch_score), moves) 
end

# general_pairwise_aligner makes a match/mismatch matrix and takes the tuple with the moves and scores, 
# making them objects of type Move or Gap and puting them into lists 
function general_pairwise_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, moves::Array{Move}) 
    # sort the moves into gaps
    # moves are all the moves
    # v - vertical, h - horizontal
    matches = Vector{Move}()
    vgaps = Vector{Gap}()
    hgaps = Vector{Gap}()

    # fill the vectors with gaps and moves
    for i ∈ moves
        if i.step[2] == 0
            push!(vgaps, Gap(i.step[1], i.score))
        elseif i.step[1] == 0
            push!(hgaps, Gap(i.step[2], i.score))
        else
            push!(matches, i)
        end
    end
    
    # make the alignment

    # Needleman-Wunsch alignment

    n, m = length(A), length(B)

    dp_matrix = fill(typemax(Float64), n + 1, m + 1)
    dp_matrix[1, 1] = .0

    # first row
    for i ∈ 2:m + 1, gap ∈ hgaps 
        if gap.step < i && dp_matrix[1, i - gap.step] + gap.score < dp_matrix[1, i]
            dp_matrix[1, i] = dp_matrix[1, i - gap.step] + gap.score 
        end
    end

    # first column
    for i ∈ 2:n + 1, gap ∈ vgaps 
        if gap.step < i && dp_matrix[i - gap.step, 1] + gap.score < dp_matrix[i, 1]
            dp_matrix[i, 1] = dp_matrix[i - gap.step, 1] + gap.score
        end
    end

    # whole dp matrix
    for i ∈ 2:n + 1
        for j ∈ 2:m + 1

            # finding the lowest score move

            # the lowest score of vertical gaps is asigned to the current position in the matrix
            for k ∈ vgaps

                # filters out moves that are too big
                if k.step < i
                    if dp_matrix[i - k.step, j] + k.score < dp_matrix[i, j]
                        dp_matrix[i, j] = dp_matrix[i - k.step, j] + k.score
                    end
                end
            end

            # the lowest score of horizontal gaps is asigned to the current position in the matrix if the score is lower than the current one
            for k ∈ hgaps

                # filters out moves that are too big
                if k.step < j
                    if dp_matrix[i, j - k.step] + k.score < dp_matrix[i, j]
                        dp_matrix[i, j] = dp_matrix[i, j - k.step] + k.score
                    end
                end
            end

            for k ∈ matches
                if k.step[1] < i && k.step[2] < j

                    # calculate mismatch score
                    s = 0
                    for l ∈ i - k.step[1]:i - 1
                        s += match_score_matrix[toInt(A[l]), toInt(B[l + j - i])]
                    end

                    # asign a new value to the matrix if the move score is lower
                    if dp_matrix[i - k.step[1], j - k.step[2]] + k.score + s < dp_matrix[i, j]
                        dp_matrix[i, j] = dp_matrix[i - k.step[1], j - k.step[2]] + k.score + s 
                    end
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
            res_A *= DNA_Gap
            push!(res_B, B[y - 1])
            y -= 1
        elseif y == 1
            push!(res_A, A[x - 1])
            push!(res_B, DNAGap)
            x -= 1
        else
            # iterate through vertical gap moves
            for k ∈ vgaps
                # check that the move isn't out of bounds
                if k.step < x
                    # opening gap
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

            # iterate through horizontal gap moves
            for k ∈ hgaps
                # check that the move isn't out of bounds
                if k.step < y
                    # opening gap
                    # check if the move lead to the current cell
                    if dp_matrix[x, y] == dp_matrix[x, y - k.step] + k.score
                        for i ∈ 1:k.step
                            res_A *= '_'
                            res_B *= string(B[y - i])
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
                if k.step[1] < x && k.step[2] < y

                    # calculate score of mismatches
                    s = 0
                    for i ∈ x - k.step[1]:x - 1
                        s += match_score_matrix[toInt(A[i]), toInt(B[i + y - x])]
                    end

                    # check if the move lead to the current cell
                    if dp_matrix[x, y] == dp_matrix[x - k.step[1], y - k.step[2]] + k.score + s
                        
                        # write the resulting sequences
                        # k.step[1] and k.step[2] is the same
                        for i ∈ 1:k.step[1]
                            push!(res_A, A[x - i])
                            push!(res_B, B[y - i])
                        end
                        x -= k.step[1]
                        y -= k.step[2]

                        # break to stop iterating through moves
                        break
                    end
                end
            end
        end
    end

    return reverse(res_A), reverse(res_B)
end

# Another method for affine gap penalties
function general_pairwise_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score::Float64, mismatch_score::Float64, moves::Array{Move}, affine_gap::Float64) 
    
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
function general_pairwise_aligner(A::LongDNA{2}, B::LongDNA{2}, match_score_matrix::Array{Float64, 2}, moves::Array{Move}, affine_gap::Float64) 

    # sort the moves into gaps 
    matches = Vector{Move}()
    vgaps = Vector{Gap}()
    hgaps = Vector{Gap}()

    for i ∈ moves
        if i.step[2] == 0
            push!(vgaps, Gap(i.step[1], i.score))
        elseif i.step[1] == 0
            push!(hgaps, Gap(i.step[2], i.score))
        else
            push!(matches, i)
        end
    end

    # Needleman-Wunsch alignment

    n, m = length(A), length(B)

    dp_matrix = fill(typemax(Float64), n + 1, m + 1)
    dp_matrix[1, 1] = .0
    vaffine_matrix = zeros(Float64, n + 1, m + 1)
    haffine_matrix = zeros(Float64, n + 1, m + 1)


    # first row
    for i ∈ 2:m + 1, gap ∈ hgaps 
        if gap.step + 1 < i && dp_matrix[1, i - gap.step] + affine_gap * gap.step < dp_matrix[1, i]
            dp_matrix[1, i] = dp_matrix[1, i - gap.step] + affine_gap * gap.step
            haffine_matrix[1, i] = dp_matrix[1, i]
            vaffine_matrix[1, i] = Inf
        elseif gap.step < i && dp_matrix[1, i - gap.step] + gap.score < dp_matrix[1, i]
            dp_matrix[1, i] = dp_matrix[1, i - gap.step] + gap.score 
            haffine_matrix[1, i] = dp_matrix[1, i]
            vaffine_matrix[1, i] = Inf
        end
    end

    # first column
    for i ∈ 2:n + 1, gap ∈ vgaps 
        if gap.step + 1 < i && dp_matrix[i - gap.step, 1] + affine_gap * gap.step < dp_matrix[i, 1]
            dp_matrix[i, 1] = dp_matrix[i - gap.step, 1] + affine_gap * gap.step
            vaffine_matrix[i, 1] = dp_matrix[i, 1]
            haffine_matrix[i, 1] = Inf
        elseif gap.step < i && dp_matrix[i - gap.step, 1] + gap.score < dp_matrix[i, 1]
            dp_matrix[i, 1] = dp_matrix[i - gap.step, 1] + gap.score
            haffine_matrix[i, 1] = Inf
        end
    end

    for i ∈ 2:n + 1
        for j ∈ 2:m + 1

            # finds the best vertical move
            vaffine_matrix[i, j] = Inf

            for k ∈ vgaps
                if k.step < i
                    # affine gap
                    if vaffine_matrix[i - k.step, j] + affine_gap * k.step < vaffine_matrix[i, j]
                        vaffine_matrix[i, j] = vaffine_matrix[i - k.step, j] + affine_gap * k.step
                    end
                    # opening gap
                    if dp_matrix[i - k.step, j] + k.score < vaffine_matrix[i, j]
                        vaffine_matrix[i, j] = dp_matrix[i - k.step, j] + k.score
                    end
                end
            end

            # finds the best horizontal move
            haffine_matrix[i, j] = Inf

            for k ∈ hgaps
                if k.step < j
                    # affine_gap
                    if haffine_matrix[i, j - k.step] + affine_gap * k.step < haffine_matrix[i, j]
                        haffine_matrix[i, j] = haffine_matrix[i, j - k.step] + affine_gap * k.step
                    end
                    # opening gap
                    if dp_matrix[i, j - k.step] + k.score < haffine_matrix[i, j]
                        haffine_matrix[i, j] = dp_matrix[i, j - k.step] + k.score
                    end
                end
            end

            # the fastest diagonal move

            dp_matrix[i, j] = Inf

            for k ∈ matches
                if k.step[1] < i && k.step[2] < j
                    s = 0
                    for l ∈ i - k.step[1]:i - 1
                        s += match_score_matrix[toInt(A[l]), toInt(B[l + j - i])]
                    end 

                    if dp_matrix[i - k.step[1], j - k.step[2]] + k.score + s < dp_matrix[i, j]
                        dp_matrix[i, j] = dp_matrix[i - k.step[1], j - k.step[2]] + k.score + s 
                    end
                end
            end

            # fastest move to dp_matrix
            dp_matrix[i, j] = min(dp_matrix[i, j], haffine_matrix[i, j], vaffine_matrix[i, j])
        end
    end

    # backtracking
    y = m + 1
    x = n + 1
    res_A = LongDNA{4}("")
    res_B = LongDNA{4}("")

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
            # iterate through vertical gap moves
            for k ∈ vgaps
                # check that the move isn't out of bounds
                if k.step < x
                    # affine gap
                    # check if the affine move lead to the current cell
                    if dp_matrix[x, y] == vaffine_matrix[x - k.step, y] + affine_gap * k.step
                        
                        # write the result
                        for i ∈ 1:k.step
                            push!(res_A, A[x - i])
                            push!(res_B, DNA_Gap)
                        end
                        x -= k.step
                        # how long the affine gap is
                        while x > k.step && vaffine_matrix[x, y] == vaffine_matrix[x - k.step, y] + affine_gap * k.step
                            for i ∈ 1:k.step
                                push!(res_A, A[x - i])
                                push!(res_B, DNA_Gap)
                            end
                            x -= k.step
                        end

                        # include the opening gap
                        for i ∈ 1:k.step
                            push!(res_A, A[x - i])
                            push!(res_B, DNA_Gap)
                        end
                        x -= k.step
                        # break to stop iterating through moves
                        break

                    # opening gap
                    # check if the move lead to the current cell
                    elseif dp_matrix[x, y] == dp_matrix[x - k.step, y] + k.score
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

            # iterate through horizontal gap moves
            for k ∈ hgaps
                # check that the move isn't out of bounds
                if k.step < y

                    # affine gap
                    # check if the affine move lead to the current cell
                    if dp_matrix[x, y] == haffine_matrix[x, y - k.step] + affine_gap * k.step
                        
                        # write the resulting sequences
                        for i ∈ 1:k.step
                            push!(res_A, DNA_Gap)
                            push!(res_B, B[y - i])
                        end
                        y -= k.step

                        # how long the affine gap is
                        while y > k.step && haffine_matrix[x, y] == haffine_matrix[x, y - k.step] + affine_gap * k.step
                            for i ∈ 1:k.step
                                push!(res_A, DNA_Gap)
                                push!(res_B, B[y - i])
                            end
                            y -= k.step
                        end
                        for i ∈ 1:k.step
                            push!(res_A, DNA_Gap)
                            push!(res_B, B[y - i])
                        end
                        y -= k.step
                        break

                    # opening gap
                    # check if the move lead to the current cell
                    elseif dp_matrix[x, y] == dp_matrix[x, y - k.step] + k.score
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
                if k.step[1] < x && k.step[2] < y

                    # calculate score of mismatches
                    s = 0
                    for i ∈ x - k.step[1]:x - 1
                        s += match_score_matrix[toInt(A[i]), toInt(B[i + y - x])]
                    end

                    # check if the move lead to the current cell
                    if dp_matrix[x, y] == dp_matrix[x - k.step[1], y - k.step[2]] + k.score + s
                        
                        # write the resulting sequences
                        # k.step[1] and k.step[2] is the same
                        for i ∈ 1:k.step[1]
                            push!(res_A, A[x - i])
                            push!(res_B, B[y - i])
                        end
                        x -= k.step[1]
                        y -= k.step[2]

                        # break to stop iterating through moves
                        break
                    end
                end
            end
        end
    end

    return reverse(res_A), reverse(res_B)
end

#print(general_pairwise_aligner(LongDNA{2}("TTCGACTG"), LongDNA{2}("TACGACGACTG"), .0, 0.5, [Move((1, 1), 0), Move((1, 0), 1), Move((0, 1), 1), Move((3, 3), 0), Move((3, 0), 2), Move((0, 3), 2)], 0.5))
