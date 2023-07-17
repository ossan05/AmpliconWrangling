nucleotides = ['A', 'C', 'G', 'T']

Dictionary = Dict(
    "GCT"=> "A", "GCC"=> "A", "GCA"=> "A", "GCG" => "A",
    "CGT"=> "R", "CGC"=> "R", "CGA"=> "R", "CGG"=> "R", "AGA"=> "R", "AGG"=> "R",
    "AAT"=> "N", "AAC"=> "N",
    "GAT"=> "D", "GAC"=> "D",
    "TGT"=> "C", "TGC"=> "C",
    "CAA"=> "Q", "CAG"=> "Q",
    "GAA"=> "E", "GAG"=> "E",
    "GGT"=> "G", "GGC"=> "G", "GGA"=> "G", "GGG"=> "G",
    "CAT"=> "H", "CAC"=> "H",
    "ATT"=> "I", "ATC"=> "I", "ATA"=> "I",
    "CTT"=> "L", "CTC"=> "L", "CTA"=> "L", "CTG"=> "L", "TTA"=> "L", "TTG"=> "L",
    "AAA"=> "K", "AAG"=> "K",
    "ATG" => "M",
    "TTT"=> "F", "TTC"=> "F",
    "CCT"=> "P", "CCC"=> "P", "CCA"=> "P", "CCG"=> "P",
    "TCT"=> "S", "TCC"=> "S", "TCA"=> "S", "TCG"=> "S", "AGT"=> "S", "AGC"=> "S",
    "ACT"=> "T", "ACC"=> "T", "ACA"=> "T", "ACG"=> "T",
    "TGG" => "W",
    "TAT"=> "Y", "TAC"=> "Y",
    "GTT"=> "V", "GTC"=> "V", "GTA"=> "V", "GTG"=> "V",
    "TAA" => ' ', "TAG" => ' ', "TGA" => ' ')

occurences1 = 0

for key in keys(Dictionary)
    for i in nucleotides
        global occurences1
        if Dictionary[key[2:3] * i] == Dictionary[key]
            occurences1 += 1
        end
    end
end

println((256 - occurences1) / (256))
println(occurences1 / (256))

occurences2 = 0

for key in keys(Dictionary)
    for i in nucleotides
        global occurences2
        if Dictionary[i * key[1:2]] == Dictionary[key]
            occurences2 += 1
        end
    end
end

println((256 - occurences2) / (256))
println(occurences2 / (256))

occurences2 = 0

for key in keys(Dictionary)
    for i in nucleotides
        global occurences2
        if Dictionary[key[1] * key[3] * i] == Dictionary[key]
            occurences2 += 1
        end
    end
end

println((256 - occurences2) / (256))
# 0.80078125

println(occurences2 / (256))







# + (length(seq) - current_indel_pos) / 3 * 0.9453125 * codon_penalty

# - (length(seq) - current_indel_pos) / 3 * 0.9453125 * codon_penalty