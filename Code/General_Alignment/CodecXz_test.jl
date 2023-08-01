using CodecXz
using BioSequences, FASTX

#Because the old one is broken by a FASTX update
function read_fasta(filename)
    stream = XzDecompressorStream(open(filename))
    decompressed = FASTA.Reader(stream)
    seqs = [String(sequence(record)) for record in decompressed]
    close(stream)
    return seqs
end

# println(read_fasta("sequences.fasta.xz"))
seqs = read_fasta("sequences.fasta.xz")
length(seqs[1])