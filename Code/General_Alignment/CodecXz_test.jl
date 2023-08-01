using CodecXz
stream = XzDecompressorStream(open("sequences.fasta.xz"))
println((stream))
println((stream))