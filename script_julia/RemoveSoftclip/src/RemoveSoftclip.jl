module RemoveSoftclip

export bam_rm_softclip
export auxtype 
export bam_tag 
export build_rm_cl_sam_record

import BioAlignments; const bioalign = BioAlignments
import XAM.BAM; const bam = BAM
import XAM.SAM; const sam = SAM
using Fmt

is_skipop(op) = bioalign.OP_SKIP == op
is_softclipop(op) = bioalign.OP_SOFT_CLIP == op

function extract_op_len(x::bam.Record, is_op::Function)::Vector{Int64}
    bam.cigar_rle(x) |> x -> x[2][is_op.(x[1])]
end

function extract_left_softclip(x::bam.Record)::Int64
    bam.cigar_rle(x) |> c -> is_softclipop(c[1][1]) ? c[2][1] : 0
end

function extract_right_softclip(x::bam.Record)::Int64
    bam.cigar_rle(x) |> c -> is_softclipop(c[1][end]) ? c[2][end] : 0
end

function extract_softclip(x::bam.Record)::Tuple{Int64,Int64}
    bam.cigar_rle(x) |>
        c -> (is_softclipop(c[1][1]) ? c[2][1] : 0, 
              is_softclipop(c[1][end]) ? c[2][end] : 0)
end

bam_qname(r::bam.Record)::String = bam.tempname(r)
bam_flag(r::bam.Record)::UInt16 = bam.flag(r)
bam_refname(r::bam.Record)::String = bam.refname(r)
bam_position(r::bam.Record)::Int64 = bam.position(r)
bam_mapq(r::bam.Record)::Int = bam.mappingquality(r)
bam_cigar(r::bam.Record)::String = bam.cigar(r)
bam_sequence(r::bam.Record)::String = String(bam.sequence(r))
bam_quality(r::bam.Record)::String = String((Char).(bam.quality(r) .+ 0x21))

"""
    auxtype(b) -> (Int8, Char, DataType)

Interpret a byte as the tag data type. Return a tuple containing "data size (byte)", "tag type character" and "`DataType` of tag".

# Arguments
- `b::UInt8`: A byte data encoding the auxiliary tag data type
"""
function auxtype(b::UInt8)::Tuple{Int8,Char,DataType}
    return (
        b == UInt8('A') ? (4, 'A', Char) :
        b == UInt8('c') ? (1, 'i', Int8) :
        b == UInt8('C') ? (1, 'i', UInt8) :
        b == UInt8('s') ? (2, 'i', Int16) :
        b == UInt8('S') ? (2, 'i', UInt16) :
        b == UInt8('i') ? (4, 'i', Int32) :
        b == UInt8('I') ? (4, 'i', UInt32) :
        b == UInt8('f') ? (4, 'f', Float32) :
        b == UInt8('Z') ? (1, 'Z', String) : 
        error("invalid type tag: '$(Char(b))'")
    )
end

"""
    bam_tag(r) -> String

Construct and return the auxiliary tag string of the bam record.

# Arguments
- `r::XAM.BAM.Record`: A BAM record
"""
function bam_tag(r::bam.Record)::String
    array = view(r.data, bam.auxdata_position(r):bam.data_size(r))
    len_array = length(array)
    ret = Vector{Char}()

    i = 0
    while i < len_array
        append!(ret, '\t', Char(array[i + 1]), Char(array[i + 2]), ':')
        s = auxtype(array[i + 3])
        append!(ret, s[2])
        append!(ret, ':')
        append!(ret, Fmt.format(f"{}", reinterpret(s[3], array[i + 4:i + 3 + s[1]])[1]))
        i += 3 + s[1]
    end
    String(ret)
end

"""
    build_rm_cl_sam_record(r) -> String

Remove soft-clips and then create sam record from bam record

# Arguments
- `r::XAM.BAM.Record`: A BAM record
"""
function build_rm_cl_sam_record(r::bam.Record)::String
    s = extract_softclip(r)
    Fmt.format(
        f"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}\tLS:i:{}\tRS:i:{}\n",
        bam.tempname(r),
        bam.flag(r),
        bam.refname(r),
        bam.position(r),
        bam.mappingquality(r),
        bam_cigar(r)[(s[1] > 0 ? Fmt.ndigits_decimal(s[1]) + 2 : 1):(s[2] > 0 ? end - (Fmt.ndigits_decimal(s[2]) + 1) : end)],
        "*",
        "0",
        "0",
        bam_sequence(r)[(s[1] + 1):(end - s[2])],
        bam_quality(r)[(s[1] + 1):(end - s[2])],
        bam_tag(r),
        s[1],
        s[2]
    )
end

"""
    bam_rm_softclip(path_in, path_out)

Remove softclips from a bam file

# Example
```
bam_rm_softclip("path/to/input.bam", "path/to/output.bam")
```

# Arguments
- `path_in::String`: a path to input bam file
- `path_out::String`: a path to output sam file
"""
function bam_rm_softclip(path_in::String, path_out::String)
    reader = open(bam.Reader, path_in)
    b_record = bam.Record()
    writer = sam.Writer(open(path_out, "w"), reader.header)
    while !eof(reader)
    # i = 0
    # @time while i < 1000000
    #     i += 1
        empty!(b_record)
        read!(reader, b_record)
        v = Vector{UInt8}(build_rm_cl_sam_record(b_record))
        unsafe_write(writer.stream, pointer(v, 1), length(v))
    end
    close(writer)
    close(reader)
end

end # module