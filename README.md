# CountESS MiniMap2 Plugin v0.0.10

This is a plugin to allow [MiniMap2](https://github.com/lh3/minimap2) 
to be used from within [CountESS](https://github.com/CountESS-Project/CountESS/)

This might be useful, but it also stands as a handy example of how to write
a CountESS plugin which wraps an external Python library.

## Parameters

### Output Column Prefix

The default output column prefix is `mm` but this can be changed.

### Ref FA / Ref MMI

A local FASTA or MMI file to look up the sequences in.
Acceptable formats are MMI or FASTA, optionally gzipped.

### Req Sequence

Alternatively, enter a single reference DNA sequence directly here.

### Preset

See [minimap2 Preset Options](https://lh3.github.io/minimap2/minimap2.html#8)

### Minimum Match Length

Reject matches with an overall length (`mm_r_en - mm_r_st`) less than this.

### Drop Unmatched

Rows with no matches will be dropped.

## Output Columns

### `mm_ctg`

The matched location.

### `mm_r_st`, `mm_r_en`

Start and end positions of the matched sequence.

### `mm_strand`

`1` if query/target on the same strand.
`-1` if opposite.
`0` if unknown.

### `mm_cigar`

CIGAR string expressing the differences within the match.  Note that not all differences are
encoded in the CIGAR string, for a more detailed description of what variations are present see
the next two fields.

### `mm_cs`

A difference string as documented in the [Minimap2 Output Format](https://lh3.github.io/minimap2/minimap2.html#10) section.

### `mm_hgvs`

The difference string translated into HGVS format relative to the matched location.

## License

BSD 3-clause.  See [LICENSE.txt](LICENSE.txt)
