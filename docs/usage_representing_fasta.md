# Representing FASTA sequences
With FastaParser, sequences can be represented by using the [`FastaSequence`](api_fastasequence.md) class.

A [`FastaSequence`](api_fastasequence.md) object takes a string of nucleotide/aminoacid letter codes as first parameter.
This would represent a sequence without a header. Optionally the id and description of the sequence can also be provided,
as the second and third parameter, respectively, as `_id` and `description`.

The type of sequence can also be defined, as the fourth parameter `sequence_type`, which is `None` by default, and can be
defined as either `'nucleotide'` or `'aminoacid'`. If this parameter is not provided, the type of sequence can be inferred
by passing the last parameter, `infer_type`, as `True`.

```python
FastaSequence(sequence, id_, description, sequence_type, infer_type)
```

[`FastaSequence`](api_fastasequence.md) objects have many features, like calculating the AC/TG ratio, generating
complementary sequences, iterating over the sequence and more
(for more details, go to the [API Specification page](api_fastasequence.md)).

# Representing letter codes
Individual letter codes can also be represented with FastaParser by using the [`LetterCode`](api_lettercode.md) class.

A [`LetterCode`](api_lettercode.md) object takes a single character string of nucleotide/aminoacid letter code as
first parameter.

The letter code type can also be defined with the second parameter `letter_type`,
which is `None` by default, and can be defined as either `'nucleotide'` or `'aminoacid'`.
It is important to define the type of letter code as this allows the LetterCode object to associate the name and/or
description of the letter code, as well as knowing whether the letter code is degenerate or not.

```python
LetterCode(letter_code, letter_type)
```
