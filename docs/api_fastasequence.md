# fastaparser.FastaSequence
Represents one single FASTA sequence.

## Parameters
The FastaSequence class can be instantiated with the following parameters
```Python
fastaparser.FastaSequence(sequence, id_='', description='', sequence_type=None, infer_type=False)
```

| Parameter | Type / Value | Default | Description|
|:---:|:---:|:---:|---|
| sequence | str | | String of characters representing a DNA, RNA or aminoacid sequence. Cannot be empty. **Must be provided** |
| id_ | str | '' | ID portion of the definition line (header). '>' and newlines will be removed, if any. Spaces will be converted to '_'. Can be an empty string. **Optional** |
| description | str | '' | Description portion of the definition line (header). Newlines will be removed, if any. Can be an empty string. **Optional** |
| sequence_type | 'nucleotide', 'aminoacid' or None | None | Indicates the sequence type. If not defined. **Optional** |
| infer_type | bool | False | Indicates if `FastaSequence` should try to infer aminoacid sequence type. If `True`, `FastaSequence` will analyse the whole sequence and, in the worst case scenario, can only identify aminoacid sequences. **Optional** |

#### Raises
**TypeError**

* If `sequence`, `id_`, `description`, `sequence_type` or `infer_type` are of the wrong type.

## Attributes
Instances of the FastaSequence class have the following attributes

| Attribute | Type / Value | Editable | Description |
|:---:|:---:|:---:|---|
| id | str | No | ID portion of the definition line (header). Can be empty |
| description | str | No | Description portion of the definition line (header). Can be empty |
| sequence | list(LetterCode) | No | Sequence |
| sequence_type | 'nucleotide', 'aminoacid' or None | Yes | Indicates the sequence type. Can be `None` if not known |
| inferred_type | bool | No | `True` if `FastaSequence` inferred the sequence type, `False` otherwise.

## Methods
Instances of the FastaSequence class have the following methods

### complement
Returns a new `FastaSequence` object containing the complementary sequence (ideally, of a nucleotide sequence).
Description is updated to mention the changes relative to the original sequence.

Non-nucleotide `LetterCodes` don't have a complement and, therefore, stay the same.
In order not to impose the setting of `sequence_type` as `'nucleotide'`, this method will work for any sequence and
`LetterCode` (as long as `sequence_type` is not `'aminoacid'`), which has the side effect of returning nonsensical
results when `LetterCodes` are not nucleotides.

Ex: For aminoacid letter codes that overlap with nucleotide letter codes, the output will be the complement of
the nucleotide represented by the same letter code, which makes no sense.

```Python
FastaSequence.complement(reverse=False)
```

| Parameter | Type / Value | Default | Description |
|:---:|:---:|:---:|---|
| reverse | bool | False | If sequence should be reversed. **Optional** |

#### Returns
**FastaSequence**

Complement of the current nucleotide `FastaSequence`. Non-nucleotide `LetterCodes` will stay the same.

#### Raises
**TypeError**

* If `sequence_type` is `'aminoacid'`.
* If `reverse` is not `bool`.

### gc_content
Calculates and returns the GC content of nucleotide sequence (as a ratio, by default).
Ignores degenerate letter codes besides S (G or C).
GC content is calculated the first time the method is called. Later calls will retrieve the same value.
GC content can also be calculated in `at_gc_ratio`.
If `sequence_type` is not `'nucleotide'` (or the sequence is not inherently a nucleotide sequence) the GC content
might be nonsensical.

```Python
FastaSequence.gc_content(as_percentage=False)
```

| Parameter | Type / Value | Default | Description |
|:---:|:---:|:---:|---|
| as_percentage | bool | False | Indicates whether the computed value should be returned as a percentage instead of the default ratio. **Optional** |

#### Returns
**float**

GC content of sequence.

#### Raises
**TypeError**

* If `sequence_type` is `'aminoacid'`.
* If `as_percentage` is not `bool`.

### at_gc_ratio
Calculates and returns the AT/GC ratio of nucleotide sequence.
Ignores degenerate letter codes besides W (A or T) and S (G or C).
AT/GC ratio is calculated the first time the method is called. Later calls will retrieve the same value.
Also uses previously calculated GC content or calculates and saves it if it hasn't been calculated yet.
If `sequence_type` is not `'nucleotide'` (or the sequence is not inherently a nucleotide sequence) the AT/GC ratio
might be nonsensical.

```Python
FastaSequence.at_gc_ratio()
```

#### Returns
**float**

AT/GC ratio of sequence.

#### Raises
**TypeError**

* If `sequence_type` is `'aminoacid'`.

### count_letter_codes
Returns a dictionary of letter code counts.
By default shows counts for all existing letter codes in the sequence,
but specific letter codes can be specified.

```Python
FastaSequence.count_letter_codes(letter_codes=None)
```

| Parameter | Type / Value | Default | Description |
|:---:|:---:|:---:|---|
| letter_codes | iterable or None | None | Iterable of all letter codes to count. **Optional** |

#### Returns
**dict**

Counts for every letter code in `letter_codes` or all letter codes in the sequence if `letter_codes` is not specified.

#### Raises
**TypeError**

* If `letter_codes` is neither an `iterable` or `None`.

### count_letter_codes_degenerate
Returns a dictionary of degenerate letter code counts. `sequence_type` must be explicitly defined.

```Python
FastaSequence.count_letter_codes_degenerate()
```

#### Returns
**dict**

Counts for every degenerate letter code in the sequence.

#### Raises
**TypeError**

* If `sequence_type` is not explicitly defined.

### formatted_definition_line
Returns a formatted FASTA definition line (header).

```Python
FastaSequence.formatted_definition_line()
```

#### Returns
**str**

FASTA definition line properly formatted.

### formatted_sequence
Formatted FASTA sequence (only the sequence, without the definition line). Lines are separated by '\n'.

```Python
FastaSequence.formatted_sequence(max_characters_per_line=70)
```

| Parameter | Type / Value | Default | Description |
|:---:|:---:|:---:|---|
| max_characters_per_line | int | 70 | Maximum number of characters per line. This value should not go above 80, as per the FASTA specification. A very low value is also not recommended. **Optional** |

#### Returns
**str**

Returns a FASTA sequence properly formatted.

#### Raises
**TypeError**

* If `max_characters_per_line` is not an `int`.

### formatted_fasta
Returns a formatted FASTA (definition line and sequence).

```Python
FastaSequence.formatted_fasta()
```

#### Returns
**str**

FASTA properly formatted.

### sequence_as_string
Returns the sequence as string. Converts the list of `LetterCode` objects to a single string.

```Python
FastaSequence.sequence_as_string()
```

#### Returns
**str**

Sequence as string.

### reverse
Iterates over the sequence in reverse order (same as calling `reversed()` on a `FastaSequence` object).
Returns a new reverse `iterator` of the sequence every time `reverse` is called.

```Python
FastaSequence.reverse()
```

#### Returns
**iterator**

Iterator over the reversed sequence.

## Class Methods
The FastaSequence class has the following class method

### from_fastasequence
Initializes with the given `FastaSequence` object (alternate `__init__` method).

```Python
FastaSequence.from_fastasequence(fastasequence)
```

| Parameter | Type / Value | Default | Description |
|:---:|:---:|:---:|---|
| fastasequence | FastaSequence | | `FastaSequence` object. **Must be provided** |

#### Returns
**FastaSequence**

Copy of `fastasequence` (`FastaSequence` object).

#### Raises
**TypeError**

* If `fastasequence` is not a `FastaSequence`.

## Special Methods
* \_\_iter__
* \_\_reversed__
* \_\_next__
* \_\_getitem__
* \_\_eq__
* \_\_len__
* \_\_repr__
* \_\_str__
