# fastaparser.Reader
Parser/Reader for a given FASTA file.
Iterates over the FASTA file using one of two parsing mechanisms:

* **rich**:
Returns `FastaSequence` objects (default). Slower, but feature rich.
* **quick**:
Generates objects containing just the FASTA `header` and `sequence` attributes
for each sequence in the FASTA file.
Parses FASTA files faster but lacks some features.

## Parameters
The Reader class can be instantiated with the following parameters

| Parameter | Type / Value | Default | Description|
|---|---|---|---|
| fasta_file | file object | | An opened file handle (for reading). **Must be provided** |
| sequences_type | 'nucleotide', 'aminoacid' or None | None | Indicates the type of sequences to expect. None if unknown. **Optional** |
| infer_type | bool | False | Indicates if Reader should try to infer aminoacid sequence type for each sequence. Can only identify aminoacid sequences. **Optional** |
| parse_method | 'rich' or 'quick' | 'rich' | Parse method to use. `'quick'` parsing method just parses the header and the sequence into individual properties, so it's much faster and less memory intensive. If selected, `sequences_type` and `infer_type` parameters are ignored. `'rich'` implements more functionality (`FastaSequence`), but is slower. **Optional** |

```Python
fastaparser.Reader(fasta_file, sequences_type=None, infer_type=False, parse_method='rich')
```

### Raises
Errors that can occur when instantiating a Reader class

**TypeError**

* If `fasta_file`, `sequences_type`, `infer_type` or `parse_method` are of the wrong type.
* If `fasta_file` is not a file object, is closed or is not readable.

## Attributes
Instances of the Reader class have the following attributes

| Attribute | Type / Value | Description |
|---|---|---|
| fasta_file | file object | The FASTA file passed as parameter |
| sequences_type | 'nucleotide', 'aminoacid' or None | Indicates the type of sequences to expect. Can be None if not known |
| infer_type | bool | True if Reader was set to infer the sequence type, False otherwise |
| parse_method | 'rich' or 'quick' | Parse method used |

## Special Methods
* \_\_iter__
* \_\_next__
* \_\_repr__


# fastaparser.Writer
Writer for the given FASTA file.
Writes `FastaSequence` objects or tuples of (`header`, `sequence`) to the given file.

## Parameters
The Writer class can be instantiated with the following parameter

| Parameter | Type / Value | Default | Description|
|---|---|---|---|
| fasta_file | file object | | An opened file handle (for writing). **Must be provided** |

```Python
fastaparser.Writer(fasta_file)
```

### Raises
Errors that can occur when instantiating a Writer class

**TypeError**

* If `fasta_file` is of the wrong type.
* If `fasta_file` is not a file object, is closed or is not writable.

## Attributes
Instances of the Writer class have the following attribute

| Attribute | Type / Value | Description |
|---|---|---|
| fasta_file | file object | The FASTA file passed as parameter |

## Methods
Instances of the Writer class have the following methods

### writefasta
Writes a single FASTA sequence to the provided file. Open the file with mode `'a'` if you want to append sequences to an existing FASTA file.

| Parameter | Type / Value | Default | Description|
|---|---|---|---|
| fasta_sequence | `FastaSequence` or (`header` : str, `sequence` : str) | | A FASTA sequence is built from the data contained in the provided `FastaSequence` object or the tuple of (`header`, `sequence`). `header` may contain or not the starting `'>'`. `header` can be an empty string. `sequence` must be a non empty string. **Must be provided** |

```Python
Writer.writefasta(fasta_sequence)
```

#### Raises
Errors that can occur when calling Writer.writefasta

**TypeError**

* If `fasta_sequence` is of the wrong type.

### writefastas
Writes multiple FASTA sequences to the provided file.
Simply calls Writer.writefasta for each object in `fasta_sequences`.
Open the file with mode `'a'` if you want to append multiple sequences to an existing FASTA file.

| Parameter | Type / Value | Default | Description|
|---|---|---|---|
| fasta_sequence | iterable of `FastaSequence` or iterable of (`header` : str, `sequence` : str) | | FASTA sequences are built from the data contained in the provided `FastaSequence` objects or the tuples of (`header`, `sequence`). `header`s may contain or not the starting `'>'`. `header`s can be empty strings. `sequence`s must be non empty strings. **Must be provided** |

```Python
Writer.writefastas(fasta_sequences)
```

#### Raises
Errors that can occur when calling Writer.writefastas

**TypeError**

* If `fasta_sequences` is not iterable.

## Special Methods
* \_\_repr__