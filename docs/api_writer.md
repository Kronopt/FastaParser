# fastaparser.Writer
Writer for the given FASTA file.
Writes [`FastaSequence`](api_fastasequence.md) objects or tuples of (`header`, `sequence`) to the given file.

## Parameters
The Writer class can be instantiated with the following parameter
```Python
fastaparser.Writer(fasta_file)
```

| Parameter | Type / Value | Default | Description|
|:---:|:---:|:---:|---|
| fasta_file | file object | | An opened file handle (for writing). **Must be provided** |

#### Raises
**TypeError**:

* If `fasta_file` is of the wrong type.
* If `fasta_file` is not a file object, is closed or is not writable.

## Attributes
Instances of the Writer class have the following attribute

| Attribute | Type / Value | Editable | Description |
|:---:|:---:|:---:|---|
| fasta_file | file object | No | The FASTA file passed as parameter |

## Methods
Instances of the Writer class have the following methods

### writefasta
Writes a single FASTA sequence to the provided file. Open the file with mode `'a'` if you want to append sequences to an existing FASTA file.

```Python
Writer.writefasta(fasta_sequence)
```

| Parameter | Type / Value | Default | Description |
|:---:|:---:|:---:|---|
| fasta_sequence | [FastaSequence](api_fastasequence.md) or (header: str, sequence: str) | | A FASTA sequence is built from the data contained in the provided [`FastaSequence`](api_fastasequence.md) object or the tuple of (`header`, `sequence`). `header` may contain or not the starting `'>'`. `header` can be an empty string. `sequence` must be a non empty string. **Must be provided** |

#### Raises
**TypeError**

* If `fasta_sequence` is of the wrong type.

### writefastas
Writes multiple FASTA sequences to the provided file.
Simply calls Writer.writefasta for each object in `fasta_sequences`.
Open the file with mode `'a'` if you want to append multiple sequences to an existing FASTA file.

```Python
Writer.writefastas(fasta_sequences)
```

| Parameter | Type / Value | Default | Description |
|:---:|:---:|:---:|---|
| fasta_sequence | iterable of [FastaSequence](api_fastasequence.md) or iterable of (header: str, sequence: str) | | FASTA sequences are built from the data contained in the provided [`FastaSequence`](api_fastasequence.md) objects or the tuples of (`header`, `sequence`). `header`s may contain or not the starting `'>'`. `header`s can be empty strings. `sequence`s must be non empty strings. **Must be provided** |

#### Raises
**TypeError**

* If `fasta_sequences` is not iterable.

## Special Methods
* \_\_repr__
