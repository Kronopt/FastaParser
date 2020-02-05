# fastaparser.Reader
Parser/Reader for a given FASTA file.
Iterates over the FASTA file using one of two parsing mechanisms:

* **rich**:
Returns `FastaSequence` objects (default). Slower, but feature rich.
* **quick**:
Generates objects containing just the FASTA header and sequence attributes
for each sequence in the FASTA file.
Parses FASTA files faster but lacks some features.

## Parameters
The Reader class can be instantiated with the following parameters
```Python
fastaparser.Reader(fasta_file, sequences_type=None, infer_type=False, parse_method='rich')
```

| Parameter | Type / Value | Default | Description|
|:---:|:---:|:---:|---|
| fasta_file | file object | | An opened file handle (for reading). **Must be provided** |
| sequences_type | 'nucleotide', 'aminoacid' or None | None | Indicates the type of sequences to expect. `None` if unknown. **Optional** |
| infer_type | bool | False | Indicates if `Reader` should try to infer aminoacid sequence type for each sequence. Can only identify aminoacid sequences. **Optional** |
| parse_method | 'rich' or 'quick' | 'rich' | Parse method to use. `'quick'` parsing method just parses the header and the sequence into individual properties, so it's much faster and less memory intensive. If selected, `sequences_type` and `infer_type` parameters are ignored. `'rich'` implements more functionality (`FastaSequence`), but is slower. **Optional** |

#### Raises
**TypeError**

* If `fasta_file`, `sequences_type`, `infer_type` or `parse_method` are of the wrong type.
* If `fasta_file` is not a file object, is closed or is not readable.

## Attributes
Instances of the Reader class have the following attributes

| Attribute | Type / Value | Editable | Description |
|:---:|:---:|:---:|---|
| fasta_file | file object | No | The FASTA file passed as parameter |
| sequences_type | 'nucleotide', 'aminoacid' or None | No | Indicates the type of sequences to expect. Can be `None` if not known |
| infer_type | bool | No | `True` if `Reader` was set to infer the sequence type, `False` otherwise |
| parse_method | 'rich' or 'quick' | No | Parse method used |

## Special Methods
* \_\_iter__
* \_\_next__
* \_\_repr__
