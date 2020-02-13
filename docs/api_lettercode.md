# fastaparser.LetterCode
Represents a single letter code.

## Parameters
The LetterCode class can be instantiated with the following parameters
```Python
fastaparser.LetterCode(letter_code, letter_type=None)
```

| Parameter | Type / Value | Default | Description|
|:---:|:---:|:---:|---|
| letter_code | str | | Letter code. **Must be provided** |
| letter_type | 'nucleotide', 'aminoacid' or None | None | Type of letter code, `None` if there is no information. **Optional** |

#### Raises
**TypeError**

* If `letter_code` or `letter_type` are of the wrong type.

## Attributes
Instances of the LetterCode class have the following attributes

| Attribute | Type / Value | Editable | Description |
|:---:|:---:|:---:|---|
| letter_code | str | No | Upper case letter code. |
| letter_type | str or None | Yes | `'nucleotide'` or `'aminoacid'`. `None` if there is no information about sequence type. |
| description | str | No | Description or nucleotide/aminoacid name of letter code (can be an empty string). |
| degenerate | bool or None | No | Indicates if a letter code is degenerate or not (can be `None` if letter code is not defined in the FASTA specification or `letter_type` is unknown). |
| supported | bool | No | Indicates if letter code is supported or not (ie, if `letter_type` is provided and letter code is defined in the FASTA specification). |
| in_fasta_spec | bool | No | Indicates if letter code is defined in the FASTA specification. |

Editable attributes can be set by standard variable assignment and deleted/reset with the del keyword:
```Python
lettercode_object.letter_type = 'nucleotide'
del lettercode_object.letter_type
```

## Methods
Instances of the LetterCode class have the following method

### complement
Complementary letter code (ideally, of a nucleotide).
If `letter_code` is not a nucleotide letter code, the complementary will be `letter_code`.

In order not to impose the setting of `letter_type` as `'nucleotide'`, this method will work for any letter code
(as long as `letter_type` is not `'aminoacid'`), which has the side effect of returning nonsensical results when
letter code is not a nucleotide.

Ex: For aminoacid letter codes that overlap with nucleotide letter codes, the output will be the complement of
the nucleotide represented by the same letter code, which makes no sense.

```Python
LetterCode.complement()
```

#### Returns
**LetterCode**

Complement of current `LetterCode`. Same `LetterCode` is returned if letter code is not a valid nucleotide.

#### Raises
**TypeError**

* If `letter_type` is `'aminoacid'`.

## Class Methods
The LetterCode class has the following class method

### from_lettercode
Initializes with the given `LetterCode` object (alternate `__init__` method).

```Python
LetterCode.from_lettercode(lettercode)
```

| Parameter | Type / Value | Default | Description |
|:---:|:---:|:---:|---|
| lettercode | LetterCode | | `LetterCode` object. **Must be provided** |

#### Returns
**LetterCode**

Copy of `lettercode` (`LetterCode` object).

#### Raises
**TypeError**

* If `lettercode` is not a `LetterCode`.

## Special Methods
* \_\_eq__
* \_\_repr__
* \_\_str__
