# Constants
Constants used throughout the whole package.
These include letter codes and their respective names for all "good" and degenerate nucleotides/aminoacids.

### NUCLEOTIDE_LETTER_CODES_GOOD
**dict**

| key | value |
|:---:|:---:|
| nucleotide letter codes | names/description |

### NUCLEOTIDE_LETTER_CODES_DEGENERATE
**dict**

| key | value |
|:---:|:---:|
| degenerate nucleotide letter codes | names/description |

### NUCLEOTIDE_LETTER_CODES_COMPLEMENT
**dict**

| key | value |
|:---:|:---:|
| nucleotide letter codes | complement letter codes |

### AMINOACID_LETTER_CODES_GOOD
**dict**

| key | value |
|:---:|:---:|
| aminoacid letter codes | names/description |

### AMINOACID_LETTER_CODES_DEGENERATE
**dict**

| key | value |
|:---:|:---:|
| degenerate aminoacid letter codes | names/description |

### LETTER_CODES
**dict**

| key | value |
|:---:|:---:|
| 'nucleotide' / 'aminoacid' | respective (LETTER_CODES_GOOD, LETTER_CODES_DEGENERATE) |


### LETTER_CODES_ALL
**set**

All valid FASTA letter codes.

### NUCLEOTIDE_LETTER_CODES_ALL
**set**

All valid FASTA nucleotide letter codes.

### AMINOACID_LETTER_CODES_ALL
**set**

All valid FASTA aminoacid letter codes.

### AMINOACIDS_NOT_IN_NUCLEOTIDES
**set**

All letter codes that only represent aminoacids (and not also nucleotides).
