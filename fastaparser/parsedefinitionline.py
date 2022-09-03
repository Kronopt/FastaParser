#!python
# coding: utf-8

"""
ParseDefinitionLine - Class intended to be extended by the Reader and Writer classes.
"""


class ParseDefinitionLine:
    """
    Implements a parser of FASTA definition lines (to be used by the Reader and Writer classes)

    Methods
    -------
    _parse_definition_line(definition_line)
        Parses FASTA definition lines

    Raises
    ------
    TypeError
        When calling _parse_definition_line, if definition_line is of the wrong type.
    """

    @staticmethod
    def _parse_definition_line(definition_line):
        """
        Parses a FASTA definition line and returns an id and a description.

        Parameters
        ----------
        definition_line : str
            FASTA sequence definition line (header). May contain or not the '>' symbol at the start.
            Can be an empty string.

        Returns
        -------
        tuple
            ID and description. Can both be empty strings.

        Raises
        ------
        TypeError
            If definition_line is of the wrong type.
        """
        if isinstance(
            definition_line, str
        ):  # '>id|more_id description ...' with or without the '>' at the start
            id_and_description = definition_line.split(
                maxsplit=1
            )  # first space separates id from description

            # both id and description can be empty
            if len(id_and_description) == 0 or (
                len(id_and_description) == 1 and id_and_description[0] == ">"
            ):
                _id = ""
                _description = ""
            else:
                if id_and_description[0].startswith(">"):
                    id_and_description[0] = id_and_description[0][1:]
                if (
                    len(id_and_description) == 1
                ):  # description can be empty (assumes only id present if len == 1)
                    _id = id_and_description[0]
                    _description = ""
                else:
                    _id, _description = id_and_description
        else:
            raise TypeError("definition_line must be str")

        return _id, _description
