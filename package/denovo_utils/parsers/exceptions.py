class DenovoEngineNotSupported(ValueError):
    """
    Exception raised when a de novo search engine is not supported.

    Attributes
    ----------
    parameter : str
        The provided engine label that caused the exception.
    allowed_values : list
        The list of allowed engine labels.
    """
    def __init__(self, parameter, allowed_values):
        """
        Initialize the exception with the provided parameter and allowed values.

        Parameters
        ----------
        parameter : str
            The provided engine label that caused the exception.
        allowed_values : list
            The list of allowed engine labels.
        """
        self.parameter = parameter
        self.allowed_values = allowed_values
        super().__init__(f"Invalid parameter value: {parameter}. Allowed values are: {allowed_values}")

class NoResultsToMergeException(ValueError):
    """
    Exception raised when no results are available to merge.

    Attributes
    ----------
    parameter : str
        The provided engine label that caused the exception.
    allowed_values : list
        The list of allowed engine labels.
    """
    def __init__(self):
        super().__init__(
            "No results in the 'results' attribute dictionary. Call the method 'read_result' first before trying to merge!"
        )

class SeparatorCharacterInTitle(ValueError):
    """
    Exception raised when the character used to split search engine names is in the title.

    Attributes
    ----------
    parameter : str
        The provided engine label that caused the exception.
    allowed_values : list
        The list of allowed engine labels.
    """
    def __init__(self, char, title):
        super().__init__(
            f"Cannot properly use '{char}' in title '{title}' to create a splittable title name. Use a different character as separator."
        )