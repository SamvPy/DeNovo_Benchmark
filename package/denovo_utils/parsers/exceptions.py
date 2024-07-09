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