class ParseInputError(Exception):
    """
    Exception occured while parsing an input file
    """
    pass

class ValidateError(Exception):
    """
    Exception used to signal a failure while validating an input/output
    """
    pass
