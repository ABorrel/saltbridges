from re import sub

def suppSpace(string):
    """ delt space in string
    in : string
    out : string without space
    """

    string = sub(" ", "", string)
    return string


def formatFloat (string):
    """format string to float
    in : string
    out : float"""

    string = suppSpace(string)
    floatOut = float (string)

    return floatOut


def formatInt (string):
    """format string to integer
    in : string
    out : integer"""

    string = suppSpace(string)
    integer = int (string)

    return integer

def date(date):
    """format date in time module
    in : date
    out : date formated"""
    
    date = sub("/", "-", date)
    date = sub(" ", "_", date)

    return date + "_"
