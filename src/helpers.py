import numpy as np

def num_str(num):
    try:
        return str(int(num)) if int(num) == num else str(num)
    # Catch if we're trying to turn a list into an int
    # also I hate weak typing
    except TypeError:
        return " ".join(map(num_str, num))