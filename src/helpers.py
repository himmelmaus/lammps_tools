import numpy as np

def num_str(num):
    try:
        return str(int(num)) if int(num) == num else str(num)
    # Catch if we're trying to turn a list into an int
    # also I hate weak typing
    except TypeError:
        return " ".join(map(num_str, num))
    
def snake_case(string):
    return string.lower().replace(" ", "_")

def readlines(file_obj, n):
    lines = [file_obj.readline() for _ in range(n)]
    return lines, file_obj