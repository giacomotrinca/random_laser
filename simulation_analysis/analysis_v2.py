import os


def list_directories(path):
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

def find_PT(string):
    if "PT" in string:
        return True
    else:
        return False


