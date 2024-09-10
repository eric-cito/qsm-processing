from typing import List, AnyStr
import os

def ListDirWithFullPaths(dir, filenameEndsWith:str = None, errorIfNotFound=True) -> List[AnyStr]:
    '''
    :returns: Sorted list of full paths in the dir provided
    '''
    if errorIfNotFound or os.path.exists(dir): 
        found = [os.path.join(dir,fn) for fn in os.listdir(dir)]
        found.sort()

        if filenameEndsWith is not None:
            found = [loc for loc in found if os.path.split(loc)[-1].endswith(filenameEndsWith)]
        return found
    else:
        return list()
    
def DeleteAllIfFound(locs:List[str]):
    for loc in locs:
        if os.path.exists(loc):
            os.remove(loc)


def GetWithDifferentSuffix(loc, newSuffix:str):
    '''Gets the file path with a new suffix'''
    if len(newSuffix) > 0 and (newSuffix[0] != '.'):
        newSuffix = '.' + newSuffix
    return GetWithoutSuffix(loc) + newSuffix


def GetParentDirectory(loc_forAFile):
    return os.path.split(os.path.abspath(loc_forAFile))[0]


def GetFilename(loc:str) -> str:
    return os.path.split(loc)[-1]


def GetWithoutSuffix(loc:str) ->str:
    loc = os.path.abspath(loc)
    return os.path.join(GetParentDirectory(loc), GetFilename(loc).split(".")[0])