from typing import List, AnyStr
import os

def ListDirWithFullPaths(dir, filenameEndsWith:str = None) -> List[AnyStr]:
    '''
    :returns: Sorted list of full paths in the dir provided
    '''
    found = [os.path.join(dir,fn) for fn in os.listdir(dir)]
    found.sort()

    if filenameEndsWith is not None:
        return [loc for loc in found if os.path.split(loc)[-1].endswith(filenameEndsWith)]