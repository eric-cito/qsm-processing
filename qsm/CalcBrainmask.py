import sys
from .IO import WriteImageIfPathProvided
from .SkullStripper import SkullStripper


if __name__ == "__main__":
    loc_orig = sys.argv[1]
    loc_to = sys.argv[2]

    ss = SkullStripper(loc_orig)
    ss.GetOrCalcBrainmask(loc_to)