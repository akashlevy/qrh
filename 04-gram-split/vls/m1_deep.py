import sys; sys.path.insert(0, '.')
from chamber_analysis import analyse
for cap, mc in [(60000, 40), (200000, 80)]:
    reps = analyse('M1', cap=cap, maxcoef=mc, verbose=True)
    print('cap', cap, 'maxcoef', mc, 'types', len(reps), flush=True)
