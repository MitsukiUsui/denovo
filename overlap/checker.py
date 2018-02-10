#!/usr/bin/env python3

import sys
import pandas as pd

sys.path.append("../helper")
from myio import *

target=sys.argv[1]
strain_lst=get_strain_lst(target)

# check event
eventFilepath="./out/{}/event.csv".format(target)
try:
    event_df=pd.read_csv(eventFilepath, dtype={"phase": object})
except FileNotFoundError:
    print("ERROR: {} dose not exist".format(eventFilepath), file = sys.stderr)
    exit(1)
except KeyError:
    print("ERROR: {} dose not column {}".format(eventFilepath), file = sys.stderr)
    exit(1)
exit(0)
    
