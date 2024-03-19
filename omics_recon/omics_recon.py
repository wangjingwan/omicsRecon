import pandas as pd

import evaluation
import plotting
import preprocess
import predict
import utils

class omicsRecon():
    def __init__(self):
            self.eva = evaluation
            self.pp = preprocess
            self.pl = plotting
            self.pred = predict
            pass