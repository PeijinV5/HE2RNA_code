"""
HE2RNA: Computation of correlations
Copyright (C) 2020  Owkin Inc.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import numpy as np
import pickle as pkl
from joblib import Parallel, delayed
from scipy import stats
from sklearn.metrics import accuracy_score

def corr(pred, label, i):
    return round(stats.pearsonr(label[:, i], pred[:, i])[0],4), round(stats.pearsonr(label[:, i], pred[:, i])[1],4)

def accuracy(pred, label, i):
    return accuracy_score(label[:,i], y_pred[:,i])

def compute_metrics(label, pred):
    res = Parallel(n_jobs=16)(
        delayed(corr)(pred, label, i) for i in range(label.shape[1])
    )
    return res

def compute_metrics_acc(pred, label):
    res = Parallel(n_jobs=16)(
        delayed(accuracy)(pred, label, i) for i in range(label.shape[1])
    )
    return res
