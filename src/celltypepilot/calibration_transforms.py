from __future__ import annotations
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from scipy.optimize import minimize
from typing import Any

from .calibration import calibration_diagnostics

class IsotonicCalibrator:
    def __init__(self):
        self.model = IsotonicRegression(out_of_bounds='clip')
        self._fitted = False
        
    def fit(self, scores: np.ndarray, y_correct: np.ndarray):
        self.model.fit(scores, y_correct)
        self._fitted = True
        return self
        
    def transform(self, scores: np.ndarray) -> np.ndarray:
        if not self._fitted:
            raise ValueError("Calibrator not fitted")
        return self.model.predict(scores)
        
    def save(self, path: str | Path):
        path = Path(path)
        data = {
            'type': 'IsotonicCalibrator',
            'X_min_': float(self.model.X_min_) if hasattr(self.model, 'X_min_') else None,
            'X_max_': float(self.model.X_max_) if hasattr(self.model, 'X_max_') else None,
            'X_thresholds_': self.model.X_thresholds_.tolist() if hasattr(self.model, 'X_thresholds_') else None,
            'y_thresholds_': self.model.y_thresholds_.tolist() if hasattr(self.model, 'y_thresholds_') else None,
        }
        path.write_text(json.dumps(data), encoding='utf-8')
        
    def load(self, path: str | Path):
        data = json.loads(Path(path).read_text(encoding='utf-8'))
        self.model.X_min_ = data.get('X_min_')
        self.model.X_max_ = data.get('X_max_')
        self.model.X_thresholds_ = np.array(data.get('X_thresholds_')) if data.get('X_thresholds_') else None
        self.model.y_thresholds_ = np.array(data.get('y_thresholds_')) if data.get('y_thresholds_') else None
        self.model.f_ = lambda x: np.interp(x, self.model.X_thresholds_, self.model.y_thresholds_)
        self._fitted = True
        
    def calibration_artifact(self) -> dict:
        return {
            "schema_version": "celltypepilot.calibration-transform.v1",
            "method": "isotonic"
        }

class PlattCalibrator:
    def __init__(self):
        self.model = LogisticRegression(solver='lbfgs')
        self._fitted = False
        
    def fit(self, scores: np.ndarray, y_correct: np.ndarray):
        X = scores.reshape(-1, 1)
        self.model.fit(X, y_correct)
        self._fitted = True
        return self
        
    def transform(self, scores: np.ndarray) -> np.ndarray:
        if not self._fitted:
            raise ValueError("Calibrator not fitted")
        X = scores.reshape(-1, 1)
        return self.model.predict_proba(X)[:, 1]
        
    def save(self, path: str | Path):
        path = Path(path)
        data = {
            'type': 'PlattCalibrator',
            'coef_': self.model.coef_.tolist() if hasattr(self.model, 'coef_') else None,
            'intercept_': self.model.intercept_.tolist() if hasattr(self.model, 'intercept_') else None,
            'classes_': self.model.classes_.tolist() if hasattr(self.model, 'classes_') else None,
        }
        path.write_text(json.dumps(data), encoding='utf-8')
        
    def load(self, path: str | Path):
        data = json.loads(Path(path).read_text(encoding='utf-8'))
        self.model.coef_ = np.array(data['coef_'])
        self.model.intercept_ = np.array(data['intercept_'])
        self.model.classes_ = np.array(data['classes_'])
        self._fitted = True
        
    def calibration_artifact(self) -> dict:
        return {
            "schema_version": "celltypepilot.calibration-transform.v1",
            "method": "platt"
        }

class TemperatureCalibrator:
    def __init__(self):
        self.T = 1.0
        self._fitted = False
        
    def fit(self, scores: np.ndarray, y_correct: np.ndarray):
        def nll(T):
            probs = 1 / (1 + np.exp(-scores / T[0]))
            probs = np.clip(probs, 1e-15, 1 - 1e-15)
            return -np.sum(y_correct * np.log(probs) + (1 - y_correct) * np.log(1 - probs))
            
        res = minimize(nll, [1.0], bounds=[(1e-3, 10.0)])
        self.T = res.x[0]
        self._fitted = True
        return self
        
    def transform(self, scores: np.ndarray) -> np.ndarray:
        if not self._fitted:
            raise ValueError("Calibrator not fitted")
        return 1 / (1 + np.exp(-scores / self.T))
        
    def save(self, path: str | Path):
        path = Path(path)
        data = {'type': 'TemperatureCalibrator', 'T': float(self.T)}
        path.write_text(json.dumps(data), encoding='utf-8')
        
    def load(self, path: str | Path):
        data = json.loads(Path(path).read_text(encoding='utf-8'))
        self.T = data['T']
        self._fitted = True
        
    def calibration_artifact(self) -> dict:
        return {
            "schema_version": "celltypepilot.calibration-transform.v1",
            "method": "temperature"
        }

def auto_select_calibrator(scores: np.ndarray, y_correct: np.ndarray, n_cv: int = 5) -> tuple[Any, dict]:
    kf = KFold(n_splits=n_cv, shuffle=True, random_state=42)
    
    models = {
        'isotonic': IsotonicCalibrator,
        'platt': PlattCalibrator,
        'temperature': TemperatureCalibrator
    }
    
    cv_eces = {k: [] for k in models}
    
    # Needs fake strings to reuse calibration_diagnostics which expects labels
    y_true_str = np.where(y_correct == 1, "A", "B")
    
    for train_idx, test_idx in kf.split(scores):
        X_train, X_test = scores[train_idx], scores[test_idx]
        y_train, y_test = y_correct[train_idx], y_correct[test_idx]
        y_true_test = y_true_str[test_idx]
        y_pred_test = np.full_like(y_true_test, "A") # Predict always A (correct if 1)
        
        for name, cls in models.items():
            model = cls().fit(X_train, y_train)
            calib_scores = model.transform(X_test)
            try:
                diag, _, _, _ = calibration_diagnostics(y_true_test, y_pred_test, calib_scores)
                ece = diag['ece']
            except Exception:
                ece = np.nan
            cv_eces[name].append(ece)
            
    mean_eces = {k: np.nanmean(v) for k, v in cv_eces.items()}
    best_name = min(mean_eces, key=lambda k: mean_eces[k] if not np.isnan(mean_eces[k]) else float('inf'))
    
    best_model = models[best_name]().fit(scores, y_correct)
    
    report = {
        "selected_method": best_name,
        "mean_eces": mean_eces
    }
    
    return best_model, report
