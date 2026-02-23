import xgboost as xgb
from sklearn.metrics import f1_score, roc_auc_score
import numpy as np
import pandas as pd
import shap
import matplotlib.pyplot as plt

from src.models.base_model import BaseModel
from src.dataset import culture_name


class XGBoost(BaseModel):
    def __init__(self, **kwargs):
        super().__init__()
        self.model = None
        self.printing_name = 'XGBoost'
        self.params = {
            'objective': 'binary:logistic',
            'eval_metric': 'logloss',
            'learning_rate': 0.005,
            'max_depth': 3,
            'subsample': 0.7
        }
        self.evals_result = {}

    def _init_model(self, custom_params=None):
        if custom_params:
            self.params.update(custom_params)
        self.model = None

    def fit(self, X, y, X_val=None, y_val=None, num_rounds=5000, **kwargs):
        # Compute scale_pos_weight for imbalanced dataset
        num_pos = np.sum(y == 1)
        num_neg = np.sum(y == 0)
        self.params['scale_pos_weight'] = num_neg / num_pos

        dtrain = xgb.DMatrix(X, label=y)
        evals = [(dtrain, 'train')]

        if X_val is not None and y_val is not None:
            dval = xgb.DMatrix(X_val, label=y_val)
            evals.append((dval, 'eval'))

        self.model = xgb.train(
            self.params,
            dtrain,
            num_boost_round=num_rounds,
            evals=evals,
            evals_result=self.evals_result,
            early_stopping_rounds=10,
            verbose_eval=0
        )

    def predict(self, X):
        dtest = xgb.DMatrix(X)
        return self.model.predict(dtest)

    def get_feature_importance(self, headers, **kwargs):
        importance_type = kwargs.get('importance_type', 'gain')

        importance = self.model.get_score(importance_type=importance_type)
        importance = {int(k.replace('f', '')): v for k, v in
                      sorted(importance.items(), key=lambda item: item[1], reverse=True)}
        importance_df = pd.DataFrame(importance.items(), columns=['Feature', 'Importance'])

        importance_df['Feature'] = headers[importance_df['Feature']]

        return importance_df

    def shapley_values(
            self,
            X: np.ndarray,
            labels: list,
            save_path: str = 'shap_summary_plot.png'
    ):
        # Create a new figure for the plot
        plt.figure()

        X = pd.DataFrame(X)
        explainer = shap.Explainer(self.model)
        shap_values = explainer(X)

        # add feature names
        shap_values.feature_names = labels
        shap.summary_plot(shap_values, X, show=False)
        plt.title(f'{culture_name} {self.printing_name} SHAP Summary Plot')
        plt.savefig(save_path, dpi=300)
        plt.close()
