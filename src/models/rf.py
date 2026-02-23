from sklearn.ensemble import RandomForestClassifier
import pandas as pd

from src.models.base_model import BaseModel


class RandomForest(BaseModel):
    def __init__(self, **kwargs):
        super().__init__()
        self.model = None
        self.printing_name = 'RandomForest'
        # Set default parameters for the Random Forest model
        self.params = {
            'n_estimators': 200,
            'max_depth': 9,
            'class_weight': 'balanced',
        }

    def _init_model(self, custom_params=None):
        # Update parameters if custom_params are provided
        if custom_params:
            self.params.update(custom_params)
        # Initialize RandomForestClassifier with the given parameters
        self.model = RandomForestClassifier(**self.params)

    def fit(self, X, y, X_val=None, y_val=None, **kwargs):
        # Ensure model is initialized before fitting
        if self.model is None:
            self._init_model()

        # Fit the model
        self.model.fit(X, y)

    def predict(self, X):
        # Predict probabilities of the positive class
        return self.model.predict_proba(X)[:, 1]

    def get_feature_importance(self, headers, **kwargs):
        # Get feature importance from the Random Forest model
        importance = self.model.feature_importances_
        importance_df = pd.DataFrame(
            list(zip(headers, importance)),
            columns=['Feature', 'Importance']
        ).sort_values(by='Importance', ascending=False)

        return importance_df

    @staticmethod
    def _roc_auc_eval(preds, dtrain):
        # Not needed for RandomForest, but kept for consistency
        pass
