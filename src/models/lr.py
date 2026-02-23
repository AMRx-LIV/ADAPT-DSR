import pandas as pd

from sklearn.linear_model import LogisticRegression as LR

from src.models.base_model import BaseModel


class LogisticRegression(BaseModel):
    def __init__(self, **kwargs):
        super().__init__()
        self.model = None
        self.printing_name = 'Logistic Regression'

    def _init_model(self, custom_params=None):
        pass

    def fit(self, X, y, X_val=None, y_val=None, num_rounds=10000, **kwargs):
        """Fit the model to the data."""
        self.model = LR(max_iter=num_rounds, class_weight='balanced', **kwargs)
        self.model.fit(X, y)

    def predict(self, X):
        """Make predictions using the trained model."""
        return self.model.predict_proba(X)[:, 1]

    def get_feature_importance(self, headers, **kwargs):
        """Logistic Regression doesn't have direct feature importance; return coefficients."""
        if self.model is None:
            raise ValueError("Model is not trained yet.")

        importance = self.model.coef_[0]
        importance = {i: v for i, v in enumerate(importance)}
        importance = {k: v for k, v in sorted(importance.items(), key=lambda item: item[1], reverse=True)}
        importance_df = pd.DataFrame(importance.items(), columns=['Feature', 'Importance'])

        importance_df['Feature'] = headers[importance_df['Feature']]

        return importance_df
