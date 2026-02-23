import pandas as pd

from sklearn.svm import SVC
from src.models.base_model import BaseModel


class SVM(BaseModel):
    def __init__(self, **kwargs):
        super().__init__()
        self.model = None
        self.printing_name = 'Support Vector Machine'

    def _init_model(self, custom_params=None):
        pass

    def fit(self, X, y, X_val=None, y_val=None, num_rounds=-1, kernel='linear', **kwargs):
        """Fit the model to the data."""
        self.model = SVC(
            probability=True,
            kernel=kernel,
            max_iter=num_rounds,
            class_weight='balanced',
            **kwargs
        )
        self.model.fit(X, y)

    def predict(self, X):
        """Make predictions using the trained model."""
        return self.model.predict_proba(X)[:, 1]

    def get_feature_importance(self, headers, **kwargs):
        """Return feature importance if available."""
        if self.model is None:
            raise ValueError("Model is not trained yet.")

        if hasattr(self.model, 'coef_'):
            importance = self.model.coef_[0]
            importance = {i: v for i, v in enumerate(importance)}
            importance = {k: v for k, v in sorted(importance.items(), key=lambda item: item[1], reverse=True)}
            importance_df = pd.DataFrame(importance.items(), columns=['Feature', 'Importance'])
            importance_df['Feature'] = headers[importance_df['Feature']]
            return importance_df
        else:
            raise ValueError("Feature importances are not available for non-linear SVM kernels.")
