class BaseModel:
    def __init__(self, **kwargs):
        self.printing_name = ''
        self.model = None

    def _init_model(self):
        raise NotImplementedError

    def fit(self, X, y, **kwargs):
        raise NotImplementedError

    def predict(self, X):
        raise NotImplementedError

    def get_feature_importance(self, headers, **kwargs):
        raise NotImplementedError