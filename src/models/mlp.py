import pandas as pd

from sklearn.utils.class_weight import compute_class_weight
import numpy as np

import torch
from torch import nn
from torch.utils.data import DataLoader

from src.models.base_model import BaseModel


class Net(nn.Module):
    def __init__(self, input_size):
        super(Net, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, 512),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 16),
            nn.ReLU(),
            nn.Linear(16, 8),
            nn.ReLU(),
            nn.Linear(8, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.model(x)


class MLP(BaseModel):
    def __init__(self, **kwargs):
        super().__init__()
        self.input_size = kwargs.get('input_size', 257)
        self.criterion = None
        self.optimizer = None
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        self.model = self._init_model()
        self.printing_name = 'MLP'

    def _init_model(self):
        return Net(self.input_size).to(self.device)

    def fit(
            self,
            X_train,
            y_train,
            X_val=None,
            y_val=None,
            num_rounds=250,
            batch_size=256,
            early_stopping_rounds=100,
            **kwargs
    ):
        class_weights = torch.FloatTensor(compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)).to(
            self.device)
        weight_ratio = class_weights[1] / class_weights[0]

        self.criterion = nn.BCELoss(weight=torch.tensor([weight_ratio]).to(self.device))
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)

        data = torch.utils.data.TensorDataset(torch.from_numpy(X_train).float().to(self.device),
                                              torch.from_numpy(y_train).float().to(self.device))
        data_loader = DataLoader(data, batch_size=batch_size, shuffle=True)

        best_val_loss = float('inf')
        epochs_no_improve = 0
        early_stopping = False
        best_model_wts = None

        for epoch in range(num_rounds):
            self.model.train()

            for inputs, labels in data_loader:
                self.optimizer.zero_grad()
                outputs = self.model(inputs).squeeze()
                loss = self.criterion(outputs, labels)
                loss.backward()
                self.optimizer.step()

            # Check validation loss for early stopping
            if X_val is not None and y_val is not None:
                self.model.eval()
                with torch.no_grad():
                    inputs = torch.from_numpy(X_val).float().to(self.device)
                    labels = torch.from_numpy(y_val).float().to(self.device)
                    outputs = self.model(inputs).squeeze()
                    val_loss = self.criterion(outputs, labels).item()

                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    epochs_no_improve = 0
                    best_model_wts = self.model.state_dict().copy()
                else:
                    epochs_no_improve += 1
                    if epochs_no_improve >= early_stopping_rounds:
                        print(f'Early stopping at epoch {epoch + 1}')
                        early_stopping = True
                        break

        # Load best model weights
        if best_model_wts is not None:
            self.model.load_state_dict(best_model_wts)

    def predict(self, X):
        self.model.eval()
        with torch.no_grad():
            inputs = torch.from_numpy(X).float().to(self.device)
            outputs = self.model(inputs).squeeze()
            probabilities = outputs.cpu().numpy()
            return probabilities

    def get_feature_importance(self, headers, **kwargs):
        importance_df = pd.DataFrame([], columns=['Feature', 'Importance'])

        return importance_df