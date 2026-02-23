import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_auc_score, confusion_matrix, f1_score, recall_score, precision_recall_curve, roc_curve, \
    precision_score, average_precision_score
import matplotlib.pyplot as plt

from src.dataset import dataset_name, culture_name


class Evaluator:
    def __init__(
            self,
            model_best_fold,
            model_bast_fold_threshold,
            model_from_x_eval,
            labels: pd.DataFrame
    ):
        self.model_best_fold = model_best_fold
        self.model_bast_fold_threshold = model_bast_fold_threshold
        self.model_from_x_eval = model_from_x_eval
        self.labels = labels

        self._rename_labels()

    def _rename_labels(self, keep=20):
        # read mapping
        mapping = pd.read_csv('variable_key.csv')
        mapping_dict = dict(zip(mapping['Long_name'], mapping['Short_name']))
        # rename columns if they are in the mapping
        self.labels = self.labels.map(lambda x: mapping_dict.get(x, x))

        # Clamp labels to keep most important ones
        # labels are dp index
        self.labels = self.labels.map(lambda x: x[:keep] + '...' if len(x) > keep else x)

    def evaluate_holdout(self, X_holdout, y_holdout, holdout_ori, phase='best'):
        if phase == 'best':
            y_holdout_pred = self.model_best_fold.predict(X_holdout)
            model_name = self.model_best_fold.printing_name
            threshold = 0.5
            print(f'using best threshold: {threshold:.2f}')
        else:
            y_holdout_pred = self.model_from_x_eval.predict(X_holdout)
            model_name = self.model_from_x_eval.printing_name
            threshold = 0.5

        # Save probabilities and original inputs to CSV
        holdout_results = holdout_ori.copy()
        holdout_results['probabilities'] = y_holdout_pred
        holdout_results.to_csv(
            f'outputs/{dataset_name}_{model_name}_holdout_prob_{phase}.csv', index=False)

        y_holdout_pred_binary = y_holdout_pred > threshold
        # print roc auc, average precision, correlation, f1, precision, recall, accuracy
        holdout_auc = roc_auc_score(y_holdout, y_holdout_pred)
        print(f'{model_name} holdout ROC AUC ({phase}): {holdout_auc:.2f}')

        # average precision
        ap = average_precision_score(y_holdout, y_holdout_pred)
        print(f'{model_name} holdout Average Precision ({phase}): {ap:.2f}')

        # F1 score
        f1 = f1_score(y_holdout, y_holdout_pred_binary)
        print(f'{model_name} holdout F1 ({phase}): {f1:.2f}')

        # correlation
        corr = np.corrcoef(y_holdout, y_holdout_pred)[0, 1]
        print(f'{model_name} holdout correlation ({phase}): {corr:.2f}')

        # precision
        precision = precision_score(y_holdout, y_holdout_pred_binary)
        print(f'{model_name} holdout precision ({phase}): {precision:.2f}')

        # recall
        recall = recall_score(y_holdout, y_holdout_pred_binary)
        print(f'{model_name} holdout recall ({phase}): {recall:.2f}')

        # accuracy
        accuracy = np.mean(y_holdout == y_holdout_pred_binary)
        print(f'{model_name} holdout accuracy ({phase}): {accuracy:.2f}')

        # Plot Confusion Matrix as a heatmap
        cm = confusion_matrix(y_holdout, y_holdout_pred_binary)
        # flip y
        cm = cm[::-1]

        plt.figure(figsize=(6, 4))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                    xticklabels=['0', '1'],
                    yticklabels=['1', '0'])
        plt.yticks(rotation=90)
        plt.title(f'{culture_name} {model_name} Confusion Matrix ({phase})')
        plt.ylabel('True Label')
        plt.xlabel('Predicted Label')
        plt.tight_layout()
        plt.savefig(f'outputs/{dataset_name}_{model_name}_confusion_matrix_{phase}.png')
        plt.close()

        # plot precision-recall curve
        precision, recall, thresholds = precision_recall_curve(y_holdout, y_holdout_pred)
        plt.figure(figsize=(6, 4))
        plt.plot(recall, precision, label=f'AP={ap:.2f}')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim(bottom=0)
        plt.title(f'{culture_name} {model_name} Precision-Recall Curve ({phase})')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'outputs/{dataset_name}_{model_name}_pr_curve_{phase}.png')
        plt.close()

        # plot roc curve
        fpr, tpr, _ = roc_curve(y_holdout, y_holdout_pred)
        plt.figure(figsize=(6, 4))
        plt.plot(fpr, tpr, label=f'AUC={holdout_auc:.2f}')
        plt.plot([0, 1], [0, 1], 'k--')  # Add diagonal line
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'{culture_name} {model_name} ROC Curve ({phase})')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'outputs/{dataset_name}_{model_name}_roc_curve_{phase}.png')
        plt.close()

    def compute_feature_importance(self):
        # Compute feature importance
        fi = self.model_from_x_eval.get_feature_importance(self.labels)
        if not fi.empty:
            fi['Abs_Importance'] = np.abs(fi['Importance'])
            fi = fi.sort_values('Abs_Importance', ascending=False)
            fi.drop(columns='Abs_Importance', inplace=True)
            fi.to_csv(f'outputs/{dataset_name}_{self.model_from_x_eval.printing_name}_fi.csv', index=False)

            # Plot Feature Importance
            plt.figure(figsize=(10, 6))
            fi.iloc[:20, :].plot(kind='bar', x='Feature', y='Importance', legend=False)
            plt.title(f'{culture_name} {self.model_from_x_eval.printing_name} Feature Importance')
            plt.tight_layout()
            plt.savefig(f'outputs/{dataset_name}_{self.model_from_x_eval.printing_name}_feature_importance.png')
            plt.close()

    def compute_shap_values(self, X_eval):
        try:
            self.model_best_fold.shapley_values(
                X_eval,
                pd.Series(self.labels),
                f'outputs/{dataset_name}_{self.model_from_x_eval.printing_name}_shap.png'
            )
        except Exception as e:
            print(f"SHAP values computation failed for {self.model_from_x_eval.printing_name}: {e}")
