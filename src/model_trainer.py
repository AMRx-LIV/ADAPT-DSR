import numpy as np
import pandas as pd
from PIL.features import features
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score, f1_score, \
    confusion_matrix

from src.dataset import dataset
from src.utils import calc_best_threshold


class ModelTrainer:
    def __init__(self, model_class, data_folded, X_eval, y_eval):
        self.model_class = model_class
        self.data_folded = data_folded
        self.X_eval = X_eval
        self.y_eval = y_eval

        self.best_model_from_folds = None
        self.bast_model_from_folds_threshold = None
        self.model_from_x_eval = None
        self.best_auc = 0
        self.aucs = []

        self.precisions = []
        self.aps = []
        self.predictions = []
        self.interp_fpr = np.linspace(0, 1, 1000)
        self.interp_tprs = []
        self.interp_recall = np.linspace(0, 1, 1000)
        self.roc_mean_tpr = None
        self.roc_mean_fpr = None
        self.pr_mean_precision = None
        self.pr_mean_recall = None
        self.roc_mean_auc = None
        self.pr_mean_ap = None
        self.pr_std_ap = None
        self.roc_std_auc = None
        self.pr_max_precision = None
        self.pr_min_precision = None
        self.roc_min_tpr = None
        self.roc_max_tpr = None
        self.f1 = None
        self.corr = None

    def train_on_folds(self):
        for fold_idx, (X_train, X_test, y_train, y_test) in enumerate(self.data_folded):
            input_size = X_train.shape[1]
            model = self.model_class(input_size=input_size)

            # Train the model
            model.fit(X_train, y_train, X_val=X_test, y_val=y_test)
            y_pred = model.predict(X_test)
            self.predictions.append(y_pred)

            # Update the best model
            auc_score = roc_auc_score(y_test, y_pred)
            print(f'{model.printing_name} Fold {fold_idx + 1} AUC: {auc_score:.2f}')

            if auc_score > self.best_auc:
                self.best_auc = auc_score
                self.best_model_from_folds = model
                self.bast_model_from_folds_threshold = calc_best_threshold(y_test, y_pred)

            # Compute ROC curve and AUC
            fpr, tpr, _ = roc_curve(y_test, y_pred)
            interp_tpr = np.interp(self.interp_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            self.interp_tprs.append(interp_tpr)
            auc_score = roc_auc_score(y_test, y_pred)
            self.aucs.append(auc_score)

            # Compute Precision-Recall curve and AP score
            precision, recall, _ = precision_recall_curve(y_test, y_pred)
            interp_precision = np.interp(self.interp_recall, recall[::-1], precision[::-1], left=1.0, right=0.0)
            self.precisions.append(interp_precision)
            ap_score = average_precision_score(y_test, y_pred)
            self.aps.append(ap_score)

    def _get_indices_for_race(self, subject_ids, raw, race_str):
        """Return the list of indices i for subjects with race == race_str."""
        column_names = raw.columns
        columns_to_check = [col for col in column_names if race_str in col.lower()]
        race_indices = []
        for i, sbj_id in enumerate(subject_ids):
            subject_data = raw[raw["subject_id"] == sbj_id].iloc[0]
            data = subject_data[columns_to_check]
            if 1 in data.values:
                race_indices.append(i)
        return race_indices

    def _get_indices_for_other_race(self, subject_ids, raw):
        """
        Return the list of indices i for subjects who do not fall into
        White, Black, Asian, or Hispanic columns.
        """
        column_names = raw.columns
        # columns that are not white/black/asian/hispanic
        columns_to_check = [
            col for col in column_names
            if (
                    'race_white' not in col.lower() and
                    'race_black' not in col.lower() and
                    'race_asian' not in col.lower() and
                    'race_hispanic' not in col.lower() and
                    'race_' in col.lower()
            )
        ]
        other_indices = []
        for i, sbj_id in enumerate(subject_ids):
            subject_data = raw[raw["subject_id"] == sbj_id].iloc[0]
            data = subject_data[columns_to_check]
            if 1 in data.values:
                other_indices.append(i)
        return other_indices

    def _get_indices_for_sex(self, subject_ids, raw, is_male=True):
        """Return the indices for subjects who have male == is_male (boolean)."""
        sex_indices = []
        for i, sbj_id in enumerate(subject_ids):
            subject_data = raw[raw["subject_id"] == sbj_id].iloc[0]
            # column is named 'male' and is boolean
            if subject_data["male"] == is_male:
                sex_indices.append(i)
        return sex_indices

    def _compute_auc_for_indices(self, indices, labels, preds):
        """Compute ROC AUC given a subset of indices."""
        if len(indices) == 0:
            return np.nan
        return float(roc_auc_score(labels[indices], preds[indices]))

    def _get_age_bins_map(self, raw, bins, labels):
        """
        Create a dictionary: subject_id -> bin_label,
        using pd.cut on the 'anchor_age' column.
        """
        # We assume anchor_age in [18..91], though we'll let pd.cut handle edge cases.
        age_series = pd.cut(raw["anchor_age"], bins=bins, right=False, labels=labels)
        age_bin_map = {}
        for sbj_id, age_label in zip(raw["subject_id"], age_series):
            age_bin_map[sbj_id] = age_label
        return age_bin_map

    def fairness_metrics(self):
        # Gather overall subject IDs and labels
        subject_ids = np.concatenate(dataset.subject_ids_folded)
        all_testing_labels = np.concatenate([y_test for _, _, _, y_test in self.data_folded])
        all_testing_preds = np.concatenate(self.predictions)

        # raw data for referencing
        raw = dataset.raw_data

        # ---------------- Race-based analysis ----------------
        white_ids = self._get_indices_for_race(subject_ids, raw, "race_white")
        black_ids = self._get_indices_for_race(subject_ids, raw, "race_black")
        asian_ids = self._get_indices_for_race(subject_ids, raw, "race_asian")
        hispanic_ids = self._get_indices_for_race(subject_ids, raw, "race_hispanic")
        other_ids = self._get_indices_for_other_race(subject_ids, raw)

        white_auc = self._compute_auc_for_indices(white_ids, all_testing_labels, all_testing_preds)
        black_auc = self._compute_auc_for_indices(black_ids, all_testing_labels, all_testing_preds)
        asian_auc = self._compute_auc_for_indices(asian_ids, all_testing_labels, all_testing_preds)
        hispanic_auc = self._compute_auc_for_indices(hispanic_ids, all_testing_labels, all_testing_preds)
        other_auc = self._compute_auc_for_indices(other_ids, all_testing_labels, all_testing_preds)

        # ---------------- Sex-based analysis (male/female) ----------------
        male_ids = self._get_indices_for_sex(subject_ids, raw, is_male=True)
        female_ids = self._get_indices_for_sex(subject_ids, raw, is_male=False)

        male_auc = self._compute_auc_for_indices(male_ids, all_testing_labels, all_testing_preds)
        female_auc = self._compute_auc_for_indices(female_ids, all_testing_labels, all_testing_preds)

        # ---------------- Age-based analysis (10–19, 20–29, etc.) ----------------
        # Bins: [10, 20, 30, 40, 50, 60, 70, 80, 90, 101]
        # => intervals: [10,20), [20,30), [30,40), ..., [90,101)
        # Labels: "10-19", "20-29", "30-39", ..., "90-99"
        bins = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        labels = [f"{bins[i]}-{bins[i + 1] - 1}" for i in range(len(bins) - 1)]
        age_bin_map = self._get_age_bins_map(raw, bins, labels)

        age_bin_aucs = {}
        for age_label in labels:
            bin_indices = []
            for i, sbj_id in enumerate(subject_ids):
                if age_bin_map[sbj_id] == age_label:
                    bin_indices.append(i)
            age_bin_aucs[age_label] = (
                {
                    "auc": self._compute_auc_for_indices(
                        bin_indices, all_testing_labels, all_testing_preds
                    ),
                    "qty": len(bin_indices),
                    "percentage": len(bin_indices) / len(subject_ids)
                }
            )

        # ---------------- Combine metrics into a dictionary ----------------
        race_results = {
            "white": {
                "auc": white_auc,
                "qty": len(white_ids),
                "percentage": len(white_ids) / len(subject_ids)
            },
            "black": {
                "auc": black_auc,
                "qty": len(black_ids),
                "percentage": len(black_ids) / len(subject_ids)
            },
            "asian": {
                "auc": asian_auc,
                "qty": len(asian_ids),
                "percentage": len(asian_ids) / len(subject_ids)
            },
            "hispanic": {
                "auc": hispanic_auc,
                "qty": len(hispanic_ids),
                "percentage": len(hispanic_ids) / len(subject_ids)
            },
            "other": {
                "auc": other_auc,
                "qty": len(other_ids),
                "percentage": len(other_ids) / len(subject_ids)
            }
        }
        sex_results = {
            "male": {
                "auc": male_auc,
                "qty": len(male_ids),
                "percentage": len(male_ids) / len(subject_ids)
            },
            "female": {
                "auc": female_auc,
                "qty": len(female_ids),
                "percentage": len(female_ids) / len(subject_ids)
            }
        }
        age_results = age_bin_aucs

        # Example: storing them or just returning them
        fairness_results = {
            "race": race_results,
            "sex": sex_results,
            "age": age_results
        }

        print(f"Race fairness: {race_results}")
        print(f"Sex fairness: {sex_results}")
        print(f"Age fairness: {age_results}")

        return fairness_results

    def train_on_full(self):
        input_size = self.X_eval.shape[1]
        model = self.model_class(input_size=input_size)
        model.fit(self.X_eval, self.y_eval)
        self.model_from_x_eval = model

    def compute_mean_metrics(self):
        y_tests = np.concatenate([y_test for _, _, _, y_test in self.data_folded])
        y_preds = np.concatenate(self.predictions)

        self.roc_mean_auc = roc_auc_score(y_tests, y_preds)
        self.roc_mean_fpr, self.roc_mean_tpr, _ = roc_curve(y_tests, y_preds)
        self.roc_std_auc = np.std(self.aucs)

        self.roc_max_tpr = np.max(self.interp_tprs, axis=0)
        self.roc_min_tpr = np.min(self.interp_tprs, axis=0)

        self.pr_mean_ap = average_precision_score(y_tests, y_preds)
        self.pr_mean_precision, self.pr_mean_recall, _ = precision_recall_curve(y_tests, y_preds)
        self.pr_std_ap = np.std(self.aps)

        self.pr_min_precision = np.min(self.precisions, axis=0)
        self.pr_max_precision = np.max(self.precisions, axis=0)

        y_predicted_labels = np.asarray(y_preds > 0.5, dtype=int)
        self.f1 = f1_score(y_tests, y_predicted_labels)
        self.corr = np.corrcoef(y_tests, y_predicted_labels)[0, 1]

        confusion = confusion_matrix(y_tests, y_predicted_labels)
        precision = confusion[1, 1] / (confusion[1, 1] + confusion[0, 1])
        recall = confusion[1, 1] / (confusion[1, 1] + confusion[1, 0])
        accuracy = (confusion[0, 0] + confusion[1, 1]) / np.sum(confusion)

        return {
            'roc_mean_fpr': self.roc_mean_fpr,
            'roc_mean_tpr': self.roc_mean_tpr,
            'roc_interp_fpr': self.interp_fpr,
            'roc_interp_tprs': self.interp_tprs,
            'roc_min_tpr': self.roc_min_tpr,
            'roc_max_tpr': self.roc_max_tpr,
            'roc_mean_auc': self.roc_mean_auc,
            'roc_std_auc': self.roc_std_auc,
            'pr_mean_recall': self.pr_mean_recall,
            'pr_mean_precision': self.pr_mean_precision,
            'pr_interp_recall': self.interp_recall,
            'pr_interp_precision': self.precisions,
            'pr_mean_ap': self.pr_mean_ap,
            'pr_std_ap': self.pr_std_ap,
            'pr_min_precision': self.pr_min_precision,
            'pr_max_precision': self.pr_max_precision,
            'f1': self.f1,
            'corr': self.corr,
            'precision': precision,
            'recall': recall,
            'accuracy': accuracy
        }
