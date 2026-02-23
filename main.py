# main.py
from src.dataset import dataset
from src.evaluator import Evaluator
from src.model_trainer import ModelTrainer
from src.models.lr import LogisticRegression
from src.models.mlp import MLP
from src.models.rf import RandomForest
from src.models.svm import SVM
from src.models.xgb import XGBoost
from src.visualizer import Visualizer


def main():
    # Load the data
    data_folded, labels, X_eval, y_eval, X_holdout, y_holdout, holdout_ori = dataset.load_classification_k_fold(10)

    # models = [RandomForest, MLP, XGBoost, LogisticRegression, SVM]
    models = [XGBoost]

    for model_class in models:
        # Train models
        trainer = ModelTrainer(model_class, data_folded, X_eval, y_eval)
        trainer.train_on_folds()
        print("-" * 10 + " " + "Fairness" + " " + "-" * 10)
        trainer.fairness_metrics()
        trainer.train_on_full()
        metrics = trainer.compute_mean_metrics()

        print(f'ROC-AUC: {metrics["roc_mean_auc"]:.2f}')
        print(f'Average Precision: {metrics["pr_mean_ap"]:.2f}')
        print(f'F1: {metrics["f1"]:.2f}')
        print(f'Correlation: {metrics["corr"]:.2f}')
        print(f'Precision: {metrics["precision"]:.2f}')
        print(f'Recall: {metrics["recall"]:.2f}')
        print(f'Accuracy: {metrics["accuracy"]:.2f}')

        # Evaluate models
        evaluator = Evaluator(
            trainer.best_model_from_folds,
            trainer.bast_model_from_folds_threshold,
            trainer.model_from_x_eval,
            labels,
        )

        evaluator.evaluate_holdout(X_holdout, y_holdout, holdout_ori, phase='best')
        evaluator.compute_feature_importance()
        evaluator.compute_shap_values(X_eval)
        evaluator.evaluate_holdout(X_holdout, y_holdout, holdout_ori, phase='whole')

        # Visualize results
        Visualizer.plot_roc_curve(
            metrics['roc_mean_fpr'],
            metrics['roc_mean_tpr'],
            metrics['roc_interp_fpr'],
            metrics['roc_interp_tprs'],
            metrics['roc_min_tpr'],
            metrics['roc_max_tpr'],
            metrics['roc_mean_auc'],
            metrics['roc_std_auc'],
            model_class().printing_name
        )

        Visualizer.plot_pr_curve(
            metrics['pr_mean_recall'],
            metrics['pr_mean_precision'],
            metrics['pr_interp_recall'],
            metrics['pr_interp_precision'],
            metrics['pr_min_precision'],
            metrics['pr_max_precision'],
            metrics['pr_mean_ap'],
            metrics['pr_std_ap'],
            model_class().printing_name
        )


if __name__ == "__main__":
    main()
