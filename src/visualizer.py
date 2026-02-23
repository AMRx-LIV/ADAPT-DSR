import matplotlib.pyplot as plt

from src.dataset import dataset_name, culture_name


class Visualizer:
    @staticmethod
    def plot_roc_curve(mean_fpr, mean_tpr, interp_fpr, interp_tprs, min_tpr, max_tpr, mean_auc, std_auc, model_name):
        plt.figure(dpi=300)
        plt.plot(mean_fpr, mean_tpr, color='b',
                 label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
                 lw=2, alpha=.8)
        _ = [plt.plot(interp_fpr, interp_tpr, color='b', alpha=.1, lw=0.5) for interp_tpr in interp_tprs]

        plt.fill_between(interp_fpr, min_tpr, max_tpr, color='grey', alpha=.15,
                         label=r'Band')

        plt.plot([0, 1], [0, 1], 'k--')  # Diagonal line
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'{culture_name} {model_name} ROC Curve')
        plt.legend(loc="lower right")
        plt.savefig(f'outputs/roc_{dataset_name}_{model_name}.png', dpi=300)
        plt.close()

    @staticmethod
    def plot_pr_curve(
            mean_recall,
            mean_precision,
            interp_recall,
            interp_precision,
            min_precision,
            max_precision,
            mean_ap,
            std_ap,
            model_name
    ):
        plt.figure(dpi=300)
        plt.plot(mean_recall, mean_precision, color='b',
                 label=r'Mean PR (AP = %0.2f $\pm$ %0.2f)' % (mean_ap, std_ap),
                 lw=2, alpha=.8)

        _ = [plt.plot(interp_recall, ip, color='b', alpha=.1, lw=0.5) for ip in interp_precision]

        plt.fill_between(interp_recall, min_precision, max_precision, color='grey', alpha=.2,
                         label=r'Band')

        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(f'{culture_name} {model_name} Precision-Recall Curve')
        plt.legend(loc="lower left")
        plt.savefig(f'outputs/pr_{dataset_name}_{model_name}.png', dpi=300)
        plt.close()
