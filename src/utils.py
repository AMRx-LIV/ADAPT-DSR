import numpy as np
from sklearn.metrics import precision_recall_curve


def calc_best_threshold(labels, predictions):
    # Find the best threshold
    precisions, recalls, thresholds = precision_recall_curve(labels, predictions)
    f1_scores = np.divide(2 * (precisions * recalls), (precisions + recalls), out=np.zeros_like(precisions),
                          where=(precisions + recalls) != 0)
    threshold = thresholds[np.argmax(f1_scores[:-1])]  # Exclude the last threshold

    return threshold
