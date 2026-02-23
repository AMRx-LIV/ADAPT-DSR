import argparse

parser = argparse.ArgumentParser(description='Train and evaluate models.')
parser.add_argument('--dataset', type=str, help='Path to the dataset file', default='preprocessed_blood2.csv')
parser.add_argument('--culture', type=str, help='Culture name, e.g., urine or blood', default='unknown')
args = parser.parse_args()