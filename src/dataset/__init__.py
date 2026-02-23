from src.args import args
from src.dataset.dataset import Dataset

dataset = Dataset(args.dataset)
dataset_name = args.dataset.replace('.csv', '')
culture_name = args.culture