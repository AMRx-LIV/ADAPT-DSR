import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, KFold


class Dataset:
    def __init__(
            self,
            path: str,
    ):
        self.max = None
        self.min = None
        self.path = path
        self.raw_data = None
        self.unique_subject_ids = None
        self._filter_data()

        self.data_folded = None
        self.labels = None
        self.X_eval = None
        self.y_eval = None
        self.X_holdout = None
        self.y_holdout = None
        self.holdout_ori = None
        self.subject_ids_folded = None

    def standardize(self, X):
        """
        Standardize the data columns to have mean 1 and std 1
        """
        mean = np.nanmean(X, axis=0)
        std = np.nanstd(X, axis=0)

        # Standardize the data
        X = (X - mean) / std

        return X

    def column_min_max_scale(
            self,
            X
    ):
        """
        Scale the data columns to have min 0 and max 1.
        """
        if self.min is None:
            self.min = np.nanmin(X, axis=0)
        if self.max is None:
            self.max = np.nanmax(X, axis=0)

        # Scale the data
        for i in range(X.shape[1]):
            if self.max[i] != self.min[i]:  # Avoid division by zero
                X[:, i] = (X[:, i] - self.min[i]) / (self.max[i] - self.min[i])
            else:
                X[:, i] = 0  # or some other appropriate value when min equals max

        return X

    def column_mean_fill(
            self,
            X
    ):
        """
        Get the mean of each column.
        """
        for i in range(X.shape[1]):
            col = X[:, i]
            mean = np.nanmean(col)
            col[np.isnan(col)] = mean
            X[:, i] = col

        return X

    def column_sample_fill(
            self,
            X
    ):
        """
        Get the mean of each column.
        """
        for i in range(X.shape[1]):
            col = X[:, i]
            # randomly sample a non-nan from the column
            sample = np.random.choice(col[~np.isnan(col)])
            col[np.isnan(col)] = sample
            X[:, i] = col

        return X

    def summarize_data(self, df):
        # Assuming df is your DataFrame
        ages = df['anchor_age']
        male = df['male']
        race = df['race']
        org_genus = df['org_genus']

        # Total number of records
        total = len(df)

        # Age summary statistics
        age_median = ages.median()
        age_q1 = ages.quantile(0.25)
        age_q3 = ages.quantile(0.75)
        print(f"Age, years, median (IQR): {age_median:.0f} ({age_q1:.0f}–{age_q3:.0f})\n")

        # Gender female counts and percentages
        num_female = (male == 0).sum()
        perc_female = num_female / total * 100
        print(f"Gender female (F), n (%): {num_female} ({perc_female:.2f}%)\n")

        # Race counts and percentages
        print("Race, n (%):")
        race_counts = race.value_counts()
        race_percent = race_counts / total * 100
        race_summary = pd.DataFrame({'Count': race_counts, 'Percentage': race_percent.round(2)})
        print(race_summary.to_string())
        print()

        # Organism grown counts and percentages
        print("Organism grown, n (%):")
        org_counts = org_genus.value_counts()
        org_percent = org_counts / total * 100
        org_summary = pd.DataFrame({'Count': org_counts, 'Percentage': org_percent.round(2)})
        print(org_summary.to_string())

    def _filter_data(
            self,
    ):
        df = pd.read_csv(self.path, low_memory=False)

        # shuffle data
        df = df.sample(frac=1, random_state=42)
        filtered_data_part = df.iloc[:, 102:]
        other_part = df.iloc[:, :102]

        # remove ur_bl_org column
        filtered_data_part = filtered_data_part.drop('ur_bl_org', axis=1)

        # summary statistics
        self.summarize_data(df)

        # find the column that contains 'UNABLE TO REPORT'
        unable_to_report = filtered_data_part.columns[filtered_data_part.isin(['Dyspnea']).any()].tolist()
        o = filtered_data_part['disposition']

        # check again if there is any 'UNABLE TO REPORT' in o
        o = o[o.isin(['ERROR'])]

        # convert admission_location to one-hot encoding
        filtered_data_part = pd.get_dummies(filtered_data_part, columns=[
            'adm_loc', 'curr_service', 'arrival_transport', 'disposition', 'chiefcomplaint', 'race'
        ])

        filtered_data_part = filtered_data_part.replace('UNABLE TO REPORT', None)
        filtered_data_part = filtered_data_part.replace('ERROR', None)
        filtered_data_part = filtered_data_part.astype(float)

        # filtered_data_part = filtered_data_part.fillna(0)

        # combine with other part
        raw = pd.concat([other_part, filtered_data_part], axis=1)

        self.raw_data = raw
        self.unique_subject_ids = raw['subject_id'].unique()

    # def load_classification(self, subset=False):
    #     # count nan values
    #     o = filtered.isnull().sum().sum()
    #
    #     X = filtered.values[:, 1:]
    #     y = filtered.values[:, 0]
    #
    #     return *train_test_split(X, y, stratify=y, test_size=0.15, random_state=42), filtered.columns[1:]

    def load_classification_k_fold(self, n_splits=5):
        filtered = self.raw_data

        eval_ids, holdout_ids = train_test_split(self.unique_subject_ids, test_size=0.2, random_state=42)
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

        data_folded = []

        # Instead of storing just "folded subject IDs" by unique ID,
        # we'll store row-level subject IDs in each fold.
        subject_ids_folded = []

        # Build your hold-out set
        X_eval = filtered[filtered['subject_id'].isin(eval_ids)].iloc[:, 103:].values
        X_eval = self.column_min_max_scale(X_eval)
        # X_eval[np.isnan(X_eval)] = 0.5
        y_eval = filtered[filtered['subject_id'].isin(eval_ids)].iloc[:, 102].values

        for train_index, test_index in kf.split(eval_ids):
            # Identify which subjects go into train vs. test
            train_subjects = eval_ids[train_index]
            test_subjects = eval_ids[test_index]

            # Filter DataFrame down to rows in each fold
            train_filter = filtered[filtered['subject_id'].isin(train_subjects)]
            test_filter = filtered[filtered['subject_id'].isin(test_subjects)]

            # Extract X and y
            X_eval_train = train_filter.iloc[:, 103:].values
            y_eval_train = train_filter.iloc[:, 102].values
            subject_ids_train = train_filter['subject_id'].values

            X_eval_test = test_filter.iloc[:, 103:].values
            y_eval_test = test_filter.iloc[:, 102].values
            subject_ids_test = test_filter['subject_id'].values

            # Scale data and fill NaNs
            X_eval_train = self.column_min_max_scale(X_eval_train)
            # X_eval_train[np.isnan(X_eval_train)] = 0.5

            X_eval_test = self.column_min_max_scale(X_eval_test)
            # X_eval_test[np.isnan(X_eval_test)] = 0.5

            # Append the data as a tuple including the subject IDs
            data_folded.append(
                (X_eval_train, X_eval_test, y_eval_train, y_eval_test)
            )

            # Append the subject IDs
            subject_ids_folded.append(subject_ids_test)

        X_holdout = filtered[filtered['subject_id'].isin(holdout_ids)].iloc[:, 103:].values
        X_holdout = self.column_min_max_scale(X_holdout)
        # X_holdout[np.isnan(X_holdout)] = 0.5
        y_holdout = filtered[filtered['subject_id'].isin(holdout_ids)].iloc[:, 102].values
        holdout_ori = filtered[filtered['subject_id'].isin(holdout_ids)]

        # Save everything to self
        self.data_folded = data_folded
        self.subject_ids_folded = subject_ids_folded  # unique IDs, optional
        self.labels = filtered.columns[103:]
        self.X_eval = X_eval
        self.y_eval = y_eval
        self.X_holdout = X_holdout
        self.y_holdout = y_holdout
        self.holdout_ori = holdout_ori

        return (data_folded,
                filtered.columns[103:],
                X_eval,
                y_eval,
                X_holdout,
                y_holdout,
                holdout_ori)


if __name__ == '__main__':
    dataset = Dataset('../preprocessed_blood.csv')
    o = dataset.load_classification()
    pass
