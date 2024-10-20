# Retrives the original diabetes dataset 
# from sklearn and saves it to a CSV file

import pandas as pd
from sklearn.datasets import load_diabetes

# Load the diabetes dataset
diabetes = load_diabetes()

# Create a DataFrame
df = pd.DataFrame(diabetes.data, columns=diabetes.feature_names)

# Add the target variable
df['target'] = diabetes.target

# Save to CSV
df.to_csv('./dat/ddat_og.csv', index=False)

print("Original diabetes dataset saved to 'diabetes_original.csv'")
print("\nDataset summary:")
print(df.describe())

print("\nFirst few rows:")
print(df.head())