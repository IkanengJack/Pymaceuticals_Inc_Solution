# Pymaceuticals Inc. Analysis

## Analysis Summary

Overall, it is clear that **Capomulin** is a viable drug regimen to reduce tumor growth.  
Capomulin had the most number of mice complete the study. With the exception of Remicane, all other regimens observed a number of mice deaths across the duration of the study.  
There is a strong correlation between mouse weight and tumor volume, indicating that mouse weight may be contributing to the effectiveness of any drug regimen.  
There was one potential outlier within the **Infubinol** regimen. While most mice showed tumor volume increase, there was one mouse that had a reduction in tumor growth in the study.

---

## Dependencies and Setup

```python
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
```

---

## Load Study Data

```python
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

study_data_complete = pd.merge(study_results, mouse_metadata, how="left", on="Mouse ID")
study_data_complete.head()
```

---

## Data Cleaning

```python
# Checking the number of mice
len(study_data_complete["Mouse ID"].unique())

# Check duplicates
duplicate_mouse_ids = study_data_complete.loc[
    study_data_complete.duplicated(subset=['Mouse ID', 'Timepoint']), 'Mouse ID'].unique()
duplicate_mouse_data = study_data_complete.loc[study_data_complete["Mouse ID"] == "g989"]

# Drop duplicates
clean_study_data_complete = study_data_complete[
    study_data_complete['Mouse ID'].isin(duplicate_mouse_ids) == False]
```

---

## Summary Statistics

```python
# Using groupby
means = clean_study_data_complete.groupby('Drug Regimen')['Tumor Volume (mm3)'].mean()
medians = clean_study_data_complete.groupby('Drug Regimen')['Tumor Volume (mm3)'].median()
variances = clean_study_data_complete.groupby('Drug Regimen')['Tumor Volume (mm3)'].var()
sds = clean_study_data_complete.groupby('Drug Regimen')['Tumor Volume (mm3)'].std()
sems = clean_study_data_complete.groupby('Drug Regimen')['Tumor Volume (mm3)'].sem()

summary_table = pd.DataFrame({
    "Mean Tumor Volume": means,
    "Median Tumor Volume": medians,
    "Tumor Volume Variance": variances,
    "Tumor Volume Std. Dev.": sds,
    "Tumor Volume Std. Err.": sems
})
```

Or using one-liner:

```python
summary_table = clean_study_data_complete.groupby("Drug Regimen").agg({
    "Tumor Volume (mm3)": ["mean", "median", "var", "std", "sem"]
})
```

---

## Bar and Pie Charts

```python
# Bar plots
counts = clean_study_data_complete['Drug Regimen'].value_counts()
counts.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.ylabel("# of Observed Mouse Timepoints")
plt.xticks(rotation=90)
plt.show()

# Pie plots (Pandas and Pyplot)
mice_df = clean_study_data_complete.loc[:, ["Mouse ID", "Sex"]].drop_duplicates()
counts = mice_df.Sex.value_counts()
counts.plot(kind="pie", autopct='%1.1f%%')
plt.show()

plt.pie(counts.values, labels=counts.index.values, autopct='%1.1f%%')
plt.ylabel("count")
plt.show()
```

---

## Quartiles, Outliers and Boxplots

```python
max_tumor = clean_study_data_complete.groupby(["Mouse ID"])['Timepoint'].max().reset_index()
merged_data = max_tumor.merge(clean_study_data_complete, on=['Mouse ID','Timepoint'], how="left")
treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

tumor_vol_list = []

for drug in treatment_list:
    final_tumor_vol = merged_data.loc[merged_data["Drug Regimen"] == drug, 'Tumor Volume (mm3)']
    tumor_vol_list.append(final_tumor_vol)
    quartiles = final_tumor_vol.quantile([.25,.5,.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq - lowerq
    lower_bound = lowerq - (1.5 * iqr)
    upper_bound = upperq + (1.5 * iqr)
    outliers = final_tumor_vol.loc[(final_tumor_vol < lower_bound) | (final_tumor_vol > upper_bound)]
    print(f"{drug}'s potential outliers: {outliers}")
```

Boxplot:

```python
orange_out = dict(markerfacecolor='red', markersize=12)
plt.boxplot(tumor_vol_list, labels=treatment_list, flierprops=orange_out)
plt.ylabel('Final Tumor Volume (mm3)')
plt.show()
```

---

## Line and Scatter Plots

```python
# Line plot
capomulin_table = clean_study_data_complete.loc[clean_study_data_complete['Drug Regimen'] == "Capomulin"]
mousedata = capomulin_table.loc[capomulin_table['Mouse ID'] == 'l509']
plt.plot(mousedata['Timepoint'], mousedata['Tumor Volume (mm3)'])
plt.xlabel('Timepoint (days)')
plt.ylabel('Tumor Volume (mm3)')
plt.title('Capomulin treatment of mouse l509')
plt.show()

# Scatter plot
capomulin_average = capomulin_table.groupby(['Mouse ID'])[['Weight (g)', 'Tumor Volume (mm3)']].mean()
plt.scatter(capomulin_average['Weight (g)'], capomulin_average['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.show()
```

---

## Correlation and Regression

```python
corr = round(st.pearsonr(
    capomulin_average['Weight (g)'],
    capomulin_average['Tumor Volume (mm3)']
)[0], 2)
print(f"The correlation between mouse weight and the average tumor volume is {corr}")

model = st.linregress(
    capomulin_average['Weight (g)'],
    capomulin_average['Tumor Volume (mm3)']
)

y_values = capomulin_average['Weight (g)'] * model[0] + model[1]
plt.scatter(capomulin_average['Weight (g)'], capomulin_average['Tumor Volume (mm3)'])
plt.plot(capomulin_average['Weight (g)'], y_values, color="red")
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.show()
```

> **The correlation between mouse weight and the average tumor volume is 0.84**
