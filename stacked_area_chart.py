import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the data from the file
file_path = 'chicago_s0058.txt'  # Update this with your actual file path
data = []

# Define lineage groups
grouping = {
    'BA.2.X': ['BA.2', 'BA.2.75'], 
    'BA.5.1.X': ['BA.5.1', 'BA.5.1.10', 'BA.5.1.13', 'BA.5.1.5', 'BA.5.1.21'],
    'BA.5.2.X': ['BA.5.2.4', 'BA.5.2.6', 'BA.5.2.8', 'BA.5.2.13', 'BA.5.2.16', 'BA.5.2.2', 'BA.5.2.21', 'BA.5.2.23', 'BA.5.2.25', 'BA.5.2.28'],
    'BA.4.X': ['BA.4.6', 'BA.4.6.1', 'BA.4.8'],
    'BQ.1.X': ['BQ.1', 'BQ.1.1', 'BQ.1.6', 'BQ.1.3', 'BQ.1.5', 'BQ.1.10', 'BQ.1.12', 'BQ.1.14'],
    'BF.X': ['BF.1', 'BF.2', 'BF.7', 'BF.10', 'BF.11', 'BF.12', 'BF.14', 'BF.26', 'BF.29']
}

with open(file_path, 'r') as file:
    lines = file.readlines()
    current_date = None
    freyja_lineages = {}
    
    for line in lines:
        if line.startswith('Site'):
            # Append the previous date's data before starting the new date
            if current_date:
                data.append({'date': current_date, **freyja_lineages})
            current_date = line.split()[2]
            freyja_lineages = {}  
        elif line.strip() and not line.startswith('Lineage'):
            parts = line.split()
            lineage = parts[0]
            # parts[1] -> freyja, parts[2] -> WEPP
            freyja_abundance = float(parts[1]) if len(parts) > 1 else 0.0
            if freyja_abundance > 0.0:
                freyja_lineages[lineage] = freyja_abundance

    # Append the last date's data after reading all lines
    if current_date:
        data.append({'date': current_date, **freyja_lineages})

# Function to group lineages and sum abundances with custom group names
def group_lineages(data, grouping):
    grouped_data = data.copy()

    for group_name, lineages in grouping.items():
        for lineage in lineages:
            if lineage in grouped_data.columns:
                # Add the lineage abundance to the group's column and then drop the individual lineage
                grouped_data[group_name] = grouped_data.get(group_name, 0) + grouped_data[lineage]
                grouped_data.drop(columns=[lineage], inplace=True)

    return grouped_data

# Convert to a DataFrame
df = pd.DataFrame(data)
df.set_index('date', inplace=True)

## Normalize the data to 100% (as per prevalence in the graph)
#for date, row in df.iterrows():
#    total = sum(row.dropna())  # Sum only non-NaN values
#    for lineage in row.dropna().index:
#        row[lineage] = row[lineage] / total if total > 0 else 0

# Expand the data into a more accessible form for plotting
expanded_data = pd.DataFrame(columns=list(df.columns))
# Create a list to hold dates for proper indexing
dates = []

for date, row in df.iterrows():
    # Convert date to datetime with a fixed year (2022 in this case)
    date = pd.to_datetime(date + "-2022", format='%b-%d-%Y')
    dates.append(date)
    
    # Add the row values to the expanded_data
    for lineage, abundance in row.items():
        expanded_data.at[date, lineage] = abundance

# Set date as the index
expanded_data.index = pd.to_datetime(dates)
# Convert the DataFrame to numeric values, coercing any errors to NaN
expanded_data = expanded_data.apply(pd.to_numeric, errors='coerce')
# Replace NaN values with 0.0
expanded_data = expanded_data.fillna(0.0)

# Apply grouping
expanded_data_grouped = group_lineages(expanded_data, grouping)

# Sort the columns alphabetically to ensure consistency in color assignment
expanded_data_grouped = expanded_data_grouped.sort_index(axis=1)

# Plot stacked area chart
#palette = sns.color_palette("husl", n_colors=len(expanded_data_grouped.columns))
palette = sns.color_palette("tab20", n_colors=len(expanded_data_grouped.columns))
expanded_data_grouped.plot(kind='area', stacked=True, figsize=(12, 8), color=palette)

# Update x-ticks to only show the actual dates from the data
plt.xticks(ticks=expanded_data_grouped.index, labels=expanded_data_grouped.index.strftime('%b-%d'), rotation=45)

plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=8)
plt.title('Freyja detection in Chicago "S0058" Over Time')
plt.ylabel('Lineage Abundance')
plt.xlabel('Date')
plt.tight_layout(pad=4.0)
plt.savefig('chicago_s0058.png')
