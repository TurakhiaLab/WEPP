import sys
import pandas as pd
from math import ceil
from collections import defaultdict

#Check arguments
if len(sys.argv) != 2:
    print("USAGE: python src/WBE/ivar_correction.py <directory>")
    sys.exit(1)

dir_name = sys.argv[1]

# Read the TSV file
variants_file = dir_name + "/cwap_variants.tsv"
depth_file = dir_name + "/cwap_depth.tsv"

df_variants = pd.read_csv(variants_file, sep='\t')
df_depth = pd.read_csv(depth_file, sep='\t', header=None, names=['Region', 'POS', 'Base', 'Value'])

# Convert the second column (POS) to the key and the last column (Value) to the value in a dictionary
site_depth_counts = dict(zip(df_depth['POS'], df_depth['Value']))

# Keep only specified columns
df_variants = df_variants[['POS', 'REF', 'ALT', 'ALT_FREQ', 'TOTAL_DP']]

# Dictionary to store site-wise deletion counts
site_del_counts = defaultdict(int)

# Find rows where 'ALT' contains "-"
rows_with_deletions = df_variants[df_variants['ALT'].str.contains("-", na=False)]
indices_to_drop = rows_with_deletions.index
df_variants = df_variants.drop(indices_to_drop)
df_variants = df_variants.reset_index(drop=True)

# Iterate through rows_with_deletions
for _, row in rows_with_deletions.iterrows():
    pos = int(row['POS'])
    alt = row['ALT'] 
    alt_freq = float(row['ALT_FREQ'])
    total_dp = int(row['TOTAL_DP'])
    
    # Calculate allele count
    allele_count = ceil(alt_freq * total_dp)
    
    # Extract the deletion sequence without the leading "-"
    deletion_sequence = alt[1:]
    
    # Iterate through the deletion sequence to calculate site positions
    for i in range(len(deletion_sequence)):
        site = pos + i + 1
        # Storing ref_nuc + site in site
        site_del_counts[deletion_sequence[i] + str(site)] += allele_count


new_rows = []  # List to store new rows
insert_positions = []  # List to store positions where new rows should be inserted
# Iterate through each site in site_del_counts and look for matching POS in df_variants
for mut, del_count in site_del_counts.items():
    # Find all rows with matching POS in df_variants
    matching_rows = df_variants[df_variants['POS'] == int(mut[1:])]
    
    allele_counts = defaultdict(int)
    total_site_count = site_depth_counts[int(mut[1:])]

    # Iterate through the matching rows to update allele counts
    if not matching_rows.empty:
        # Iterate through the matching rows
        for _, row in matching_rows.iterrows():
            alt = row['ALT']
            alt_freq = float(row['ALT_FREQ'])
            total_dp = int(row['TOTAL_DP'])

            # Update allele counts
            allele_counts[alt] = ceil(alt_freq * total_dp)
    
        # Update the matching rows with new ALT_FREQ values and TOTAL_DP
        for index, row in matching_rows.iterrows():
            # Update ALT_FREQ
            row['ALT_FREQ'] = min(allele_counts[row['ALT']] / total_site_count, 1.0)
            # Update TOTAL_DP
            row['TOTAL_DP'] = total_site_count

            # Apply the changes to the dataframe
            df_variants.loc[index, 'ALT_FREQ'] = row['ALT_FREQ']  # Update the ALT_FREQ in df_variants
            df_variants.loc[index, 'TOTAL_DP'] = row['TOTAL_DP']  # Update the TOTAL_DP in df_variants
    
    # Add a row for the deletion at that site with ALT = "-"
    new_row = {
        'POS': mut[1:],
        'REF': mut[0],
        'ALT': '-',
        'ALT_FREQ': min(del_count / total_site_count, 1.0) if total_site_count > 0 else 0.0,
        'TOTAL_DP': total_site_count
    }
    
    # Check if there are any matching rows
    if matching_rows.empty:
        # If no matching rows, determine the insertion position
        smaller_pos_rows = df_variants[df_variants['POS'] < int(mut[1:])]
        insert_position = 0 if smaller_pos_rows.empty else df_variants.index.get_loc(smaller_pos_rows.index[-1]) + 1
    else:
        # If matching_rows is not empty, insert after the last matching row
        insert_position = df_variants.index.get_loc(matching_rows.index[-1]) + 1

    # Collect new row and insert position
    new_rows.append(new_row)
    insert_positions.append(insert_position)

# Sort insert positions and corresponding new rows
insert_positions_sorted = sorted(zip(insert_positions, new_rows), key=lambda x: x[0])

# Insert rows at desired positions
offset = 0  # Offset adjusts positions as rows are inserted
for position, new_row in insert_positions_sorted:
    actual_position = position + offset
    top = df_variants.iloc[:actual_position]
    bottom = df_variants.iloc[actual_position:]
    new_row_df = pd.DataFrame([new_row])
    df_variants = pd.concat([top, new_row_df, bottom], ignore_index=True)
    offset += 1  # Increment offset as a new row is added

# Save the updated df_variants to a file
df_variants.to_csv(dir_name + "/corrected_cwap_variants.tsv", sep='\t', index=False)