import pandas as pd
import sys
import re

def format_species_name(species_name):
    # Replace any non-alphanumeric characters (excluding underscore) with a single underscore
    species_name = re.sub(r'[^A-Za-z0-9_]', '_', species_name)
    # Remove leading and trailing underscores
    species_name = species_name.strip('_')
    parts = species_name.split("_")
    if len(parts) >= 2:
        return f"{parts[0].capitalize()}_{'_'.join(parts[1:]).lower()}"
    return species_name.capitalize()

def remove_underscore(species_name):
    return species_name.replace("_", " ")

def process_files(file1, file2, output_file, output_species_file):
    # Read the input CSV files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Combine the dataframes
    combined_df = pd.concat([df1, df2])

    # Format the species names
    combined_df['Species'] = combined_df['Species'].apply(format_species_name)

    # Create a pivot table to handle BOTH cases
    combined_df = combined_df.pivot_table(index='Species', values='Index', aggfunc=lambda x: 'BOTH' if 'CAN' in x.values and 'ALT' in x.values else x.values[0]).reset_index()

    # Save the combined dataframe to a TSV file
    combined_df.to_csv(output_file, sep='\t', index=False)
    print("Combined DataFrame:\n", combined_df)
    
    # Prepare species names without underscores for the second TSV file
    combined_df['Species_No_Underscore'] = combined_df['Species'].apply(remove_underscore)
    
    # Save the species column without underscores to a separate TSV file
    combined_df[['Species_No_Underscore']].to_csv(output_species_file, sep='\t', index=False, header=['Species'])
    print("Species DataFrame without underscores:\n", combined_df[['Species_No_Underscore']])

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python process_species_files.py <file1.csv> <file2.csv> <output_file.tsv> <output_species_file.tsv>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]
    output_species_file = sys.argv[4]

    process_files(file1, file2, output_file, output_species_file)
    print(f"Combined species information saved to {output_file}")
    print(f"Species list saved to {output_species_file}")
