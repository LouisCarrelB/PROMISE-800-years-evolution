import json

def count_species_in_json(file_path):
    """
    Load a JSON file and count the number of key-value pairs in it.
    Args:
        file_path (str): Path to the JSON file.
    Returns:
        int: Number of key-value pairs in the JSON dictionary.
    """
    try:
        with open(file_path, 'r') as json_file:
            data = json.load(json_file)
            return len(data)
    except FileNotFoundError:
        print(f"File not found: {file_email}")
        return 0
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the file: {file_path}")
        return 0

# Example usage:
file_path = 'DATA/ENSG00000107643/Dict_species.json'  # Ensure this is the correct relative or absolute path
print(count_species_in_json(file_path))  # Explicitly print the result
