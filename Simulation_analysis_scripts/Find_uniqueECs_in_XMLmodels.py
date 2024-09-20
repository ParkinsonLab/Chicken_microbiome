import os
import xml.etree.ElementTree as ET
import pandas as pd
import re
import numpy

# Function to extract EC numbers from an XML file
def extract_ec_numbers(xml_file):
    ec_numbers = set()
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Search for the pattern "ec-code:" followed by digits and dots in the XML file
    ec_pattern = re.compile(r"ec-code:(\d+\.\d+\.\d+\.\d+)")
    for match in ec_pattern.findall(ET.tostring(root).decode()):
        ec_numbers.add(match)

    return ec_numbers

# Define the folder containing XML files
xml_folder_models = "models_sbml"

# Read the table into a pandas DataFrame
table_abund = pd.read_table("Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt", delimiter='\t')

ec_numbers_models = {}

# Iterate over rows in the DataFrame
for model in table_abund["Specie"]:
        xml_file = os.path.join(xml_folder_models,
                                model + "_gapFilled-adapt_EXadded.xml")

        # Check if the XML file exists
        if os.path.exists(xml_file):
            ec_numbers_in_xml = extract_ec_numbers(xml_file)

            ec_numbers_models[model] = ec_numbers_in_xml

# Save the ec numbers for each model to a  file
with open("Table_ECnumbers_inEachModel_gapfilled_corn_absorbed_2fibre_May2024.txt", 'w') as file:
    # Write the header
    file.write("SpecieName\tEC_numbers\n")

    # Iterate over each item in the dictionary
    for species, ec_numbers in ec_numbers_models.items():
        # Join the EC numbers into a comma-separated string
        ec_numbers_string = ', '.join(sorted(ec_numbers))
        # Write the species name and EC numbers to the file
        file.write(f"{species}\t{ec_numbers_string}\n")