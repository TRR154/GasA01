import os

# === AUTO-DETECTED CURRENT DIRECTORY ===
folder_path = os.getcwd()
output_file = "output.txt"

# === STRING PREFIXES ===
prefix1 = "./testdata/gaslib582-v2.net"
prefix2 = "./testdata/" 
prefix3 = "GL582/"

# Get list of files in the folder
files = sorted(os.listdir(folder_path))

# Write formatted output to the text file
with open(output_file, "w") as f:
   for filename in files:
      if os.path.isfile(os.path.join(folder_path, filename)):
         line = f"{prefix1}\t{prefix2} {prefix3}/{filename}\n"
         f.write(line)

print(f"Output written to {output_file} in {folder_path}")
