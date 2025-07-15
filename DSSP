import os
import gzip
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def replace_char_at_index(original_string, index, new_char):
    if not (0 <= index < len(original_string)):
        return original_string
    if len(new_char) != 1:
        return original_string
    return original_string[:index] + new_char + original_string[index+1:]

ss_types = {
    "BEND": 'S',
    "HELX_LH_PP": 'P',
    "HELX_RH_3T": 'G',
    "HELX_RH_AL": 'H',
    "STRN": 'E',
    "TURN_TY1": 'T'
}

folder_path = r"C:/Users/gayat/OneDrive/Desktop/structure_data"
excel_path = r"C:/Users/gayat/OneDrive/Desktop/Sequences.xlsx"

df = pd.read_excel(excel_path)

df["Predicted_SS"] = ""

for idx, row in df.iterrows():
    uniprot_id = row["UNIPROT_ID"]
    cif_path = os.path.join(folder_path, f"AF-{uniprot_id}-F1-model_v4.cif.gz")

    if not os.path.isfile(cif_path):
        print(f"File not found: {cif_path}")
        continue

    with gzip.open(cif_path, 'rt') as handle:
        mmcif_dict = MMCIF2Dict(handle)

    try:
        seq = mmcif_dict["_entity_poly.pdbx_seq_one_letter_code"][0].replace("\n", "")
        ss = 'C' * len(seq)

        for i in range(len(mmcif_dict["_struct_conf.conf_type_id"])):
            raw_type = mmcif_dict["_struct_conf.conf_type_id"][i]
            norm_type = raw_type.split('_P')[0]  # Remove _P suffix like _P1
            conf_type = ss_types.get(norm_type, 'C')

            start = int(mmcif_dict["_struct_conf.beg_auth_seq_id"][i]) - 1
            end = int(mmcif_dict["_struct_conf.end_auth_seq_id"][i])

            for j in range(start, end):
                ss = replace_char_at_index(ss, j, conf_type)

        df.at[idx, "Predicted_SS"] = ss

    except Exception as e:
        print(f"Error processing {uniprot_id}: {e}")
        df.at[idx, "Predicted_SS"] = "ERROR"

output_path = r"C:/Users/gayat/OneDrive/Desktop/Sequences_with_SS.xlsx"
df.to_excel(output_path, index=False)
print(f"Done! Results saved to: {output_path}")
