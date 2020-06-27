import yaml
import sys

def make_yaml(definition_file):
    YAML = {
        "GENE": "", 
        "TRANSLATION_TABLE": "", 
        "VCF": "",
        "SAMPLES": []}
    with open(definition_file, "r") as def_in:
        for rawline in def_in.readlines():
            line = rawline.split("\t")
            alleles_raw = line[1].strip()
            if "x" in alleles_raw or "+" in alleles_raw:
                continue
            YAML["SAMPLES"].append({
                "ID": line[0],
                "ALLELES": [sys.argv[2] + a.replace("*", "(star)") for a in alleles_raw.split("/")],
                "NOTES": alleles_raw
            })
    print(yaml.dump(YAML))
    



if __name__ == "__main__":
    make_yaml(sys.argv[1])
