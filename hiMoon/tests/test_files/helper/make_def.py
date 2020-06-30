#    Copyright 2020 Solomon M. Adams

#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at

#        http://www.apache.org/licenses/LICENSE-2.0

#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

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
