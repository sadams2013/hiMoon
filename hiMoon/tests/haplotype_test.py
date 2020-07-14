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
import glob
import os
from parameterized import parameterized

from hiMoon import gene, vcf, subject, himoon, get_config

PATH = os.path.dirname(os.path.abspath(__file__))

CONFIG = get_config()

def load_definition_file(file_path: str) -> dict:
    """
    Wrapper for yaml.load

    Args:
        file_path (str): path/to/yml/file

    Returns:
        dict: dict representation of yml
    """
    with open(file_path, "r") as sample_allele_data:
        try:
            return(yaml.load(sample_allele_data, Loader = yaml.CLoader))
        except AttributeError: # Happens if libyaml headers are not available
            return(yaml.load(sample_allele_data, Loader = yaml.Loader))

def prep_samples() -> []:
    gene_samples = []
    definition_files = glob.glob(PATH + "/test_files/*.yaml")
    for definition_file in definition_files:
        sample_alleles = load_definition_file(definition_file)
        vcf_file = vcf.VarFile(PATH + "/test_files/" + sample_alleles["VCF"])
        gene_obj = gene.AbstractGene(PATH + "/test_files/" + sample_alleles["TRANSLATION_TABLE"], vcf = vcf_file, config = CONFIG)
        gene_samples += [(subject.Subject(sample["ID"], genes = [gene_obj], config = CONFIG), gene_obj, sample) for sample in sample_alleles["SAMPLES"]]
    return(gene_samples)

def rm_sub_allele(allele: str) -> str:
    """
    Remove suballele from a called named allele.
    (star)2.001 -> (star)2
    (star)2.001_x2 -> (star)2x2

    Args:
        allele (str): [description]

    Returns:
        str: [description]
    """
    sv = allele.split("_")[-1] if "_" in allele else None
    main_allele = allele.split(".")[0]
    return "".join([main_allele, sv]) if sv else main_allele

@parameterized.expand(prep_samples())
def test_sample_haplotypes(subj, gene_obj, sample):
    possible_haplotypes = []
    haps = subj.called_haplotypes[str(gene_obj)]["HAPS"][0]
    for hap in haps:
        possible_haplotypes.append(sorted([rm_sub_allele(a) for a in hap]))
    print(f"{sample['ID']} @ {str(gene_obj)} should be {'/'.join(sample['ALLELES'])}. Found {possible_haplotypes}")
    assert sorted(sample["ALLELES"]) in possible_haplotypes
    