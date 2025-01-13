from os import listdir 
from os.path import join

build_inps = []
for root, _, files in os.walk("src"):
    for f in files:
        build_inps.append(join(root, f))

rule install:
    input:
        "CMakeLists.txt"
    output:
        "build/Makefile"
    shell:
        "./workflow/scripts/install.sh"

rule build_wbe:
    input:
        "build/Makefile",
        build_inps
    output:
        "build/wbe"
    shell:
        "cd build && make -j"

rule freyja_preprocess:
    input:
        "data/{dataset}/{file_prefix}_alignment.sam"
    output:
        "data/{dataset}/{file_prefix}_variants.tsv",
        "data/{dataset}/{file_prefix}_depth.tsv"
    shell:
        '''
            # some conda bug cant run in strict mode?
            # https://github.com/conda/conda/issues/9966
            set +u
            source ~/miniconda3/etc/profile.d/conda.sh
            conda activate freyja-env
            cp {input} ./Freyja/my_vcf_alignment.sam
            sed -i 's/NC_045512v2/NC_045512.2/g' ./Freyja/my_vcf_alignment.sam
            samtools view -Sb ./Freyja/my_vcf_alignment.sam > ./Freyja/my_vcf_alignment.bam
            samtools sort ./Freyja/my_vcf_alignment.bam -o ./Freyja/resorted.bam 
            freyja variants ./Freyja/resorted.bam --variants ./data/{wildcards.dataset}/{wildcards.file_prefix}_variants.tsv --depths ./data/{wildcards.dataset}/{wildcards.file_prefix}_depth.tsv
        '''
