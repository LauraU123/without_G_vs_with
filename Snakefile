configfile: "config/configfile.yaml"

build_dir = 'results'
auspice_dir = 'auspice'

rule all:
    input:
        expand("auspice/rsv_{subtype}_{build}_{with_or_without}.json",
               subtype = config.get("subtypes",['a']),
               build = config.get("buildstorun", ['genome']),
               with_or_without = config["with_or_without"])

include: "workflow/snakemake_rules/core.smk"

include: "workflow/snakemake_rules/export.smk"

#include: "workflow/snakemake_rules/download.smk"

#include: "workflow/snakemake_rules/glycosylation.smk"

#include: "workflow/snakemake_rules/clades.smk"

rule clean:
    params:
        targets = ["auspice", "results"]
    shell:
        """
        rm -rf {params.targets}
        """

rule clobber:
    params:
        targets = ["data", "auspice", "results"]
    shell:
        """
        rm -rf {params.targets}
        """
