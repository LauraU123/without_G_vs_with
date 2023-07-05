'''
This part of the workflow expects input files
            sequences = "data/sequences.fasta"
            metadata = "data/metadata.tsv"
'''
with_or_without = ["with", "without"]
rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta"
    output:
        sequence_index = build_dir + "/{a_or_b}/{build_name}/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule newreference:
    message:
        """
        Making new reference
        """
    input:
        oldreference = "config/{a_or_b}reference.gbk"
    output:
        newreferencegbk = build_dir + "/{a_or_b}/{build_name}/newreference.gbk",
        newreferencefasta = build_dir + "/{a_or_b}/{build_name}/newreference.fasta",
        greference = build_dir + "/{a_or_b}/{build_name}/greference.fasta"
    params:
        gene = lambda w: w.build_name,
        newreference = build_dir + "/{a_or_b}/{build_name}/newreference",
        oldreference = 'config/{a_or_b}reference',
        greference = build_dir + "/{a_or_b}/{build_name}/greference"
    shell:
        """
        python scripts/newreference.py \
            --greference {params.greference} \
            --reference {params.oldreference} \
            --output {params.newreference} \
            --gene {params.gene}
        """


rule filter:
    message:
        """
        Aligning sequences to {input.reference}
            - gaps relative to reference are considered real
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta",
        reference = "config/{a_or_b}reference.gbk",
        metadata = "data/{a_or_b}/metadata.tsv",
        sequence_index = rules.index_sequences.output
    output:
    	sequences = build_dir + "/{a_or_b}/{build_name}/filtered.fasta"
    params:
    	group_by = config["filter"]["group_by"],
    	min_length = lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
    	subsample_max_sequences = config["filter"]["subsample_max_sequences"],
    	strains = "config/dropped_strains.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --exclude {params.strains} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = "config/{a_or_b}reference.fasta"
    output:
        alignment = build_dir + "/{a_or_b}/{build_name}/sequences.aligned.fasta",
        insertionsfile = build_dir + "/{a_or_b}/{build_name}/insertions.csv"
    threads: 4
    shell:
        """
        nextalign run -j {threads}\
            --reference {input.reference} \
            --output-fasta {output.alignment} \
            --output-insertions {output.insertionsfile} \
            {input.sequences}
        """

rule cut_G_out:
    input:
        oldalignment = rules.align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        newalignment = build_dir + "/{a_or_b}/{build_name}/aligned_without_G.fasta"
    shell:
        """
        python scripts/cut_G_out.py \
            --oldalignment {input.oldalignment} \
            --newalignment {output.newalignment} \
            --reference {input.reference} \
        """
        
rule cut:
    input:
        oldalignment = rules.align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        newalignment = build_dir + "/{a_or_b}/{build_name}/newalignment.fasta"
    shell:
        """
        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --newalignment {output.newalignment} \
            --reference {input.reference} \
        """

rule realign:
    input:
        newalignment = rules.cut.output.newalignment,
        reference = rules.newreference.output.greference
    output:
        realigned = build_dir + "/{a_or_b}/{build_name}/realigned.fasta"
    threads: 4
    shell:
        """
        augur align --nthreads {threads} \
            --sequences {input.newalignment} \
            --reference-sequence {input.reference} \
            --output {output.realigned}
        """

rule alignment_for_tree:
    input:
        realigned = rules.realign.output.realigned,
        original = rules.align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        aligned_for_tree = build_dir + "/{a_or_b}/{build_name}/alignment_for_tree.fasta"
    params:
        gene = lambda w: w.build_name
    shell:
        """
        python scripts/align_for_tree.py \
            --realign {input.realigned} \
            --original {input.original} \
            --reference {input.reference} \
            --output {output.aligned_for_tree} \
            --gene {params.gene}
        """

rule tree:
    message: "Building tree without G gene"
    input:
        alignment = rules.cut_G_out.output.newalignment
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree_raw_without.nwk"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
        """

rule tree_with_G:
    message: "Building tree with G included"
    input:
        alignment = rules.alignment_for_tree.output
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree_raw_with.nwk"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
        """


rule tree_untangle:
    message:
        """untangling tree branches"""
    input:
        tree_with_G = build_dir + "/{a_or_b}/{build_name}/tree_raw_with.nwk",
        tree_without_G = build_dir + "/{a_or_b}/{build_name}/tree_raw_without.nwk"
    output:
        untangled_with_G = build_dir + "/{a_or_b}/{build_name}/tree_with_resolve.nwk",
        untangled_without_G = build_dir + "/{a_or_b}/{build_name}/tree_without_resolve.nwk"
    shell:
        """
        julia "scripts/treeknit.jl" \
        --treel {input.tree_without_G} \
        --treerest {input.tree_with_G} \
        --outputl {output.untangled_without_G} \
        --outputrest {output.untangled_with_G}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
        """
    input:
        tree = build_dir + "/{a_or_b}/{build_name}/tree_{with_or_without}_resolve.nwk",
        alignment = rules.tree.input.alignment,
        metadata = rules.filter.input.metadata
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree_{with_or_without}_refined.nwk",
        node_data = build_dir + "/{a_or_b}/{build_name}/branch_lengths_{with_or_without}.json"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --no-covariance \
            --keep-polytomies \
            --keep-root 
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree = build_dir + "/{a_or_b}/{build_name}/tree_{with_or_without}_refined.nwk",
        alignment = rules.alignment_for_tree.output.aligned_for_tree
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/{with_or_without}nt_muts.json"
    params:
    	inference = config["ancestral"]["inference"]
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --keep-ambiguous
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree =  build_dir + "/{a_or_b}/{build_name}/tree_{with_or_without}_refined.nwk",
        node_data = rules.ancestral.output.node_data,
        reference = rules.newreference.output.newreferencegbk
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/{with_or_without}aa_muts.json",
        aa_data = build_dir + "/{a_or_b}/{build_name}/{with_or_without}alignedG.fasta"
    params:
    	alignment_file_mask = build_dir + "/{a_or_b}/{build_name}/{with_or_without}aligned%GENE.fasta"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
            --alignment-output {params.alignment_file_mask}
        """

rule traits:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.filter.input.metadata
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/{with_or_without}traits.json"
    log:
        "logs/{a_or_b}/traits_{build_name}_{with_or_without}rsv.txt"
    params:
    	columns = config["traits"]["columns"]
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "config/clades_G_{a_or_b}.tsv"
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/{with_or_without}clades_G.json"
    log:
        "logs/{a_or_b}/clades_{build_name}{with_or_without}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """
