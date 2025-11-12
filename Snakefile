import pandas as pd
import os
import glob
import re

configfile: "config.yaml"

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
MANIFEST = config.get("MANIFEST", "manifest.tab")
REF_SITE = config.get("REF_SITE", f"{SNAKEMAKE_DIR}/db_source/human_sites_n10.fa")
SITE_CENTER = config.get("SITE_CENTER", f"{SNAKEMAKE_DIR}/db_source/human_sites_center.txt")
ROTATION_MATRIX = config.get("ROTATION_MATRIX", f"{SNAKEMAKE_DIR}/db_source/human_sites_rotationalMatrix.tsv")
EXTERNAL_COUNTS_DIR = config.get("EXTERNAL_COUNTS_DIR", "")
COUNT_FILE_EXP = config.get("COUNT_FILE_EXP", "count")
COMPARE_ONLY_EXT = config.get("COMPARE_ONLY_EXT", False)
ALN_PARAMS = config.get('ALN_PARAMS', '')
ALN_REF = "/net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa" # fixed
SVDPREFIX = config.get("SVDPREFIX", "/net/eichler/vol26/7200/software/pipelines/vbi-smk/resource_files/finalout.vcf.gz")

command_dict = {}
command_dict['ONT'] = 'minimap2 -ax map-ont -I 8G -t'
command_dict['Illumina'] = 'bwa mem -Y -K 100000000 -t'
command_dict['PacBio'] = 'minimap2 -ax map-hifi -I 8G -t'

df = pd.read_csv(MANIFEST, sep="\t", header=0).set_index("ID", drop=True)

def find_fofn(wildcards):
    return df.loc[wildcards.id, "FOFN"]

def get_external_counts(wildcards):
    return sorted(list(glob.iglob("ntsm/count_files/external/*.count")))

def find_command(wildcards):
    return command_dict[df.loc[wildcards.id, 'TYPE']]

def find_aln_params(wildcards):
    data_type = df.loc[wildcards.id, 'TYPE']
    return ALN_PARAMS+" "+f"-R '@RG\\tID:{wildcards.id}\\tSM:{wildcards.id}\\tPL:{data_type}'"

def get_bed(wildcards):
    return f"{SVDPREFIX}.bed"


def check_svd_files(wildcards):
    return [f"{SVDPREFIX}.{x}" for x in ["bed", "mu", "UD", "V"]]


wildcard_constraints:
    id = r"[^/.]+"

localrules:
    all,
    link_external_count,
    check_external_links,
    vbi_summary,

rule all:
    input:
        "summary/matched_result.tsv"

rule get_plot:
    input:
        "summary/all_pairwise.pdf"

rule get_all_count_files:
    input:
        expand("ntsm/count_files/{id}.count",
            id = df.index,
        ),

rule link_external_count:
    output:
        link_external_done = temp("ntsm/.external_link_done")
    threads: 1,
    run:
        external_counts = sorted(glob.glob(os.path.join(EXTERNAL_COUNTS_DIR, f"*.{COUNT_FILE_EXP}")))
        external_counts_dir = "ntsm/count_files/external"
        os.makedirs(external_counts_dir, exist_ok=True)
        for external_count in external_counts:
            link_name = os.path.basename(external_count).replace(f".{COUNT_FILE_EXP}","-EXT.count")
            link_dest_path = f"ntsm/count_files/external/{link_name}"
            if not os.path.lexists(link_dest_path):
                os.symlink(os.path.abspath(external_count), link_dest_path)
        open(output.link_external_done, "w").close()
        
rule make_count:
    input:
        fofn = find_fofn,
    output:
        count_file = "ntsm/count_files/{id}.count",
    params:
        ref_site = REF_SITE,
    threads: 8,
    resources:
        mem=lambda wildcards, attempt: 4 * attempt,
        hrs=12,
    singularity:
        "docker://eichlerlab/ntsm:1.2.1",
    shell: """
        ntsmCount -t {threads} -s {params.ref_site} $(cat {input.fofn}) > {output.count_file}
        """

rule check_external_links:
    input:
        rules.link_external_count.output.link_external_done
    output:
        check_link_done = temp(".check_external_links_done")
    threads: 1,
    shell: """ 
        touch {output.check_link_done}
        """

rule get_ntsm_quick_summary:
    input:
        all_counts = expand("ntsm/count_files/{id}.count",
                id = df.index,
            ),
        check_external_links = rules.check_external_links.output.check_link_done,
        external_all_counts = get_external_counts,
    output:
        summary = "summary/ntsm/ntsm_summary.tsv",
    threads: 4,
    params:
        site_center = SITE_CENTER,
        rotation_matrix = ROTATION_MATRIX,
    resources:
        mem=lambda wildcards, attempt: 4 * attempt,
        hrs=24,
    singularity:
        "docker://eichlerlab/ntsm:1.2.1",
    shell: """
        ntsmEval -t {threads} -n {params.site_center} -p {params.rotation_matrix} {input.all_counts} {input.external_all_counts} | sed -e 's/external\///g' -e 's/ntsm\/count_files\///g' -e 's/\.count//g' > {output.summary}
        """    


rule get_ntsm_all_pairwise_summary:
    input:
        all_counts = expand("ntsm/count_files/{id}.count",
                id = df.index,
            ),
        check_external_links = rules.check_external_links.output.check_link_done,
        external_all_counts = get_external_counts,
    output:
        summary = "summary/ntsm/all_pairwise.tsv",
    threads: 24,
    resources:
        mem=lambda wildcards, attempt: 2 * attempt,
        hrs=24,
    singularity:
        "docker://eichlerlab/ntsm:1.2.1",
    shell: """
        ntsmEval -a -t {threads} {input.all_counts} {input.external_all_counts} | sed -e 's/external\///g' -e 's/ntsm\/count_files\///g' -e 's/\.count//g' > {output.summary}
        """    

rule get_matched_summary:
    input:
        ntsm_quick_summary = "summary/ntsm/ntsm_summary.tsv",
        ntsm_summary = "summary/ntsm/all_pairwise.tsv",
        all_vbi_summary = "summary/vbi/vbi_summary.tsv"
    output:
        matched_summary = "summary/matched_result.tsv"
    threads: 1,
    params:
        compare_only_ext = COMPARE_ONLY_EXT
    resources:
        mem=16,
        hrs=4,
    run:
        compare_only_ext = params.compare_only_ext
        manifest_df = pd.read_csv(MANIFEST, sep="\t", header=0).set_index("ID", drop=True)
        ntsm_summary_df = pd.read_csv(input.ntsm_summary, sep="\t")
        vbi_summary_df = pd.read_csv(input.all_vbi_summary, sep="\t", usecols=["#SEQ_ID","AVG_DP","FREEMIX"])
        
        matched_data = []
        index_order = manifest_df.index.tolist()
        index_order
        for sample in index_order:
            vbi_subset = vbi_summary_df[vbi_summary_df["#SEQ_ID"] == sample]
            matched_samples = []
            matched_distance = []
            matched_relate = []
            subset_df = ntsm_summary_df[((ntsm_summary_df["sample1"] == sample) | (ntsm_summary_df["sample2"] == sample)) & (ntsm_summary_df["same"]>0)]
            if not subset_df.empty:
                for idx,row in subset_df.iterrows():
                    if row["sample1"] == sample:
                        matched_sample = row["sample2"]
                        # matched_samples.append(row["sample2"])
                    elif row["sample2"] == sample:
                        matched_sample = row["sample1"]
                        # matched_samples.append(row["sample1"])
                    if (compare_only_ext) and (not "-EXT" in matched_sample):
                        continue
                    matched_samples.append(matched_sample)
                    matched_distance.append(format(row["score"],".4f"))
                    matched_relate.append(format(row["relate"],".4f"))
            else: # matched sample not found
                subset_df = ntsm_summary_df[((ntsm_summary_df["sample1"] == sample) | (ntsm_summary_df["sample2"] == sample))] # subset without matching.
                print (sample, subset_df)
                closest_row = subset_df.loc[subset_df["relate"].idxmax()] # The highiest relate score result.
                if closest_row["sample1"] == sample:
                    closest_sample = closest_row["sample2"]
                else:
                    closest_sample = closest_row["sample1"]
                matched_samples = [f"NotFound. Closest:{closest_sample}"]
                matched_distance = [f"NA. Closest:{closest_row['score']}"]
                matched_relate = [f"NA. Closest:{closest_row['relate']}"]
            if vbi_subset.empty:
                vbi_dp = "NA"
                vbi_freemix = "NA"
                vbi_warning = "NA"
            else:
                vbi_dp = vbi_subset.iloc[0]["AVG_DP"]
                vbi_freemix = vbi_subset.iloc[0]["FREEMIX"]
                if vbi_freemix >= 0.01:
                    vbi_warning = "WARNING"
                else:
                    vbi_warning = ""

            matched_data.append([sample, ",".join(matched_samples), ",".join(matched_distance), ",".join(matched_relate), vbi_dp, vbi_freemix, vbi_warning])
        matched_df = pd.DataFrame(matched_data, columns = ["ID","MATCHED_SAMPLE","SCORE","RELATE","VBI_AVG_DP","VBI_FREEMIX","VBI_WARNING"])
        matched_df.to_csv(output.matched_summary, sep="\t", index=False)

rule plot_ntsm_summary:
    input:
        summary = rules.get_ntsm_all_pairwise_summary.output.summary,
    output:
        summary_plot = "summary/ntsm/all_pairwise.pdf",
    threads: 1,
    resources:
        mem=16,
        hrs=4,
    run:
        import numpy as np
        import seaborn as sns
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap, Normalize

        df = pd.read_csv(input.summary, sep="\t", header=0, usecols=["sample1","sample2","score"])
        
        samples = sorted(set(df["sample1"]).union(df["sample2"]))
        
        distance_matrix = pd.DataFrame(np.inf, index=samples, columns=samples)


        for idx, row in df.iterrows():
            distance_matrix.loc[row['sample1'], row['sample2']] = row['score']
            distance_matrix.loc[row['sample2'], row['sample1']] = row['score']

        np.fill_diagonal(distance_matrix.values, 0)

        data_max = df["score"].max()
        vmax = max(4.0, data_max)

        ## color code
        cool_blue = "#3b4cc0"
        middle_color = "#f7f7f7"
        warm_orange = "#f07c52"
        warm_red= "#b40426"
        
        colors = [
            (0.0, cool_blue),
            (0.5, middle_color),
            (1.0, warm_orange),
            (max(data_max, 4.0), warm_red)
        ]


        color_stops = [(min(v / vmax, 1.0), color) for v, color in sorted(colors, key=lambda x: x[0])]

        custom_cmap = LinearSegmentedColormap.from_list("my_map", color_stops)
        norm = Normalize(vmin=0.0, vmax=vmax)

        plt.figure(figsize=(30, 24))

        sns.heatmap(distance_matrix, annot=True, fmt=".2f", cmap=custom_cmap, annot_kws={"size":15}, cbar=True, norm=norm)

        plt.title(f"PCA Distance Heatmap", fontsize=30, pad=20)
        plt.xlabel("")
        plt.ylabel("")
        plt.xticks(rotation=90, fontsize=15)
        plt.yticks(rotation=0, fontsize=15)
        
        plt.savefig(output.summary_plot, format="pdf", bbox_inches="tight")
        plt.close()

rule map_reads:
    input:
        fofn = find_fofn,
    output:
        bam = "alignment/GRCh38/{id}.bam",
        bai = "alignment/GRCh38/{id}.bam.bai"
    resources:
        mem = 12,
        smem = 4,
        hrs = 96
    threads: 8
    params:
        ref = ALN_REF,
        command = find_command,
        aln_params = find_aln_params
    singularity:
        "docker://eichlerlab/align-basics:0.2",
    shell: """
        {params.command} {threads} {params.aln_params} {params.ref} $(cat {input.fofn}) | samtools view -b - | sambamba sort -t {threads} -o {output.bam} -m {resources.mem}G /dev/stdin
        samtools index {output.bam}
        """

rule run_pileup:
    input:
        ref = ALN_REF,
        bed = get_bed,
        bam = rules.map_reads.output.bam,
        bai = rules.map_reads.output.bai,
    output:
        pileup="vbi/pileup_files/{id}.pile",
    threads: 1
    params:
        ref = ALN_REF,
    resources:
        mem=16,
        hrs=12,
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    shell: """
        samtools mpileup -B -f {input.ref} -l {input.bed} {input.bam} -o {output.pileup}
        """

rule run_vbi:
    input:
        pileup_files=rules.run_pileup.output.pileup,
        svdprefix=check_svd_files,
    output:
        raw_selfsm=temp("vbi/results/{id}.raw.selfSM"),
    resources:
        mem=8,
        hrs=12,
    threads: 8
    params:
        ref = ALN_REF,
        svdprefix = SVDPREFIX,
    singularity:
        "docker://eichlerlab/vbi:2.0.1"
    shell: """
        AVG_DP=$(awk -F '\\t' '{{sum += $4}} END {{if (NR>0) print sum/NR; else print 0}}' {input.pileup_files})
        if awk -v a="$AVG_DP" 'BEGIN {{exit (a >= 3) ? 0 : 1}}'; then
            VerifyBamID \
            --DisableSanityCheck \
            --PileupFile {input.pileup_files} \
            --SVDPrefix {params.svdprefix} \
            --NumThread {threads} \
            --Reference {params.ref} \
            --Output $(echo {output.raw_selfsm} | sed 's/.selfSM//')
        else
            echo -e "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\nDefaultSampleName\tNA\tNA\tNA\tNA\t$AVG_DP\tNA" > {output.raw_selfsm}
        fi
        """

rule add_name_to_report:
    input:
        raw_selfsm = rules.run_vbi.output.raw_selfsm
    output:
        selfsm = "vbi/results/{id}.selfSM",
    threads: 1
    resources:
        mem=4,
        hrs=1,
    shell: """
        sed "s/DefaultSampleName/{wildcards.id}/" {input.raw_selfsm} > {output.selfsm}
        """

rule vbi_summary:
    input:
        sm_all=expand("vbi/results/{id}.selfSM", id=df.index),
    output:
        vbi_summary="summary/vbi/vbi_summary.tsv",
    threads: 1        
    shell: """
        for file in $( echo {input.sm_all} ); do cat $file; done | awk '!seek[$0]++' > {output.vbi_summary}
        """
