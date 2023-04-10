version 1.0

workflow filtervcf {
	input {
		String outputFileNamePrefix
		String tumorSampleName
		File tumorvcf
		File tumorvcfindex
		String controlFileList = "/.mounts/labs/gsi/src/pwgs_hbc/1.0/HBC.bam.list"
	}

	parameter_meta {
		tumorvcf: "tumor vcf file, bgzip"
		tumorvcfindex: "tumor vcf index file"
		outputFileNamePrefix: "Prefix for output file"
		tumorSampleName: "ID for WGS tumor sample, must match .vcf header"
		controlFileList: "tab seperated list of bam and bai files for healthy blood controls"
	}

	call filterVCF {
		input:
		tumorvcf = tumorvcf,
		tumorvcfindex = tumorvcfindex,
		tumorSampleName = tumorSampleName
	}

	meta {
		author: "Felix Beaudry"
		email: "fbeaudry@oicr.on.ca"
		description: "Workflow for pre-filtering VCF for MRDetect or fragmentomics analysis"
		dependencies:
		[
			{
				name: "bcftools/1.9",
				url: "https://github.com/samtools/bcftools"
			}
		]
		output_meta: {
			snvDetectionResult: "Result from SNV detection incl sample HBCs",
			final_call: "Final file of mrdetect results",
			snvDetectionVAF: "VAF from SNV detection for sample",
			snpcount: "number of SNPs in vcf after filtering"
		}
	}
	output {
		File snpcount = filterVCF.snpcount
		File filteredvcf = filterVCF.filteredvcf
	}
}

task filterVCF {
	input {
		File tumorvcf
		File tumorvcfindex
		String tumorSampleName
		String tumorVCFfilter = "FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"
		String tumorVAF = "0.1"
		String genome = "$HG38_ROOT/hg38_random.fa"
		String difficultRegions = "--regions-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"
		String modules = "bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		tumorvcf: "tumor vcf file, bgzip"
		tumorvcfindex: "tumor vcf index file"
		tumorSampleName: "ID for WGS tumor sample"
		tumorVCFfilter: "set of filter calls to incl. in tumor VCF (any line with these flags will be included"
		tumorVAF: "Variant Allele Frequency for tumor VCF"
		genome: "Path to loaded genome .fa"
		difficultRegions: "Path to .bed excluding difficult regions, string must include the flag --regions-file "
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"

	}

	command <<<
		set -euo pipefail

		$BCFTOOLS_ROOT/bin/bcftools view -s ~{tumorSampleName} ~{difficultRegions} ~{tumorvcf} |\
		$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome} |\
		$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='snps'" |\
		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{tumorVCFfilter}" |\
		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{tumorVAF}" > ~{tumorSampleName}.SNP.vcf

		awk '$1 !~ "#" {print}' ~{tumorSampleName}.SNP.vcf | wc -l >SNP.count.txt

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File filteredvcf = "~{tumorSampleName}.SNP.vcf"
		File snpcount = "SNP.count.txt"
	}

	meta {
		output_meta: {
			filteredvcf: "Filtered vcf",
			snpcount: "Number of SNPs after filtering"
		}
	}
}
