workflow referencepanel_workflow {
	call referencepanel {
	}
}

task referencepanel {
	File thousg_vcf
	File thousg_vcf_idx
	File? exclude_samples
	String chr
	Int diskSpaceGb
	Int memoryGb
	Int preemptible


	command <<<

		## keep only SNPs and remove multiallelic records
		bcftools view -m 2 -M 2 -v snps ${"-s " + exclude_samples} --threads 4 ${thousg_vcf} -Ob -o thousGP.chr${chr}.subset.bcf

		## index the vcf output
		bcftools index -f thousGP.chr${chr}.subset.bcf


	>>>

	output {
		File refpanel_curated = "$thousGP.chr${chr}.subset.bcf"
		File refpanel_curated_index = "thousGP.chr${chr}.subset.bcf.csi"
	}

	runtime {
		docker: "quay.io/biocontainers/bcftools:1.9--ha228f0b_3"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}