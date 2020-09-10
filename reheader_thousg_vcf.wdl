workflow reheader_workflow {
	call reheader {
	}
}

task reheader {
	File thousg_vcf
	File thousg_vcf_idx
	File GTline
	String name = basename(thousg_vcf, "_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
	Int diskSpaceGb
	Int memoryGb
	Int preemptible


	command <<<

		## add line to specify genotype field 
		bcftools annotate --header-lines ${GTline} ${thousg_vcf} -Oz -o ${name}.reheader.genotypes.vcf.gz
		tabix ${name}.reheader.genotypes.vcf.gz


	>>>

	output {
		File referencepanel_reheader = "${name}.reheader.genotypes.vcf.gz"
		File referencepanel_reheader_idx = "${name}.reheader.genotypes.vcf.gz.tbi"
	}

	runtime {
		docker: "vanallenlab/samtools_bcftools_htslib:1.9"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}