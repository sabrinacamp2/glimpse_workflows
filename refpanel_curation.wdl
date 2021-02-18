workflow referencepanel_workflow {
	call referencepanel {
	}
}

task referencepanel {
	File thousg_vcf
	File thousg_vcf_idx
	File? keep_samples
	String name = basename(thousg_vcf, ".genotypes.vcf.gz")
	Int diskSpaceGb
	Int memoryGb
	Int cpu = 4
	Int preemptible


	command <<<

		## keep only SNPs and remove multiallelic records
		bcftools view -m 2 -M 2 -v snps ${"-S " + keep_samples} --threads 4 ${thousg_vcf} -Ob -o thousGP.${name}.subset.bcf

		## index the vcf output
		bcftools index -f thousGP.${name}.subset.bcf


	>>>

	output {
		File refpanel_curated = "thousGP.${name}.subset.bcf"
		File refpanel_curated_index = "thousGP.${name}.subset.bcf.csi"
	}

	runtime {
		docker: "quay.io/biocontainers/bcftools:1.9--ha228f0b_3"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		cpu: "${cpu}"
		preemptible: "${preemptible}"
	}
}