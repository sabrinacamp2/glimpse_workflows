workflow extractsites_workflow {
	call extractsites {
	}
}

task extractsites {
	File refpanel_curated
	File refpanel_curated_index
	String name = basename(refpanel_curated, ".subset.bcf")
	Int diskSpaceGb
	Int memoryGb
	Int preemptible


	command <<<

		## get a snp only vcf
		bcftools view -G -m 2 -M 2 -v snps ${refpanel_curated} -Oz -o thousGP.${name}.sites.vcf.gz
		bcftools index -f thousGP.${name}.sites.vcf.gz

		## get the snp only vcf in a tsv format 
		bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' thousGP.${name}.sites.vcf.gz | bgzip -c > thousGP.${name}.sites.tsv.gz
		tabix -s1 -b2 -e2 thousGP.${name}.sites.tsv.gz


	>>>

	output {
		File refpanel_sitesvcf = "thousGP.${name}.sites.vcf.gz"
		File refpanel_sitesvcf_index = "thousGP.${name}.sites.vcf.gz.tbi"
		File refpanel_sitestsv = "thousGP.${name}.sites.tsv.gz"
		File refpanel_sitestsv_index = "thousGP.${name}.sites.tsv.gz.tbi"
	}

	runtime {
		docker: "quay.io/biocontainers/bcftools:1.9--ha228f0b_3"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}