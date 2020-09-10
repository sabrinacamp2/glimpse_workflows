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
		bcftools view -G -m 2 -M 2 -v snps ${refpanel_curated} -Oz -o ${name}.sites.vcf.gz
		bcftools index -f ${name}.sites.vcf.gz

		## get the snp only vcf in a tsv format 
		bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${name}.sites.vcf.gz | bgzip -c > ${name}.sites.tsv.gz
		tabix -s1 -b2 -e2 ${name}.sites.tsv.gz
        
        ls -lh


	>>>

	output {
		File refpanel_sitesvcf = "${name}.sites.vcf.gz"
		File refpanel_sitesvcf_index = "${name}.sites.vcf.gz.csi"
		File refpanel_sitestsv = "${name}.sites.tsv.gz"
		File refpanel_sitestsv_index = "${name}.sites.tsv.gz.tbi"
	}

	runtime {
		docker: "vanallenlab/samtools_bcftools_htslib:1.9"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}