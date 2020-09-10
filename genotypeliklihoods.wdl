workflow gtlikelihoods_workflow {
	Array[File] refpanel_sitestsv = ${refpanel_sitestsv}
	Array[File] refpanel_sitesvcf = ${refpanel_sitesvcf}
	Array[String] chr = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX']

	scatter(pair in zip(zip(refpanel_sitestsv, refpanel_sitesvcf), chr)){
		call gtlikelihoods {
			input: sitesvcf = pair.left.right,
					sitestsv = pair.left.left,
					chr_names = pair.right
	}
	}
	
}

task gtlikelihoods {
	
	File samplebam
	File sitesvcf
	File sitestsv
	String chr_names
	File reffasta
	File refdict
	File refidx
	String sampleid
	Int diskSpaceGb
	Int memoryGb
	Int preemptible


	command <<<

		## compute genotype likelihoods
		bcftools mpileup -f ${reffasta} -I -E -a 'FORMAT/DP' -T ${sitesvcf} -r ${chr_names} ${samplebam} -Ou | bcftools call -Aim -C alleles -T ${sitestsv} -Oz -o ${sampleid}.${chr_name}.vcf.gz
		bcftools index -f ${sampleid}.${chr_name}.vcf.gz

	>>>

	output {
		File genotypelikelihood = "${sampleid}.${chr_name}.vcf.gz"
		File genotypelikelihood_index = "thousGP.${name}.sites.vcf.gz.csi"
	}

	runtime {
		docker: "vanallenlab/samtools_bcftools_htslib:1.9"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}