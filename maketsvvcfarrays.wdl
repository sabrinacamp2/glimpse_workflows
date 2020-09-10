workflow makearrays_workflow {
	call makearrays {
	}
}

task makearrays {
	Array[File] refpanel_sitestsv
	Array[File] refpanel_sitestsv_index
	Array[File] refpanel_sitesvcf
	Array[File] refpanel_sitesvcf_index
	Int diskSpaceGb
	Int memoryGb
	Int preemptible


	command <<<

		echo "please work"


	>>>

	output {
		Array[File] refpanel_sitestsv = glob("*.tsv.gz")
		Array[File] refpanel_sitestsv_index = glob("*.tsv.gz.tbi")
		Array[File] refpanel_sitesvcf = glob(".vcf.gz")
		Array[File] refpanel_sitesvcf_index = glob(".vcf.gz.tbi")
	}

	runtime {
		docker: "vanallenlab/samtools_bcftools_htslib:1.9"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}