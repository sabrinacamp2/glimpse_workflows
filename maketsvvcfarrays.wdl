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
        mv ${sep = ' ' refpanel_sitestsv} /cromwell_root/
        mv ${sep = ' ' refpanel_sitestsv_index} /cromwell_root/
        mv ${sep = ' ' refpanel_sitesvcf} /cromwell_root/
        mv ${sep = ' ' refpanel_sitesvcf_index} /cromwell_root/

	>>>

	output {
		Array[File] refpanel_sitestsv_array = glob("*.tsv.gz")
		Array[File] refpanel_sitestsv_index_array = glob("*.tsv.gz.tbi")
		Array[File] refpanel_sitesvcf_array = glob("*.vcf.gz")
		Array[File] refpanel_sitesvcf_index_array = glob("*.vcf.gz.csi")
	}

	runtime {
		docker: "vanallenlab/samtools_bcftools_htslib:1.9"
		memory: "${memoryGb} GB"
		disks: "local-disk ${diskSpaceGb} HDD"
		preemptible: "${preemptible}"
	}
}