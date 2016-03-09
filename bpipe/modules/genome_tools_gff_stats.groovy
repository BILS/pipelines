
genome_tools_gff_stats = {

  doc title: "Generate statistics from a gff3 file using the genometools package",

  desc: "Genometools compiles feature counts from a gff3-formatted annotation file",

  author: "marc.hoeppner@bils.se"

  var sample_dir : false
  if (branch.sample_dir) { sample_dir = true }
  
  requires GENOME_TOOLS : "Must provide path to genometools (GENOME_TOOLS)"

  if (sample_dir) {
    output.dir = branch.outdir
  }
    
  transform(".gff") to (".gtcounts") {
    exec "gt stat $input > $output 2>/dev/null"
  }
  
}
