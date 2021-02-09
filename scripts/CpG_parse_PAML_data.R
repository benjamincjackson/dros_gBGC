setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source('../scripts/PAML_parse_output_file.R')

CpG.pe <- readOutput("../CpG/PAML/CpG/pointestimate/run_1/output")
nonCpG.pe <- readOutput("../CpG/PAML/nonCpG/pointestimate/run_1/output")
CpGprone.pe <- readOutput("../CpG/PAML/CpGprone/pointestimate/run_1/output")
nonCpGprone.pe <- readOutput("../CpG/PAML/nonCpGprone/pointestimate/run_1/output")

CpG.pe$sim_blength + CpG.pe$mel_blength
nonCpG.pe$sim_blength + nonCpG.pe$mel_blength
CpGprone.pe$sim_blength + CpGprone.pe$mel_blength
nonCpGprone.pe$sim_blength + nonCpGprone.pe$mel_blength


((CpG.pe$sim_blength + CpG.pe$mel_blength) * CpG.pe$length + (nonCpG.pe$sim_blength + nonCpG.pe$mel_blength) * nonCpG.pe$length) / (CpG.pe$length + nonCpG.pe$length)
