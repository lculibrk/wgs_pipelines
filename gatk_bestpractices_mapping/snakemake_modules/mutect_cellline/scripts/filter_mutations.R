#!/usr/bin/env Rscript
library(getopt)
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help'    , 'h', 0, "logical",
	'daughter', 'd', 1, "character",
	'parent'  , 'p', 1, "character",
	'threads' , 't', 2, "integer"
), byrow=TRUE, ncol=4)


args = getopt(spec)
if(!all(c("daughter", "parent") %in% names(args))){
	print(getopt(spec, usage = T))
	stop()
}

if(!'threads' %in% names(args)){
	threads = 1
}else{
	threads = args$threads
}

library(data.table)
message(paste0("Using ", threads, " threads for data.table"))
setDTthreads(threads)
library(stringr)

parent_table = fread(args$parent, header = F)
daughters_table = fread(args$daughter, header = F)
names(parent_table) = c("path", "name", "line")
names(daughters_table) = c("path", "germ", "name", "parent", "line")



## function to fix GT columns
fix_gt = function(ref, gt_string){
	if(grepl("^[A-Z](/|\\|)[A-Z]$", gt_string)){
		return(gsub("[A-Z]", "1", gsub(ref, "0", gt_string)))
	}
	slash = grepl("/", gt_string)
	if(slash){
		s = "/"
		os = s
	}else{
		s = "\\|"
		os = "|"
	}
	gts = unlist(strsplit(gt_string, split = s))
	gts = gts[nchar(gts) == 1][1:2]
	fix_gt(ref, paste0(gts, collapse = os))
}
## Split analysis by cell line, so that we limit resource usage

lines = unique(daughters_table$line)
big_qc_list = list()
for(current_line in lines){
	message(paste0("Processing line ", current_line))
	## Subset tables
	daughter_subtable = daughters_table[line == current_line]
	parent_subtable = parent_table[line == current_line]
	## Process parents for this line
	## We want to get a table of one line per mutation, and the n of parents with the mutation
	message("Reading parent variants...")
	parents = lapply(1:nrow(parent_subtable), function(x){
		message(paste0("Reading parent ", parent_subtable[x]$name[1]))
		tbl = fread(parent_subtable[x]$path)
		message("Table read successfully")
		tbl$name = parent_subtable[x]$name
		samplenamestring = gsub("(.*\\.)F2R1", "\\1", names(tbl)[grepl("\\.F2R1", names(tbl))])
		names(tbl) = gsub(samplenamestring, "", names(tbl))
		tbl = tbl[!grepl(",0", AD)]
		return(tbl)
	})
	parents = rbindlist(parents)
	parents = parents[nchar(REF) == 1 & nchar(ALT) == 1]
	parents[,n:=.N, by=c("CHROM", "POS", "REF", "ALT")]
	## Need to fix multiallelic sites. - remember to use DP instead of a+r
	parents[,als:=str_count(AD, ",")]
	parents[als>1,id:=1:.N]
	parents[als > 1,AD:=paste0(unlist(strsplit(AD, ","))[c(1,which(ALT == unlist(strsplit(GT, "(/|\\|)"))))], collapse = ","), by = id]
	parents[,c("r", "a"):=lapply(tstrsplit(AD, ","), as.numeric)]
	qc_list = parents[n == 1 & a/(r+a) > 0.25 & DP >= 15 & FILTER == "PASS"]
	big_qc_list = c(big_qc_list, list(qc_list))
	## Now load daughters and start filtering
	message("Reading daughter variants...")
	daughters = lapply(1:nrow(daughter_subtable), function(x){
		tbl_std = fread(daughter_subtable[x]$path)[nchar(REF) == 1 & nchar(ALT) == 1]
		tbl_comp = fread(daughter_subtable[x]$germ)[nchar(REF) == 1 & nchar(ALT) == 1]
		tbl_std$type = "std"
		tbl_comp$type = "comp"
		tbl = rbind(tbl_std, tbl_comp)
		tbl[,n:=.N, by= c("CHROM", "POS", "REF", "ALT")]
		tbl = tbl[n == 1 | type == "std"]
		tbl$n = NULL
		tbl$type = NULL
		tbl$name = daughter_subtable[x]$name
		tbl$parent = daughter_subtable[x]$parent
		rm(tbl_std)
		rm(tbl_comp)
		return(tbl)
	})
	daughters = rbindlist(daughters)
	daughters[,n_lineage:=.N, by = c("CHROM", "POS", "REF", "ALT", "parent")]
	daughters[,n:=.N, by = c("CHROM", "POS", "REF", "ALT")]
	
	message("Filtering based on parental and daughter sequence content...")
	m2_pass = daughters[FILTER == "PASS"]
	setkey(m2_pass, CHROM, POS, REF, ALT)
	parent_filtset = unique(parents[n > nrow(parent_subtable)/2][,c("CHROM", "POS", "REF", "ALT")])
	setkey(parent_filtset, CHROM, POS, REF, ALT)
	m2_pass[,iden:=1:.N]
	p_remove = m2_pass[parent_filtset][!is.na(name)]$iden
	m2_pass$iden = NULL
	m2_pass[n == 1,tumor.LAB:="unq"]
	m2_pass[n_lineage > 1,tumor.LAB:="shared_lineage"]
	m2_pass[n > 1 & n_lineage != n,tumor.LAB:="shared_external"]
	m2_pass[p_remove, tumor.LAB:="shared_parental"]
	m2_pass[normal.DP < 15, tumor.LAB:="depth"]
	m2_pass[tumor.LAB %in% c("shared_external", "shared_parental", "depth"),FILTER:=tumor.LAB]
	daughters[FILTER != "PASS",tumor.LAB:="mutect_filtered"]
	daughters = rbind(daughters[FILTER != "PASS"], m2_pass)

	## Filtering done
	## Now we perform the cell-lineage-of-origin QC for each line
	message("Performing parent-of-origin QC...")
	setkey(daughters, CHROM, POS, REF, ALT)
	setkey(qc_list, CHROM, POS, REF, ALT)
	overlaps = daughters[qc_list[,c("CHROM", "POS", "REF", "ALT", "name")]]
	combins = data.table(expand.grid(sort(unique(daughter_subtable$name)), sort(unique(qc_list$name))))
	for(i in 1:nrow(combins)){
		daughter = combins[i]$Var1
		test_parent = combins[i]$Var2
		ols = nrow(overlaps[name == daughter & i.name == test_parent])/nrow(qc_list[name == test_parent])
		combins[i,sim:=ols]
	}
	combins = combins[order(Var1, Var2)]
	
	info_columns = c("AS_FilterStatus", "AS_SB_TABLE", "AS_UNIQ_ALT_READ_COUNT", "CONTQ", "DP", "ECNT", "GERMQ", "MBQ", "MFRL", "MMQ", "MPOS", "NALOD", "NCount", "NLOD", "OCM", "PON", "POPAF", "ROQ", "RPA", "RU", "SEQQ", "STR", "STRANDQ", "STRQ", "TLOD")
	
	format_phase   = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "PGT", "PID", "PS", "SB", "LAB")
	format_nophase = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD",                     "SB", "LAB")

	message("Writing QC file and final filtered VCFs...")
	for(i in 1:nrow(daughter_subtable)){
		message(paste0("Converting ", daughter_subtable[i]$name, " table to VCF and writing..."))
		## First get the path we'll dump files to
		out_path = gsub("(.*/)[^/]*", "\\1", gsub("std", "proc", daughter_subtable[i]$path))
		## QC results
		out_qc = combins[Var1 == daughter_subtable[i]$name]
		if(!exists(out_path)){
			dir.create(out_path)
		}
		
		fwrite(out_qc, paste0(out_path, "parent_of_origin.txt"), sep = "\t")
		
		## Filtered variants, convert to VCF without header
		out_daughters = daughters[name == daughter_subtable[i]$name]
		out_daughters[,c("n", "name", "parent", "n_lineage"):=NULL]
		
		for(col in info_columns){
			out_daughters[,(col):=paste0(col, "=", out_daughters[[col]])]	
		}
		out_daughters[,(info_columns):=lapply(.SD, function(x)gsub(".*=NA", "", x)), .SDcols = info_columns]
		out_daughters[,INFO:=do.call(paste, c(.SD, sep = ";")), .SDcols = info_columns]
		out_daughters[,INFO:=gsub(";;+", ";", INFO)]
		out_daughters[,(info_columns):=NULL]
		out_daughters[,normal.LAB:="norm"]
		## Now handle the format column
		out_daughters[ is.na(normal.PGT),FORMAT:=paste0(format_nophase, collapse=":")]
		out_daughters[!is.na(normal.PGT),FORMAT:=paste0(format_phase,   collapse=":")]
	
		
		out_daughters[,normal.GT:=fix_gt(REF, normal.GT), by = 1:nrow(out_daughters)]
		out_daughters[,tumor.GT:=fix_gt(REF, tumor.GT), by = 1:nrow(out_daughters)]

		
		out_daughters[is.na(normal.PGT),normal:=do.call(paste, c(.SD, sep = ":")), .SDcols = paste0("normal.", format_nophase)]
		out_daughters[!is.na(normal.PGT),normal:=do.call(paste, c(.SD, sep = ":")), .SDcols = paste0("normal.", format_phase)]
		
		out_daughters[is.na(normal.PGT),tumor:=do.call(paste, c(.SD, sep = ":")), .SDcols = paste0("tumor.", format_nophase)]
		out_daughters[!is.na(normal.PGT),tumor:=do.call(paste, c(.SD, sep = ":")), .SDcols = paste0("tumor.", format_phase)]
		
		out_daughters[,FILTER:=gsub(",", ";", FILTER)]
		out_vcf = out_daughters[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "tumor", "normal")]
		names(out_vcf)[1] = "#CHROM"
		
		fwrite(out_vcf, paste0(out_path, "variants.vcf"), sep = "\t")
		
		}
}

