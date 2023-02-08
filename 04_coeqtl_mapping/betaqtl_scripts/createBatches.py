import gzip
import sys
import os
import glob

if len(sys.argv) < 6:
	print("Usage: createbatches.py expfile.txt.gz gte.txt genotype.vcf.gz genelist.txt.gz annotation.txt.gz template.sh nrmaxgenesperbatch outdir")
	sys.exit(0)

expfile = sys.argv[1]
gte = sys.argv[2]
genotype = sys.argv[3]
genelist = sys.argv[4]
annotation = sys.argv[5]
template = sys.argv[6]
nrgenes = int(sys.argv[7])
out = sys.argv[8]
condition = sys.argv[9]
celltype = sys.argv[10]

if not out.endswith("/"):
	out = out + "/"

def writeJob(exp, gte, gt, template, batchfile, jobfile, outprefix, logprefix, chr, condition, celltype):
	print("Writing job: "+jobfile)
	fh = open(template,'r')
	lines = fh.readlines()
	fh.close()
	fho = open(jobfile,'w')
	for line in lines:
		line = line.replace("GENOTYPE",gt)
		line = line.replace("GTE",gte)
		line = line.replace("EXPRESSION",exp)
		line = line.replace("CHROM",str(chr))
		line = line.replace("BATCHFILE",batchfile)
		line = line.replace("OUTPREFIX",outprefix)
		line = line.replace("LOGPREFIX",logprefix)
		line = line.replace("CONDITION", condition)
		line = line.replace("CELLTYPE", celltype)
		fho.write(line)
	fho.close()

def checkDir(path):
	if os.path.exists(path):
		# delete contents
		files = glob.glob(path+"*")
		for file in files:
			print("Removing: "+file)
			os.remove(file)
	else:
		print("Creating dir: "+path)
		os.mkdir(path)

abspath = os.path.abspath(out)
checkDir(abspath+"/batches/")
checkDir(abspath+"/output/")
checkDir(abspath+"/jobs/")
checkDir(abspath+"/logs/")

# read expression file
fh = None
genesinfile=genelist
print("Reading: "+genesinfile)
if genesinfile.endswith(".txt.gz"):
	fh = gzip.open(genesinfile,'rt')
else:
	fh = open(genesinfile,'r')
genesInExp = set()
fh.readline()
for line in fh:
#	gene = line.split("\t", maxsplit=1)[0]
	gene = line.strip()
	genesInExp.add(gene)
#	print(gene)

fh.close()
print("{} genes in {}".format(len(genesInExp),expfile))

# read gene set
geneset = set()
fh = None
print("Genelist: "+genelist)
if genelist.endswith(".txt.gz"):
	fh = gzip.open(genelist,'rt')
else:
	fh = open(genelist,'r')
for line in fh:
	gene  = line.strip()
	if gene in genesInExp:
		geneset.add(line.strip())
fh.close()
print("Genes in genelist: {}".format(len(geneset)))

# read annotation
print("Annotation: "+annotation)
fh = None
if annotation.endswith(".txt.gz"):
	fh = gzip.open(annotation,'rt')
else:
	fh = open(annotation,'r')
fh.readline()
genesPerChr = {}
annotread = 0
for line in fh:
	elems = line.strip().split("\t")
	gene = elems[1]
	if gene in geneset:
		chr = -1
		try:
			chr = int(elems[3])
		except:
			print(gene+" has non-numeric chromosome: "+elems[3])
		if chr < 23 and chr > 0:
			pos = int(elems[4])
			chrgenes = genesPerChr.get(chr)
			if chrgenes is None:
				chrgenes = []
			chrgenes.append(gene)
			genesPerChr[chr] = chrgenes
			annotread = annotread + 1
fh.close()
print("Annotation read for {} genes".format(annotread))

# create batches
for chr in genesPerChr.keys():
	bctr = 1
	chrgenes = genesPerChr.get(chr)
	gctr = 0
	bgctr = 0
	batchname = "chr"+str(chr)+"-batch-"+str(bctr)
	# write job script for first batch
	batchfile = abspath+"/batches/"+batchname+".txt"
	print("Writing batch: "+batchfile)
	jobfile = abspath+"/jobs/"+batchname+".sh"
	outprefix = abspath+"/output/"+batchname
	logprefix = abspath+"/logs/"+batchname
	print()
	print("Writing job: "+template+"\n"+batchfile+"\n"+jobfile+"\n"+outprefix+"\n"+str(chr))
	# exp, gte, gt, template, batchfile, jobfile, outprefix, chr
	chrgenotype = genotype.replace("CHR",str(chr))
	writeJob(expfile, gte, chrgenotype, template, batchfile, jobfile, outprefix, logprefix, chr, condition, celltype)
	bgout = open(batchfile,'w')
	while gctr < len(chrgenes):
		bgout.write(chrgenes[gctr]+"\n")
		bgctr = bgctr + 1
		if bgctr == nrgenes:
			# start new batch
			bgout.close()
			bctr = bctr + 1
			# write job script for new batch
			batchname = "chr"+str(chr)+"-batch-"+str(bctr)
			batchfile = abspath+"/batches/"+batchname+".txt"
			print("Writing batch: "+batchfile)
			jobfile = abspath+"/jobs/"+batchname+".sh"
			outprefix = abspath+"/output/"+batchname
			logprefix = abspath+"/logs/"+batchname
#			writeJob(template, batchfile, jobfile, outprefix, chr)
			writeJob(expfile, gte, chrgenotype, template, batchfile, jobfile, outprefix, logprefix, chr, condition, celltype)

			bgout = open(batchfile,'w')
			bgctr = 0
		gctr = gctr + 1
	# if there are any genes left, close batch
	if bgctr > 0:
		bgout.close()

