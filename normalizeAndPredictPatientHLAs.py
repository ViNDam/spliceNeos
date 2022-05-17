#!/user/bin/python
############################################################################
# This codes process and normalize gene expression data, splicing data, and mutation data in parallel.. It also predict patients' HLA types.
# Usage:
# python quantifyAndNormalizeData.py -e sample_gene_expression.tsv -m TCGA.PAAD.maf -f fastqFolder -source 1 -norm True


from __future__ import print_function
from collections import defaultdict
import pickle
from config import *
import matplotlib.pyplot as plt
from Bio import SeqIO
import os
import stat
import pandas as pd
import argparse
import pymp
import subprocess
import glob
import numpy as np
import extract9mersMissenseMutation
from k9merCodes import get9merAbundanceTCGA
import shutil


############################################################################

def processGeneExpression(geArg, sourceArg, normArg):
    """ Process or/and normalize gene expression

    Keyword arguments:
    geArg -- gene expression data
    sourceArg -- 1 if data is tcga , 0: others
    normArg -- True if user want to normalize gene expression data
    """
    ge = pd.read_csv(geArg, delimiter="\t",
                     header=0, index_col=0)  # gene expression data
    gc = pd.read_csv(gencode, delimiter="\t", header=0)  # gencode data
    # if normalize
    if normArg == "True":
        tpms = normalizeGeneExpressionData(ge, gc, int(sourceArg))
    else:
        tpms = ge
    tpms.index = [ENST.split(".")[0] for ENST in tpms.index]
    tpms.index.names = ['ENSG']
    tpms.to_csv(out + "/Norm_tpms.tsv",
                header=True, index=True, sep="\t")


def processSplicing(fastqFoldArg, sourceArg):
    """ Quantify and normalize splicing expression

    Keyword arguments:
    fastqFoldArg -- directory of folder contains fastqs
    sourceArg -- 1 if data is tcga , 0: others
    """
    print("processSplicing")
    spl_ENST_Expr_dict = quantifySplFastq(fastqFoldArg, sourceArg)
    HLAA_dict = extractHLAA_ENSTs_toDict(spl_ENST_Expr_dict)
    # -------- convert ENST_splicing data to 9MERS_splicing data --------
    os.chdir(codes + "/k9merCodes")
    get9merAbundanceTCGA.main(
        kmerDict, out + "/splicing/allSplicingTPM.csv", out + "/splicing/allSplicingTPM_9mers.tsv")

    # quantile normalize splicing data
    spl_qr = quantileNormalizeSpl(sourceArg)
    spl_qr.to_csv(out + "/splicing/allSplicingTPM_9mers_qr.tsv", sep="\t")
    splQrDict = toDictAndSave(spl_qr, "spl_qr_dict.pkl")
    print("Done processing splicing.")


def toDictAndSave(df, name):
    if not os.path.exists(wdir + "/bin/"):
        os.makedirs(wdir + "/bin/")
    df_dict = df.to_dict(orient='index')  # 11,167,266 peptides
    f = open(wdir + "/bin/" + name, "wb")
    pickle.dump(df_dict, f)
    f.close()
    return(df_dict)


def extractHLAA_ENSTs_toDict(splQrDict):
    print("extractHLAA_ENSTs_toDict")
    with open(ENSG_Gene_dict, "rb") as f:
        ENSGGeneDict = pickle.load(f)
    with open(ENST_ENSG_dict, "rb") as f:
        ENSTENSGDict = pickle.load(f)
    # with open(spl_qr_dict, "rb") as f:
    #    splQrDict = pickle.load(f)
    HLAA_ENSGs = [ENST for ENST, symb in ENSGGeneDict.items()
                  if symb == "HLA-A"]
    HLAA_ENSTs = [ENST.split(".")[0] for ENST,
                  ENSG in ENSTENSGDict.items() if ENSG in HLAA_ENSGs]
    HLAA_dict = defaultdict(dict)
    for ENST in HLAA_ENSTs:
        id_expr = {ID: splQrDict[ENST][ID]
                   for ID in splQrDict[ENST].keys()}
        for ID in id_expr.keys():
            HLAA_dict.setdefault(ID, [])
            HLAA_dict[ID].append(id_expr[ID])
    f = open(root + "/bin/HLAA_dict.pkl", "wb")
    pickle.dump(HLAA_dict, f)
    f.close()
    return(HLAA_dict)


def normalizeGeneExpressionData(ge, gc, sourceArg):
    """ Normalize gene expression data. If data is TCGA type, filter TCGA patient IDs before normalize.

    Keyword arguments:
    ge -- gene expression dataframe (columns: patient IDs, rows : gene ID ENSG)
    gc -- gene code info dataframe (columns: genes' info, rows : gene ID ENSG)
    sourceArg -- 1 = tcga type; 0 = other type
    """
    print("normalizeGeneExpressionData")
    if (int(float(sourceArg)) == 1):
        ge = filterCancerIDsOnly(ge)
        tpms = normalize(ge, gc)
    else:
        tpms = normalize(ge, gc)
    print("Done processing gene expression with normalization.")
    return(tpms)


def filterCancerIDsOnly(ge):
    """ Filter TCGA patient IDs columns that are cancer only

    Keyword arguments:
    ge -- gene expression dataframe (columns: patient IDs, rows : gene ID ENSG)
    """
    ind = [ID[13:15] for ID in ge.columns]  # remove normal samples in ge
    ind = [i for i, v in enumerate(ind) if int(float(v)) < 11]
    ge = ge.iloc[:, ind]
    ge.columns = [ID[0:12] for ID in ge.columns]
    return(ge)


def normalize(ge, gc):
    """ Transcript Per Million (TPM) Normalization

    Keyword arguments:
    ge -- gene expression dataframe (columns: patient IDs, rows : gene ID ENSG)
    gc -- gene code info dataframe (columns: genes' info, rows : gene ID ENSG)
    """
    gc.index = gc.ENSG
    gc = gc.reindex(index=ge.index)
    tmp = ge.divide(gc['Gene_length'], axis=0)
    sums = tmp.sum(axis=0)
    tpm = tmp / sums * 1e6
    return(tpm)


def quantifySplFastq(fastqFolder, sourceArg):
    """quantify fastq files and convert their transcript's expression to 9_mers expression

    Keyword arguments:
    fastqFolder -- folder contains all fastq files ending with .fastq
    """
    print("quantifySplFastq")
    outdir = out + "/salmon_quant103/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(out+"/splicing/"):
        os.makedirs(out+"/splicing/")
    salmonQuantification(fastqFolder, outdir)
    # -------- combine all quantified splicing files --------
    os.chdir(outdir)
    dirs = [f for f in os.listdir(outdir) if f.endswith(".quant")]
    splDF = combineQuantFiles(dirs, outdir, sourceArg)
    # -------- save as dictionary for later use
    spl_ENST_Expr_dict = splDF.to_dict(orient='index')  # 180,869 ENSTs
    f = open(root + "/bin/spl_ENST_Expr_dict.pkl", "wb")
    pickle.dump(spl_ENST_Expr_dict, f)
    f.close()
    return(spl_ENST_Expr_dict)


def salmonQuantification(fastqFolder, outdir):
    """ Parallel quantify splicing expression from fastq using salmon

    Keyword arguments:
    fastqFolder -- folder contains fastq files
    files -- list of files
    """
    print("salmonQuantification")
    files = glob.glob(fastqFolder + "/*_1.fastq")
    if len(files) == 0:
        files = glob.glob(fastqFolder + "/*.fastq")  # single-end reads
        print(files)
        with pymp.Parallel(ncore) as p1:
            for index in p1.range(0, len(files)):
                #for index in range(0, len(files)):
                outf = outdir + os.path.basename(files[index]) + ".quant"
                # print('salmon quant -i ' + data + '/salmon_index103 --libType A -p 12 --validateMappings -o '
                #                         + outf + ' -r ' + files[index])
                # salmon quant -i /data/vdam/proj/immgenGithub/data/salmon_index103 --libType A -p 12 --validateMappings -o /data/vdam/proj/immgenGithub/outfiles/salmon_quant103/TCGA-2J-AAB9-01A.fastq.bam -r /data/vdam/proj/immgenGithub/inputs/fastqFolder/TCGA-2J-AAB9-01A.fastq
                print('docker run --rm -v ' + wdir + ':' + wdir + ' -t combinelab/salmon salmon quant -i '
                      + root + '/data/salmon_index103/ --libType A -p 12 --validateMappings -r ' + files[index] + ' -o ' + outf)
                subprocess.call(['docker run --rm -v ' + wdir + ':' + wdir + ' -t combinelab/salmon salmon quant -i '
                                 + root + '/data/salmon_index103/ --libType A -p 12 --validateMappings -r ' + files[index] + ' -o ' + outf], shell=True)
                subprocess.call(['docker run --rm -v ' + wdir + ':'
                                 + wdir + ' -t combinelab/salmon chmod -R 777 ' + outf], shell=True)
                # docker run --rm -v /data/vdam/proj/immgenGithub:/tmp -t combinelab/salmon salmon quant -i /data/vdam/proj/immgenGithub/data/salmon_index103/ --libType A -p 12 --validateMappings -r /data/vdam/proj/immgenGithub/inputs/fastqFolder/sample1.fastq -o /mnt/fast/dario_and_vi/proj/neoantigens_local/immgenGithub/outfiles/salmon_quant103/sample1.fastq.quant
    else:  # paired-end reads
        with pymp.Parallel(ncore) as p1:
            for index in p1.range(0, len(files)):
                #for index in range(0, len(files)):
                f2 = files[index].replace("_1.fastq", "_2.fastq")
                if os.path.exists(f2):
                    outf = outdir + os.path.basename(files[index]).replace(
                        "_1.fastq", ".fastq") + ".quant"
                    subprocess.call(['docker run --rm -v ' + wdir + ':' + wdir + ' -t combinelab/salmon salmon quant -i '
                                     + root + '/data/salmon_index103/ --libType A -p 12 --validateMappings -o '
                                     + outf + ' -1 ' + files[index] + ' -2 ' + f2], shell=True)
                    subprocess.call(['docker run --rm -v ' + wdir + ':'
                                     + wdir + ' -t combinelab/salmon chmod -R 777 ' + files[index]], shell=True)
                    subprocess.call(['docker run --rm -v ' + wdir + ':'
                                     + wdir + ' -t combinelab/salmon chmod -R 777 ' + f2], shell=True)
    print("DONE salmonQuantification")


def combineQuantFiles(dirs, outdir, sourceArg):
    """ Combine all quantified expression into table

    Keyword arguments:
    dirs -- directory of all folders contained quantified expressions
    """
    n = 1
    splDF = pd.DataFrame()
    for d in dirs:
        if os.path.exists(outdir + "/" + d + "/quant.sf"):
            ID = d.partition(".fastq.quant")[0]
            if (sourceArg == 1):
                ID = ID[0:12]
            df = pd.read_csv(outdir + "/" + d + "/quant.sf", delimiter="\t",
                             header=0)
            if n == 1:
                splDF["Name"] = df["Name"]
            splDF[ID] = df.TPM
            n += 1
    splDF["Name"] = [ENST.split(".")[0] for ENST in splDF["Name"]]
    splDF = splDF.set_index('Name', drop=True)
    splDF.to_csv(out + "/splicing/allSplicingTPM.csv",
                 sep=",", header=True, index=True)
    return(splDF)


def quantileNormalizeSpl(sourceArg):
    """ Quantile normalize splicing expression

    Keyword arguments:
    sourceArg -- 1:TCGA type, 0:others
    """
    spl9mers = pd.read_csv(out + "/splicing/allSplicingTPM_9mers.tsv",
                           delimiter="\t", header=0, index_col=0)  # rows ENSTs
    spl9mers = spl9mers.drop(["MEAN", "MEDIAN", "SUM", "STD"], axis=1)
    if sourceArg == 1:
        spl9mers = filterSplicingCancerIDsOnly(spl9mers)
        spl_qr = quantile_normalize(spl9mers)  # about 30 mins
    else:
        spl_qr = quantile_normalize(spl9mers)  # about 30 mins
    return(spl_qr)


def filterSplicingCancerIDsOnly(splDF):
    """ Filter patient IDs columns that are cancer only

    Keyword arguments:
    splDF -- splicing expression dataframe with numerical values (columns: patient IDs, rows : transcript ID ENST)
    """
    ind = [ID[13:15] for ID in splDF.columns]  # remove normal samples in ge
    ind = [i for i, v in enumerate(ind) if int(v) < 6]
    splDF = splDF.iloc[:, ind]
    splDF.columns = splDF.columns.str[:12]
    return(splDF)


def quantile_normalize(splDF):
    """ Quantile normalize splicing expression

    Keyword arguments:
    splDF: gene expression dataframe with numerical values (columns: patient IDs, rows : transcript ID ENST)
    """
    df_sorted = pd.DataFrame(np.sort(splDF.values, axis=0),
                             index=splDF.index,
                             columns=splDF.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn = splDF.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)


def formatHLAsToMhcFlurry(hlas, ID):
    '''
    Convert list of HLAs format to MHCFlurry format. i.e. HLA-A*02:01 to HLA-A01:01

    hlas -- list of Optitype's predicted HLA
    ID -- patient ID
    '''
    print("FormatHLAsToMhcFlurry")
    hlas = " ".join([h for h in hlas])  # HLA-A*02:01 HLA-A*02:02
    hlas = hlas.replace("*", "")  # HLA-A02:01 HLA-A02:02
    with open(out + "/optitype_result/" + ID + "_optitype.txt", "w") as textFile:
        textFile.write(hlas)


def optitype(source, fastqFolder):
    '''
    Check if fastq is single-end or pair-end reads.
    Write HLA types file of predicted patients. Else, predict HLA types of new patients.

    hlas -- list of Optitype's predicted HLA
    ID -- patient ID
    '''
    print("optitype")
    files = glob.glob(fastqFolder + "/*_1.fastq")
    preHLAsDF = pd.read_csv(optitypeHLAs, delimiter="\t",
                            header=0, index_col=0)
    if len(files) != 0:  # paired-end reads
        allHLAs = []
        for f1 in files:
            fbase = f.replace(fastqFolder + "/",
                              "").replace("_1", "")  # sample1.fastq
            # if source is TCGA, get 12 digit
            if (source == 1):
                ID = f1[0:12]
            else:
                ID = fbase.replace(".fastq", "")  # sample1
            # if ID in pre-predicted TCGA HLAs
            if ID in list(preHLAsDF.sample_id):
                hlas = preHLAsDF.iloc[preHLAsDF['sample_id']
                                      == ID, 'allele']
            else:
                f2 = f1.replace("_1.fastq", "_2.fastq")
                deleteAndCreateOptitypeFolder(fbase)
                callOptitypeDocker(f1, f2, fbase)
                allHLAs = formatOptitypeOutput(fbase, ID, allHLAs)
        printAllHLAsToFile(allHLAs)
    else:  # single-end reads
        files = glob.glob(fastqFolder + "/*.fastq")
        #print(files)
        allHLAs = []
        for f1 in files:
            print(f1)  # wdir + sample1.fastq
            f2 = ""
            fbase = f1.replace(fastqFolder + "/", "")  # sample1.fastq
            if (source == 1):
                ID = fbase[0:12]
            else:
                ID = fbase.replace(".fastq", "")  # sample1
            if ID in list(preHLAsDF.sample_id):
                print(ID)
                hlas = preHLAsDF.iloc[preHLAsDF['sample_id']
                                      == ID, 'allele']
                allHLAs.extend([h for h in hlas])
                allHLAs = list(set(allHLAs))
                formatHLAsToMhcFlurry(hlas, ID)
            else:
                deleteAndCreateOptitypeFolder(fbase)
                callOptitypeDocker(f1, f2, fbase)
                allHLAs = formatOptitypeOutput(fbase, ID, allHLAs)
        printAllHLAsToFile(allHLAs)
    print("DONE ALL Optitype")


def deleteAndCreateOptitypeFolder(fname):
    tmp = out + '/optitype_result/' + fname + '_optitype'
    if os.path.exists(tmp):
        try:
            shutil.rmtree(tmp)
        except OSError as e:
            print("Error: %s : %s" % (tmp, e.strerror))
    os.makedirs(tmp)
    os.chmod(tmp, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)


def callOptitypeDocker(f1, f2, outF):
    print("Start Optitype")
    #
    tmp = 'docker run --rm -v ' + wdir + ':' + wdir
    tmpOut = out + '/optitype_result/' + outF + '_optitype'
    if f2 != "":
        subprocessCmd = tmp + ' -t fred2/optitype --input ' + \
            f1 + ' ' + f2 + ' --rna --outdir ' + tmpOut
        subprocess.call([subprocessCmd], shell=True)
        subprocessCmd = tmp + ' -t combinelab/salmon chmod -R 777 ' + tmpOut + '*'
        subprocess.call([subprocessCmd], shell=True)
    else:
        subprocessCmd = tmp + ' -t fred2/optitype --input ' + \
            f1 + ' --rna --outdir ' + tmpOut
        #print(subprocessCmd)
        subprocess.call([subprocessCmd], shell=True)
        #
        subprocessCmd = tmp + ' -t combinelab/salmon chmod -R 777 ' + tmpOut + '*'
        subprocess.call([subprocessCmd], shell=True)
    print("Done Optitype")


def formatOptitypeOutput(f, ID, allHLAs):
    for rt, dirs, files in os.walk(out + '/optitype_result/' + f + "_optitype"):
        for file in files:
            if file.endswith("result.tsv"):
                tmp = os.path.join(rt, file)
    #print(tmp)
    hlas = pd.read_csv(tmp, delimiter="\t",
                       header=0, index_col=0)
    hlas = [h for h in list(hlas.iloc[0, 0:5])
            if pd.isnull(h) == False]
    allHLAs.extend(hlas)
    formatHLAsToMhcFlurry(hlas, ID)
    return(allHLAs)


def printAllHLAsToFile(allHLAs):
    allHLAs = " ".join(set([h.replace("*", "") for h in allHLAs]))
    print(allHLAs)
    #allHLAs.replace("*", "")
    outf = open(out + "/optitype_result/all_optitype.txt", 'w')
    outf.write(allHLAs)
    outf.close()


# #####################################   MAIN  ###############################


def main():
    """
    Normalize gene expression data if True.
    Extract missense mutation peptides.
    Quantify and noramlize splicing expression.
    Predict patient's HLA types
    """
    # ######  parse parameters  ###############################
    parser = argparse.ArgumentParser(description='Prep data')

    parser.add_argument("-e", help="gene expression data")
    parser.add_argument("-m", help="mutation annotation data")
    parser.add_argument("-f", help="fastq folder that contains fastq files")
    parser.add_argument("-source", required=False, default=1,
                        help="1 if data source is tcga , else 0")
    parser.add_argument("-norm", required=False, default=False,
                        help="if normalize data")
    parser.add_argument("-ncore", required=False, default=3,
                        help="core number for parallel processing")

    args = parser.parse_args()

    if args.e is None and args.m is None and args.f is None:
        parser.error('No action requested, add -e, -m or -f')

    if args.f is not None:
        if len(glob.glob(args.f + "/*.fastq")) == 0:
            parser.error("No fastq file is found in this directory.")

    args = parser.parse_args()

    # ######  process data  ###############################
    os.chdir(wdir)  # user's current working directory

    arr = []
    if args.e is not None:
        arr.append("ge")
    if args.m is not None:
        arr.append("mut")
    if args.f is not None:
        arr.append("spl")
        # arr.append("opt")
    print(arr)

    pymp.config.nested = True
    global ncore
    ncore = int(args.ncore)
    # parallel processing each dataType
    with pymp.Parallel(ncore) as p:
        for ind in p.range(0, len(arr)):
            if arr[ind] == "ge":
                # print(args.e)
                # print(args.source)
                # print(args.norm)
                processGeneExpression(args.e, args.source, args.norm)
                print("DoneGeneExpression")
            elif arr[ind] == "mut":
                # print(args.m)
                # print(ncore)
                extract9mersMissenseMutation.main(args.m, ncore)
                print("DoneMut")
            elif arr[ind] == "spl":
                # print(args.f)
                # print(args.source)
                with pymp.Parallel(2) as p2:
                    for ind in p2.range(0, 2):
                        if ind == 1:
                            processSplicing(args.f, args.source)
                        else:
                            optitype(args.source,  args.f)
                        print("DoneFastq")

if __name__ == '__main__':
    main()
