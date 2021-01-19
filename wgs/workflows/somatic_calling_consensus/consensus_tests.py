import random
import vcf
import gzip
import os
import consensus
import pypeliner
from shutil import copyfile
import csv
from collections import namedtuple

def get_testdata(testdir):

    museq = os.path.join(testdir, "museq_paired_germline.vcf.gz")
    mutect = os.path.join(testdir,  "mutect_germline.vcf.gz")
    strelka_indel = os.path.join(testdir, "strelka_indel_germline.vcf.gz")
    strelka_snv = os.path.join(testdir, "strelka_snv_germline.vcf.gz")


    pypeliner.commandline.execute("wgt", "-O", museq, "https://wgs-test-sets.s3.amazonaws.com/consensus_test_vcfs/somatic/museq_paired_germline.vcf.gz")
    pypeliner.commandline.execute("wgt",  "-O", museq+".tbi","https://wgs-test-sets.s3.amazonaws.com/consensus_test_vcfs/somatic/museq_paired_germline.vcf.gz.tbi")
    pypeliner.commandline.execute("wgt",  "-O", mutect,"https://wgs-test-sets.s3.amazonaws.com/consensus_test_vcfs/somatic/mutect_germline.vcf.gz")
    pypeliner.commandline.execute("wgt",  "-O", mutect+".tbi","https://wgs-test-sets.s3.amazonaws.com/consensus_test_vcfs/somatic/mutect_germline.vcf.gz.tbi")
    pypeliner.commandline.execute("wgt",  "-O", strelka_indel,"https://wgs-test-sets.s3.amazonaws.com/consensus_test_vcfs/somatic/strelka_indel_germline.vcf.gz")
    pypeliner.commandline.execute("wgt",  "-O", strelka_indel+".tbi","https://wgs-test-sets.s3.amazonaws.com/consensus_test_vcfs/somatic/strelka_indel_germline.vcf.gz.tbi")
    pypeliner.commandline.execute("wgt",  "-O", strelka_snv,"https://wgs-test-sets.s3.amazonaws.com/consensus_test_vcfs/somatic/strelka_snv_germline.vcf.gz")
    pypeliner.commandline.execute("wgt",  "-O", strelka_snv+".tbbi","https://wgs-test-sets.s3.amazonaws.com/consensus_test_vcfs/somatic/strelka_snv_germline.vcf.gz.tbi")

    return {"museq": museq, "mutect": mutect, "strelka_indel": strelka_indel, "strelka_snv": strelka_snv}, "SA609N", "SA609T"


def get_random_record(reader):
    chrom = random.choice(list(map(str, list(range(1, 23)) + ["X"])))
    lines1 = reader.fetch(chrom)
    return random.choice([next(lines1) for _ in range(10000)])

def get_test_line(v1, v2, n=10000):
    #choose a random chrom
    not_matched = False
    reader1 = vcf.Reader(filename=v1)
    reader2 = vcf.Reader(filename=v2)
    line = get_random_record(reader1)
    while not_matched == False:
        try:
            next(reader2.fetch(line.CHROM, start=line.POS-1, end=line.POS))
            line = get_random_record(reader1)

        except StopIteration:
            not_matched=True
            return line


def record_to_str(record, f, tmpdir):

    dummyfile = os.path.join(tmpdir, "dummy")
    pypeliner.commandline.execute("scp", f, dummyfile)

    r = vcf.Reader(open(dummyfile, "rt"))
    w = vcf.Writer(open(dummyfile, "at"), r)  

    template_record = next(r)

    ffs = w._map(str, [record.CHROM, record.POS, template_record.ID, record.REF] \
            + [w._format_alt(record.ALT), template_record.QUAL or '.', w._format_filter(template_record.FILTER),
                w._format_info(template_record.INFO)])

    if template_record.FORMAT:
        ffs.append(template_record.FORMAT)

    samples = [w._format_sample(template_record.FORMAT, sample)
        for sample in template_record.samples]   

    w.close()

    #delete dummy file
    pypeliner.commandline.execute("rm", dummyfile)

    output = "\t".join(ffs+samples)
    return output


def add_entry_to_vcf(record, v, tmpdir):
    '''
    you CANNOT append to a file wiht pyvcf, as it automatically writes out a vcf header with the first call to
    write_record. the work around is the make a writer on a dummy tmp vcf, use that writer to  get the formatted
    str version of thevariant, than normal io append that string to an unzipped v, then delete the dummy and recipe and index v
    '''
    #unzip
    pypeliner.commandline.execute("gunzip", v)
    gunzipped = v.replace(".gz", "")

    recordstr = record_to_str(record, gunzipped, tmpdir)

    #write
    with open(gunzipped, "at") as f:
        f.write(recordstr + "\n")
    f.close()

    #rezip
    s = gunzipped+"_s"
    pypeliner.commandline.execute("vcf-sort", gunzipped ,">", s)
    pypeliner.commandline.execute("mv", s, gunzipped)
    pypeliner.commandline.execute("bgzip", gunzipped)

    #now reindex
    pypeliner.commandline.execute("tabix", "-f", "-p", "vcf", v)



def run_consensus_test(vcfs, test1, test2, tmpdir):

    testfile1 = os.path.join(tmpdir, os.path.basename(vcfs[test1]))
    testfile2 = os.path.join(tmpdir, os.path.basename(vcfs[test2]))

    copyfile(vcfs[test1], testfile1)
    copyfile(vcfs[test1] + ".tbi", testfile1 + ".tbi")
    copyfile(vcfs[test2], testfile2)

    copyfile(vcfs[test2] + ".tbi", testfile2 + ".tbi")

    vcfs[test1] = testfile1
    vcfs[test2] = testfile2

    #make test
    record = get_test_line(testfile1, testfile2)
    
    add_entry_to_vcf(record, testfile2, tmpdir) 

    reader1 = vcf.Reader(filename=testfile1)
    assert next(reader1.fetch(record.CHROM, start=record.POS-1, end=record.POS))
    reader2 = vcf.Reader(filename=testfile2)
    assert next(reader2.fetch(record.CHROM, start=record.POS-1, end=record.POS))

    #run consensus
    chroms = list(map(str, list(range(1, 23)) + ["X"]))
    consensus_vcf = os.path.join(tmpdir, "consensus.vcf")
    counts = os.path.join(tmpdir, "counts.tsv")

    consensus.main(vcfs["museq"], 
        vcfs["strelka"], 
        vcfs["mutect"], 
        vcfs["strelka_indel"], 
        consensus_vcf, 
        counts, 
        chroms    
    )

    #read in consensus vcf and make sure the variant is included
    consensus_reader = pd.read_csv(consensus_vcf, sep="\t")
    assert not consensus_reader[(consensus_reader.POS==record.POS) & consensus_reader["#CHROM"]==record.CHROM].empty


def run_normalize_test(vcfs, test1, test2, a1, r1, a2, r2, tmpdir):

    testfile1 = os.path.join(tmpdir, os.path.basename(vcfs[test1]))
    testfile2 = os.path.join(tmpdir, os.path.basename(vcfs[test2]))

    copyfile(vcfs[test1], testfile1)
    copyfile(vcfs[test1] + ".tbi", testfile1 + ".tbi")
    copyfile(vcfs[test2], testfile2)

    copyfile(vcfs[test2] + ".tbi", testfile2 + ".tbi")

    vcfs[test1] = testfile1
    vcfs[test2] = testfile2

    #make test
    record = get_test_line(testfile1, testfile2)

    record.REF = r1
    record.ALT = a1
    add_entry_to_vcf(record, testfile1, tmpdir) 

    record.REF = r2
    record.ALT = a2
    add_entry_to_vcf(record, testfile2, tmpdir) 

    reader1 = vcf.Reader(filename=testfile1)
    assert next(reader1.fetch(record.CHROM, start=record.POS-1, end=record.POS))
    reader2 = vcf.Reader(filename=testfile2)
    assert next(reader2.fetch(record.CHROM, start=record.POS-1, end=record.POS))

    #run consensus
    chroms = list(map(str, list(range(1, 23)) + ["X"]))
    consensus_vcf = os.path.join(tmpdir, "consensus.vcf")
    counts = os.path.join(tmpdir, "counts.tsv")

    consensus.main(vcfs["museq"], 
        vcfs["strelka"], 
        vcfs["mutect"], 
        vcfs["strelka_indel"], 
        consensus_vcf, 
        counts, 
        chroms    
    )

    #read in consensus vcf and make sure the variant is included
    consensus_reader = pd.read_csv(consensus_vcf, sep="\t")
    assert not consensus_reader[(consensus_reader.POS==record.POS) & consensus_reader["#CHROM"]==record.CHROM].empty
    


def test_museq_strelka(testdir):
    vcfs, _, _  =  get_testdata(testdir)
    run_consensus_test(vcfs, "museq", "strelka", testdir)

def test_museq_strelka_indel(testdir):
    vcfs, _, _ =  get_testdata(testdir)
    run_consensus_test(vcfs, "museq", "strelka_indel", testdir)

def test_museq_mutect(testdir):
    vcfs, _, _  =  get_testdata(testdir)
    run_consensus_test(vcfs, "museq", "mutect", testdir)

def test_strelka_mutect(testdir):
    vcfs, _, _  =  get_testdata(testdir)
    run_consensus_test(vcfs, "strelka", "mutect", testdir)

def test_strelka_strelka_indel(testdir):
    vcfs, _, _  =  get_testdata(testdir)
    run_consensus_test(vcfs, "strelka", "strelka_indel", testdir)

def test_mutect_strelka_indel(testdir):
    vcfs, _, _ =  get_testdata(testdir)
    run_consensus_test(vcfs, "mutect", "strelka_indel", testdir)


def test_normalization_1(testdir):
    # A-> AT and AT->ATT
    vcfs, _, _  =  get_testdata(testdir)
    r1 =  "A"
    a1 = "AT"
    r2 = "AT"
    a2 = "ATT"
    consensus_variant =  run_normalize_test(vcfs, "strelka", "museq",  a1, r1, a2, r2, testdir)

    assert consensus_variant.REF == "A"
    assert consensus_variant.ALT == "AT"

def test_normalization_2( testdir):
    # A-> AT and AT->ATT
    vcfs, _, _  =  get_testdata(testdir)
    r1 =  "GTA"
    a1 = "G"
    r2 = "GTATA"
    a2 = "GTA"
    consensus_variant =  run_normalize_test(vcfs, "strelka", "museq",  a1, r1, a2, r2, testdir)

    assert consensus_variant.REF == "GTA"
    assert consensus_variant.ALT == "G"

def test_normalization_3(testdir):
    # A-> AT and AT->ATT
    vcfs, _, _  =  get_testdata(testdir)
    r1 =  "CA"
    a1 = "CTA"
    r2 = "C"
    a2 = "CT"
    consensus_variant =  run_normalize_test(vcfs, "strelka", "museq",  a1, r1, a2, r2, testdir)

    assert consensus_variant.REF == "C"
    assert consensus_variant.ALT == "CT"

def test_normalization_4( testdir):
    # A-> AT and AT->ATT
    vcfs, _, _ =  get_testdata(testdir)
    r1 =  "T"
    a1 = "TAC"
    r2 = "T"
    a2 = "TACACAC"
    consensus_variant =  run_normalize_test(vcfs, "strelka", "museq",  a1, r1, a2, r2, testdir)

    assert consensus_variant.REF == "T"
    assert consensus_variant.ALT == "TAC"

def test_normalization():
    refs =  ["T", "T"]
    alts = ["TAC", "TACACAC"]

    consensus_ref, consensus_alt = consensus.normalize(refs, alts)
    assert consensus_ref == "T"
    assert consensus_alt == "TAC"

    refs =  ["CA", "C"]
    alts = ["CTA", "CT"]

    consensus_ref, consensus_alt = consensus.normalize(refs, alts)
    assert consensus_ref == "C"
    assert consensus_alt == "CT"
 
    refs =  ["GTA", "GTATA"]
    alts = ["G", "GTA"]

    consensus_ref, consensus_alt = consensus.normalize(refs, alts)
    assert consensus_ref == "GTA"
    assert consensus_alt == "G"
 
    refs =  ["A", "AT"]
    alts = ["AT", "ATT"]

    consensus_ref, consensus_alt = consensus.normalize(refs, alts)
    assert consensus_ref == "A"
    assert consensus_alt == "AT"


def test_get_counts():

    #make example  call
    Call = namedtuple("Call", "DP TAR TIR RC AC AD AU TU")
    normal_data = vcf.model._Call(None, "normal", Call(1,[1],[1],1,1, [1,1], [1], [1]) )
    tumor_data = vcf.model._Call(None, "tumor", Call(1,[1],[1],1,1, [1,1], [1], [1] ))

    test = vcf.model._Record("1", 10, 1, "A", [None], 1, 1, None, None, None, [normal_data, tumor_data])

    assert consensus.get_counts(test, "museq_snv", "tumor", "normal", "A", "T") == (1, [1], 1, 1, [1], 1)
    assert consensus.get_counts(test, "strelka_indel",  "tumor", "normal","A", "T") == (1, [1], 1, 1, [1], 1)
    assert consensus.get_counts(test, "strelka_snv",  "tumor", "normal","A", "T") == (1, [1], 1, 1, [1], 1)
    assert consensus.get_counts(test, "mutect",  "tumor", "normal","A", "T") == (1, [1], 1, 1, [1], 1)



museq = "/juno/work/shah/isabl_data_lake/analyses/78/36/7836/results/somatic/SA609T/SA609T_museq_paired_annotated.vcf.gz" 
strelka = "/juno/work/shah/isabl_data_lake/analyses/78/36/7836/results/somatic/SA609T/SA609T_strelka_snv_annotated.vcf.gz"
strelka_indel = "/juno/work/shah/isabl_data_lake/analyses/78/36/7836/results/somatic/SA609T/SA609T_strelka_indel_annotated.vcf.gz"
mutect = "/juno/work/shah/isabl_data_lake/analyses/78/36/7836/results/somatic/SA609T/SA609T_mutect.vcf.gz"


test_get_counts()
# make_test("/juno/work/shah/abramsd/CODE/TEST/SA1181_N_freebayes_germline.vcf.gz",
# "/juno/work/shah/abramsd/CODE/TEST/SA1181_N_museq_single_annotated.vcf.gz")

