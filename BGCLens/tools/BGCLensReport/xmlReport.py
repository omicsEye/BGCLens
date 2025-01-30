#!/usr/bin/python
# Classes to represent genome identification XML report

import xml.etree.ElementTree as ET
import xml.dom.minidom
# from pathoscope.pathodb import dbUtils


class Gene:
    value = ""
    protein_id = ""
    locus_tag = ""
    product = ""
    ref_name = ""

    def __init__(self, value):
        self.value = value

    def buildElement(self):
        root = ET.Element("gene")
        if self.protein_id:
            root.set("protein_id", self.protein_id)
        if self.locus_tag:
            root.set("locus_tag", self.locus_tag)
        if self.product:
            root.set("product", self.product)
        if self.ref_name:
            root.set("ref_name", self.ref_name)
        root.text = self.value
        return root


class Read:
    readName = ""
    readSequence = ""

    def __init__(self, readName):
        self.readName = readName

    def buildElement(self):
        root = ET.Element("read_sequence")
        if self.readName:
            root.set("read_name", self.readName)
        root.text = self.readSequence
        return root


class Contig:
    value = ""
    contig = ""
    length = 0
    ref_name = ""

    def __init__(self, value):
        self.value = value

    def buildElement(self):
        root = ET.Element("ti_contig")
        if self.ref_name:
            root.set("ref_name", self.ref_name)
        if self.length:
            root.set("length", self.length)
        if self.contig:
            root.set("contig", self.contig)
        root.text = self.value
        return root


class RelativeAmount:
    value = 100
    count = 0

    def __init__(self, value):
        self.value = value

    def buildElement(self):
        root = ET.Element("relativeAmount")
        root.set("count", str(self.count))
        root.text = str(self.value)
        return root


class Taxonomy:
    value = ""
    taxon_id = 0

    def __init__(self, value):
        self.value = value

    def buildElement(self):
        root = ET.Element("taxonomy")
        root.set("taxon_id", str(self.taxon_id))
        root.text = self.value
        return root


class Organism:
    name = ""
    relativeAmount = None
    taxonomy = None
    genus = ""
    species = ""
    strain = ""
    nearestNeighbor = ""
    genes = []
    reads = []
    contigs = []

    def __init__(self, name):
        self.name = name

    def setRelativeAmount(self, relativeAmount):
        self.relativeAmount = relativeAmount

    def setTaxonomy(self, taxonomy):
        self.taxonomy = taxonomy

    def addGene(self, gene):
        self.genes.append(gene)

    def addContig(self, contig):
        self.contigs.append(contig)

    def buildElement(self):
        root = ET.Element("organism")
        if self.relativeAmount:
            root.append(self.relativeAmount.buildElement())
        if self.taxonomy:
            root.append(self.taxonomy.buildElement())
        organismName = ET.SubElement(root, "organismName")
        organismName.text = self.name
        if self.genus:
            genusElem = ET.SubElement(root, "genus")
            genusElem.text = self.genus
        if self.species:
            speciesElem = ET.SubElement(root, "species")
            speciesElem.text = self.species
        if self.strain:
            strainElem = ET.SubElement(root, "strain")
            strainElem.text = self.strain
        if self.nearestNeighbor:
            nearestNeighborElem = ET.SubElement(root, "nearestNeighbor")
            nearestNeighborElem.text = self.nearestNeighbor
        if self.genes:
            genesElem = ET.SubElement(root, "genes")
            for gene in self.genes:
                genesElem.append(gene.buildElement())
        if self.contigs:
            contigsElem = ET.SubElement(root, "contigs")
            for contig in self.contigs:
                contigsElem.append(contig.buildElement())
        if self.reads:
            readsElem = ET.SubElement(root, "reads")
            for read in self.reads:
                readsElem.append(read.buildElement())
        return root


class Organisms:
    organisms = []
    numAlignedReads = 0
    numMappedGenomes = 0

    def __init__(self):
        self.organisms = []

    def buildElement(self):
        root = ET.Element("organisms")
        root.set("num_aligned_reads", str(self.numAlignedReads))
        root.set("num_mapped_genomes", str(self.numMappedGenomes))
        for organism in self.organisms:
            root.append(organism.buildElement())
        return root


def writeElementXML(element, elementXMLFile):
    xmlString = ET.tostring(element, encoding="UTF-8", method="xml")
    xml1 = xml.dom.minidom.parseString(xmlString)
    prettyXmlString = xml1.toprettyxml(encoding="UTF-8")
    # print prettyXmlString #debug
    with open(elementXMLFile, 'w') as f:
        f.write(prettyXmlString)


# tree = ET.ElementTree(element)
# ET.dump(tree)
# tree.write("element.xml")

# =======================================
def writePathoXML(run_parameters, outputFileName, h_annoT, h_ti_contig,
                  h_refRead, h_refScore, h_gisPerTi, h_tiRef, reads, h_readSequence, samFile, mySqlConf):
    root = ET.Element("root")
    program = ET.SubElement(root, "program")
    name = ET.SubElement(program, "name")
    name.text = "pathoscope2.0"
    options = ET.SubElement(program, "options")
    options.text = samFile + " " + run_parameters  # TODO, split this long string to sub fields
    hostTaxon = ""
    organismsElement = buildOrganismsElement(h_annoT, h_ti_contig, hostTaxon,
                                             h_refRead, h_refScore, h_gisPerTi, h_tiRef, reads, h_readSequence, samFile,
                                             mySqlConf)
    root.append(organismsElement)
    writeElementXML(root, outputFileName)


# =======================================
def buildOrganismsElement(h_annoT, h_ti_contig, hostTaxon, h_refRead, h_refScore,
                          h_gisPerTi, h_tiRef, reads, h_readSequence, samFile, mySqlConf):
    NAs = 'X'
    useMysql = True
    con = None
    # (hostname,port,user,passwd,defaultDb)=range(5)
    (_, _, _, passwd, _) = range(5)
    if mySqlConf[passwd] == NAs:  # then, we do not use mysql
        useMysql = False
    organismsObj = Organisms()

    readCnt = len(reads)
    hostScore = 0
    if len(hostTaxon) > 0:
        try:
            hostScore = h_refScore[hostTaxon]
        except:
            hostScore = 0
    numTargetReads = readCnt - hostScore
    organismsObj.numAlignedReads = numTargetReads
    organismsObj.numMappedGenomes = len(h_gisPerTi)
    # if useMysql:
    #     con = dbUtils.init_mysql_innocentive(mySqlConf, 0)
    for ti in h_gisPerTi:
        refIdName = h_tiRef.get(ti, [ti, ti])
        refId = refIdName[0]
        score = h_refScore.get(refId, 0)

        organismName = refIdName[1]
        lineage = ''
        # if taxonomyLevelF:
        # if useMysql:
        #     organismName, lineage = dbUtils.findOrganismLineage(con, ti)
        organism = Organism(organismName)
        if useMysql:
            words = organismName.split()
            length = len(words)
            if length > 0:
                organism.genus = words[0]
            if length > 1:
                organism.species = words[1]
            if length > 2:
                organism.strain = words[2]

        organism.relativeAmount = RelativeAmount(score)
        organism.relativeAmount.count = len(h_refRead.get(refId, [-1]))
        organism.taxonomy = Taxonomy(lineage)
        organism.taxonomy.taxon_id = ti
        genes = []
        if h_annoT.get(ti, -1) != -1:
            for giList in h_annoT[ti]:
                gene = Gene(giList[1])
                if giList[2] and giList[2] != NAs:
                    gene.locus_tag = giList[2]
                if giList[3] and giList[3] != NAs:
                    gene.protein_id = giList[3]
                if giList[4] and giList[4] != NAs:
                    gene.ref_name = giList[4]
                if giList[5] and giList[5] != NAs:
                    gene.product = giList[5]
                genes.append(gene)
            organism.genes = genes

        # add contig information
        contigs = []
        j = 0
        ctgs = h_ti_contig.get(ti, [])
        for c in ctgs:
            ti_contig = ti + '_ctg_' + str(j)
            contig2 = Contig(ti_contig)
            contig2.ref_name = c[0]  # make sure that all string format only available in xml
            contig2.length = str(c[1])
            contig2.contig = c[2]
            j += 1
            contigs.append(contig2)
        organism.contigs = contigs

        # add read information
        reads = []
        readnames = h_refRead.get(refId, [])
        for readname in readnames:
            read = Read(readname)
            read.readSequence = h_readSequence[readname]
            reads.append(read)
        organism.reads = reads

        organismsObj.organisms.append(organism)

    organismsObj.organisms = sorted(organismsObj.organisms,
                                    key=lambda x: x.relativeAmount.value, reverse=True)
    organismsElement = organismsObj.buildElement()

    # if con:
    #     dbUtils.mysql_close(con)
    return organismsElement


# =======================================
def buildOrganismsElementExample():
    organismsObj = Organisms()
    organism = Organism("Host")
    organism.genus = "Homo"
    organism.species = "Sapiens"
    # organism.strain = "American"
    organism.relativeAmount = RelativeAmount(97.482)
    organism.relativeAmount.count = 1332476
    organism.taxonomy = Taxonomy("Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; "
                                 "Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo.")
    organism.taxonomy.taxon_id = 9606
    genes = []
    gene = Gene("xylB3")
    gene.protein_id = "AAS63914.1"
    gene.locus_tag = "YP_3766"
    gene.product = "putative carbohydrate kinase"
    gene.ref_name = "AAS63914.1"
    genes.append(gene)
    gene = Gene("firA")
    gene.protein_id = "CAL19719.1"
    gene.locus_tag = "YPO1054"
    gene.product = "UDP-3-o-[3-hydroxymyristoyl] glucosamineN-acyltransferase"
    gene.ref_name = "CAL19719.1"
    genes.append(gene)
    organism.genes = genes
    organism.reads = ["HT8EEZ101A215J", "HT8EEZ101EPOHR"]
    organismsObj.organisms.append(organism)
    organismsElement = organismsObj.buildElement()
    return organismsElement


# =======================================
def buildReportXMLExample(descriptionFile):
    string = ""
    with open(descriptionFile, 'r') as infile:
        for line in infile:
            string = string + line.strip()
    root = ET.fromstring(string)
    organismsElement = buildOrganismsElementExample()
    root.append(organismsElement)
    writeElementXML(root, 'organismreport.xml')


# =======================================
def buildOrganismsXMLExample():
    organismsElement = buildOrganismsElementExample()
    writeElementXML(organismsElement, 'organisms.xml')
