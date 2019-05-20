import mysql.connector
import os
import pysam
import pybedtools

__author__ = 'Glueck Tobias'

class Assignment1:
    
    def __init__(self):
        ## Your gene of interest
        self.gene = "CLDN14"
        self.geneList = []
        self.lineList = []
        self.genedict = self.download_gene_coordinates("hg38", "fetched_genes")
        self.bamfile = os.path.join(os.getcwd(), "chr21.bam")
        self.baifile = os.path.join(os.getcwd(), "chr21.bam.bai")
        self.samfile = pysam.AlignmentFile(self.bamfile, "rb")
        self.reads = list(self.samfile.fetch(int(self.genedict["txStart"]), int(self.genedict["txEnd"])))

    def download_gene_coordinates(self, genome_reference, file_name):
        ## TODO concept
        self.file = file_name
        
        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        
        ## Execute query
        cursor.execute(query)
        
        ## Write to file
        ## TODO this may need some work
        line = 0
        with open("output.txt", "w") as fh:
            for row in cursor:
            	line = line + 1
                fh.write(str(row) + "\n")
                if (row[0] == self.gene):
                	self.lineList.append(line)
                	self.geneList.append(row) ## all lines relevant to the gene are in this list now
        cursor.close()
        cnx.close()
        
        print("Done fetching data")
        
    def get_coordinates_of_gene(self):
        ## Use UCSC file
        strCoordinates = str(self.geneList[0][3]) + " - " + str(self.geneList[0][4])
        return(strCoordinates)
        
    def get_gene_symbol(self):
        return(self.geneList[0][0])
                        
    def get_sam_header(self):
        self.alignfile = pysam.AlignmentFile(self.file, "rb")
        strHeader = ""
        for i, j in self.alignfile.header["HD"].items():
            if i == "SO":
                strHeader += "\tSO (Sorting order of Alignments): " + str(j)
            if i == "VN":
                strHeader += "\tVN (Format version): " + str(j)
            if i == "GO":
                strHeader += "\tGO: (Grouping of alignments): " + str(j)
        return(strHeader)
        
    def get_properly_paired_reads_of_gene(self):
        strPaired = ""
        i = 0
        for read in self.reads:
            if read.is_proper_pair:
                i += 1
        if i==0:
            strPaired = "No properly paired read"
        else:
            strPaired = "Number of properly paired reads: " + str(i)
        return(strPaired)
        
    def get_gene_reads_with_indels(self):
        strIndels = ""
        i=0
        for read in self.reads:
            if not read.is_unmapped:
                cigar = read.cigar
                for (type, length) in cigar:
                    if (type ==1) or (type == 2):
                        i+=1
        if i == 0:
            strIndels = "No gene reads with indels"
        else:
            strIndels= "Number of gene reads with indels: " + str(i)
        return(strIndels)
        
    def calculate_total_average_coverage(self):
        bam = pybedtools.BedTool(self.bamfile)
        cov = bam.genome_coverage(bg=True, genome="hg38")
        i = 0
        sum = 0
        for line in cov:
            number = float(line[3])
            sum += number
            i+=1
        return(sum/i)
        
    def calculate_gene_average_coverage(self):
        a = pybedtools.BedTool(self.bamfile)
        b = a.genome_coverage(bg=True)
        sum = 0
        i=0

        for line in b:
            number = float(line[3])
            beg = int(line[1])
            if beg > self.genedict["tx:Start"]:
                if int(line[2]) <= self.genedict["txEnd"]:
                    sum += number
                    i+=1
        return(sum/i)
        
    def get_number_mapped_reads(self):
        i=0
        for read in self.reads:
            if not read.is_unmapped:
                i+=1
        return(i)

    def get_region_of_gene(self):
        return(self.geneList[0][2])
        
    def get_number_of_exons(self):
        return(self.geneList[0][6])
    
    
    def print_summary(self):
        print("Coordinateds of gene: " + str(self.get_coordinates_of_gene()))
        print("Gene symbol: " + str(self.get_gene_symbol()))
        print("Sam Header: " + str(self.get_sam_header()))
        print("Properly paired reads: " + str(self.get_properly_paired_reads_of_gene()))
        print("Reads with indels: " + str(self.get_gene_reads_with_indels()))
        print("Total average coverage: " + str(self.calculate_total_average_coverage()))
        print("Gene average coverage: " + str(self.calculate_gene_average_coverage()))
        print("Number of mapped reads: " + str(self.get_number_mapped_reads()))
        print("Region of gene: " + str(self.get_region_of_gene()))
        print("Number of exons: " + str(self.get_number_of_exons()))
    
    
def main():
    print("Assignment 1")
    assignment1 = Assignment1()
    assignment1.download_gene_coordinates("hg38", "chr21.bam")
    assignment1.print_summary()
    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    main()
    
    
