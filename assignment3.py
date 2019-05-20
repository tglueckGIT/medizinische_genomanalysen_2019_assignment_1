#! /usr/bin/env python3

import vcf
import httplib2
import json

__author__ = 'Glueck Tobias'


##
##
## Aim of this assignment is to annotate the variants with various attributes
## We will use the API provided by "myvariant.info" - more information here: https://docs.myvariant.info
## NOTE NOTE! - check here for hg38 - https://myvariant.info/faq
## 1) Annotate the first 900 variants in the VCF file
## 2) Store the result in a data structure (not in a database)
## 3) Use the data structure to answer the questions
##
## 4) View the VCF in a browser
##

class Assignment3:
    
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
        
        ## Call annotate_vcf_file here
        self.vcf_path = "chr16.vcf"  # TODO

    def annotate_vcf_file(self):
        '''
        - Annotate the VCF file using the following example code (for 1 variant)
        - Iterate of the variants (use first 900)
        - Store the result in a data structure
        :return:
        '''    
        print("TODO")
                
        ##
        ## Example loop
        ##
        
        ## Build the connection
        h = httplib2.Http()
        headers = {'content-type': 'application/x-www-form-urlencoded'}
                
        params_pos = []  # List of variant positions
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for counter, record in enumerate(vcf_reader):
                params_pos.append(record.CHROM + ":g." + str(record.POS) + record.REF + ">" + str(record.ALT[0]))
                
                if counter >= 899:
                    break
        
        ## Build the parameters using the list we just built
        params = 'ids=' + ",".join(params_pos) + '&hg38=true'
        
        ## Perform annotation
        res, con = h.request('http://myvariant.info/v1/variant', 'POST', params, headers=headers)
        annotation_result = con.decode('utf-8')
        
        self.annotation_json = json.loads(annotation_result)
        ## TODO now do something with the 'annotation_result'
        
        ##
        ## End example code
        ##
        
        return self.annotation_json  ## return the data structure here
    
    
    def get_list_of_genes(self):
        for g in self.annotation_json:
            if "cadd" in g:
                if "genename" in g["cadd"]["gene"]:
                    return(g["cadd"]["gene"]["genename"])
    
    
    def get_num_variants_modifier(self):
        c = 0
        for i in self.annotation_json:
            if "snpeff" in i:
                for j in i["snpeff"]["ann"]:
                    if j["putative_impact"] == "MODIFIER":
                        c += 1
        return c


        
    
    def get_num_variants_with_mutationtaster_annotation(self):
        c = 0
        for i in self.annotation_json:
            if "dbnsfp" in i:
                if "mutationtaster"] in i["dbnsfp"]:
                    c += 1
        return c
        
    
    def get_num_variants_non_synonymous(self):
        c = 0
        for i in self.annotation_json:
            if 'cadd' in i:
                if i['cadd']['consequence'] == 'NON_SYNONYMOUS':
                    c += 1
        return c
        
    
    def view_vcf_in_browser(self):
        '''
        - Open a browser and go to https://vcf.iobio.io/
        - Upload the VCF file and investigate the details
        :return:
        '''
   
        ## Document the final URL here
        print("TODO")
            
    
    def print_summary(self):
        self.annotate_vcf_file()
        print("List of genes: " + self.get_list_of_genes())
        print("Number of variants: " + self.get_num_variants_modifier())
        print("Number of variants with mutations: " + self.get_num_variants_with_mutationtaster_annotation())
        print("Number of non synonymous variants: " + self.get_num_variants_non_synonymous())

    
    
def main():
    print("Assignment 3")
    assignment3 = Assignment3()
    assignment3.print_summary()
    print("Done with assignment 3")
        
        
if __name__ == '__main__':
    main()
   
    



