import copy
import numpy as np
import pandas as pd


class Seqlib:
	"""
	A seqlib object for PDSB class.
	"""
    def __init__(self, ninds, nsites):
    	## generate the full sequence array
        self.ninds = ninds
        self.nsites = nsites
        self.seqs = self._simulate()

        ## store maf of the full seq array
        self.maf = self._get_maf()


    ## private functions only used during the init ------
    def _mutate(self, base):
        "converts a base to another base"
        diff = set("ACTG") - set(base) #Create a set of bases without the input base
        return np.random.choice(list(diff)) #Choose a random base from the subset "diff"
    
    def _simulate(self, ninds, nsites):
        "returns a random array of DNA bases with mutations"
        oseq = np.random.choice(list("ACGT"), size=self.nsites) #Creates a random array of As, Ts, Cs and Gs of length nsites 
        arr = np.array([oseq for i in range(self.ninds)]) #Creates ninds number of arrays of length nsites using the code above.
        muts = np.random.binomial(1, 0.1, (self.ninds, self.nsites)) #Creates a mask with a 10% success rate with dimensions of ninds x nsites 
        for col in range(self.nsites):
            newbase = mutate(self.arr[0, col]) #Selecting a random new base other than the base present for each position in the sequences.
            mask = muts[:, col].astype(bool) #creates a mask of bools for each position. 
            self.arr[:, col][mask] = newbase #mutate the base with the newbase if it passes the mask filter for mutations.
        missing = np.random.binomial(1, 0.1, (self.ninds, self.nsites)) #creates another mask with a 10% success rate
        arr[missing.astype(bool)] = "N" #replaces the base with an N if it passes the mask filter for missing bases.
        return arr

    def _get_maf(sefl):
        "returns the maf of the full seqarray while not counting Ns"
        ## init an array to fill and iterate over the columns
        maf = np.zeros(self.nsites)
        for col in range(self.nsites):
            ##select this column of bases
            thiscol = self.seqs[:, col]
            
            ## mask "N" bases and get the new length
            nmask = thiscol != "N"
            no_n_len = np.sum(nmask)

            ## mask "N" bases and get the first base
            first_non_n_base = thiscol[nmask][0]
            
            ## calculate maf of "N" masked bases
            freq = np.sum(thiscol[nmask] != first_non_n_base) / no_n_len
            	if freq > 0.5: 
               	 	maf[col] = 1 - freq
            	else: 
                	maf[col] = freq
        return maf
        
    ## Private functions that are called within other functions 
    def _filter_missing(self, maxmissing):
        "Returns a booleon filter True for columns with Ns > maxmissing"
        #This function identifies the seq positions with missing bases, removes those columns from all sequences and returns the filtered array.
        freqmissing = np.sum(self.arr == "N", axis=0) / self.arr.shape[0] #Sums # of times N is found in a column and divides it by the number of sequences (i.e.rows.
        return freqmissing > maxmissing #selects those columns that have N frequencies less than or equal to the given maxfreq.

    def _filter_maf(self, minmaf): #Calculates the minor allele freq and returns an array of sequences with elements/positions that surpass the maf. 
        "returns a boolean filter True for columns with a freq less than the minimum allele freq"
        return self.maf < minmaf
    
    ## public funtions ------
    def filter(self, minfreq, maxmissing):
        """
        Applies both the maf and missing filter functions to return a view of the filtered sequence array"
        """
        filter1 = self._filter_maf(minmaf)
        filter2 = self._filter_missing(maxmissing)
        fullfilter = filter1 + filter2
        return self.seqs[:, np.invert(fullfilter)]
    
    def filter_seqlib(self, minmaf, maxmissing): 
        """
        Applies maf and missing filters to the array and returns a copy of the seqlib object where the .seqs array has been filtered
        """        
        ## apply filters to get new array size
        newseqs = self.filter(minmaf, maxmissing)
        
        ## make new copy of the seqlib object
        newself = copy.deepcopy(self)
        newself.__init__(newseqs.shape[0], newseqs.shape[1])
        
        ##store the array (overwrite it)
        newself.seqs = newseqs
        
        ## call the _get_maf to match new array
        newself._get_maf()
        return newself
    
    def calculcate_statistics(self):
        """ 
        Returns a dataframe of statistics on the seqs array."""
        if self.seqs.size:
            nd = np.var(self.seqs == self.seqs[0], axis=0).mean() #Calculates the mean of the variance across the sequences -> mean nucleotide diversity
            mf = np.mean(
                np.sum(self.seqs != self.seqs[0], axis=0) / self.seqs.shape[0]) #Calculates average of invariant sites frequency -> the minor allele frequency  
            inv = np.all(self.seqs == self.seqs[0], axis=0).sum() #Sums the number of invariant sites (mutations)
            var = self.arr.shape[1] - inv #Subtracts the number of invariant sites from the number of sites -> variable sites
            return pd.Series(
                {"mean nucleotide diversity": nd,
                 "mean minor allele frequency": mf,
                 "invariant sites": inv,
                 "variable sites": var,
                })
        else: 
            print("seqs array is empty")