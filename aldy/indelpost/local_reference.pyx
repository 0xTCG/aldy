from pysam.libcfaidx cimport FastaFile


cdef class UnsplicedLocalReference:
    
    def __init__(self, str chrom, int pos, int ref_len, int window, FastaFile reference):
            
        cdef:
            int local_ref_start = max(0, pos - window * 10) # 0-based coordinate
            str unspliced_local_reference = reference.fetch(chrom, local_ref_start, min(pos + window * 10, ref_len))
        
        self.chrom = chrom
        self.pos = pos
        self.ref_len = ref_len
        self.window = window
        self.local_ref_start = local_ref_start
        self.unspliced_local_reference = unspliced_local_reference

    def fetch_ref_seq(self, target_pos, window):
        #pos_diff = target_pos - self.pos
        #lt_end = self.window + pos_diff
        #local_ref_len = min(target_pos + self.window * 3, self.ref_len) - max(0, target_pos - self.window * 3)
        
        self.left_len = target_pos - max(0, target_pos - window * 3)
        return self.get_ref_seq(max(0, target_pos - window * 3), min(target_pos + window * 3, self.ref_len))
        

    def get_ref_seq(self, start, end):
        # start, end: 0-based
        start_idx = start - self.local_ref_start
        return self.unspliced_local_reference[start_idx : start_idx + (end - start)]
        

        
