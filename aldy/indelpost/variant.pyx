#cython: profile=False

from .utilities import *
from .localn import make_aligner, align, findall_indels
from pysam.libcfaidx cimport FastaFile
from pysam.libcbcf cimport VariantFile


cdef class NullVariant:
    """This class is returned when :class:`~indelpost.VariantAlignment` cannot find the target indel in the BAM file. Boolean expression evaluates to `False <https://docs.python.org/3/library/stdtypes.html#boolean-values>`__. Ref and Alt alleles are the 
    reference base at the locus.
    
    Parameters
    ----------
    chrom : string
        contig name

    pos : integer
        1-based genomic position

    reference : pysam.FastaFile
        reference FASTA file supplied as
        `pysam.FastaFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.FastaFile>`__ object.
    """
    def __init__(self, str chrom, int pos, FastaFile reference):
        self.chrom = chrom
        self.pos = pos
        self.ref = reference.fetch(chrom, pos - 1, pos)
        self.alt = self.ref
        self.reference = reference

    def __getstate__(self):
        return (self.chrom, self.pos, self.ref, self.alt, self.reference.filename)

    def __setstate__(self, state):
        self.chrom = state[0]
        self.pos = state[1]
        self.ref = state[2]
        self.alt = state[3]

        self.reference = FastaFile(state[4])

    def __bool__(self):
        return False
    
    def __eq__(self, other):
        
        if isinstance(other, Variant):
            return False
       
        chrom_equal = (self.chrom == other.chrom)
        pos_equal = (self.pos == other.pos)
        ref_equal = (self.ref == other.ref)
        alt_equal = (self.alt == other.alt)
        return all([chrom_equal, pos_equal, ref_equal, alt_equal])

    def __hash__(self):
        hashable = (self.chrom, self.pos, self.ref, self.alt)
        return hash(hashable)


cdef class Variant:
    """This class accepts a VCF-style variant representation as input. 
    Equality holds between :class:`~indelpost.Variant` objects 
    if they are indentical in the normalized form.  

    Parameters
    ----------
    chrom : string
        contig name. 
    
    pos : integer
        1-based genomic position.
    
    ref : string
        VCF-style reference allele. 
    
    alt : string
        VCF-style alternative allele.
    
    reference : pysam.FastaFile
        reference FASTA file supplied as 
        `pysam.FastaFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.FastaFile>`__ object.
    
    Raises
    ------
    ValueError
        if the input locus is not defined in the reference or the input alleles contain letters 
        other than A/a, C/c, G/g, T/t, and N/n. 
        
    """
    def __init__(self, str chrom, int pos, str ref, str alt, FastaFile reference, skip_validation=False):
        self._chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.reference = reference

        
        if not skip_validation:
            self.chrom = self.__format_chrom_name(self._chrom, reference=reference)
            self.__validate()
        else:
            self.chrom = chrom
        
    def __getstate__(self):
        return (self.chrom, self.pos, self.ref, self.alt, self.reference.filename) 
        

    def __setstate__(self, state):        
        self.chrom = state[0]
        self.pos = state[1]
        self.ref = state[2]
        self.alt = state[3]
        
        self.reference = FastaFile(state[4])
    
    
    def __format_chrom_name(self, chrom, **kwargs):
        if kwargs.get("vcf", False):
            chrom_names = list(kwargs["vcf"].header.contigs)
        elif kwargs.get("reference", False):
            chrom_names = kwargs["reference"].references

        is_prefixed = True if chrom_names[0].startswith("chr") else False
        is_mt = True if "chrMT" in chrom_names or "MT" in chrom_names else False

        chrom = chrom.replace("chr", "")
        if chrom == "M" and is_mt:
            chrom = "MT"
        elif chrom == "MT" and not is_mt:
            chrom = "M"

        if is_prefixed:
            chrom = "chr" + chrom

        return chrom


    cpdef __validate(self):
        if not self.ref or not self.alt:
            raise ValueError("Allele may not be empty")

        if self.ref == self.alt:
            raise ValueError(
                "Not a variant: reference allele and alternate allele may not be identical"
            )

        bases = {"A", "C", "T", "G", "N", "a", "t", "c", "g", "n"}
        ref_lst, alt_lst = list(self.ref), list(self.alt)
        if not set(ref_lst) <= bases or not set(alt_lst) <= bases:
            self.ref = "".join([base if base in bases else "N" for base in ref_lst])
            self.alt = "".join([base if base in bases else "N" for base in alt_lst])

        # check if contig is valid
        try:
            if not self.reference.fetch(self.chrom, self.pos - 1 , self.pos):
                raise ValueError("The locus is not defined in the reference")
        except:
            raise ValueError("The locus is not defined in the reference")
    
    @property
    def variant_type(self):
        """ returns "I" if the net allele-length change is gain, "D" if loss, 
            "S" if signle-nucleotide substitution (zero net change), "M" if multi-nucleotide substitution 
            (zero net change). 
        """ 
        
        cdef int r_len, a_len
        cdef str var_type
            
        r_len, a_len = len(self.ref), len(self.alt)
        if r_len < a_len:
            var_type = "I"
        elif r_len > a_len:
            var_type = "D"
        elif a_len == 1:
            var_type = "S"
        else:
            var_type = "M"

        return var_type


    @property
    def is_del(self):
        """evaluates if :attr:`~indelpost.Variant.variant_type` is "D".
        """
        return self.variant_type == "D"


    @property
    def is_ins(self):
        """evaluates if :attr:`~indelpost.Variant.variant_type` is "I".
        """
        return self.variant_type == "I"

        
    @property
    def is_indel(self):
        """returns True if :attr:`~indelpost.Variant.variant_type` is "I" or "D".
        """
        return self.is_ins or self.is_del


    @property
    def indel_seq(self):
        """returns the inserted/deleted sequence for non-complex indels. None for substitutions.
        """
        if self.is_ins:
            return self.alt[len(self.ref) :]
        elif self.is_del:
            return self.ref[len(self.alt) :]
        else:
            return ""


    def __eq__(self, other):
        
        if isinstance(other, NullVariant):
            return False
        
        cdef Variant i, j
        i, j = self.normalize(), other.normalize()

        chrom_eq = (
            (i.chrom.replace("chr", "") == j.chrom.replace("chr", ""))
            or
            (i._chrom.replace("chr", "") == j._chrom.replace("chr", ""))
        )
             
        equivalent = (
                chrom_eq
                and i.pos == j.pos
                and j.ref.upper() == i.ref.upper()
                and i.alt.upper() == j.alt.upper()
        )
        
        return equivalent


    def __hash__(self):
        cdef Variant i = self.normalize() if self.is_indel else self
        hashable = (i._chrom, i.pos, i.ref, i.alt)

        return hash(hashable)


    def __dealloc__(self):
        pass


    @property
    def is_leftaligned(self):
        """returns True if the indel position is the smallest possible value (left-aligned). 
        """
        if self.ref[-1].upper() != self.alt[-1].upper():
            return True
        elif "N" in self.ref.upper() or "N" in self.alt.upper():
            return True

    
    @property
    def is_normalized(self):
        """returns True if left-aligned and the allele representations are minimal. 
        """
        if self.is_leftaligned:
            if len(self.ref) > 1 and len(self.alt) and (self.ref[0].upper() == self.alt[0].upper()):
                return False
            else:
                return True
        else:
            return False


    def normalize(self, inplace=False):
        """normalizes :class:`~indelpost.Variant` object.
        
        Parameters
        ----------

        inplace : bool
            returns None and normalizes this :class:`~indelpost.Variant` object if True. 
            Otherwise, returns a normalized copy of this object. Default to False.

        """
        cdef Variant i
        cdef str ihs
        cdef int n = 0, lhs_len = 0


        if inplace:
            i = self
        else:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference, skip_validation=True)
        
        condition_1 = i.ref[-1].upper() == i.alt[-1].upper() != "N" 
        lhs = i.reference.fetch(i.chrom, max(0, i.pos - 1 - 300), i.pos - 1)[::-1]
        lhs_len = len(lhs)
        while condition_1 and n < lhs_len:
            
            left_base = lhs[n]

            i.ref = left_base + i.ref[:-1]
            i.alt = left_base + i.alt[:-1]
            i.pos -= 1
            
            condition_1 = i.ref[-1].upper() == i.alt[-1].upper() != "N"

            n += 1 

        condition_2 = i.ref[0].upper() == i.alt[0].upper()
        condition_3 = len(i.ref) > 1 and len(i.alt) > 1
        while condition_2 and condition_3:
            i.ref = i.ref[1:]
            i.alt = i.alt[1:]
            i.pos += 1
            condition_2 = i.ref[0].upper() == i.alt[0].upper()
            condition_3 = len(i.ref) > 1 and len(i.alt) > 1

        if inplace:
            return None
        else:
            return i
    
    
    def generate_equivalents(self):
        """generates non left-aligned copies of :class:`~indelpost.Variant` object.
        """ 
        
        i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference, skip_validation=True).normalize()
        
        pos, ref, alt = i.pos, i.ref, i.alt
        is_ins = i.is_ins
        
        res = [i]
        if not i.is_indel:
            return res
        
        
        window = 300
        ref_lim = i.reference.get_reference_length(i.chrom)
        if i.is_non_complex_indel() and i.variant_type == "I":
            rt_flank = i.reference.fetch(i.chrom, i.pos, min(i.pos + window, ref_lim))
        else:
            if i.is_non_complex_indel() and i.variant_type == "D":
                event_len = len(i.indel_seq)
            else:
                event_len = len(i.ref) - 1  
            rt_flank = i.reference.fetch(i.chrom, i.pos + event_len, min(i.pos + event_len + window, ref_lim))

        n = 0
        while self == i and n < window:
            right_base = rt_flank[n]
            if is_ins:
                ref = alt[1]
                alt = alt[1:] + right_base
            else:
                alt = ref[1]
                ref = ref[1:] + right_base
            
            pos += 1 
            
            i = Variant(self.chrom, pos, ref, alt, self.reference, skip_validation=True)
            
            if self == i: 
                res.append(i)
            
            n += 1

        return res


    def _generate_equivalents_private(self):
        
        if self.is_non_complex_indel():
            return self.generate_equivalents()
        else:
            # define complex indel at the start and end of the deleted sequence 
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference, skip_validation=True)
            j = Variant(self.chrom, self.pos + len(self.ref), self.ref, self.alt, self.reference, skip_validation=True)
            return [i, j]
    

    def _get_indel_seq(self, how=None):
        if self.is_non_complex_indel():
            return self.indel_seq
        else:
            if how == "I":
                return self.alt[1:]
            elif how == "D":
                return self.ref[1:]
    
    def _reduce_complex_indel(self, to=None):
        if self.is_non_complex_indel():
            return NullVariant(self.chrom, self.pos, self.reference)
        else:
            if to == "I":
                return Variant(self.chrom, self.pos, self.alt[0], self.alt, self.reference, skip_validation=True)
            elif to == "D":
                return Variant(self.chrom, self.pos, self.ref, self.ref[0], self.reference, skip_validation=True)
         
    
    def query_vcf(self, VariantFile vcf, matchby="normalization", window=50, indel_only=True, as_dict=True):
        """returns a `list <https://docs.python.org/3/library/stdtypes.html#list>`__ of VCF records matching 
        this :class:`~indelpost.Variant` object as `dictionary <https://docs.python.org/3/library/stdtypes.html#dict>`__ (default).

        Parameters
        ----------
        vcf : pysam.VariantFile
            VCF file to be queried. 
            Supply as 
            `pysam.VariantFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantFile>`__ object.

        matchby : string
            - "normalization" (default) matches by normalizing VCF records.
            - "locus" matches by the normalized genomic locus.  
            - "exact" finds the exact matche without normalization. 
                
        window : integer
            searches the VCF records in indel position +/- window.
            
        indel_only : bool
            returns matched indel and SNV records if False. Meaningful when matchby is "locus".
               
        as_dict : bool
            returns `pysam.VariantRecord <https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantRecord>`__ if False.
        """
        matchbys = ["normalization", "locus", "exact"]
        if not matchby in matchbys:
            raise ValueError("match by one of: %s" % matchbys)
        
        if self.variant_type == "S":
            leftaligned_pos, window = self.pos, 1
        else:
            leftaligned_pos = self.normalize().pos
        
        chrom = self.__format_chrom_name(self.chrom, vcf=vcf)
        
        searchable = vcf.fetch(chrom, leftaligned_pos - 1, leftaligned_pos - 1 + window)
        
        if not searchable:
            return []

        records = to_flat_list(
            [
                to_flat_vcf_records(rec)
                for rec in searchable
            ]
        )
        
        hits = [
                record.orig
                for record in records
                if match_indels(
                    Variant(self.chrom, record.pos, record.ref, record.alt, self.reference), 
                    self, 
                    matchby,
                    indel_only,
                )
        ]
        
        if as_dict:
            hits = [
                {
                    "CHROM": hit.chrom,
                    "POS": hit.pos,
                    "ID": hit.id,
                    "REF": hit.ref,
                    "ALT": ",".join(list(hit.alts)),
                    "QUAL": hit.qual,
                    "FILTER": to_dict(hit.filter),
                    "INFO": to_dict(hit.info),
                    "FORMAT": to_dict(hit.format),
                    "SAMPLES": to_dict(hit.samples),
                }
                for hit in hits
            ]
        
        return hits


    def left_flank(self, window=50, normalize=False):
        """extracts the left-flanking reference sequence. See also :meth:`~indelpost.Variant.right_flank`.

        Parameters
        ----------
        window : integer
            extract the reference sequence [variant_pos - window, variant_pos]. 
        normalize : bool
            if True, the normalized indel position is used as the end of the flanking sequence.
        """ 
        if normalize:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference, skip_validation=True)
        else:
            i = self

        if i.is_non_complex_indel():
            pos = i.pos
        else:
            pos = i.pos - 1
        
        lt_flank = i.reference.fetch(i.chrom, max(0, pos - window), pos)
        
        return lt_flank

    
    def right_flank(self, window=50, normalize=False):
        """extracts the right-flanking reference sequence. See also :meth:`~indelpost.Variant.left_flank`. 

        Parameters
        ----------
        window : integer
            extract the reference sequence [variant_end_pos, variant_end_pos + window].
        normalize : bool
            if True, the normalized indel position is used as the start of the flanking sequence.
        """
        if normalize:
            i = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference, skip_validation=True)
        else:
            i = self

        ref_lim = i.reference.get_reference_length(i.chrom)
        if i.is_non_complex_indel() and i.variant_type == "I":
            rt_flank = i.reference.fetch(i.chrom, i.pos, min(i.pos + window, ref_lim))
        else:
            if i.is_non_complex_indel() and i.variant_type == "D":
                event_len = len(i.indel_seq)
            else:
                event_len = len(i.ref) - 1  
            rt_flank = i.reference.fetch(i.chrom, i.pos + event_len, min(i.pos + event_len + window, ref_lim))

        return rt_flank

    
    def count_repeats(self, by_repeat_unit=True):
        """counts indel repeats in the flanking reference sequences. The search window is
        defined by :meth:`~indelpost.Variant.left_flank` and :meth:`~indelpost.Variant.right_flank`.
        
        Parameters
        ----------
        by_repeat_unit : bool
            count by the smallest tandem repeat unit. For example, the indel sequence "ATATATAT" has
            tandem units "ATAT" and "AT". The occurrence of "AT" will be counted if True (default).
        """ 

        if self.is_non_complex_indel():
            seq = self.indel_seq
        else:
            seq = self.alt
        
        if by_repeat_unit:
            seq = to_minimal_repeat_unit(seq)

        lt_flank = self.left_flank()
        lr_repeat = repeat_counter(seq, lt_flank[::-1])
        rt_flank = self.right_flank()
        rt_repeat = repeat_counter(seq, rt_flank)
        
        return lr_repeat + rt_repeat
    

    def is_non_complex_indel(self):
        """returns True only if non-complex indel (False if complex indel or substitution).
        """
        i = self.normalize()
        ref, alt = i.ref, i.alt
        if len(ref) == len(alt):
            return False
        
        if ref[0] != alt[0]:
            return False
        
        the_shorter = ref if i.is_ins else alt
        if len(the_shorter) > 1:
            return False
        
        return True


    def decompose_complex_variant(self, match_score=3, mismatch_penalty=2, gap_open_penalty=4, gap_extension_penalty=0):
        """returns a `list <https://docs.python.org/3/library/stdtypes.html#list>`__ of 
        non-complex :class:`~indelpost.Variant` objects decomposed by the Smith-Waterman local alignment 
        with a given set of score/penalty.

        Parameters
        ----------
        match_score : integer
            default to 3.
        mismatch_penalty : integer
            default to 2.
        gap_open_penalty : integer
            default to 4.
        gap_extension_penalty : integer
            default to 0.    
        """
        if self.is_non_complex_indel():
            return [self]
        
        var = Variant(self.chrom, self.pos, self.ref, self.alt, self.reference, skip_validation=True).normalize()
        
        lt_pos = var.pos - 1
        rt_pos = var.pos - 1 + len(var.ref)
        
        window = 100
        mut_seq = self.reference.fetch(var.chrom, lt_pos - window, lt_pos) + var.alt + self.reference.fetch(var.chrom, rt_pos, rt_pos + window)
        ref_seq = self.reference.fetch(var.chrom, lt_pos - window, lt_pos + len(var.ref) + window)
        
        aln = align(make_aligner(ref_seq, match_score, mismatch_penalty), mut_seq, gap_open_penalty, gap_extension_penalty) 
        
        genome_aln_pos = lt_pos + 1 - window + aln.reference_start
        
        indels, snvs = findall_indels(aln, genome_aln_pos, ref_seq, mut_seq, report_snvs=True)
        
        variants = []
        if indels:
            for idl in indels:
                padding_base = idl["lt_ref"][-1]
                if idl["indel_type"] == "D":
                    ref = padding_base + idl["del_seq"]
                    alt = padding_base
                else:
                    ref = padding_base
                    alt = padding_base + idl["indel_seq"]

                variants.append(Variant(self.chrom, idl["pos"], ref, alt, self.reference, skip_validation=True))
        
        if snvs:
            for snv in snvs:
                variants.append(Variant(self.chrom, snv["pos"], snv["ref"], snv["alt"], self.reference, skip_validation=True))
             
        return variants
