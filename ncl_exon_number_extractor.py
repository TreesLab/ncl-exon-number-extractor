#! /usr/bin/env python
from __future__ import print_function

import sys
import re
import gzip
from collections import namedtuple
from itertools import groupby
from operator import itemgetter


JuncSite = namedtuple("JuncSite", ['chr', 'pos', 'strand', 'gene'])

class NCLevent:
    def __init__(self, raw_data):
        self.raw_data = raw_data[:]
        self.donor = JuncSite(raw_data[0], raw_data[1], raw_data[2], raw_data[6])
        self.accepter = JuncSite(raw_data[3], raw_data[4], raw_data[5], raw_data[7])
        self.intragenic = int(raw_data[8])
    
    
class ExonExtractor:
    def __init__(self):
        self.tid_pattern = re.compile('transcript_id "([^;]*)";')

    def _basic_info_getter(self, data):
        return [data[0], data[3], data[4], data[6]]
        
    def extract(self, data_iter):
        self.transcripts = {}
        self.exons = []
        
        exon_count = 0
        for line in data_iter:
            
            if line.startswith('#'):
                continue

            data = line.rstrip('\n').split('\t')
            region_type = data[2]

            if region_type == 'transcript':
                tid = re.search(self.tid_pattern, data[8]).group(1)
                transcript_len = int(data[4]) - int(data[3]) + 1
                self.transcripts[tid] = transcript_len
                
                exon_count = 0

            elif region_type == 'exon':
                exon_count += 1

                basic_info = self._basic_info_getter(data)
                exon_tid = re.search(self.tid_pattern, data[8]).group(1)
                assert(exon_tid == tid)

                self.exons.append(basic_info + [exon_tid, exon_count])

                
class JunctionSitesDB:
    def __init__(self):
        self.ncl_donor = {}
        self.ncl_accepter = {}
    
    def _get_donor_accepter(self, exons, di, ai):
        def group_exons(exs):
            for k, gp in groupby(sorted(exs, key=itemgetter(0)), key=itemgetter(0)):
                gp = [data[1:] for data in gp]
                yield [k, gp]

        donor = map(itemgetter(di, 4, 5), exons)
        accepter = map(itemgetter(ai, 4, 5), exons)

        grouped_donor = dict(group_exons(donor))
        grouped_accepter = dict(group_exons(accepter))

        return grouped_donor, grouped_accepter
    
    def generate_db(self, exon_data):
        for chrm, chrm_exons in groupby(exon_data, key=itemgetter(0)):
            chrm_exons = list(chrm_exons)
            plus_part = list(filter(lambda data: data[3]=='+', chrm_exons))
            minus_part = list(filter(lambda data: data[3]=='-', chrm_exons))

            plus_donor, plus_accepter = self._get_donor_accepter(plus_part, 2, 1)
            minus_donor, minus_accepter = self._get_donor_accepter(minus_part, 1, 2)

            self.ncl_donor[chrm] = {'+': plus_donor, '-': minus_donor}
            self.ncl_accepter[chrm] = {'+': plus_accepter, '-': minus_accepter}
            
    def get_junc_site(self, junc_site, donor_accepter):
        if donor_accepter == "donor":
            db = self.ncl_donor
        elif donor_accepter == "accepter":
            db = self.ncl_accepter
        else:
            raise ValueError("Only 'donor' or 'accepter' is valid!")
        
        try:
            res = db[junc_site.chr][junc_site.strand][junc_site.pos]
        except KeyError:
            print("No such site: {}".format(junc_site))
        else:
            return res
        
        
def get_longest_tid(tids, tid_len_dict):
    the_longest = max(tids, key=lambda tid: (tid_len_dict[tid], tid))
    return the_longest

def get_longest_common_tid(donor_tids, accepter_tids, tid_len_dict):
    common = set(donor_tids) & set(accepter_tids)
    if len(common) > 0:
        the_longest = get_longest_tid(common, tid_len_dict)
        return the_longest

def get_tid_exon_number(tid_data, tid):
    return dict(tid_data).get(tid)


def print_usage():
    usage_msg = \
        "Usage:\n"\
        "  cat test_NCLscan.result "\
        "| ./ncl_exon_number_extractor.py gencode.v28.annotation.gtf.gz "\
        "> test_NCLscan.exon_numbers.result"
    print(usage_msg, file=sys.stderr)


def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt')
    else:
        return open(filename)



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print_usage()
        exit(1)

    # generate db
    exon_extractor = ExonExtractor()
    with open_file(sys.argv[1]) as anno_data:
        exon_extractor.extract(anno_data)

    junc_sites_db = JunctionSitesDB()
    junc_sites_db.generate_db(exon_extractor.exons)

    # get exon numbers
    for line in sys.stdin:
        data = line.rstrip('\n').split('\t')
        ncl_event = NCLevent(data)
        donor = junc_sites_db.get_junc_site(ncl_event.donor, 'donor')
        accepter = junc_sites_db.get_junc_site(ncl_event.accepter, 'accepter')
        donor_tids = list(map(itemgetter(0), donor))
        accepter_tids = list(map(itemgetter(0), accepter))

        if ncl_event.intragenic:
            the_longest_common_tid = get_longest_common_tid(donor_tids, accepter_tids, \
                                                            exon_extractor.transcripts)
            if the_longest_common_tid:
                res_data = ncl_event.raw_data \
                            + [the_longest_common_tid, \
                                str(get_tid_exon_number(donor, the_longest_common_tid)), \
                                str(get_tid_exon_number(accepter, the_longest_common_tid))]
                print(*res_data, sep='\t')
                
            else:
                if (',' in ncl_event.donor.gene) or (',' in ncl_event.accepter.gene):
                    ncl_event.intragenic = 0
                else:
                    if len(sys.argv) >= 3:
                        if sys.argv[2] == 'all':
                            ncl_event.intragenic = 0
                    else:
                        res_data = ncl_event.raw_data + ['', '', '']
                        print(*res_data, sep='\t')
        
        if not ncl_event.intragenic:
            the_longest_tid_donor = get_longest_tid(donor_tids, exon_extractor.transcripts)
            the_longest_tid_accepter = get_longest_tid(accepter_tids, \
                                                        exon_extractor.transcripts)
            res_data = ncl_event.raw_data \
                        + ["{}|{}".format(the_longest_tid_donor, the_longest_tid_accepter), \
                            str(get_tid_exon_number(donor, the_longest_tid_donor)), \
                            str(get_tid_exon_number(accepter, the_longest_tid_accepter))]
            print(*res_data, sep='\t')
