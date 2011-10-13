import stats

from collections import namedtuple

FIELDS = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()
FILTER = FIELDS.index('FILTER')
Vcf = namedtuple('Vcf', FIELDS)


CALL_DELIM = ':'
FIELD_DELIM = '\t'
INFO_DELIM = ';'


TYPES = { 
    'RA': int, 
    'AA': int,
    'GQ': float
}

def parse_call(format, call):
    dat = dict([(f, TYPES.get(f, str)(v)) 
        for f, v in zip(format.split(CALL_DELIM), call.split(CALL_DELIM))])
    if 'GT' in dat and dat['GT'] != '.':
        dat['GT'] = int(dat['GT'][0]), int(dat['GT'][-1])
    return dat
    
def parse_info(info):
    info = info.split(';')
    dat = dict((x.split('=') for x in info if '=' in x))
    for x in info:
        if '=' not in x:
            dat[x] = True
    return dat
    
def parse_calls(fields, line): 
    toks = line.split()
    if ',' in fields.ALT: return [] 
    format = fields.FORMAT
    calls = [parse_call(format, x) for x in toks[9:]]
    return [x for x in calls if x['GT'] != '.']

def parse_fields(line):
    toks = line.split(FIELD_DELIM)[:len(FIELDS)]
    toks[7] = parse_info(toks[7])
    vcf = Vcf(*toks)
    return vcf
    
def parse_repeats(rep, SEP='|', EQ=':'):
    toks = rep.split(SEP)
    rep = {}
    for dat in toks: 
        bases, count = dat.split(EQ)
        rep[bases] = int(count)
    return rep 
    
def add_filter(line, extra_filter):
    toks = line.split(FIELD_DELIM)
    filt = toks[FILTER] if toks[FILTER] != '.' else ''
    filt = filt.split(INFO_DELIM) + [extra_filter]
    filt = INFO_DELIM.join(filt).lstrip(INFO_DELIM)
    toks[FILTER] = filt
    return FIELD_DELIM.join(toks)
    
def is_homo_ins(field):
    ref = field.REF
    alt = field.ALT
    c = alt.count(ref)
    return c > 1 and alt == (ref * c)

def is_homo_del(field):
    ref = field.REF
    alt = field.ALT
    c = ref.count(alt)
    return c > 1 and ref == (alt * c)


def apply_filters(line, error_bias=True, homopolymer_indel=True):      
    
    max_indel_repeat = 3
    min_genotype_quality = 30
            
    fields = parse_fields(line)
    calls = parse_calls(fields, line)
    
    
    # cannot handle composite classes
    if ',' in fields.ALT: 
        line = add_filter(line, 'composite')
        return line
    
    # biased site filter: tests for all sites coming from an error probability
    if error_bias:
        if calls: 
            passed, tv, ab = stats.bias_test(calls)
            if not passed: 
                line = add_filter(line, 'errbias=%(tv).2f,%(ab).2f' % locals())
    
    if homopolymer_indel: 
        if 'REPEAT' in fields.INFO: 
            repeats = parse_repeats(fields.INFO['REPEAT'])
            is_ins, is_del = is_homo_ins(fields), is_homo_del(fields) 
            if is_del:
                ref = fields.REF
                if ref == ref[0] * len(ref):
                    ref = ref[0]
                ref_count = repeats.get(ref, 0)
            else: 
                ref_count = repeats.get(fields.REF, 0)
            # print is_ins, is_del, ref_count
            if (is_ins or is_del) and ref_count > max_indel_repeat:
                line = add_filter(line, 'homo_indel=%s' % ref_count)
    
    if min_genotype_quality and calls: 
        quals = [x['GQ'] for x in calls if x['GT'] != 0]
        # print quals, any([x>min_genotype_quality for x in quals])
        if not any([x>min_genotype_quality for x in quals]):
            line = add_filter(line, 'max_GQ=%s' % max(quals))
            
        
    return line