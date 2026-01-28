import pysam
import sys
import collections
import math
import argparse

# Script: annotate_vcf_with_pop_stats.py
# Description: 
#   Calculates population statistics for a VCF file, stratified by Sex and Country.
#   It updates the VCF INFO field with these statistics.
#
# Generated Tags (where {C} is the country code, e.g., 'PT'):
#   - Base Metrics: AC, AN, AF, AC_Het, AC_Hom, AC_Hemi, MAF, HWE, NS, F_MISSING
#   - Stratifications: Country, Sex, Country+Sex for all base metrics.
#   - If --raw is specified, base metrics become AC_raw, AN_raw, etc.

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description="Annotate VCF with population statistics")
parser.add_argument("input_vcf", help="Input VCF file")
parser.add_argument("output_vcf", help="Output VCF file")
parser.add_argument("sample_info", help="Sample info file (TSV: ID, Sex)")
parser.add_argument("country_code", help="Country code for stratification")
parser.add_argument("--raw", action="store_true", help="Suffix stats with _raw (e.g. AC_raw)")

args = parser.parse_args()

# --- HWE Calculation Helper ---
def hwe_test(obs_hets, obs_hom1, obs_hom2):
    """
    Calculates the Hardy-Weinberg Equilibrium p-value using a chi-squared test.
    For bi-allelic sites only.
    """
    total_samples = obs_hets + obs_hom1 + obs_hom2
    if total_samples == 0:
        return 1.0 # Cannot calculate

    # Allele counts
    p_allele_count = 2 * obs_hom1 + obs_hets
    q_allele_count = 2 * obs_hom2 + obs_hets
    total_alleles = p_allele_count + q_allele_count
    
    if total_alleles == 0:
        return 1.0

    # Allele frequencies
    p_freq = p_allele_count / total_alleles
    q_freq = q_allele_count / total_alleles

    # Expected genotype counts
    exp_hom1 = p_freq * p_freq * total_samples
    exp_hom2 = q_freq * q_freq * total_samples
    exp_hets = 2 * p_freq * q_freq * total_samples

    # Chi-squared test
    chi_sq = 0.0
    if exp_hom1 > 0: chi_sq += (obs_hom1 - exp_hom1)**2 / exp_hom1
    if exp_hom2 > 0: chi_sq += (obs_hom2 - exp_hom2)**2 / exp_hom2
    if exp_hets > 0: chi_sq += (obs_hets - exp_hets)**2 / exp_hets
    
    # p-value from chi-squared distribution with 1 degree of freedom
    # This is a simplification; a proper implementation would use a lookup table or a library function.
    # Using a common approximation for p-value from chi-squared value.
    # For 1 d.f., p-value can be calculated from the normal distribution CDF.
    # P(X^2 > k) = 2 * P(Z > sqrt(k)) where Z is standard normal.
    z = math.sqrt(chi_sq)
    p_value = math.erfc(z / math.sqrt(2)) # More stable than 1 - erf(...)

    return p_value

input_vcf_path = args.input_vcf
sample_info_path = args.sample_info
output_vcf_path = args.output_vcf
country_code = args.country_code
is_raw = args.raw

# --- 1. Load Sample Metadata ---
sample_sex = {}
print(f"Loading sample info from {sample_info_path}...")
with open(sample_info_path, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        parts = line.strip().split('\t')
        if len(parts) < 2: continue
        sid, raw_sex = parts[0], parts[1].strip().upper()
        if raw_sex in ['M', 'MALE', 'XY', '1']: sample_sex[sid] = 'M'
        elif raw_sex in ['F', 'FEMALE', 'XX', '2']: sample_sex[sid] = 'F'
        else: sample_sex[sid] = None

# --- 2. Open VCF & Setup Groups ---
vcf_in = pysam.VariantFile(input_vcf_path)
vcf_samples = list(vcf_in.header.samples)

indices_map = {'ALL': [], 'M': [], 'F': []}
for i, s in enumerate(vcf_samples):
    indices_map['ALL'].append(i)
    s_sex = sample_sex.get(s)
    if s_sex == 'M': indices_map['M'].append(i)
    elif s_sex == 'F': indices_map['F'].append(i)

print(f"Stats: Total={len(indices_map['ALL'])}, M={len(indices_map['M'])}, F={len(indices_map['F'])}")

# --- 3. Prepare Header ---
base_metrics = ['AC', 'AN', 'AF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'MAF', 'HWE', 'NS', 'F_MISSING']
metrics = []

# Modify metrics name if raw
for m in base_metrics:
    if is_raw:
        metrics.append(f"{m}_raw")
    else:
        metrics.append(m)

header = vcf_in.header

def add_header_line(tag_id, metric_base, desc):
    if tag_id in header.info: return
    
    number, dtype = 'A', 'Integer'
    if metric_base in ['AF', 'MAF', 'HWE', 'F_MISSING']: dtype = 'Float'
    if metric_base in ['AN', 'NS']: number = '1'
    
    header.info.add(tag_id, number=number, type=dtype, description=desc)

for m_name, m_base in zip(metrics, base_metrics):
    add_header_line(f"{m_name}", m_base, f"{m_name} for all samples")
    add_header_line(f"{m_name}_{country_code}", m_base, f"{m_name} for all samples in {country_code}")
    for sex in ['M', 'F']:
        add_header_line(f"{m_name}_{sex}", m_base, f"{m_name} for {sex} samples")
        add_header_line(f"{m_name}_{country_code}_{sex}", m_base, f"{m_name} for {sex} samples in {country_code}")

vcf_out = pysam.VariantFile(output_vcf_path, 'w', header=header)

# --- 4. Calculation Logic ---
def calc_stats(record, indices):
    if not record.alts: return None
    
    num_alts = len(record.alts)
    ac = [0] * num_alts
    an = 0
    n_het = [0] * num_alts
    n_hom = [0] * num_alts
    n_hemi = [0] * num_alts
    n_missing = 0
    obs_hom_ref = 0 # For HWE

    for idx in indices:
        sample = record.samples[vcf_samples[idx]]
        gt = sample.get('GT')
        
        if gt is None or None in gt:
            n_missing += 1
            continue
            
        ploidy = len(gt)
        an += ploidy
        
        # Count Alleles (AC)
        for allele in gt:
            if allele > 0: ac[allele - 1] += 1
                
        # Count Genotypes (Het, Hom, Hemi) and HomRef for HWE
        if ploidy == 1:
            if gt[0] > 0: n_hemi[gt[0] - 1] += 1
        elif ploidy == 2:
            a1, a2 = gt[0], gt[1]
            if a1 == a2:
                if a1 == 0: obs_hom_ref += 1
                elif a1 > 0: n_hom[a1 - 1] += 1
            else:
                if a1 > 0: n_het[a1 - 1] += 1
                if a2 > 0: n_het[a2 - 1] += 1
                    
    # --- Derived Statistics ---
    ns = len(indices) - n_missing
    f_missing = n_missing / len(indices) if len(indices) > 0 else 0.0
    
    af = [count / an for count in ac] if an > 0 else [0.0] * num_alts
    
    # MAF: Frequency of the second most common allele.
    # Ref allele freq + all alt freqs should sum to 1.
    ref_af = 1.0 - sum(af)
    all_freqs = sorted([ref_af] + af, reverse=True)
    maf = all_freqs[1] if len(all_freqs) > 1 else 0.0

    # HWE: only for bi-allelic sites
    hwe = 1.0
    if num_alts == 1:
        hwe = hwe_test(n_het[0], obs_hom_ref, n_hom[0])

    # Return dict with keys matching base_metrics
    return {
        'AC': ac, 'AN': an, 'AF': af, 'AC_Het': n_het, 'AC_Hom': n_hom, 
        'AC_Hemi': n_hemi, 'NS': ns, 'F_MISSING': f_missing, 'MAF': maf, 'HWE': hwe
    }

# --- 5. Iterate and Process ---
for record in vcf_in:
    stats_all = calc_stats(record, indices_map['ALL'])
    stats_m = calc_stats(record, indices_map['M'])
    stats_f = calc_stats(record, indices_map['F'])
    
    if stats_all is None: 
        vcf_out.write(record)
        continue
        
    def set_info(tag_suffix, stats):
        if not stats: return
        for m_name, m_base in zip(metrics, base_metrics):
            # Construct final tag: e.g. AC_raw + _PT_M
            final_tag = f"{m_name}{tag_suffix}"
            
            val = stats[m_base]
            if m_base == 'HWE':
                val = [val] 
            elif m_base == 'MAF':
                 val = [val]
            try:
                record.info[final_tag] = val
            except Exception as e:
                print(f"Error setting {final_tag} with value {val}: {e}")

    # tag_suffix logic needs to align with header creation
    # Header: {m_name}, {m_name}_{country}, {m_name}_{sex}, {m_name}_{country}_{sex}
    # set_info suffix should be "", "_{country}", "_{sex}", "_{country}_{sex}"
    
    set_info("", stats_all)
    set_info(f"_{country_code}", stats_all)
    set_info("_M", stats_m)
    set_info(f"_{country_code}_M", stats_m)
    set_info("_F", stats_f)
    set_info(f"_{country_code}_F", stats_f)

    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
