import csv
import os

def parse_sample_info(sample_info_path):
    """
    Parses a sample info file (TSV) which may or may not have a header.
    Assumes:
    - If header exists: looks for 'sample'/'id' and 'sex' (case-insensitive).
    - If no header: assumes col 0 is ID, col 1 is Sex.
    
    Returns a dictionary: {sample_id: sex_code}
    Sex codes are normalized to:
    - '1' (Male): M, MALE, 1, XY
    - '2' (Female): F, FEMALE, 2, XX
    - '0' (Unknown): anything else
    """
    # Normalization map
    norm = {
        'M': '1', 'MALE': '1', '1': '1', 'XY': '1',
        'F': '2', 'FEMALE': '2', '2': '2', 'XX': '2'
    }
    
    samples = {}
    
    if not os.path.exists(sample_info_path):
        raise FileNotFoundError(f"Sample info file not found: {sample_info_path}")

    with open(sample_info_path, 'r') as f:
        # Read a few lines to guess format and detect header
        lines = [f.readline() for _ in range(5)]
        f.seek(0)
        
        if not lines:
            return {}

        # Heuristic for header detection
        first_line = lines[0].lower()
        has_header = any(keyword in first_line for keyword in ["sample", "sex", "id"])
            
        reader = csv.reader(f, delimiter='\t')
        
        if has_header:
            header = next(reader)
            header_low = [x.lower() for x in header]
            
            # Find indices
            try:
                id_idx = next(i for i, x in enumerate(header_low) if any(k in x for k in ["sampleid", "sample", "id"]))
                sex_idx = next(i for i, x in enumerate(header_low) if "sex" in x)
            except StopIteration:
                # Fallback if keywords not explicitly found despite heuristic
                id_idx, sex_idx = 0, 1
        else:
            id_idx, sex_idx = 0, 1
            
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) <= max(id_idx, sex_idx):
                continue
            
            sid = row[id_idx].strip()
            raw_sex = row[sex_idx].strip().upper()
            
            samples[sid] = norm.get(raw_sex, '0')
            
    return samples

def get_nextflow_metadata(sample_info_path, country_code):
    """
    Generates rows for the Nextflow metadata CSV.
    Returns a list of lists: [['SAMPLE', 'SEX', 'ANCESTRY'], [id, sex, country], ...]
    """
    samples = parse_sample_info(sample_info_path)
    rows = [['SAMPLE', 'SEX', 'ANCESTRY']]
    for sid, sex in samples.items():
        rows.append([sid, sex, country_code])
    return rows
