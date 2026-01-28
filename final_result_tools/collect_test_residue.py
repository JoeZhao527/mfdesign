"""
Load each sample from samples_gt_tagged.jsonl. Extract a mapping of entry -> CDR sequences.
Output format:
{
    "7lxn_H_L_A": {
        "HCDR1": "...",
        "HCDR2": "...",
        "HCDR3": "...",
        "LCDR1": "...",
        "LCDR2": "...",
        "LCDR3": "...",
    },
    "7noz_R__CD": {
        "HCDR1": "...",
        "HCDR2": "...",
        "HCDR3": "...",
        "LCDR1": null,
        "LCDR2": null,
        "LCDR3": null
    }
}
"""
import json
import re
import sys

def extract_cdr_sequence(seq_string):
    """Extract the actual sequence from a string like '<HCDR1>XXXXXXX</HCDR1>'."""
    # Remove XML-like tags and extract the sequence
    match = re.search(r'<[^>]+>(.*?)</[^>]+>', seq_string)
    if match:
        return match.group(1)
    return seq_string

def process_jsonl(input_file, output_file, is_mid_train: bool = True):
    """Process JSONL file and create CDR sequence mapping."""
    entry_to_cdrs = {}
    
    # Mapping from internal CDR names to output names
    cdr_name_mapping = {
        'H1': 'HCDR1',
        'H2': 'HCDR2',
        'H3': 'HCDR3',
        'L1': 'LCDR1',
        'L2': 'LCDR2',
        'L3': 'LCDR3'
    }
    
    print(f"Processing {input_file}...")
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            data = json.loads(line)
            
            # Extract entry name from meta
            if is_mid_train:
                entry = data.get('meta', {}).get('entry')
                sample_id = data.get('sample_id', None)
                if not entry:
                    continue
                
                if "dsasa" not in sample_id:
                    continue
            else:
                entry = data.get('sample_id').split('_')[0]
            
            # Initialize entry if not seen before
            if entry not in entry_to_cdrs:
                entry_to_cdrs[entry] = {
                    'HCDR1': None,
                    'HCDR2': None,
                    'HCDR3': None,
                    'LCDR1': None,
                    'LCDR2': None,
                    'LCDR3': None
                }
            
            # Parse target_text to get CDR sequences
            target_text = data.get('target_text', '')
            if target_text:
                target_data = json.loads(target_text)
                cdr_sequences = target_data.get('answer', {}).get('cdr_sequences', {})
                
                # Extract sequences for each CDR
                for cdr_key, cdr_info in cdr_sequences.items():
                    if cdr_key in cdr_name_mapping:
                        output_key = cdr_name_mapping[cdr_key]
                        seq_string = cdr_info.get('seq', '')
                        if seq_string:
                            sequence = extract_cdr_sequence(seq_string)
                            # Check for duplicates and assert they match
                            if entry_to_cdrs[entry][output_key] is None:
                                entry_to_cdrs[entry][output_key] = sequence
                            else:
                                # Duplicate entry - assert the sequences match
                                existing_seq = entry_to_cdrs[entry][output_key]
                                if existing_seq != sequence:
                                    raise AssertionError(
                                        f"Mismatch for entry '{entry}' {output_key}: "
                                        f"existing='{existing_seq}' vs new='{sequence}'"
                                    )
    
    print(f"Processed {line_num} lines, found {len(entry_to_cdrs)} unique entries")
    
    # Write output JSON file
    with open(output_file, 'w') as f:
        json.dump(entry_to_cdrs, f, indent=2)
    
    print(f"Output written to {output_file}")
    
    # Print summary
    entries_with_all_heavy = sum(1 for cdrs in entry_to_cdrs.values() 
                                 if all(cdrs[f'HCDR{i}'] is not None for i in [1, 2, 3]))
    entries_with_all_light = sum(1 for cdrs in entry_to_cdrs.values() 
                                  if all(cdrs[f'LCDR{i}'] is not None for i in [1, 2, 3]))
    entries_with_both = sum(1 for cdrs in entry_to_cdrs.values() 
                           if all(cdrs[f'HCDR{i}'] is not None for i in [1, 2, 3]) and
                              all(cdrs[f'LCDR{i}'] is not None for i in [1, 2, 3]))
    
    print(f"\nSummary:")
    print(f"  Entries with all heavy chain CDRs: {entries_with_all_heavy}")
    print(f"  Entries with all light chain CDRs: {entries_with_all_light}")
    print(f"  Entries with both heavy and light: {entries_with_both}")
    
    return entry_to_cdrs

if __name__ == '__main__':
    input_file = 'stage3_test_20260124_234002.jsonl'
    output_file = 'cdr_sequences_mapping.json'
    
    final_res_files = [
        dict(
            input_file = 'final_result_cdr/trail_0.jsonl',
            output_file = 'final_result_cdr/trail_0.json'
        ),
        dict(
            input_file = 'final_result_cdr/trail_1.jsonl',
            output_file = 'final_result_cdr/trail_1.json'
        ),
    ]
    
    for i_o in final_res_files:
        process_jsonl(i_o['input_file'], i_o['output_file'], False)
    
    
    
    
    
    # if len(sys.argv) > 1:
    #     input_file = sys.argv[1]
    # if len(sys.argv) > 2:
    #     output_file = sys.argv[2]
    
    # entry_to_cdrs = process_jsonl(input_file, output_file)