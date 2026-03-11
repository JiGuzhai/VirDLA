from alphagenome.models import dna_client
import numpy as np
import argparse
import json
import os
config_file = 'test_box.json'
config = {}
if os.path.exists(config_file):
    with open(config_file, 'r', encoding='utf-8') as f:
        config = json.load(f)
parser = argparse.ArgumentParser(description='Virana')
parser.add_argument('-ot', type=str, default=config.get('ot', ''), help='')
parser.add_argument('-p1', type=str, default=config.get('p1', ''), help='')
parser.add_argument('-p2', type=str, default=config.get('p2', ''), help='')
parser.add_argument('-CDS1', type=str, default=config.get('CDS1', ''), help='')
parser.add_argument('-CDS2', type=str, default=config.get('CDS2', ''), help='')
parser.add_argument('-insulator', type=str, default=config.get('insulator', ''), help='')
parser.add_argument('-t', type=str, default=config.get('t', ''), help='')
args = parser.parse_args(args=[])
print("Config loaded from:", config_file)
model = dna_client.create("your_api_key")
def filter_noise(data, threshold_multiplier):
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    threshold = median + threshold_multiplier * mad
    filtered_data = data[data > threshold]
    if len(filtered_data) == 0:
        return np.mean(data)
    return np.mean(filtered_data)

def analyze_sequence(model, args):
    full_sequence = args.p1 + args.CDS1 + args.insulator + args.p2+ args.CDS2
    padded_sequence = full_sequence.upper().center(dna_client.SEQUENCE_LENGTH_16KB, 'N')
    print("API loaded")
    
    group1_output = None
    if args.ot:
        group1_output = model.predict_sequence(
            sequence=padded_sequence,
            requested_outputs=[dna_client.OutputType.RNA_SEQ, dna_client.OutputType.ATAC],
            ontology_terms=[args.ot]
        )
    
    def calculate_coordinates(padded_sequence, full_sequence):
        origin = (len(padded_sequence) - len(full_sequence)) // 2
        p1_s , p1_e = origin, origin + len(args.p1)
        CDS1_s, CDS1_e = p1_e, p1_e + len(args.CDS1)
        insulator_s, insulator_e = CDS1_e, CDS1_e + len(args.insulator)
        p2_s, p2_e = insulator_e, insulator_e + len(args.p2)
        CDS2_s, CDS2_e = p2_e, p2_e + len(args.CDS2)
        return origin, p1_s, p1_e, CDS1_s, CDS1_e, insulator_s, insulator_e, p2_s, p2_e, CDS2_s, CDS2_e

    origin, p1_s, p1_e, CDS1_s, CDS1_e, insulator_s, insulator_e, p2_s, p2_e, CDS2_s, CDS2_e = calculate_coordinates(padded_sequence, full_sequence)

    def analyze_region(output, group_name, region_name, start, end, noise_filtering=True, threshold_multiplier=-1):
        RNA_auc = 0.0
        try:
            RNAseq = output.rna_seq.values[start:end]
            if RNAseq.ndim == 2:
                RNAseq = RNAseq[:, 0]
            if RNAseq.size > 0:
                if np.any(np.isinf(RNAseq)) or np.any(np.isnan(RNAseq)):
                    RNAseq = RNAseq[~np.isinf(RNAseq) & ~np.isnan(RNAseq)]
                    if RNAseq.size > 0:
                        if noise_filtering:
                            filtered_data = filter_noise(RNAseq, threshold_multiplier)
                            filtered_data_array = np.full_like(RNAseq, filtered_data)
                            RNA_auc = np.trapezoid(filtered_data_array)
                        else:
                            RNA_auc = np.trapezoid(RNAseq)
                    else:
                        RNA_auc = 0.0
                        print(f"{group_name} {region_name} region RNAseq data all infinite or NaN")
                else:
                    if noise_filtering:
                        filtered_data = filter_noise(RNAseq, threshold_multiplier)
                        filtered_data_array = np.full_like(RNAseq, filtered_data)
                        RNA_auc = np.trapezoid(filtered_data_array)
                    else:
                        RNA_auc = np.trapezoid(RNAseq)
                print(f"{group_name} {region_name} region RNAseq AUC: {RNA_auc:.4f}")
            else:
                print(f"{group_name} {region_name} region RNAseq data empty")
        except (AttributeError, ValueError, IndexError) as e:
            print(f"{group_name} {region_name} region RNAseq data unavailable: {str(e)}")

        ATAC_auc = 0.0
        try:
            ATAC = output.atac.values[start:end]
            if ATAC.size > 0:
                if ATAC.ndim == 2:
                    if ATAC.shape[1] > 0:
                        ATAC = ATAC[:, 0]
                    else:
                        raise ValueError("ATAC data dimension incorrect")
                if np.any(np.isinf(ATAC)) or np.any(np.isnan(ATAC)):
                    ATAC = ATAC[~np.isinf(ATAC) & ~np.isnan(ATAC)]
                    if ATAC.size > 0:
                        if noise_filtering:
                            filtered_data = filter_noise(ATAC, threshold_multiplier)
                            filtered_data_array = np.full_like(ATAC, filtered_data)
                            ATAC_auc = np.trapezoid(filtered_data_array)
                        else:
                            ATAC_auc = np.trapezoid(ATAC)
                    else:
                        ATAC_auc = 0.0
                        print(f"{group_name} {region_name} region ATAC data all infinite or NaN")
                else:
                    if noise_filtering:
                        filtered_data = filter_noise(ATAC, threshold_multiplier)
                        filtered_data_array = np.full_like(ATAC, filtered_data)
                        ATAC_auc = np.trapezoid(filtered_data_array)
                    else:
                        ATAC_auc = np.trapezoid(ATAC)
                print(f"{group_name} {region_name} region ATAC AUC: {ATAC_auc:.4f}")
            else:
                print(f"{group_name} {region_name} region ATAC data empty")
        except (AttributeError, ValueError, IndexError) as e:
            print(f"{group_name} {region_name} region ATAC data unavailable: {str(e)}")
        return RNA_auc, ATAC_auc

    z1 = []
    z1_RNA = []
    z1_ATAC = []

    regions = [
        ("p1", p1_s, p1_e),
        ("CDS1", CDS1_s, CDS1_e),
        ("insulator", insulator_s, insulator_e),
        ("p2", p2_s, p2_e),
        ("CDS2", CDS2_s, CDS2_e)
    ]

    for region_name, start, end in regions:
        RNA_auc, ATAC_auc = analyze_region(group1_output, "Group1", region_name, start, end)
        z1.append(RNA_auc)
        z1.append(ATAC_auc)
        z1_RNA.append(RNA_auc)
        z1_ATAC.append(ATAC_auc)
    
    return group1_output, z1, z1_RNA, z1_ATAC


def analyze_empty_sequence(model, args, noise_filtering=True, threshold_multiplier=-1):
    full_sequence = args.CDS1 + args.insulator + args.p2 + args.CDS2
    padded_sequence = full_sequence.upper().center(dna_client.SEQUENCE_LENGTH_16KB, 'N')
    
    group1_output = None
    if args.ot:
        group1_output = model.predict_sequence(
            sequence=padded_sequence,
            requested_outputs=[dna_client.OutputType.RNA_SEQ],
            ontology_terms=[args.ot]
        )
    
    origin = (len(padded_sequence) - len(full_sequence)) // 2
    CDS1_s, CDS1_e = origin, origin + len(args.CDS1)
    
    group1_CDS1_RNA = 0.0
    try:
        RNAseq = group1_output.rna_seq.values[CDS1_s:CDS1_e]
        if RNAseq.ndim == 2:
            RNAseq = RNAseq[:, 0]
        if RNAseq.size > 0:
            RNAseq = RNAseq[~np.isinf(RNAseq) & ~np.isnan(RNAseq)]
            if RNAseq.size > 0:
                if noise_filtering:
                    filtered_data = filter_noise(RNAseq, threshold_multiplier)
                    filtered_data_array = np.full_like(RNAseq, filtered_data)
                    group1_CDS1_RNA = np.trapezoid(filtered_data_array)
                else:
                    group1_CDS1_RNA = np.trapezoid(RNAseq)
            else:
                group1_CDS1_RNA = 0.0
                print("Group1 CDS1 region RNAseq data all infinite or NaN")
        else:
            print("Group1 CDS1 region RNAseq data empty")
    except (AttributeError, ValueError, IndexError) as e:
        print(f"Group1 CDS1 region RNAseq data unavailable: {str(e)}")
    
    return group1_CDS1_RNA

def VirReport(z1_RNA, z1_ATAC, empty_parameter):
    try:
        group1_CDS1_RNA = z1_RNA[1]
        group1_CDS2_RNA = z1_RNA[4]
        group1_p1_ATAC = z1_ATAC[0]
        group1_p2_ATAC = z1_ATAC[3]

        second_parameter = 0.5*(group1_CDS1_RNA-empty_parameter)/(group1_CDS2_RNA)+0.5*(group1_p1_ATAC)/(group1_p2_ATAC)
        print(f"Virana result: {second_parameter:.4f}")
        return second_parameter
    except Exception as e:
        print(f"Calculation failed: {str(e)}")
        return None

empty_parameter = analyze_empty_sequence(model, args)
print(f"Empty parameter: {empty_parameter:.4f}")
group1_output, z1, z1_RNA, z1_ATAC = analyze_sequence(model, args)
VirReport(z1_RNA, z1_ATAC, empty_parameter)
