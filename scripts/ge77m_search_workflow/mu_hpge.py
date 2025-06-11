import os 

def process_mu_hpge_coinc(input, output):
    print("Processing mu_hpge_coinc...")
    print("input:", len(input))
    print("output:", len(output))    
    os.system(f"touch {output}")