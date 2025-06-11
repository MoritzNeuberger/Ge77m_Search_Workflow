import os

def process_delayed_coinc(input, output):
    print("Processing delayed_coinc...")
    print("input:", len(input))
    print("output:", len(output))    
    os.system(f"touch {output}")