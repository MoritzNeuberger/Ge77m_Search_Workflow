import os

def process_delayed_coinc(input, output):
    print("Processing delayed_coinc...")
    print("input:", input)
    print("output:", output)
    os.system(f"touch {output}")