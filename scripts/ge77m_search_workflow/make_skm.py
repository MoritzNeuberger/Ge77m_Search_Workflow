import os

def merge_lh5_files(input, output):
    print("Making skm...")
    print("input:", len(input))
    print("output:", len(output))    
    os.system("touch " + output)