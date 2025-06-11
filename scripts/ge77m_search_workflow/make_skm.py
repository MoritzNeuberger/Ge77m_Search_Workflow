import os

def merge_lh5_files(input, output):
    print("Making skm...")
    print("input:", input)
    print("output:", output)    
    os.system("touch " + output)