import pyhepmc
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path

def find_tail_cutoff(arr, percentage):
    """
    Sorts an array in increasing values and finds the index where the sum 
    of the remaining entries (from that index to the end) is less than a 
    fixed percentage of the total sum.
    
    Parameters:
    arr (array-like): The input array (assumed to be non-negative numbers)
    percentage (float): The threshold percentage (e.g., 0.15 for 15%)
    
    Returns:
    tuple: (sorted_array, cutoff_index, tail_sum)
    """
    arr = np.asarray(arr)
    
    # 1. Sort the array in increasing (ascending) order
    sorted_arr = np.sort(arr)
    
    # 2. Calculate the total sum and the target threshold
    total_sum = np.sum(sorted_arr)
    threshold = total_sum * percentage
    
    # 3. Calculate suffix sums (sum of elements from index i to the end)
    # Reversing the array, taking the cumulative sum, and reversing it back 
    # gives us an array where index i is the sum of sorted_arr[i:]
    suffix_sums = np.cumsum(sorted_arr[::-1])[::-1]
    
    # 4. Find the first index where the tail sum is strictly less than the threshold
    condition_met = suffix_sums < threshold
    
    if not np.any(condition_met):
        # If no tail sum is small enough, the cutoff is the end of the array (empty tail)
        cutoff_index = len(sorted_arr)
    else:
        # np.argmax returns the first index where the boolean condition is True
        cutoff_index = np.argmax(condition_met)
        
    return sorted_arr[cutoff_index]

parser = argparse.ArgumentParser(description="A weight plotter")
parser.add_argument("filename", type=str, help="Event file to plot")
parser.add_argument("--cutoff", type=float, default=0.1, help="Percent cutoff")
args = parser.parse_args()

stem = Path(args.filename).stem

weights = []
with pyhepmc.open(args.filename) as f:
    for event in f:
        weights.append(event.weight("CV"))

wgt = np.array(weights)
mean = np.mean(wgt)
target_percentage = args.cutoff
max_wgt = find_tail_cutoff(wgt, target_percentage)
percentile = np.percentile(wgt, 99)
eff = mean/max_wgt*100

fig, ax = plt.subplots()
plt.hist(np.log10(wgt), bins=500, histtype='step')
plt.axvline(x=np.log10(mean), color='red', ls='--')
plt.axvline(x=np.log10(max_wgt), color='green', ls='--')
plt.axvline(x=np.log10(percentile), color='blue', ls='--')
plt.text(0.05, 0.95, f"Eff = {eff:.2f}%", transform=ax.transAxes)
plt.ylabel('# Events')
plt.xlabel(r'$log_{10}(wgt)$')
plt.yscale('log')
plt.savefig(f'weights_{stem}.png')
plt.show()
