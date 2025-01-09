#!/usr/bin/env python3

import sys
import math
import re
from typing import Tuple, Optional
import argparse


EXPECTED = {
    "run_electron_12C_SF_500":  [1.75967e5, 7.88041e1],
    "run_electron_12C_SF_1000": [3.68880e4, 1.80475e1],
    "run_electron_12C_SF_2000": [5.61866e3, 2.45236e0],
    "run_nue_12C_SF_500":  [4.05979e-5, 1.90439e-8],
    "run_nue_12C_SF_1000": [5.72067e-5, 2.40873e-8],
    "run_nue_12C_SF_2000": [6.10140e-5, 2.73874e-8],
    "run_nuebar_12C_SF_500":  [1.06739e-5, 5.26641e-9],
    "run_nuebar_12C_SF_1000": [2.13327e-5, 1.00923e-8],
    "run_nuebar_12C_SF_2000": [3.47148e-5, 1.69763e-8],
}


def extract_xsec(filename: str) -> Optional[Tuple[float, float]]:
    """
    Extract cross section and error from the log file.

    Args:
        filename: Path to the log file

    Returns:
        Tuple of (cross_section, error) if found, None otherwise
    """
    try:
        with open(filename, 'r') as f:
            content = f.read()

        # Find the line containing "Total xsec:"
        number_re = r"([-+]?\d*\.?\d+([eE][-+]?\d+)?)"
        pattern = rf"Total xsec:\s+{number_re}\s+\+/-\s+{number_re}"
        match = re.search(pattern, content)

        if not match:
            print("Error: Could not find cross section values in log file")
            return None

        return float(match.group(1)), float(match.group(3))

    except FileNotFoundError:
        print(f"Error: File {filename} not found")
        return None
    except Exception as e:
        print(f"Error: Failed to process file: {str(e)}")
        return None


def calculate_significance(measured: float, measured_err: float,
                           expected: float, expected_err: float) -> float:
    """
    Calculate the statistical significance between measured and expected values.

    Args:
        measured: Measured value
        measured_err: Error on measured value
        expected: Expected value
        expected_err: Error on expected value

    Returns:
        Number of standard deviations between the measurements
    """
    diff = abs(measured - expected)
    combined_err = math.sqrt(measured_err**2 + expected_err**2)
    return diff / combined_err


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Validate cross section measurements.')
    parser.add_argument('logfile', help='Path to the log file')
    parser.add_argument('name', type=str,
                        help='Name of the example ran')
    parser.add_argument('--threshold', type=float, default=3.0,
                        help='Significance threshold in sigma (default: 3.0)')

    # Parse arguments
    args = parser.parse_args()

    # Extract values from log file
    result = extract_xsec(args.logfile)
    if result is None:
        sys.exit(1)

    measured_value, measured_error = result
    expected_value, expected_error = EXPECTED[args.name]

    # Calculate significance
    n_sigma = calculate_significance(measured_value, measured_error,
                                     expected_value, expected_error)

    # Print results
    print(f"\nMeasured Cross Section: {measured_value:.6e} ± {measured_error:.6e}")
    print(f"Expected Cross Section: {expected_value:.6e} ± {expected_error:.6e}")
    print(f"Statistical Significance: {n_sigma:.3f} sigma")

    # Evaluate compatibility
    if n_sigma < args.threshold:
        print(f"\nPASS: Values are compatible (less than {args.threshold} sigma difference)")
        sys.exit(0)
    else:
        print(f"\nFAIL: Values differ by {n_sigma:.3f} sigma (threshold: {args.threshold} sigma)")
        sys.exit(1)


if __name__ == "__main__":
    main()
