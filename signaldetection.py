import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.signal import find_peaks
import csv
import time
import os

# Constants for pcm detection
DISTANCE_THRESHOLD = 20
POWER_THRESHOLD = 0.3
CENTRAL_DIFFERENCE_THRESHOLD = 1.1
MAX_PCM_SIGNALS = 5
KHZ_THRESHOLD = 100

# Constants for csv output
MAX_NON_PCM_SIGNALS = 20
OUTPUT_CSV_FILENAME = 'output.csv'

# for interactive plots
plt.switch_backend('TkAgg')

def read_dat_file(filename, num_bins):
    data = np.fromfile(filename, dtype=np.float32)
    num_vectors = data.size // num_bins
    data = data[:num_vectors * num_bins]
    reshaped_data = data.reshape((num_vectors, num_bins))
    return reshaped_data

\def read_csv_file(filename):
    with open(filename, mode='r') as file:
        reader = csv.reader(file)
    data = [list(map(float, row)) for row in reader]
    return np.array(data)

def create_empty_output_csv(csv_filename):
    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        headers = ['Graph', 'Noise Floor', 'Signal Detected', 'PCM Detected'] + \
                  [f'PCM Frequency {i+1}' for i in range(MAX_PCM_SIGNALS)] + \
                  [f'Subcarrier Bandwidth {i+1}' for i in range(MAX_PCM_SIGNALS)] + \
                  [f'PCM 99% BW Signal Power {i+1}' for i in range(MAX_PCM_SIGNALS)] + \
                  [f'PCM OBW {i+1}' for i in range(MAX_PCM_SIGNALS)] + \
                  [f'PCM Noise Power within OBW {i+1}' for i in range(MAX_PCM_SIGNALS)] + \
                  [f'PCM SNR {i+1}' for i in range(MAX_PCM_SIGNALS)] + \
                  [f'Central Frequency {i+1}' for i in range(MAX_NON_PCM_SIGNALS)] + \
                  [f'Central Power {i+1}' for i in range(MAX_NON_PCM_SIGNALS)] + \
                  [f'3dB Bandwidth {i+1}' for i in range(MAX_NON_PCM_SIGNALS)] + \
                  [f'99% BW Signal Power {i+1}' for i in range(MAX_NON_PCM_SIGNALS)] + \
                  [f'OBW {i+1}' for i in range(MAX_NON_PCM_SIGNALS)] + \
                  [f'Noise Power within OBW {i+1}' for i in range(MAX_NON_PCM_SIGNALS)] + \
                  [f'SNR {i+1}' for i in range(MAX_NON_PCM_SIGNALS)] + \
                  ['Noise Floor Density']
        writer.writerow(headers)
    return csv_filename

# moving average of values across time
def moving_average(data, window_size):
    cumulative_sum = np.cumsum(np.insert(data, 0, 0, axis=0), axis=0)
    averaged_data = (cumulative_sum[window_size:] - cumulative_sum[:-window_size]) / window_size
    return averaged_data

# noise floor (mean with margin)
def calculate_noise_floor(data, margin):
    mean_noise = np.mean(data)
    noise_floor = mean_noise + margin
    return noise_floor

# noise floor density
def calculate_noise_floor_density(noise_floor, recording_bandwidth, num_bins):
    hz_per_bin = recording_bandwidth / num_bins
    noise_floor_linear = 10 ** (noise_floor / 10)
    noise_floor_density = noise_floor_linear / hz_per_bin
    noise_floor_density_db = 10 * np.log10(noise_floor_density)
    return noise_floor_density_db

def bin_to_frequency(bin_index, central_bandwidth_frequency, num_bins, hz_per_bin):
    central_bin = num_bins // 2
    frequency = central_bandwidth_frequency + (bin_index - central_bin) * hz_per_bin
    return frequency

# 99% BW Real Signal Power for non-PCM signals (uses 3*3db BW)
def calculate_99_bw_signal_power(data, central_freq, bandwidth):
    span = int(bandwidth * 3)
    start_idx = max(0, central_freq - span // 2)
    end_idx = min(len(data), central_freq + span // 2)
    span_data = data[start_idx:end_idx]

    linear_sum = np.sum(10 ** (span_data / 10))

    real_signal_power = 10 * np.log10(linear_sum * 0.99)

    return real_signal_power, linear_sum

# 99% BW Real Signal Power for PCM signals (uses 3*subcarrier BW)
def calculate_pcm_99_bw_signal_power(data, central_freq, subcarrier_bandwidth):
    span = int(subcarrier_bandwidth * 3)
    start_idx = max(0, central_freq - span // 2)
    end_idx = min(len(data), central_freq + span // 2)
    span_data = data[start_idx:end_idx]

    linear_sum = np.sum(10 ** (span_data / 10))

    real_signal_power = 10 * np.log10(linear_sum * 0.99)

    return real_signal_power, linear_sum

# OBW for non-PCM signals
def calculate_obw_and_noise_power(data, central_freq, linear_sum_99bw, noise_floor_density, recording_bandwidth, num_bins):
    hz_per_bin = recording_bandwidth / num_bins
    total_linear_sum = 0
    left_idx = central_freq
    right_idx = central_freq
    while total_linear_sum < linear_sum_99bw * 0.99:
        if left_idx > 0:
            left_idx -= 1
            total_linear_sum += 10 ** (data[left_idx] / 10)
        if right_idx < len(data) - 1:
            right_idx += 1
            total_linear_sum += 10 ** (data[right_idx] / 10)

    obw_bins = right_idx - left_idx
    obw_hz = obw_bins * hz_per_bin

    noise_power_linear = 10 ** (noise_floor_density / 10) * obw_hz
    noise_power_dbm = 10 * np.log10(noise_power_linear)

    return obw_hz, noise_power_dbm

# OBW for PCM signals
def calculate_pcm_obw_and_noise_power(data, central_freq, linear_sum_99bw, noise_floor_density, subcarrier_bandwidth, recording_bandwidth, num_bins):
    hz_per_bin = recording_bandwidth / num_bins
    total_linear_sum = 0
    left_idx = central_freq
    right_idx = central_freq
    while total_linear_sum < linear_sum_99bw * 0.99:
        if left_idx > 0:
            left_idx -= 1
            total_linear_sum += 10 ** (data[left_idx] / 10)
        if right_idx < len(data) - 1:
            right_idx += 1
            total_linear_sum += 10 ** (data[right_idx] / 10)

    obw_bins = right_idx - left_idx
    obw_hz = obw_bins * hz_per_bin

    noise_power_linear = 10 ** (noise_floor_density / 10) * obw_hz
    noise_power_dbm = 10 * np.log10(noise_power_linear)

    return obw_hz, noise_power_dbm

def calculate_snr(signal_power_dbm, noise_power_dbm):
    snr_db = signal_power_dbm - noise_power_dbm
    return snr_db

# compare peaks to get subcarrier frequency
def find_subcarrier_frequency(peaks, peak_powers, central_freq, power_interval=0.3):
    subcarrier_pairs = []
    peak_diff_matrix = np.abs(np.expand_dims(peak_powers, axis=1) - peak_powers)
    close_peaks = np.where(peak_diff_matrix <= power_interval)
    subcarrier_pairs = [(i, j) for i, j in zip(*close_peaks) if i < j]

    closest_subcarrier_bandwidth = None
    min_distance_to_central = float('inf')

    for pair in subcarrier_pairs:
        freq_diff = abs(peaks[pair[0]] - peaks[pair[1]])
        center_diff = min(abs(peaks[pair[0]] - central_freq), abs(peaks[pair[1]] - central_freq))
        if freq_diff <= 100 and center_diff < min_distance_to_central:
            min_distance_to_central = center_diff
            closest_subcarrier_bandwidth = freq_diff

    return closest_subcarrier_bandwidth

# remove SDR artifact by finding the center of the signal using median of peaks
# @ALB:
def remove_SDR(peaks, peak_powers, lo_offset):
    if len(peaks) == 0:
        return peaks, peak_powers

    center_freq = int(np.median(peaks))
    start_idx = max(0, center_freq - lo_offset)
    end_idx = min(peaks[-1], center_freq + lo_offset)

    center_peaks = [i for i in range(len(peaks)) if start_idx <= peaks[i] <= end_idx]
    if not center_peaks:
        return peaks, peak_powers
    # remove sdr artifact
    sdr_peak_index = center_peaks[np.argmax(peak_powers[center_peaks])]
    peaks_filtered = np.delete(peaks, sdr_peak_index)
    peak_powers_filtered = np.delete(peak_powers, sdr_peak_index)

    return peaks_filtered, peak_powers_filtered

# calculate 3db bandwidth
def calculate_3db_bandwidth(data, peaks, peak_powers):
    if len(peak_powers) == 0:
        return None, None, None, None, None

    highest_peak_index = np.argmax(peak_powers)
    central_freq = peaks[highest_peak_index]
    central_power = peak_powers[highest_peak_index]

    # calculate 3 dB bandwidth
    three_db_level = central_power - 3

    left_base = np.where(data[:central_freq] < three_db_level)[0]
    right_base = np.where(data[central_freq:] < three_db_level)[0]

    left_base = left_base[-1] if len(left_base) > 0 else 0 # TODO: @ALB: is it an issue if len = 0?
    right_base = right_base[0] + central_freq if len(right_base) > 0 else len(data) - 1

    bandwidth = right_base - left_base

    return central_freq, central_power, bandwidth, left_base, right_base

# detects multiple PCM signals in the same timestamp, MAX_PCM_SIGNALS can be removed if needed.
def detect_pcm_psk_pm(peaks, peak_powers, DISTANCE_THRESHOLD, POWER_THRESHOLD, CENTRAL_DIFFERENCE_THRESHOLD, MAX_PCM_SIGNALS):
    pcm_signals = []

    while len(pcm_signals) < MAX_PCM_SIGNALS:
        if len(peaks) == 0:
            break

        central_peak_idx = np.argmax(peak_powers)
        central_peak = peaks[central_peak_idx]
        central_power = peak_powers[central_peak_idx]

        left_subcarrier = None
        right_subcarrier = None

        for j in range(len(peaks)):
            if j == central_peak_idx:
                continue
            freq_diff = abs(peaks[j] - central_peak)

            if freq_diff <= DISTANCE_THRESHOLD:
                if peaks[j] < central_peak:
                    left_subcarrier = j
                elif peaks[j] > central_peak:
                    right_subcarrier = j

        if left_subcarrier is not None and right_subcarrier is not None:
            left_power = peak_powers[left_subcarrier]
            right_power = peak_powers[right_subcarrier]
            left_distance = abs(peaks[left_subcarrier] - central_peak)
            right_distance = abs(peaks[right_subcarrier] - central_peak)

            if abs(left_power - right_power) <= POWER_THRESHOLD:
                if left_power < central_power and right_power < central_power:
                    if abs(left_distance - right_distance) <= CENTRAL_DIFFERENCE_THRESHOLD:
                        pcm_signals.append((central_peak, central_power))
                        peaks = np.delete(peaks, central_peak_idx)
                        peak_powers = np.delete(peak_powers, central_peak_idx)
                        continue

        break

    return pcm_signals

# remove peaks around PCM signals, tolerance is kHz
def remove_peaks_around_pcm(peaks, peak_powers, pcm_frequencies, tolerance=KHZ_THRESHOLD):
    to_remove = []
    for pcm_freq in pcm_frequencies:
        to_remove += [i for i, peak in enumerate(peaks) if abs(peak - pcm_freq) <= tolerance]
    peaks_filtered = np.delete(peaks, to_remove)
    peak_powers_filtered = np.delete(peak_powers, to_remove)
    return peaks_filtered, peak_powers_filtered

# detect remaining non pcm signals
def detect_additional_signals(data, peaks, peak_powers, peak_threshold, MAX_NON_PCM_SIGNALS):
    additional_signals = []
    while len(peaks) > 0:
        peak_idx = np.argmax(peak_powers)
        central_freq = peaks[peak_idx]
        central_power = peak_powers[peak_idx]

        left = central_freq
        right = central_freq

        # group peaks to left and right until signal falls 1 dB below peak threshold
        while left > 0 and data[left] >= (peak_threshold - 1):
            left -= 1
        while right < len(data) - 1 and data[right] >= (peak_threshold - 1):
            right += 1

        left_base = max(left, 0)
        right_base = min(right, len(data) - 1)

        signal_peaks = [i for i in range(len(peaks)) if left_base <= peaks[i] <= right_base]

        bandwidth = right_base - left_base

        peaks = np.delete(peaks, signal_peaks)
        peak_powers = np.delete(peak_powers, signal_peaks)

        additional_signals.append((central_freq, central_power, bandwidth))
        if len(additional_signals) >= MAX_NON_PCM_SIGNALS:
            break

    return additional_signals

def main():
    parser = argparse.ArgumentParser(description="Signal Detection")
    parser.add_argument('filename', type=str, help="Path to the file")
    parser.add_argument('file_type', type=str, choices=['dat', 'csv'], help="Type of the file: 'dat' or 'csv'")
    parser.add_argument('num_bins', type=int, help="Number of bins (required for .dat files)")
    parser.add_argument('num_iterations', help="Number of iterations or 'end' to process until end")
    parser.add_argument('num_moving_avg', type=int, help="Number for moving average (rec. 100 for short signals or 1000 for longer)") # TODO: @ALB: what unit is this? Seconds? [number of iterations, out of num_iterations]
    parser.add_argument('num_peak_baseline', type=int, help="How much dB above noise floor to detect peaks (rec. 3)")
    parser.add_argument('lo_offset', type=int, help="Enter lo_offset for SDR artifact removal (rec. 50) (0 if no artifact)") # @ALB: unit = number of bins
    parser.add_argument('recording_bandwidth', type=float, help="Recording bandwidth in Hz")
    parser.add_argument('central_bandwidth_frequency', type=float, help="Central frequency of the recording bandwidth in Hz")
    parser.add_argument('plot_graph', type=str, help="Plot the graph (y/n) or a specific number")
    parser.add_argument('output_type', type=str, choices=['gui', 'img'], help="Graph output: 'gui' for gui or 'img' to save as image")

    args = parser.parse_args()

    filename = args.filename
    file_type = args.file_type
    num_bins = args.num_bins
    num_iterations = int(args.num_iterations) if args.num_iterations != 'end' else 'end'
    num_moving_avg = args.num_moving_avg
    num_peak_baseline = args.num_peak_baseline
    lo_offset = args.lo_offset
    recording_bandwidth = args.recording_bandwidth
    central_bandwidth_frequency = args.central_bandwidth_frequency
    plot_graph = args.plot_graph.lower()
    output_type = args.output_type.lower()

    if file_type == 'dat':
        data = read_dat_file(filename, num_bins)
    elif file_type == 'csv':
        data = read_csv_file(filename)

    averaged_data = moving_average(data, num_moving_avg)

    csv_filename = create_empty_output_csv(OUTPUT_CSV_FILENAME)

    num_vectors = averaged_data.shape[0]
    end_iteration = num_vectors if num_iterations == 'end' else min(num_iterations, num_vectors)

    hz_per_bin = recording_bandwidth / num_bins

    start_time = time.time()

    def plot_graph_func(i):
        plt.figure(figsize=(12, 10))

        plt.subplot(3, 1, 1)
        plt.plot(data[i], label=f'Original Data ({i+1})')
        plt.xlabel('Frequency')
        plt.ylabel('Power')
        plt.title(f'Original Data ({i+1})')
        plt.legend()

        plt.subplot(3, 1, 2)
        plt.plot(averaged_slice, label=f'Moving Average Data ({i+1}-{i+num_moving_avg})')
        if len(peaks) > 0:
            plt.plot(peaks_filtered, averaged_slice[peaks_filtered], "x", label='Peaks')
        plt.axhline(y=noise_floor, color='r', linestyle='--', label='Noise Floor')

        if central_freq is not None:
            plt.plot(central_freq, central_power, "ro", label='Central Frequency')
            plt.text(central_freq, central_power, f'({central_freq}, {central_power})', fontsize=9, ha='right')
            plt.hlines(y=central_power - 3, xmin=left_base, xmax=right_base, color='g', linestyle='--', label='3 dB Bandwidth')

        plt.xlabel('Frequency')
        plt.ylabel('Power')
        plt.title(f'Moving Average Data ({i+1}-{i+num_moving_avg})')
        plt.legend()

        plt.tight_layout()

        if output_type == 'gui':
            plt.show()
        elif output_type == 'img':
            img_filename = f"graph_{i+1}.png"
            plt.savefig(os.path.join(os.getcwd(), img_filename))
            plt.close()

    for i in range(end_iteration):
        averaged_slice = averaged_data[i]

        noise_floor = calculate_noise_floor(averaged_slice, margin=0)
        noise_floor_density = calculate_noise_floor_density(noise_floor, recording_bandwidth, num_bins)
        peak_baseline = calculate_noise_floor(averaged_slice, margin=num_peak_baseline)

        peaks, properties = find_peaks(averaged_slice, height=peak_baseline)
        peak_powers = averaged_slice[peaks]

        if len(peaks) > 0:
            # 1. Remove SDR artifacts if required
            if lo_offset != 0:
                peaks_filtered, peak_powers_filtered = remove_SDR(peaks, peak_powers, lo_offset)
            else:
                peaks_filtered, peak_powers_filtered = peaks, peak_powers

            # 2. Calculate the 3DB bandwidth
            central_freq, central_power, bandwidth, left_base, right_base = calculate_3db_bandwidth(averaged_slice, peaks_filtered, peak_powers_filtered)

            # 3. Extract all PCM signals. If any are found, process them.
            if central_freq is not None:
                pcm_signals = detect_pcm_psk_pm(peaks_filtered, peak_powers_filtered, DISTANCE_THRESHOLD, POWER_THRESHOLD, CENTRAL_DIFFERENCE_THRESHOLD, MAX_PCM_SIGNALS)
                graph_label = f'{i+1}-{i+num_moving_avg}'
                central_freq_detected = 'Yes'
                pcm_detected = len(pcm_signals)
                pcm_frequencies = [signal[0] for signal in pcm_signals]

                # search for subcarrier bandwidth within KHZ_THRESHOLD kHz of PCM signal closest to central freq
                subcarrier_bandwidths = []
                pcm_99bw_powers = []
                pcm_obws = []
                pcm_noise_powers_within_obw = []
                pcm_snrs = []
                for pcm_freq in pcm_frequencies:
                    pcm_peaks = [peak for peak in peaks_filtered if abs(peak - pcm_freq) <= KHZ_THRESHOLD]
                    pcm_peak_powers = [peak_powers_filtered[i] for i, peak in enumerate(peaks_filtered) if abs(peak - pcm_freq) <= KHZ_THRESHOLD]
                    subcarrier_bandwidth = find_subcarrier_frequency(np.array(pcm_peaks), np.array(pcm_peak_powers), pcm_freq)
                    if subcarrier_bandwidth:
                        subcarrier_bandwidths.append(subcarrier_bandwidth / 2)
                        pcm_99bw_power, pcm_linear_sum_99bw = calculate_pcm_99_bw_signal_power(averaged_slice, pcm_freq, subcarrier_bandwidth / 2)
                        pcm_99bw_powers.append(pcm_99bw_power)
                        pcm_obw, pcm_noise_power_within_obw = calculate_pcm_obw_and_noise_power(
                            averaged_slice, pcm_freq, pcm_linear_sum_99bw, noise_floor_density, subcarrier_bandwidth / 2, recording_bandwidth, num_bins
                        )
                        pcm_obws.append(pcm_obw)
                        pcm_noise_powers_within_obw.append(pcm_noise_power_within_obw)
                        pcm_snr = calculate_snr(pcm_99bw_power, pcm_noise_power_within_obw)
                        pcm_snrs.append(pcm_snr)

                peaks_filtered, peak_powers_filtered = remove_peaks_around_pcm(peaks_filtered, peak_powers_filtered, pcm_frequencies)
                additional_signals = detect_additional_signals(averaged_slice, peaks_filtered, peak_powers_filtered, noise_floor + 3, MAX_NON_PCM_SIGNALS)

                additional_freqs = [sig[0] for sig in additional_signals]
                additional_powers = [sig[1] for sig in additional_signals]
                additional_bandwidths = [sig[2] for sig in additional_signals]
                additional_99bw_powers = []
                additional_obws = []
                additional_noise_powers_within_obw = []
                additional_snrs = []

                for sig in additional_signals:
                    additional_99bw_power, additional_linear_sum_99bw = calculate_99_bw_signal_power(averaged_slice, sig[0], sig[2])
                    additional_99bw_powers.append(additional_99bw_power)
                    additional_obw, additional_noise_power_within_obw = calculate_obw_and_noise_power(
                        averaged_slice, sig[0], additional_linear_sum_99bw, noise_floor_density, recording_bandwidth, num_bins
                    )
                    additional_obws.append(additional_obw)
                    additional_noise_powers_within_obw.append(additional_noise_power_within_obw)
                    additional_snr = calculate_snr(additional_99bw_power, additional_noise_power_within_obw)
                    additional_snrs.append(additional_snr)

                # convert from bins to frequency
                pcm_frequencies_hz = [bin_to_frequency(f, central_bandwidth_frequency, num_bins, hz_per_bin) for f in pcm_frequencies]
                subcarrier_bandwidths_hz = [bw * hz_per_bin for bw in subcarrier_bandwidths]
                pcm_obws_hz = [obw for obw in pcm_obws]

                additional_freqs_hz = [bin_to_frequency(f, central_bandwidth_frequency, num_bins, hz_per_bin) for f in additional_freqs]
                additional_bandwidths_hz = [bw * hz_per_bin for bw in additional_bandwidths]
                additional_obws_hz = [obw for obw in additional_obws]

                row = [graph_label, noise_floor, central_freq_detected, pcm_detected]
                row.extend(pcm_frequencies_hz + ['N/A'] * (MAX_PCM_SIGNALS - len(pcm_frequencies_hz)))
                row.extend(subcarrier_bandwidths_hz + ['N/A'] * (MAX_PCM_SIGNALS - len(subcarrier_bandwidths_hz)))
                row.extend(pcm_99bw_powers + ['N/A'] * (MAX_PCM_SIGNALS - len(pcm_99bw_powers)))
                row.extend(pcm_obws_hz + ['N/A'] * (MAX_PCM_SIGNALS - len(pcm_obws_hz)))
                row.extend(pcm_noise_powers_within_obw + ['N/A'] * (MAX_PCM_SIGNALS - len(pcm_noise_powers_within_obw)))
                row.extend(pcm_snrs + ['N/A'] * (MAX_PCM_SIGNALS - len(pcm_snrs)))

                row.extend(additional_freqs_hz + ['N/A'] * (MAX_NON_PCM_SIGNALS - len(additional_freqs_hz)))
                row.extend(additional_powers + ['N/A'] * (MAX_NON_PCM_SIGNALS - len(additional_powers)))
                row.extend(additional_bandwidths_hz + ['N/A'] * (MAX_NON_PCM_SIGNALS - len(additional_bandwidths_hz)))
                row.extend(additional_99bw_powers + ['N/A'] * (MAX_NON_PCM_SIGNALS - len(additional_99bw_powers)))
                row.extend(additional_obws_hz + ['N/A'] * (MAX_NON_PCM_SIGNALS - len(additional_obws_hz)))
                row.extend(additional_noise_powers_within_obw + ['N/A'] * (MAX_NON_PCM_SIGNALS - len(additional_noise_powers_within_obw)))
                row.extend(additional_snrs + ['N/A'] * (MAX_NON_PCM_SIGNALS - len(additional_snrs)))

                row.append(noise_floor_density)

                with open(csv_filename, mode='a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(row)

            else:
                graph_label = f'{i+1}-{i+num_moving_avg}'
                row = [graph_label, noise_floor, 'No', 'N/A'] + ['N/A'] * MAX_PCM_SIGNALS + ['N/A'] * MAX_PCM_SIGNALS + ['N/A'] * (MAX_PCM_SIGNALS * 5 + MAX_NON_PCM_SIGNALS * 7) + [noise_floor_density]

                with open(csv_filename, mode='a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(row)

        else:
            # no signal edge case
            graph_label = f'{i+1}-{i+num_moving_avg}'
            row = [graph_label, noise_floor, 'No', 'N/A'] + ['N/A'] * MAX_PCM_SIGNALS + ['N/A'] * MAX_PCM_SIGNALS + ['N/A'] * (MAX_PCM_SIGNALS * 5 + MAX_NON_PCM_SIGNALS * 7) + [noise_floor_density]

            with open(csv_filename, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row)

        if plot_graph == 'y':
            plot_graph_func(i)
        elif plot_graph.isdigit() and int(plot_graph) == i:
            plot_graph_func(i)

    print(f"Time taken to process whole file: {time.time() - start_time:.3f} seconds")
    # plt.pause(378)

# interactive graphs
plt.ion()

if __name__ == '__main__':
    main()
