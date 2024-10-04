# Summer 2024 Infostellar Internship Project - Signal Detection Automation 

## Brian Li

### signaldetection.py

**Description:**
Run using the following command line arguments:
```
python signaldetection.py filename [file type (ex: dat/csv)] [number of bins] [number of iterations or 'end' to process whole file] [number for moving average] [how many db above noise floor to detect peaks] [lo_offset (0 if no SDR artifact)] [graph (y to plot all, n to not plot anything, or the interation number of the plot you want to graph)]
```
- Example: 
```
python signaldetection.py GL20000.dat dat 20000 end 1000 3 50 n
```

**Outputs:**
The code will detect and output the following to a CSV file:
- Graph Number with Moving Average
- Noise Floor
- Signal Detected
- Number of PCM-PSK-PM Signals Detected
- PCM-PSK-PM Signal Frequencies
- PCM-PSK-PM Subcarrier Bandwidths
- Signal Central Frequency 
- Signal Central Power
- Signal 3dB Bandwidth
