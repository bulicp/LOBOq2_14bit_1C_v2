set_units -time ns -resistance kOhm -capacitance pF -voltage V -current mA
create_clock [get_ports clk] -period 10  -waveform {0 5}
set_load 0.01 [get_ports p]
