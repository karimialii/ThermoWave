# -*- coding: utf-8 -*-

from tespy.networks import Network
from tespy.components import (
    Sink, Source, Turbine, HeatExchanger, CycleCloser, Compressor,
    Motor, PowerBus, PowerSource
)
from tespy.connections import Connection, PowerConnection
import pandas as pd
import numpy as np
from exerpy import ExergyAnalysis
# from exerpy.sankey import create_sankey

import json


# (format string, column header with unit) — all values passed in display units
_COL_FMT = {
    'm':       ('{:.3f}', 'm in kg/s'),
    'T':       ('{:.1f}', 'T in °C'),
    'p':       ('{:.2f}', 'p in bar'),
    'h':       ('{:.1f}', 'h in kJ/kg'),
    'e_T':     ('{:.1f}', 'e_T in kJ/kg'),
    'e_M':     ('{:.1f}', 'e_M in kJ/kg'),
    'e_PH':    ('{:.1f}', 'e_PH in kJ/kg'),
    'E_T':     ('{:.2f}', 'E_T in kW'),
    'E_M':     ('{:.2f}', 'E_M in kW'),
    'E_PH':    ('{:.2f}', 'E_PH in kW'),
    'E_F':     ('{:.2f}', 'E_F in kW'),
    'E_P':     ('{:.2f}', 'E_P in kW'),
    'E_D':     ('{:.2f}', 'E_D in kW'),
    'E_L':     ('{:.2f}', 'E_L in kW'),
    'epsilon': ('{:.1f}', 'ε in %'),
    'y_Dk':    ('{:.1f}', 'y_Dk in %'),
    'y*_Dk':   ('{:.1f}', 'y*_Dk in %'),
}


def result_to_markdown(df, filename, prefix=''):
    df = df.copy()
    rename = {}
    for col in df.columns:
        fmt, header = _COL_FMT[col]
        if prefix == 'δ ':
            df[col] = (df[col] * 100).apply(fmt.format)
            rename[col] = prefix + col + ' in %'
        else:
            df[col] = df[col].apply(fmt.format)
            rename[col] = prefix + header
    df.rename(columns=rename, inplace=True)
    df.to_markdown(filename, disable_numparse=True,
                   colalign=['left'] + ['right'] * len(df.columns))


# specification of ambient state
pamb = 1
Tamb = 25

# ambient state in SI units for exerpy
pamb_Pa = pamb * 1e5
Tamb_K = Tamb + 273.15

# setting up network
nw = Network()
nw.set_attr(
    T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s',
    s_unit="kJ / kgK"
)

# components definition
water_in = Source('Water source')
water_out = Sink('Water sink')

air_in = Source('Air source')
air_out = Sink('Air sink')

closer = CycleCloser('Cycle closer')

cp = Compressor('Compressor')
turb = Turbine('Turbine')

cold = HeatExchanger('Cooling heat exchanger')
hot = HeatExchanger('Heat sink heat exchanger')

# connections definition
# refrigerant cycle
c0 = Connection(cold, 'out2', closer, 'in1', label='00')
c1 = Connection(closer, 'out1', cp, 'in1', label='01')
c2 = Connection(cp, 'out1', hot, 'in1', label='02')
c3 = Connection(hot, 'out1', turb, 'in1', label='03')
c4 = Connection(turb, 'out1', cold, 'in2', label='04')

c11 = Connection(air_in, 'out1', cold, 'in1', label='11')
c12 = Connection(cold, 'out1', air_out, 'in1', label='12')

c21 = Connection(water_in, 'out1', hot, 'in2', label='21')
c22 = Connection(hot, 'out2', water_out, 'in1', label='22')

nw.add_conns(c0, c1, c2, c3, c4, c11, c12, c21, c22)

# power network:
# turbine and compressor sit on the same shaft; the net power deficit is
# supplied through an inverter/motor, so losses apply only to the difference
grid = PowerSource('grid')
shaft_bus = PowerBus('shaft', num_in=2, num_out=1)
inverter = Motor('Inverter')

e_grid = PowerConnection(grid, 'power', inverter, 'power_in', label='e_grid')
e_inv = PowerConnection(inverter, 'power_out', shaft_bus, 'power_in1', label='e_inv')
e_turb = PowerConnection(turb, 'power', shaft_bus, 'power_in2', label='e_turb')
e_comp = PowerConnection(shaft_bus, 'power_out1', cp, 'power', label='e_comp')

nw.add_conns(e_grid, e_inv, e_turb, e_comp)

# connection parameters
c0.set_attr(T=-30, p=1, fluid={'Air': 1, 'water': 0})
c2.set_attr(p=5.25)
c3.set_attr(p=5, T=35)
c4.set_attr(p=1.05)

c11.set_attr(fluid={'Air': 1, 'water': 0}, T=-10, p=1)
c12.set_attr(p=1, T=-20)

c21.set_attr(fluid={'Air': 0, 'water': 1}, T=25, p=1.5)
c22.set_attr(p=1.5, T=40)

# component parameters
turb.set_attr(eta_s=0.8)
cp.set_attr(eta_s=0.8)
cold.set_attr(Q=-100e3)

# inverter/motor efficiency applied only to the net power difference
inverter.set_attr(eta=0.9)

nw.solve(mode='design')
nw.print_results()

# carry out exergy analysis
ean = ExergyAnalysis.from_tespy(nw, Tamb_K, pamb_Pa, split_physical_exergy=True)

fuel = {"inputs": ["e_grid"], "outputs": []}
product = {"inputs": ["12"], "outputs": ["11"]}
loss = {"inputs": ["22"], "outputs": ["21"]}

ean.analyse(E_F=fuel, E_P=product, E_L=loss)
df_components, df_connections, df_power = ean.exergy_results()

# generate Grassmann diagram
ean.export_to_json('refrigeration_result.json')
# with open('refrigeration_result.json') as f:
#     sankey_data = json.load(f)
# fig = create_sankey(sankey_data, mode="E")
# fig.show()

# validation (connections)

df_original_data = pd.read_csv(
    'connection_validation.csv', sep=';', decimal=',', index_col='label'
)

df_tespy = df_connections.set_index('Connection').rename(columns={
    'm [kg/s]': 'm', 'T [°C]': 'T', 'p [bar]': 'p', 'h [kJ/kg]': 'h',
    'e^T [kJ/kg]': 'e_T', 'e^M [kJ/kg]': 'e_M',
})
df_tespy.index = pd.to_numeric(df_tespy.index, errors='coerce')

# zero point of enthalpy differs from original data; normalise using first index
air_idx = [1, 2, 3, 4, 11, 12]
water_idx = [21, 22]
df_tespy.loc[air_idx, 'h'] -= df_tespy.loc[air_idx[0], 'h']
df_tespy.loc[water_idx, 'h'] -= df_tespy.loc[water_idx[0], 'h']
df_original_data.loc[air_idx, 'h'] -= df_original_data.loc[air_idx[0], 'h']
df_original_data.loc[water_idx, 'h'] -= df_original_data.loc[water_idx[0], 'h']

idx = np.intersect1d(df_tespy.index, df_original_data.index)
df_tespy = df_tespy.loc[idx, df_original_data.columns]

df_diff_abs = df_tespy - df_original_data
df_diff_rel = (df_tespy - df_original_data) / df_original_data

result_to_markdown(df_diff_abs, 'connections_delta_absolute', 'Δ ')
result_to_markdown(df_diff_rel, 'connections_delta_relative', 'δ ')

# validation (components)

df_original_data = pd.read_csv(
    'component_validation.csv', sep=';', decimal=',', index_col='label'
)

df_tespy = df_components.set_index('Component').rename(columns={
    'E_F [kW]': 'E_F', 'E_P [kW]': 'E_P', 'E_D [kW]': 'E_D',
})

cols = ['E_F', 'E_P', 'E_D']
idx = np.intersect1d(df_tespy.index, df_original_data.index)
df_tespy = df_tespy.loc[idx, cols]
df_original_data = df_original_data.loc[idx, cols]

df_diff_abs = (df_tespy - df_original_data).dropna()
df_diff_rel = ((df_tespy - df_original_data) / df_original_data).dropna()

result_to_markdown(df_diff_abs, 'components_delta_absolute', 'Δ ')
result_to_markdown(df_diff_rel, 'components_delta_relative', 'δ ')

# export results

# connections: specific exergy (kJ/kg) and exergy flow (kW)
df_conn_exp = df_connections.set_index('Connection').copy()
df_conn_exp.index = pd.to_numeric(df_conn_exp.index, errors='coerce')
df_conn_exp['e_PH'] = df_conn_exp['e^PH [kJ/kg]']
df_conn_exp['e_T'] = df_conn_exp['e^T [kJ/kg]']
df_conn_exp['e_M'] = df_conn_exp['e^M [kJ/kg]']
df_conn_exp['E_PH'] = df_conn_exp['m [kg/s]'] * df_conn_exp['e_PH']
df_conn_exp['E_T'] = df_conn_exp['m [kg/s]'] * df_conn_exp['e_T']
df_conn_exp['E_M'] = df_conn_exp['m [kg/s]'] * df_conn_exp['e_M']
result_to_markdown(
    df_conn_exp[['e_PH', 'e_T', 'e_M', 'E_PH', 'E_T', 'E_M']],
    'connections_result'
)

# components: E in kW, ε and y already in % from exerpy
df_comp_exp = df_components.set_index('Component').rename(columns={
    'E_F [kW]': 'E_F', 'E_P [kW]': 'E_P', 'E_D [kW]': 'E_D',
    'epsilon [%]': 'epsilon', 'y [%]': 'y_Dk', 'y* [%]': 'y*_Dk',
}).drop(index='TOT', errors='ignore')
result_to_markdown(
    df_comp_exp[['E_F', 'E_P', 'E_D', 'epsilon', 'y_Dk', 'y*_Dk']],
    'components_result'
)

# network: E_F/P/D/L in W from ean → convert to kW; epsilon is a fraction → convert to %
network_result = pd.DataFrame({
    'E_F': [ean.E_F / 1e3],
    'E_P': [ean.E_P / 1e3],
    'E_D': [ean.E_D / 1e3],
    'E_L': [ean.E_L / 1e3],
    'epsilon': [ean.epsilon * 100],
})
result_to_markdown(network_result, 'network_result')
