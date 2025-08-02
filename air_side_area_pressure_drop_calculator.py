
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI

def calculate_air_side_area(
    tube_od_mm,
    tube_thickness_mm,
    triangular_pitch_mm,
    fin_thickness_mm,
    fpi,
    num_rows,
    face_width_m,
    face_height_m,
    air_flow_cmh,
    air_temp_C
):
    # Convert inputs
    tube_od_m = tube_od_mm / 1000
    triangular_pitch_m = triangular_pitch_mm / 1000
    fin_thickness_m = fin_thickness_mm / 1000
    fins_per_m = fpi * 39.3701
    frontal_area_m2 = face_width_m * face_height_m
    fin_depth_m = num_rows * triangular_pitch_m

    # Tubes
    tubes_per_row = math.floor(face_width_m / triangular_pitch_m)
    total_tubes = tubes_per_row * num_rows
    tube_ext_area = total_tubes * (math.pi * tube_od_m)

    # Fin areas
    fin_area_per_fin = 2 * face_width_m * fin_depth_m
    total_gross_fin_area = fin_area_per_fin * fins_per_m
    hole_area_per_tube = (math.pi / 4) * tube_od_m**2
    total_hole_area = hole_area_per_tube * total_tubes * fins_per_m
    net_fin_area = total_gross_fin_area - total_hole_area
    total_air_side_area = (tube_ext_area + net_fin_area)* face_height_m 

    # Free flow area
    fin_spacing_m = 1 / fins_per_m
    open_area_per_gap = face_width_m * (fin_spacing_m - fin_thickness_m)
    total_open_area = open_area_per_gap * fins_per_m * face_height_m
    frontal_tube_blockage = tubes_per_row * (math.pi / 4) * tube_od_m**2 * face_height_m
    net_free_flow_area = total_open_area - frontal_tube_blockage
    percent_free_area = 100 * net_free_flow_area / frontal_area_m2

    # Air flow rate and velocity
    air_flow_m3s = air_flow_cmh / 3600
    air_velocity_ms = air_flow_m3s / net_free_flow_area if net_free_flow_area > 0 else 0

    # Air properties from CoolProp
    T_K = air_temp_C + 273.15
    rho = PropsSI('D', 'T', T_K, 'P', 101325, 'Air')
    mu = PropsSI('V', 'T', T_K, 'P', 101325, 'Air')

    # Reynolds number and friction factor
    Re = rho * air_velocity_ms * tube_od_m / mu
    f = 0.25 * Re ** -0.25 if Re > 0 else 0

    # Pressure drop
    dP = (f * num_rows * rho * air_velocity_ms**2) / 2

    return {
        "Tubes per row": tubes_per_row,
        "Total tubes": total_tubes,
        "Fin depth (m)": fin_depth_m,
        "Tube external area (m²)": tube_ext_area,
        "Net fin area (m²)": net_fin_area,
        "Total air side area (m²)": total_air_side_area,
        "Free flow area (m²)": net_free_flow_area,
        "Free flow area (%)": percent_free_area,
        "Air velocity (m/s)": air_velocity_ms,
        "Air density (kg/m³)": rho,
        "Air viscosity (Pa·s)": mu,
        "Reynolds number": Re,
        "Friction factor": f,
        "Air-side Pressure Drop (Pa)": dP
    }

# Streamlit interface
st.title("Air-Side Area, Flow, and Pressure Drop Calculator")

tube_od_mm = st.number_input("Tube Outer Diameter (mm)", value=9.525)
tube_thickness_mm = st.number_input("Tube Wall Thickness (mm)", value=0.35)
triangular_pitch_mm = st.number_input("Triangular Pitch (mm)", value=25.4)
fin_thickness_mm = st.number_input("Fin Thickness (mm)", value=0.12)
fpi = st.number_input("Fins per Inch (FPI)", value=12, step=1)
num_rows = st.number_input("Number of Rows", value=4, step=1)
face_width_m = st.number_input("Coil Face Width (m)", value=1.0, step=0.0254)
face_height_m = st.number_input("Coil Face Height (m)", value=1.0, step=0.0127)
air_flow_cmh = st.number_input("Air Flow Rate (m³/h)", value=10000, step=50)
air_temp_C = st.number_input("Air Temperature (°C)", value=35.0, step=0.5)

if st.button("Calculate"):
    results = calculate_air_side_area(
        tube_od_mm, tube_thickness_mm, triangular_pitch_mm,
        fin_thickness_mm, fpi, num_rows, face_width_m,
        face_height_m, air_flow_cmh, air_temp_C
    )
    df = pd.DataFrame(list(results.items()), columns=["Parameter", "Value"])
    st.table(df)
