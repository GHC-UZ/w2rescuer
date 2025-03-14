import re
import numpy as np
import pandas as pd
import argparse

def convert_tecplot_to_vtu(input_file, output_file):
    with open(input_file, "r") as file:
        lines = file.readlines()
    
    # Extract variable names
    variable_names = lines[1].strip().split("=")[1].split()
    
    # Extract node and element counts
    zone_line = lines[2]
    num_nodes_match = re.search(r"N\s*=\s*(\d+)", zone_line)
    num_elements_match = re.search(r"E\s*=\s*(\d+)", zone_line)
    
    if not num_nodes_match or not num_elements_match:
        raise ValueError("Could not extract node and element counts from the ZONE line.")
    
    num_nodes = int(num_nodes_match.group(1))
    num_elements = int(num_elements_match.group(1))
    
    # Read node data
    node_data_start = 4  # First data line after headers
    node_data_end = node_data_start + num_nodes
    node_data = pd.read_csv(
        input_file,
        skiprows=node_data_start,
        nrows=num_nodes,
        delim_whitespace=True,
        names=variable_names,
    )
    
    # Read element connectivity
    element_data_start = node_data_end
    elements = np.loadtxt(input_file, skiprows=element_data_start, max_rows=num_elements, dtype=int)
    elements -= 1  # Convert to zero-based indexing for VTK
    
    # Write VTU file
    with open(output_file, "w") as f:
        f.write('<?xml version="1.0" ?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{num_nodes}" NumberOfCells="{num_elements}">\n')
        
        # Points
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for i in range(num_nodes):
            f.write(f'          {node_data["X"][i]} {node_data["Y"][i]} 0.0\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')
        
        # Cells
        f.write('      <Cells>\n')
        f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for element in elements:
            f.write(f'          {element[0]} {element[1]} {element[2]}\n')
        f.write('        </DataArray>\n')
        
        f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        for i in range(1, num_elements + 1):
            f.write(f'          {i * 3}\n')
        f.write('        </DataArray>\n')
        
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        f.write('          ' + ' '.join(['5'] * num_elements) + '\n')
        f.write('        </DataArray>\n')
        f.write('      </Cells>\n')
        
        # Point Data
        f.write('      <PointData>\n')
        for col in node_data.columns:
            if col not in ["X", "Y"]:
                f.write(f'        <DataArray type="Float32" Name="{col}" format="ascii">\n')
                for val in node_data[col]:
                    f.write(f'          {val}\n')
                f.write('        </DataArray>\n')
        f.write('      </PointData>\n')
        
        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')
    
    print(f"Conversion complete: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Tecplot .dat files to ParaView .vtu format.")
    parser.add_argument("input_file", help="Path to the input Tecplot .dat file")
    parser.add_argument("output_file", help="Path to the output ParaView .vtu file")
    args = parser.parse_args()
    
    convert_tecplot_to_vtu(args.input_file, args.output_file)
