# Crystal-script
Collection of script useful when working with Crystal17.

**basis_gautocry.py** converts basis set in Gaussian format as downloaded from the Basis Set Exchange (basis.gbs) to Crystal17 format (basis.d12). To run the script correctly just type the following command: `python basis_gautocry.py basis.gbs` (tested and run with python3.8).

**find_disp.py** take a Crystal17 output file (file.out) read geometries from optimization and evaluate atom displacements between the last and first step. Can be useful when an output file.gui or file.xyz is not generated and you want to use the ATOMDISP keyword for your calculations. To run the script simply run `python find_disp.py file.out` (tested and run with python3.8).

