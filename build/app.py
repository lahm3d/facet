import os
import subprocess

import streamlit as st

TAUDEM_DIR = "/dependencies/taudem/bin"
TAUDEM_DIR_OPTIMIZED = "/dependencies/taudem_accelerated_flowDirections/taudem/build/bin"
NUM_PROC = 4


def run_taudem_command(cmd):
    """execute tauDEM commands"""
    try:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        output, err = p.communicate()
        # Get some feedback from the process to print out:
        if err is None:
            text = output.decode()
            st.write("\n", text, "\n")
        else:
            st.error(err)
    except subprocess.CalledProcessError as e:
        st.error(f"failed to return code: {e}")
    except OSError as e:
        st.error(f"failed to execute shell: {e}")
    except IOError as e:
        st.error(f"failed to read file(s): {e}")


def file_selector(folder_path='.', target="background"):

    filenames = [f for f in os.listdir(folder_path) if
                 not f[0] == "."]  # get file names from dir excluding hidden files
    selected_filename = st.selectbox(f'Select a {target}', filenames)
    abs_path = os.path.join(folder_path, selected_filename)
    if os.path.isdir(abs_path):
        return file_selector(abs_path, target)
    return os.path.join(folder_path, selected_filename)


def setup_app():

    st.title("FACET User Interface")

    """ _Pre-processing_ """
    selected_path = file_selector("/data")
    st.write(f"Selected file: {selected_path}")

    """ Run TauDEM functions """
    run_taudem = st.button("Run TauDEM D8FlowDir")
    d8_flowdir_grid = "/data/d8_flowdir.tif"
    d8_slope_grid = "/data/d8_slope.tif"
    cmd = f'mpiexec -n {NUM_PROC} {TAUDEM_DIR_OPTIMIZED}/d8flowdir -fel "{selected_path}" -p "{d8_flowdir_grid}" -sd8 "{d8_slope_grid}"'

    if run_taudem:
        st.write("Running...!")
        run_taudem_command(cmd=cmd)


if __name__ == '__main__':
    setup_app()
